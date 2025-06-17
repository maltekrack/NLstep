%========================================================================
% DESCRIPTION:
% A cantelevered beam is considered, which is subjected to dry Coulomb
% friction (constant limit force) at a certain location and driven by 
% harmonic forcing near primary resonance with the fundamental bending
% mode. Euler-Bernoulli theory is used and a finite element model with two 
% degrees of freedom per node (1 translation, 1 rotation) is used as point 
% of departure. The MacNeal method is then applied to reduce the model 
% order and obtain a massless contact boundary coordinate.
% The following simulations are carried out: First, the transient frequency 
% sweep through resonance is simulated using time step integration. Second, 
% time step integration is done at a fixed near-resonant frequency (dwell).
% Finally, Harmonic Balance is applied to compute the periodic steady-state
% frequency response, using a Dynamic Lagrangian approach to treat the
% set-valued frictional contact law.
%========================================================================
% This file is part of NLstep.
% 
% If you use NLstep, please refer to the article:
% "A massless boundary component mode synthesis method for elastodynamic
% contact problems", C.D. Monjaraz-Tec, J. Gross, M. Krack:
% https://doi.org/10.1016/j.compstruc.2021.106698
%
% COPYRIGHT AND LICENSING: 
% NLstep Copyright (C) 2025  
%   Malte Krack  (malte.krack@ila.uni-stuttgart.de) 
%   Johann Gross (johann.gross@ila.uni-stuttgart.de)
%   University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLstep is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
clearvars;
close all;
clc;
% Add SRC and subfolders to search path. Require NLvib.
addpath(genpath(['..',filesep,'..',filesep,'SRC']));
requireNLvib;
%% Set up FE model

% Geometry and material parameters; intended unit system is N-mm-t-sec
len = 2e3;              % length
height = .05*len;       % height in the bending direction
thickness = 3*height;   % thickness in the third dimension
E = 185e3;              % Young's modulus
rho = 7.83e-9;          % mass density

% Construct FE model object and account for clamped-free boundary
% conditions
numberOfNodes = 9;
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,'clamped-free',...
    numberOfNodes);

% Apply friction element at node 4 in translational direction, with a limit
% force of 1500 N
iNodeNL = 4;
add_nonlinear_attachment(beam,iNodeNL,'trans','friction',...
    'friction_limit_force',1500);
% Store index of coordinate where nonlinear element is applied
INL = find_coordinate(beam,iNodeNL,'trans');

% Apply forcing to free end of beam in translational direction, with
% magnitude 200 N
iNodeFex = numberOfNodes;
add_forcing(beam,iNodeFex,'trans',200);

% Store index of response coordinate (translation at drive point)
IRESP = find_coordinate(beam,iNodeFex,'trans');

%% Derive massless-boundary reduced-order model; specify response 
% coordinate, damping and forcing properties

% Apply MacNeal method with nm dynamic modes, using INL as boundary
nm = 5;
ROM = CMS_ROM(beam,INL,nm,'MCN');

% Define response vector (associated column of matrix of component modes)
ROM.Tresp = ROM.T(IRESP,:);

% Specify mass proportional damping so that the fundamental fixed-interface
% (sticking) mode has 5% damping ratio
II = (ROM.nb+1):ROM.n;
omFixed = sort(sqrt(eigs(ROM.K(II,II),ROM.M(II,II),1,'sm')));
ROM.D = 2*5e-2*omFixed*ROM.M;

% Specify start and end frequency of sweep
omFree = sort(sqrt(eigs(ROM.K,ROM.M,1,'sm')));
OmS = 1.6*omFixed;
OmE = .9*omFree;

%% Determine steady-state response in the linear case of sticking contact
OmLIN = linspace(OmS,OmE,1e2);
QLIN = zeros(ROM.n,length(OmLIN));
for iOm=1:length(OmLIN)
    QLIN(II,iOm) = ( -OmLIN(iOm)^2*ROM.M(II,II) + ...
        1i*OmLIN(iOm)*ROM.D(II,II) + ROM.K(II,II)  )\...
        ROM.Fex1(II);
end
%% Simulate nonlinear transient frequency sweep

% Specify simulation duration and time discretization
period_s = 2*pi/((OmS+OmE)/2);
numberOfPeriods = 1500;
numberOfTimeStepsPerPeriod = 2^9;
t = linspace(0,numberOfPeriods*period_s,...
    numberOfPeriods*numberOfTimeStepsPerPeriod+1);

% Specify forcing as function handle; sweep frequency at constant rate
Om = @(t) OmS + (OmE-OmS)/(numberOfPeriods*period_s)*t;
th = @(t) OmS*t + (OmE-OmS)/(numberOfPeriods*period_s)*t.^2/2;
fex = @(t) real(ROM.Fex1*exp(1i*th(t)));

% Start from fixed-interface steady-state response
q0 = real(QLIN(:,1));
u0 = real(1i*OmLIN(1)*QLIN(:,1));

% Call integrator
[QSWEEP,USWEEP] = sim_friction1D(t,q0,u0,ROM,[],fex);
OmSWEEP = Om(t);

%% Compute periodic steady-state frequency response using HB

% Set harmonic order (H), harmonic set, and number of samples per period 
% within alternating frequency-time scheme (N)
H = 50;
Hset = 1:2:H; % only odd harmonics are considered in this example
numH = numel(Hset);
N = 2^11;

% Specify the Dynamic Lagrangian parameter epsDL. This should not affect the
% harmonically converged results. We set epsDL to the smallest eigenvalue of
% the linear dynamic stiffness matrix, which we expected to provide a good
% tradeoff between convergence speed and robustness.
ROM.nonlinear_elements{1}.epsDL = abs(eigs(...
    -omFixed^2*ROM.M(II,II) + 1i*omFixed*ROM.D(II,II) + ROM.K(II,II), ...
    1, 'sm' ));

% Solve HB equations and continue (sequentially with fixed step size) 
% w.r.t. Om, using linear initial guess
x0 = [real(QLIN(:,1));-imag(QLIN(:,1));zeros(ROM.n*2*(numH-1),1)];
ds = abs(OmS-OmE)/30;
Sopt = struct('flag',0,'stepadapt',0,...
    'Dscale',[max(abs(QLIN(:,1)))*ones(size(x0));(OmS+OmE)/2]);
X = solve_and_continue(x0,...
    @(X) HB_residual_friction1D(X,ROM,H,N),...
    OmS,OmE,ds,Sopt);

% Interpret solver output
OmHB = X(end,:);
n = ROM.n;
ID = repmat(1:n,1,numH)+n*kron(Hset,ones(1,n));
IC = repmat(1:n,1,numH)+n*kron(0:2:2*(numH-1),ones(1,n)); IS = IC+n;
QHB = zeros(n*(H+1),size(X,2));
QHB(ID,:) = X(IC,:)-1i*X(IS,:);

% Set up inverse discrete Fourier transform and evaluate amplitude measure
tau = (0:2*pi/N:2*pi-2*pi/N)';
H_iDFT = exp(1i*tau*(0:H));
QrespHB = kron(eye(H+1),ROM.Tresp)*QHB;
arespHB = max(real(H_iDFT*QrespHB));

%% Simulate nonlinear response at fixed frequency (amplitude peak)

% Determine frequency at amplitude maximum
[~,iMAX] = max(arespHB);
OmMAX = OmHB(iMAX);

% Redefine excitation to have fixed frequency
[~,itmp] = min(abs(Om(t)-OmMAX));
thMAX = th(t(itmp));
th = @(t) thMAX + OmMAX*t;
fex = @(t) real(ROM.Fex1*exp(1i*th(t)));

% Use corresponding sweep response point for initial values
q0 = QSWEEP(:,itmp);
u0 = USWEEP(:,itmp);

% Redefine time vector and call integrator for 'numberOfPeriods' periods
numberOfPeriods = 30;
t = linspace(0,numberOfPeriods*2*pi/OmMAX,...
    numberOfPeriods*numberOfTimeStepsPerPeriod+1);
[QDWELL,UDWELL] = sim_friction1D(t,q0,u0,ROM,[],fex);

%% Illustrate results

% Amplitude-frequency curve
figure; hold on;
plot(OmSWEEP,ROM.Tresp*QSWEEP,'k');
plot(OmHB,arespHB,'g-','linewidth',2);
plot(OmLIN,abs(ROM.Tresp*QLIN),'-','color',.75*[1 1 1]);
legend('time integration','HBDL','linear');
xlabel('excitation frequency in rad/s');
ylabel('tip displacement amplitude in mm');

% Phase projection
figure; hold on;
qHB = real(H_iDFT*QrespHB(:,iMAX));
uHB = real(H_iDFT*(1i*OmMAX*(0:H)'.*QrespHB(:,iMAX)));
plot(qHB,uHB,'g-','linewidth',6);
indLastTwoPeriods = t>(t(end)-2*2*pi/OmMAX); % show only last 2 periods
plot(ROM.Tresp*QDWELL(:,indLastTwoPeriods),...
    ROM.Tresp*UDWELL(:,indLastTwoPeriods),'k--');
legend('HBDL','time integration');
xlabel('tip displacement in mm'); 
ylabel('tip velocity in mm/s');