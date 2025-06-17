%========================================================================
% DESCRIPTION: 
% A model of the test rig from [1] is considered. It consists of a 
% cantilevered beam with a curved plate bolted to its free end, which is
% pressed against a flat plate. The beam is harmonically driven by a
% shaker. This leads to periodic stick-slip, mimicking the generation of
% squeak noise relevant e.g. in the automotive industry.
%
% A massless Craig-Bampton ROM is provided, which was derived from a 3D FE
% model (provided by authors of [1], not publicly available). The ROM has 5
% normal modes and 6 static constraint modes associated with the relative 
% displacements at the 2 vertex nodes of two quadrilateral surface elements
% (linear shape functions) highlighted in [1], Fig. 4b. The ROM also 
% contains the amplitude vector of the excitation force 'Fex1', and the row
% vector 'Tresp' which recovers the response in the y-direction at the
% laser vibrometer's measurement point ([1], Fig. 3). The unit system is 
% N-mm-t-sec.
% 
% REFERENCES
% [1] Utzig, L.; Weisheit, K.; Sepahvand, K.; Marburg, S.: "Innovative 
%   squeak noise prediction: An approach using the harmonic balance method 
%   and a variable normal contact force". Journal of Sound and 
%   Vibration 501 (2021): 116077. https://doi.org/10.1016/j.jsv.2021.116077
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
%% Load massless-boundary reduced-order model; determine number of contacts
load('ROM.mat','ROM');
nc = ROM.nb/3;
%% Define 3D contact elements
mu = 0.36;          % friction coefficient
lam0_N = 13.6/2;    % normal preload at each node in N
ROM.nonlinear_elements = cell(nc,1);
for inl=1:nc
    ROM.nonlinear_elements{inl} = struct('type','3D',...
        'preload',lam0_N,'frictionCoefficient',mu);
end
% After manual definition, a check is generally advisable
check_nonlinearities(ROM);
%% Define damping, forcing, and initial conditions

% Specify 1 % damping ratio for the fixed-interface normal modes. The type
% of ROM (massless Craig-Bampton) and knowledge of the corresponding 
% structure of 'ROM.K' is exploited here (it has the angular frequencies of
% the fixed-interface normal modes on the diagonal associated with the
% inner coordinates).
II = (ROM.nb+1):ROM.n;
om = sqrt(diag(ROM.K(II,II)));
ROM.D = blkdiag(zeros(ROM.nb),diag(2*1e-2*om));

% The excitation frequency is 7 Hz. The excitation force is 7.8 N, which is
% considered ROM.Fex1. The force is ramped using a half-cosine ramp of 
% 0.5 periods duration.
Om_rad_s = 7*2*pi;
period_s = 2*pi/Om_rad_s; 
fex = @(t) (( 1 - cos(pi*t/(0.5*period_s)) )/2.*double(t<0.5*period_s)+...
    double(t>=0.5*period_s)) .* ...
    real(ROM.Fex1*exp(1i*Om_rad_s*t));

% Specify homogeneous initial conditions
q0 = zeros(ROM.n,1); u0 = q0;

%% Simulate treating contact as non-regularized ('rigid')

% Specify simulation duration and time discretization
numberOfPeriods = 5;
numberOfTimeStepsPerPeriod = 2^10;
t_s = linspace(0,numberOfPeriods*period_s,...
    numberOfPeriods*numberOfTimeStepsPerPeriod+1);

% Call integrator
Sopt.eliminatePresumedStickingContacts = true;
[Q,U,LAM] = sim_contact3D(t_s,q0,u0,ROM,Sopt,fex);

%% Simulate treating contact via penalty regularization

% Set parameters of time step integration scheme.
% The asymptotic spectral radius is set = 0 here, leading to the most
% severe numerical dissipation. This is needed for the high penalty
% coefficient!
Sopt.scheme = struct('name','generalizedAlpha','rhoinfty',0,...
    'verbosity',3);

% Study the effect of the penalty stiffness, considering two different
% values
penaltyStiffness_N_mm = [2e2 2e4];
Uregular = cell(length(penaltyStiffness_N_mm),1);
for ii=1:length(Uregular)

    % Set penalty stiffnes in each contact element
    for inl=1:length(ROM.nonlinear_elements)
        ROM.nonlinear_elements{inl}.stiffness = penaltyStiffness_N_mm(ii);
    end
    
    % Call integrator for regularized contact
    [~,Uregular{ii}] = sim_regularNL(t_s,q0,u0,ROM,Sopt,fex);

end

%% Illustrate results
figure; hold on;
plot(Om_rad_s*t_s/(2*pi),ROM.Tresp*U,'k-');
plot(Om_rad_s*t_s/(2*pi),ROM.Tresp*Uregular{1},'b');
plot(Om_rad_s*t_s/(2*pi),ROM.Tresp*Uregular{2},'g--');
legend('reference','penalty, soft','penalty, stiff');
xlabel('excitation period');
ylabel('response velocity in mm/s');
set(gca,'Xlim',[2 5]);