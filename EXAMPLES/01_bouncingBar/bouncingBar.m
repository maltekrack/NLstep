%========================================================================
% DESCRIPTION: 
% The common bouncing bar problem is considered, where an elastic bar drops
% from a certain height onto the rigid ground. For the specified
% parameters, the bar bounces periodically.
% 
% A 1D finite element model with equidistant nodes and linear shape 
% functions is used as point of departure. The MacNeal method is applied, 
% and the reduced model is integrated in time using a leapfrog 
% scheme, where the interaction with the rigid ground is treated as 
% unilateral constraint enforced on displacement level. Similar simulations
% are shown in [1]. As reference, the exact solution of the spatially  
% continuus problem from [2] is used.
% 
% REFERENCES
% [1] Monjaraz-Tec, C.; Gross, J.; Krack, M.: "A massless boundary 
%       component mode synthesis method for elastodynamic contact 
%       problems." Computers & Structures 260 (2022): 106698.
%       https://doi.org/10.1016/j.compstruc.2021.106698
% [2] Doyen, David, Alexandre Ern, and Serge Piperno. "Time-integration 
%       schemes for the finite element dynamic Signorini problem." SIAM 
%       Journal on Scientific Computing 33.1 (2011): 223-249.
%       https://doi.org/10.1137/100791440
%========================================================================
% This file is part of NLstep.
% 
% If you use NLstep, please refer to the article [1].
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
clc;
close all;
% Add SRC and subfolders to search path. Require NLvib.
addpath(genpath(['..',filesep,'..',filesep,'SRC']));
requireNLvib;
%% Parameters of the problem setting
% The parameters are adopted from [2], the unit system is unknown.
ell = 10;   % length of bar
E = 900;    % Young's modulus of bar
rho = 1;    % mass density of bar
ag = 10;    % imposed acceleration (gravity)
% NOTE: The initial height is set below along with the exact solution.
%% Parameters of the simulation
numberOfPeriods = 10;                   % number of periods to simulate
dx = .5e-2;                             % element length
CourantNumber = 3;                      % Courant number
dt = CourantNumber*dx/(sqrt(E/rho));    % time step
%% Exact solution of continuous problem [2]

% Account for different naming convention:
g0 = ag; L = ell;

% Determine wave speed and set initial height in such a way that periodic
% behavior is obtained.
c0 = sqrt(E/rho); % wave speed
tauw = L/c0; % time it takes a wave to cross the bar
% tauf = sqrt(2*h0/ag); % time it takes the bar to reach the ground
% The initial height is set in such a way that one has periodic behavior.
% This is achieved here by setting the time it takes to reach the ground
% 3x the time it takes the wave to cross the bar
tauf = 3*tauw;
% From this condition, one can follow for the initial height:
h0 = tauf^2*ag/2;
T = 16*tauw; % period length

% Determine number of steps per period and time vector. NOTE: As the period
% 'T' is needed for this, it cannot be done before. The same time vector is
% adopted for the numerical simulation.
numberOfTimeStepsPerPeriod = round(T/dt);
t = linspace(0,numberOfPeriods*T,...
    numberOfPeriods*numberOfTimeStepsPerPeriod+1);

% Set up exact solution, following closely the notation in [2], 
% Section 3.1-3.2, Eq. (28)-(31).
% Fundamental parameters of the solution:
n = (1:1000)'; % modal truncation
lamn = n.*pi/L;
nun = (n-0.5).*pi/L;
an = -2*g0./(c0^2*L*((n-0.5).*pi/L).^3);
bn = 4*g0./(c0^2*lamn.^2);
vf = sqrt(2*h0*g0);
% Expressions for piecewise solution at the bottom of the bar (x=0):
x = 0;
P = @(t) h0 - 1/2*g0*(t-tauf).^2;
S1 = @(t) (an.*sin(nun*x))'*(1-cos(c0*nun*t));
S2 = @(t) -((2*g0*L^2)/(3*c0^2)) +...
    (bn.*cos(lamn*x))'*cos(c0*lamn*t);
H = @(t) -vf.*min(x/c0,tauw-abs(t-tauw));
% Transition points between pieces/phases within first period:
qbref = zeros(length(t),1);
t1 = tauf;
t2 = t1+2*tauw;
t3 = t2+2*tauf;
t4 = t3+2*tauw;
% Phase 1: free fall of undeformed bar, zero initial velocity
qbref(t<t1) = P(t(t<t1)+tauf);
% Phase 2: first impact, propagating wave
qbref(t1<=t&t<t2) = H(t(t1<=t&t<t2)-t1) + S1(t(t1<=t&t<t2)-t1);
% Phase 3: take-off; free flight + propagating wave
qbref(t2<=t&t<t3) = P(t(t2<=t&t<t3)-t2) + S2(t(t2<=t&t<t3)-t2);
% Phase 4: second impact, superimposing waves cancel at the end
qbref(t3<=t&t<t4) = H(t4-t(t3<=t&t<t4)) + S1(t4-t(t3<=t&t<t4));
% Phase 5: free flight as rigid body
qbref(t4<=t&t<=(t4+t1)) = P(t(t4<=t&t<=(t4+t1))-t4);
% Repeat periodically
qbref = [repmat(qbref(0<=t&t<T),numberOfPeriods,1);qbref(t==T)];
%% Set up FE model

% Specify number of nodes and use NLvib's Mechanical System class to set up
% the model.
numberOfNodes = round(ell/dx);
barFEmodel = FE_ElasticRod(ell,1,E,rho,'free-free',numberOfNodes);

% Apply unilateral contact element at node 1:
inode = 1;
add_nonlinear_attachment(barFEmodel,inode,'unilateral','imposedGap',h0);

% Apply gravity load:
barFEmodel.Fex1 = -barFEmodel.M*ones(size(barFEmodel.M,1),1)*ag;
%% Derive massless-boundary reduced-order model

% Specify number of dynamic modes and extract boundary coordinate:
Nm = 15;
ib = find_coordinate(barFEmodel,inode);

% Apply MacNeal's component mode synthesis method (leading to a massless
% boundary)
typeOfROM = 'MCN';
barROM = CMS_ROM(barFEmodel,ib,Nm,typeOfROM);
%% Numerical simulation

% Specify homogeneous initial conditions
q0 = zeros(length(barROM.M),1);
u0 = q0;

% Specify (time-constant) gravity load:
funcExc = @(t) barROM.Fex1;

% Run simulation
[Q,U,FNL,CSTATE] = sim_unilateralContactDisplacementLevel(t,q0,u0,barROM,[],...
    funcExc);

% Shift the boundary displacement by the initial height (as in the exact
% solution):
qb = Q(1,:) + h0;
%% Illustrate results

% FULL SIMULATED TIME SPAN
figure; hold on;
plot(t/T,qbref/h0,'k-');
plot(t/T,qb/h0,'g--');
legend('reference',typeOfROM);
xlabel('number of periods');
ylabel('qb/q0');
set(gca,'XLim',[0,numberOfPeriods],'YLim',[-.1,1.19])

% INITIAL TWO PERIODS
figure; hold on;
plot(t/T,qbref/h0,'k-');
plot(t/T,qb/h0,'g--');
legend('reference',typeOfROM);
xlabel('number of periods');
ylabel('qb/q0');
set(gca,'XLim',[0,2],'YLim',[-.1,1.19]);