%========================================================================
% DESCRIPTION:
% A 3D FE model of a cantilever with rectangular cross section is
% considered, which has a kind of miter / is inclined at its free end, see 
% FEmodel/FEmodel.png. A certain part of this free end is preloaded against 
% a rigid wall. The system is subjected to harmonic concentrated force in 
% bending direction, so that it undergoes recurrent frictional 
% stick-slip-liftoff interactions.
% 
% An important distinguishing feature of this example is that the 3D
% contact is distributed over an area (as opposed to a line/curve). This
% area is highlighted in FEmodel/FEmodel_contact.png and spans 12
% quadrilateral element faces. It is shown how the transform from the
% global cartesian to the contact coordinate system can be implemented in
% an ad-hoc way. Further, it is shown how a node-based quadrature of the
% contact stress can be implemented for this simple case (C3D8 elements
% having linear shape functions with uniform mesh at the contact). Finally,
% it is shown how the clamping of the cantilever can be considered within 
% MATLAB (as opposed to in the FE tool) using single-point constraints.
% 
% The required data of the FE model is in the FEmodel sub-folder. The unit
% system is N-mm-t-sec.
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

% Import FE model from files
FEMODEL = FEmodel(['FEmodel',filesep,'FEmodel']);

% Read node sets for contact, clamping and forcing,
% see FEmodel/nodeSets.png and FEmodel/nodeSet_contact.png. Only the top 2
% rows in FEmodel/nodeSet_contact.png are considered here as contact node
% set, as these are the nodes associated with the contact area highlighted
% in FEmodel/FEmodel_contact.png.
nodeSet_contact = [4 12 37 169 184:194 1029:1039];
nodeSet_constraint = load(['FEmodel',filesep,'nodeSet_constraint.nam']);
nodeSet_forcing = 428;

% Sort contact boundary DOFs to the top and apply single-point constraints 
% (remove DOFs associated with 'nodeSet_constraint'). NOTE: Of course, one 
% could avoid direct modification of properties of the FEMODEL object by 
% devising corresponding methods within the 'FEmodel' class.
IB = getIndex2NodalDOF(FEMODEL,nodeSet_contact);
ISPC = getIndex2NodalDOF(FEMODEL,nodeSet_constraint);
IREST = setdiff((1:FEMODEL.n)',[IB;ISPC]);
Ts = speye(FEMODEL.n);
Ts = Ts(:,[IB;IREST]);
FEMODEL.DOFmap = FEMODEL.DOFmap([IB;IREST],:);
FEMODEL.M = Ts'*FEMODEL.M*Ts;
FEMODEL.D = Ts'*FEMODEL.D*Ts;
FEMODEL.K = Ts'*FEMODEL.K*Ts;
FEMODEL.Fex1 = Ts'*FEMODEL.Fex1;
FEMODEL.n = size(FEMODEL.M,1);
IB = (1:length(IB))';
%% Derive massless-boundary reduced-order model; specify damping, forcing 
% properties, and response coordinate

% Apply massless Hurty-/Craig-Bampton method with nm dynamic modes
nm = 5;
typeOfROM = 'CBMB';
ROM = CMS_ROM(FEMODEL,IB,nm,typeOfROM);

% Specify mass proportional damping so that the second fixed-interface 
% mode has 5% damping ratio.
II = (ROM.nb+1):ROM.n;
om = sort(sqrt(eigs(ROM.K(II,II),ROM.M(II,II),2,'sm')));
ROM.D = 2*5e-2*om(2)*ROM.M;

% Forcing of 25 N is applied to the y-direction of node 428, at the
% frequency of the second fixed-interface mode
IFEX = getIndex2NodalDOF(FEMODEL,nodeSet_forcing,2);
FEMODEL.Fex1 = zeros(FEMODEL.n,1); FEMODEL.Fex1(IFEX) = 25;
ROM.Fex1 = ROM.T'*FEMODEL.Fex1;
Om = om(2);

% Response coordinate is that of the forcing. 'Tresp' is the corresponding 
% row of the matrix of component modes
ROM.Tresp = ROM.T(IFEX,:);

%% Define 3D contact elements with normal preload, and account for specific 
% areas associated with the nodes for node-based quadrature of the contact
% stress

% Determine number of contact nodes and check consistency with IB size
nc = length(nodeSet_contact);
if length(IB)~=3*nc
    error(['Number of boundary coordinates must be equal to 3 x number' ...
        ' of contact nodes for 3D contact.']);
end

% Determine area per contact element from a specific one (uniform mesh 
% presumed here), and also determine unit normal and tangential directions.
pos12 = getNodePosition(FEMODEL,12)';
pos169 = getNodePosition(FEMODEL,169)';
pos194 = getNodePosition(FEMODEL,194)';
rt1 = pos12-pos169;
rt2 = pos194-pos12;
rn = cross(rt1,rt2);
areaPerElement_mm2 = norm(rn);
en = rn/norm(rn); et1 = rt1/norm(rt1); et2 = rt2/norm(rt2);

% For mapping between consistent nodal forces and contact stress, in a
% node-based quadrature scheme, one needs to distinguish corner from edge
% nodes
nodeSet_contactCorner = [4 12 37 169];
% nodeSet_contactEdge = setdiff(nodeSet_contact,nodeSet_contactCorner);

% Transform to local contact coordinate system (assuming all elements in
% the same plane orthogonal to 'en')
R = [en et1 et2];
R = blkdiag(kron(eye(nc),R),eye(ROM.n-ROM.nb));
apply_transform(ROM,R);

% Set up contact elements
mu = 0.36;          % friction coefficient
pN0_MPa = 0.4;      % homogeneous initial pressure in MPa=N/mm^2
ROM.nonlinear_elements = cell(nc,1);
for inl=1:nc
    
    % Determine contact area associated to respective node
    if ismember(nodeSet_contact(inl),nodeSet_contactCorner)
        % Corner nodes get 1/4 of the contact element area
        areaPerNode_mm2 = areaPerElement_mm2/4;
    else
        % Edge nodes get 1/2 of the contact element area
        areaPerNode_mm2 = areaPerElement_mm2/2;
    end

    % Store collected information
    ROM.nonlinear_elements{inl} = struct('type', '3D',...
        'preload',pN0_MPa*areaPerNode_mm2,...
        'frictionCoefficient',mu,'area',areaPerNode_mm2);
end
% After manual definition, a check is generally advisable
check_nonlinearities(ROM);
%% Simulate

% Specify simulation duration and time discretization
period_s = 2*pi/Om;
numberOfPeriods = 10;
numberOfTimeStepsPerPeriod = 2^9;
t = linspace(0,numberOfPeriods*period_s,...
    numberOfPeriods*numberOfTimeStepsPerPeriod+1);

% Specify forcing as function handle
fex = @(t) real(ROM.Fex1*exp(1i*Om*t));

% Start from fixed-interface steady-state response
Qc = [zeros(ROM.nb,1); ...
    (-Om^2*ROM.M(II,II)+1i*Om*ROM.D(II,II)+ROM.K(II,II))\ROM.Fex1(II)];
q0 = real(Qc);
u0 = real(1i*Om*Qc);

% Call integrator with default settings
Qref = sim_contact3D(t,q0,u0,ROM,[],fex);

% Activate elimination of presumed sticking contacts, and rerun integrator
Sopt.eliminatePresumedStickingContacts = true;
Q = sim_contact3D(t,q0,u0,ROM,Sopt,fex);

%% Illustrate results

% Extract displacement of response coordinate
qresp = ROM.Tresp*Q;
qrespref = ROM.Tresp*Qref;

figure; hold on;
plot(t*Om/(2*pi),qresp,'g');
plot(t*Om/(2*pi),qrespref,'k--');
legend('presumed sticking contacts eliminated','reference');
xlabel('excitation period');
ylabel('response displacement in mm');