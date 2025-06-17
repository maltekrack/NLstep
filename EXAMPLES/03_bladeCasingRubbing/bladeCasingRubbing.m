%========================================================================
% DESCRIPTION: 
% A 3D model of a rotating compressor blade (NASA Rotor 37) subjected to 
% frictional impacts with a rigid oval casing is considered, as in 
% Section 4.3 of [1]. No abradable coating is considered, i.e., direct
% blade/casing contact is considered. Centrifugal effects are ignored.
% 
% The required data of the FE model is in the FEmodel sub-folder. The unit
% system is N-mm-t-sec.
% 
% REMARK 1: There is an error in [1] with regard to the reported
% number of time levels per revolution. The number should actually be
% multiplied by the number of revolutions (50), i.e., the largest stable
% time step actually corresponds to ca. 1000 instead of the ca. 20 time 
% levels per period reported in [1]. It is important to stress that this
% error affects both the massless and the mass-carrying method, so that
% effectively the x-axis in Fig. 11 would have to be scaled, without
% affecting the reported computation effort on the y-axis, and obviously
% without affecting the findings on the relative speedup of the computation
% effort.
% 
% REMARK 2: Two configurations with respect to the contact definition are
% given below. 'PAPER' yields results (up to numerical precision) identical
% to those in the paper. 'PLAUSIBLE' sets the contact definition up in
% a way that seems more plausible for the considered technical application
% to blade-casing rubbing (w.r.t. individual radius and consistent contact
% coordinate system, phasing of imposed gap among contact nodes, direction 
% of rotation, ...), and follows more closely [2]. 
% 
% REFERENCES
% [1] Monjaraz-Tec, C.; Gross, J.; Krack, M.: "A massless boundary 
%       component mode synthesis method for elastodynamic contact 
%       problems". Computers & Structures 260 (2022): 106698.
%       https://doi.org/10.1016/j.compstruc.2021.106698
% [2] Piollet, E.; Nyssen, F.; Batailly, A.: "Blade/casing rubbing 
%       interactions in aircraft engines: Numerical benchmark and design 
%       guidelines based on NASA rotor 37". Journal of Sound and Vibration 
%       460 (2019): 114878. https://doi.org/10.1016/j.jsv.2019.114878
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
close all;
clc;
% Add SRC and subfolders to search path. Require NLvib.
addpath(genpath(['..',filesep,'..',filesep,'SRC']));
requireNLvib;
%% Set up FE model

% Import FE model from files
bladeFEmodel = FEmodel(['FEmodel',filesep,'blade']);
% NOTE: The root of the blade is constrained, and this boundary condition 
% is already considered in the provided files; i.e., the corresponding 
% nodal degrees of freedom have been removed from mass and stiffness
% matrices, and the DOF map.

% Set and number of contact nodes (see FEmodel/nodeSet_contact.png)
nodeSet_contact = [188,190,192,194,196,198,200,203,253,255,...
    257,259,261,263,265,267]; 
nc = length(nodeSet_contact);

% Get index to contact boundary degrees of freedom
IB = getIndex2NodalDOF(bladeFEmodel,nodeSet_contact);
if length(IB)~=3*nc
    error(['Number of boundary coordinates must be equal to 3 x number' ...
        ' of contact nodes for 3D contact.']);
end
%% Derive massless-boundary reduced-order model; specify response 
% coordinate, and damping

% Apply MacNeal's component mode synthesis method with nm dynamic modes
nm = 50;
typeOfROM = 'MCN';
bladeROM = CMS_ROM(bladeFEmodel,IB,nm,typeOfROM);

% Response coordinate is y-translation of node 1352
isens = getIndex2NodalDOF(bladeFEmodel,1353,2);
% 'Tresp' is the corresponding row of the matrix of component modes
% It should be treated as property of the CMS_ROM object to account for
% possible coordinate transforms applied later.
bladeROM.Tresp = bladeROM.T(isens,:);

% Free-interface modes receive 0.5% damping ratio
modalDampingRatio = 0.5e-2;
om = sort(sqrt(eigs(bladeROM.K,bladeROM.M,nm,'sm')));
bladeROM.D = blkdiag(zeros(size(IB,1)),diag(2*modalDampingRatio*om));

%% Define 3D contact elements with imposed gap and velocity
mu = 0.15;              % friction coefficient
gmean_mm = 0.356;       % imposed mean gap
gampl_mm = 0.37;        % imposed gap amplitude
CONFIG = 'PLAUSIBLE';   % ['PAPER'|'PLAUSIBLE']; see comment in header
R = cell(nc,1);         % rotation matrix local to global
bladeROM.nonlinear_elements = cell(nc,1);
for inl=1:nc
    switch CONFIG
        case 'PAPER'
            % REMARK 3: The results in [1] correspond to the angular 
            % velocity specified below. The value specified in the text 
            % of [1] is by factor 2 different.
            Omrot_rad_s = 2*pi*312/4;

            % Radius and phasing as in the paper
            pos = getNodePosition(bladeFEmodel,nodeSet_contact(1));
            radius_mm = pos(1);
            phi_rad = 2*pi/50*(inl-1)/(nc-1);

            % In this case, no coordinate transform was applied for
            % simplicity
            R{inl} = eye(3);

            % Define imposed (normal) gap as scalar
            imposedGap = @(t) gmean_mm + gampl_mm*...
                cos(2*(Omrot_rad_s*t(:).'+phi_rad));
            % Define imposed gap velocity as vector
            imposedGapVelocity = @(t) [ ...
                -2*Omrot_rad_s*gampl_mm*...
                sin(2*(Omrot_rad_s*t(:).'+phi_rad)); ...
                radius_mm*Omrot_rad_s*ones(1,length(t)); ...
                zeros(1,length(t))];
        case 'PLAUSIBLE'
            % For the considered casing geometry, this rotor speed should
            % lead to a gap oscillating with fundamental frequency close
            % to that of the first blade mode.
            Omrot_rad_s = 1344;

            % Extract cartesian coordinates of the blade. The third
            % coordinate is aligned with the rotation axis. Determine 
            % radius and angle for the imposed gap and velocity.
            pos = getNodePosition(bladeFEmodel,nodeSet_contact(inl));
            r = pos(1)+1i*pos(2);
            radius_mm = abs(r);
            phi_rad = angle(r);

            % Determine unit normal and tangential directions (positive
            % displacement in normal direction increases gap in this tool),
            % and set up transformation matrix from global cartesian to
            % local contact coordinate system
            en = [-pos(1);-pos(2);0]; en = en/norm(en);
            et1 = [0;0;1];
            et2 = cross(en,et1);
            R{inl} = [en et1 et2];

            % Define imposed (normal) gap as scalar
            % For the specific function used here, cf. [2], Eqs. (3)-(4).
            theta = @(t) Omrot_rad_s*t(:).' + phi_rad;
            imposedGap = @(t) gmean_mm*(1-...
                2*exp(-((mod(theta(t),pi)/pi-.5)/.15).^2));
            % Define imposed gap velocity as vector
            imposedGapVelocity = @(t) [ ...
                4*Omrot_rad_s*gmean_mm/(pi*.15)*(mod(theta(t),pi)/pi-.5)/.15.*...
                exp(-((mod(theta(t),pi)/pi-.5)/.15).^2); ...
                zeros(1,length(t)); ...
                -radius_mm*Omrot_rad_s*ones(1,length(t))];
    end

    % Add nonlinear contact element
    bladeROM.nonlinear_elements{inl} = struct('type', '3D',...
        'imposedGap',imposedGap,...
        'imposedGapVelocity',imposedGapVelocity,...
        'frictionCoefficient',mu);
end
% After manual definition, a check is generally advisable
check_nonlinearities(bladeROM);

% Apply coordinate transform to nb=3*nc boundary coordinates
R = blkdiag(R{:},eye(bladeROM.n-3*nc));
apply_transform(bladeROM,R);

%% Simulate

% Specify homogeneous initial conditions
q0 = zeros(bladeROM.n,1); u0 = q0;

% Specify simulation duration and time discretization
numberOfRevolutions = 50;
numberOfTimeStepsPerRevolution = 1e3;
timePerRevolution_s = 2*pi/Omrot_rad_s;
t_s = linspace(0,numberOfRevolutions*timePerRevolution_s,...
    numberOfRevolutions*numberOfTimeStepsPerRevolution+1);

% Call integrator
[Q,~,LAM] = sim_contact3D(t_s,q0,u0,bladeROM);

%% Illustrate results similar to Figure 10 in [1]
revolutionNumber = t_s/timePerRevolution_s;
qR = bladeROM.Tresp*Q; % REMARK 4: -qR is actually plotted in the paper.

% OVERVIEW
figure; hold on;
plot(revolutionNumber,qR);
xlabel('revolution number'); ylabel('qR in mm');

% ZOOM INTO INITIAL PHASE
figure; hold on;
plot(revolutionNumber,qR);
xlabel('revolution number'); ylabel('qR in mm');
set(gca,'Xlim',[0 5]);

% ZOOM INTO STEADY-STATE PHASE
figure; hold on;
plot(revolutionNumber,qR);
xlabel('revolution number'); ylabel('qR in mm');
set(gca,'Xlim',[45 50]);

%% Export animation

% Truncate to last half-revolution and down-sample to reduce file size
idxANIMATE = find(revolutionNumber>revolutionNumber(end)-.5);
idxANIMATE = idxANIMATE(1:10:end);

% Expand from ROM to FE model
qFE = bladeROM.T*Q(:,idxANIMATE);

% Animate FE model
animate(bladeFEmodel,qFE);
% You can open the animation file with CalculiX, which is an open source FE
% tool available at https://www.dhondt.de/.