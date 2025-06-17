%========================================================================
% Matlab function for time step integration of massless-bounday models 
% subjected to frictional-unilateral contact (set-valued Coulomb-Signorini 
% conditions) and possible forcing with explicit time dependence, from a 
% given set of initial values of displacement and velocities, in 
% accordance with [1]. The model must be an object of the class CMS_ROM 
% (part of the tool NLvib [3]), which provides appropriate component mode 
% synthesis methods.
% 
% DYNAMIC FORCE BALANCE
% The following differential-algebraic equation system is treated:
%
%   K_bb * q_b + K_bi * q_i - lambda = fex_b(t),
%   M_ii * \dot u_i + D_ii * u_i + K_ii * q_i + K_ib * q_b = fex_i(t),
%   u_b = \dot q_b; u_i = \dot q_i.
%
% The first line is the static force balance at the boundary, the second
% line is the dynamic internal force balance, and the last line just gives
% the coordinate-velocity relation. Herein, \dot denotes derivative with 
% respect to time t, q_b and q_i are the vectors of contact boundary and 
% remaining coordinates, respectively. M_ii, K_ii, K_ib, K_bi, K_bb, D_ii 
% are the respective partitions of the symmetric mass, damping and 
% stiffness matrices. fex_b, fex_i are vectors of external forces with 
% known explicit time dependence. lambda is the vector of contact forces. 
% More specifically, it is the dynamic part of the contact force, while the 
% total one is this lambda plus a possible preload. Such a preload may 
% arise when the above CMS_ROM was obtained after applying a (possibly 
% nonlinear) static analysis step and subsequent linearization. Then the 
% coordinates also count from this reference configuration (static 
% equilibrium).
% 
% CONTACT KINEMATICS, NODE-BASED QUADRATURE OF CONTACT STRESS
% It is presumed that the model is formulated in the local relative contact 
% coordinates, so that the gap and gap velocity vectors are simply
%
%   g = q_b + g_imp(t),
%   gamma = u_b + gamma_imp(t),
%
% where g_imp/gamma_imp are imposed gap/velocity with known explicit time
% dependence. The vectors are grouped in such a way that the normal 
% coordinate is followed by the two orthogonal tangential ones, and nc 
% of such triplets are stacked, where nc denotes the number of contacts. 
% WHATCH OUT: In the normal direction, an increasing coordinate means 
% opening contact! 
% If contact is to be modeled between two finite element interfaces, it is 
% presumed that an appropriate transformation to relative displacements has 
% been applied (which implies a small sliding assumption). To this end, one 
% side of the contact interface is to be selected as dependent and the 
% other as independent. The (absolute) displacements of the independent 
% nodes are included into the set of remaining coordinates. Interpolation 
% is necessary on the independent side to treat non-conforming meshes. This 
% corresponds to a node-to-surface formulation. 
% A consistent node-based integration of the contact stress can be done. 
% Here, the contact stress is evaluated at each dependent node, and the 
% consistent nodal force is obtained by multiplying this with the area 
% associated to the respective node. If you do not specify this area, a 
% unity area will be used.
% 
% TIME STEP INTEGRATION, CONTACT TREATMENT
% An event-capturing approach in conjuction with the leapfrog integration
% scheme is implemented. Here, coordinate and velocity grid are shifted by
% a half-step. The dynamic internal force balance is stepped in an explicit
% way, while the static force balance at the boundary is treated in an
% implicit way. To this end, contact kinematics, integration
% scheme, and force balance at the boundary, evaluated at a certain time
% level, are solved for gamma as function of lambda, which yields
%
%       gamma = G*lambda + c,
%
% where G is denoted Delassus matrix.
% The Coulomb-Signorini conditions are expressed as projective non-smooth 
% equation system, relating g, gamma and lambda,
%
%       lambda - projC(lambda - epsAL.*gamma) = 0,
%
% where projC denotes the projection onto the admissible set of contact
% forces, and epsAL>0 can be interpreted as Augmented Lagrangian parameter.
% More specifically, this projective equation system is to be satisfied for
% the closed contacts. Together with the linear relation between gamma and
% lambda, one obtains an implicit equation system for lambda. This is
% iteratively solved using projected Jacobi over-relaxation. Notation and
% algorithmic treatment of the contact follow closely [2], including the
% selection of the important parameter epsAL.
% 
% SETTING THE (CONTACT) PARAMETERS
% By requiring that the input model is of massless-boundary type within
% the NLvib's class CMS_ROM, some checks have already been made and some
% properties of the model are already known. In particular, this applies to
% the mass, damping and stiffness matrices (MBmodel.M, MBmodel.D,
% MBmodel.K). The contact elements are specified manually on the other hand
% and the associated parameters are checked in this function. Recall that 
% contact is modeled in terms of dry Coulomb friction in the tangential
% contact plane, and unilateral interaction in the normal contact
% direction. The parameters are specified in MBmodel.nonlinear_elements,
% which is a cell array of structures. The number of elements must equal
% to nc, the number of contacts, and the sorting must correspond to the
% sorting of the coefficient matrices (there is no resorting or coordinate
% transform done in this function!).
% For contact element nl, the contact parameters are specified as fields of
% MBmodel.nonlinear_elements{nl}, with the names type, frictionCoefficient, 
% imposedGap, imposedGapVelocity, preload, and area. The names should be
% self-explanatory. The field type must be '3D' for this function. 
% Otherwise, the friction coefficient is the only mandatory field,
% while the rest is set to zero if not specified, and area=1 is used on
% default. Imposed gap and imposed gap velocity can have explicit time 
% dependence (function handle). In the present version of the tool, only a 
% normal preload is handled, whereas tangential preload is not allowed. For
% a given contact element, you have to specfiy the imposedGap or the
% preload (exactly ONE of those).
% NOTE: The excitation force vector within the model MBmodel is ignored in
% this function. Nonzero external forcing must be specified as function 
% handle (func_fex=varargin{2} input).
% 
% ELIMINATION OF PRESUMED STICKING CONTACTS
% As stated above, the projective non-smooth equation system holds for 
% the non-open contacts only, and must be restricted accordingly. Sticking 
% contacts have known gap velocity zero and can hence be eliminated from 
% the equation system. Of course, it is not known exactly a priori whether
% a contact is sticking. But depending on kinematic and kinetic
% considerations, a set of presumed sticking contacts can be identified.
% The benefit of eliminating those from the equation system is a smaller 
% problem dimension for the projected Jacobi over-relaxation. The drawback 
% is that the iteration matrix must be recomputed every time step where 
% the set of presumed sticking and/or active contacts changes, and the
% expression for this matrix involves a Schur complement which may
% be expensive to compute if the number of active contacts is
% large. The elimination of presumed sticking contacts is switched
% on by specifying Sopt.eliminatePresumedStickingContacts = true.
%
% REFERENCES
% [1] Monjaraz-Tec, C.; Gross, J.; Krack, M.: "A massless boundary
%       component mode synthesis method for elastodynamic contact
%       problems." Computers & Structures 260 (2022): 106698.
%       https://doi.org/10.1016/j.compstruc.2021.106698
% [2] Studer, Christian. Numerics of unilateral contacts and friction:
%       modeling and numerical time integration in non-smooth dynamics.
%       Vol. 47. Springer Science & Business Media, 2009.
%       https://doi.org/10.1007/978-3-642-01100-9
% [3] https://www.ila.uni-stuttgart.de/nlvib/
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
%
% Mandatory input variables:
%   t                       time vector with equidistant time samples; can
%                           also be specified as two-element vector with
%                           t(1) = dt, t(2) = tend, so that the actual time
%                           vector is set up in this function as 
%                           t = 0:dt:tend
%   q0                      initial coordinate vector
%   u0                      initial velocity vector
%   MBmodel                 massless boundary model (class CMS_ROM)
%
% Optional input variables:
%   varargin{1} = Sopt      Options structure; only specified options 
%                           are changed from their default values
%       eliminatePresumedStickingContacts
%                               flag whether presumed sticking contacts are
%                               to be eliminated from the projective
%                               non-smooth equation system (default is
%                               false)
%       tol                     scalar convergence tolerance; residuum 
%                               is measured in terms of the L-infinity 
%                               norm of the difference between the vectors 
%                               of contact forces ('lam') between the last 
%                               two iterations; this is required to be
%                               below 'tol' times the L-infinity norm 
%                               of 'lam'
%       epsALscl                scaling to increase/decrease Augmented 
%                               Lagrangian parameter from the one 
%                               recommended in [2]
%       itermax                 maximum number of projected Jacobi
%                               over-relaxation iterations per time step
%       dataStorageFrequency    set this to n>1 to store data only every n
%                               time levels
%       cmdWindowOutputFrequency    ... (as data storage but for command
%                                   window output)
%       verbosity               0: none; 1: basic; 2: iter; 3: detailed
%   varargin{2} = func_fex  function handle to external forcing with known
%                           explicit time dependence, fex(t)
%
% Output variables:
%   Q, U, LAM               matrices whose colums are the coordinates,
%                           velocities, and contact forces at the time
%                           instants defined in t
%   CSTATE                  matrix containing status of each contact at the
%                           time instants defined in t
%                           0: separated; 1: sticking; 2: slipping
%   Solinfo                 Structure of solver info for each time instant
%       NIT                     number of iterations
%       RESIDUUM                residuum norm
%       IEx                     1: converged; 0: not converged
%========================================================================
function [Q,U,LAM,CSTATE,Solinfo] = ...
    sim_contact3D(t,q0,u0,MBmodel,varargin)
disp('=================================================================');
disp('NLstep, Copyright (C) 2025 Malte Krack, Johann Gross');
disp('This program comes with ABSOLUTELY NO WARRANTY.');
disp('This is free software, and you are welcome to redistribute');
disp('it under certain conditions, see gpl-3.0.txt.');
disp('=================================================================');

%% SET AND CHECK OPTIONS

% Set default values
Sopt = struct( ...
    'eliminatePresumedStickingContacts',false, ...
    'tol',1e-6, ...
    'epsALscl',1,...
    'itermax',1e4,...
    'dataStorageFrequency',1,...
    'cmdWindowOutputFrequency',10,...
    'verbosity',1);

% Overwrite by user-specified values
if nargin>4 && isstruct(varargin{1})
    tmp = varargin{1};
    new_fieldnames = fieldnames(tmp);
    for ij=1:length(new_fieldnames)
        Sopt.(new_fieldnames{ij}) = ...
            tmp.(new_fieldnames{ij});
    end
end

% Check for reasonable values
if numel(Sopt.tol)~=1 || ~isreal(Sopt.tol) || Sopt.tol<=0
    error('Sopt.tol must be positive real scalar.');
end
if numel(Sopt.epsALscl)~=1 || ~isreal(Sopt.epsALscl) || Sopt.epsALscl<=0
    error('Sopt.epsALscl must be positive real scalar.');
end
if numel(Sopt.itermax)~=1 || ~isreal(Sopt.itermax) || Sopt.itermax<1
    error('Sopt.itermax must be real scalar => 1.');
end
if numel(Sopt.dataStorageFrequency)~=1 || ...
        ~isreal(Sopt.dataStorageFrequency) || Sopt.dataStorageFrequency<1
    error('Sopt.dataStorageFrequency must be real scalar => 1.');
end
if numel(Sopt.cmdWindowOutputFrequency)~=1 || ...
        ~isreal(Sopt.cmdWindowOutputFrequency) || Sopt.cmdWindowOutputFrequency<1
    error('Sopt.cmdWindowOutputFrequency must be real scalar => 1.');
end
if numel(Sopt.verbosity)~=1 || ~isreal(Sopt.verbosity)
    error(['Sopt.verbosity must be real scalar (values 0, 1, 2, 3 ' ...
    'are reasonable).']);
end

%% HANDLE INPUT

% Set up /check time vector
if length(t)==2
    % Two-element vector specified, interpret first element as time step,
    % second as end time.
    dt = t(1);
    tend = t(2);
    if dt<=0 || tend<=0 || dt>tend
        error(['Time step and end time must be positive, ' ...
            'and time step must be larger than end time.']);
    else
        t = 0:t(1):t(2);
    end
elseif length(t)>2
    t = t(:).'; % enforce row vector
    dt = t(2)-t(1);
    if max(abs((diff(t)/dt-1)))>1e-6
        error('Time vector must be equidistant.');
    end
else
    error('Time vector must contain at least two elements.');
end
nt = length(t);

% Check model
if ~isa(MBmodel,'CMS_ROM') || ~isMasslessBoundaryROM(MBmodel)
    error('This function requires a massless-boundary ROM.');
end

% Set up M,D,K matrices and its partitions
M = MBmodel.M;      % reduced mass matrix
D = MBmodel.D;      % reduced damping matrix
K = MBmodel.K;      % reduced stiffness matrix
n = size(M,1);      % total number of generalized coordinates 
nb = MBmodel.nb;    % number of contact boundary coordinates
IB = 1:nb;          % index vector to contact boundary coordinates
II = (nb+1):n;      % index vector to remaining coordinates
% Extract partitions of M,D, and K. 
% NOTE: Mbb, Mib, Mbi are zero by construction of the massless-boundary 
% ROM. D is obtained by projection, in contrast, so that it should be 
% checked if it has indeed zero partitions associated with the boundary.
Mii = M(II,II);
Kii = K(II,II);
Kbb = K(IB,IB);
Kbi = K(IB,II);
Kib = K(II,IB);
if any(any(D(IB,IB))) || any(any(D(IB,II))) || any(any(D(II,IB)))
    warning(['Damping matrix has non-zero partitions related to ' ...
        'contact boundary coordinates. These are ignored in this ' ...
        'function.']);
end
Dii = D(II,II);

% Check initial value dimensions and reasonable numbers
q0 = q0(:); u0 = u0(:); % enforce column vectors
if size(q0,1)~=n || size(u0,1)~=n
    error(['Initial coordinate and velocity must be provided as ' ...
        'vector of length matching the provided model.']);
elseif any(any(isnan([q0 u0]))) || any(any(isinf([q0 u0])))
    error(['Initial coordinate and velocity values must not be ' ...
        'NaN or INF.']);
end

% Check for 3D contact and determine parameters
nl = MBmodel.nonlinear_elements;
% Check model dimension
if nb~=3*length(nl)
    error(['This function is for 3D contact only. ' ...
        'Every nonlinear element must be assigned to exactly three '...
        'contact boundary coordinate.']);
else
    % In the 3D case, number of contacts equals number of boundary
    % coordinates divided by 3
    nc = nb/3;
    % Define indices to normal and orthogonal tangential directions
    iN = (1:3:nb)';
    iT1 = iN+1;
    iT2 = iN+2;
end
% Determine imposed gap, velocity, preload, friction coefficient, and area 
% associated with contact node for quadrature from stress
gimp = cell(1,nc);
gamimp = cell(1,nc);
lam0 = cell(1,nc);
mu = zeros(nc,1);
area = zeros(nc,1);
for inl = 1:length(nl)
    % Check 3D contact type
    if ~strcmpi(nl{inl}.type,'3D')
        error('This function is for 3D contact only.');
    end
    % Handle imposed gap
    if isfield(nl{inl},'imposedGap')
        if isa(nl{inl}.imposedGap,'function_handle')
            gimp{inl} = feval(nl{inl}.imposedGap,t);
            if ~(size(gimp{inl},1)==1 || size(gimp{inl},1)==3) || ...
                    size(gimp{inl},2)~=nt
                error(['Imposed gap function handle must return either '...
                    '1D information (normal gap) or 3D information as ' ...
                    'column vector, and return an array whose second ' ....
                    'dimension matches the input time vector.']);
            end
            % If 1D information is given, interpret as imposed normal gap,
            % and set imposed tangential gaps to zero.
            if size(gimp{inl},1)==1
                gimp{inl} = [gimp{inl};zeros(2,nt)];
            end
        else
            if isscalar(nl{inl}.imposedGap)
                gimp{inl} = [repmat(nl{inl}.imposedGap,1,nt);zeros(2,nt)];
            elseif numel(nl{inl}.imposedGap)==3
                gimp{inl} = repmat(nl{inl}.imposedGap,1,nt);
            else
                error(['Imposed gap specified per contact element must' ...
                    ' be 1D or 3D.']);
            end
        end
    else
        gimp{inl} = spalloc(3,nt,0);
    end
    % Handle imposed gap velocity
    if isfield(nl{inl},'imposedGapVelocity')
        if isa(nl{inl}.imposedGapVelocity,'function_handle')
            gamimp{inl} = feval(nl{inl}.imposedGapVelocity,t);
            if ~(size(gamimp{inl},1)==1 || size(gamimp{inl},1)==3) || ...
                    size(gamimp{inl},2)~=nt
                error(['Imposed gap velocity function handle must return either '...
                    '1D information (normal gap) or 3D information as ' ...
                    'column vector, and return an array whose second ' ....
                    'dimension matches the input time vector.']);
            end
            % If 1D information is given, interpret as imposed normal gap,
            % and set imposed tangential gaps to zero.
            if size(gamimp{inl},1)==1
                gamimp{inl} = [gamimp{inl};zeros(2,nt)];
            end
        else
            if isscalar(nl{inl}.imposedGapVelocity)
                gamimp{inl} = [repmat(nl{inl}.imposedGapVelocity,1,nt);...
                    zeros(2,nt)];
            elseif numel(nl{inl}.imposedGapVelocity)==3
                gamimp{inl} = repmat(nl{inl}.imposedGapVelocity,1,nt);
            else
                error(['Imposed gap velocity specified per contact element must' ...
                    ' be 1D or 3D.']);
            end
        end
    else
        gamimp{inl} = spalloc(3,nt,0);
    end
    % Handle preload
    if isfield(nl{inl},'preload')
        if isa(nl{inl}.preload,'function_handle')
            error('Preload must not be function handle.');
        else
            if isscalar(nl{inl}.preload)
                lam0{inl} = [nl{inl}.preload;zeros(2,1)];
            else
                error(['(Normal) preload specified per contact element must' ...
                    ' be scalar.']);
            end
        end
    else
        lam0{inl} = zeros(3,1);
    end
    % Handle friction coefficient
    if isfield(nl{inl},'frictionCoefficient')
        if isa(nl{inl}.frictionCoefficient,'function_handle')
            error('Friction coefficient must not be function handle.');
        else
            if isscalar(nl{inl}.frictionCoefficient)
                mu(inl) = nl{inl}.frictionCoefficient;
            else
                error(['Friction coefficient specified per contact element must' ...
                    ' be scalar.']);
            end
        end
    else
        error(['A friction coefficient must be specified. Set it to ' ...
            'zero if you want frictionless contact.']);
    end
    % Handle area associated with contact node for quadrature from stress
    if isfield(nl{inl},'area')
        if isa(nl{inl}.area,'function_handle')
            error('Area associated with node must not be function handle.');
        else
            if isscalar(nl{inl}.area)
                area(inl) = nl{inl}.area;
            else
                error(['Area associated with node must' ...
                    ' be scalar.']);
            end
        end
    else
        area(inl) = 1;
    end
end
gimp = vertcat(gimp{:}); % nb x nt matrix ~ 3D gaps
gamimp = vertcat(gamimp{:}); % nb x nt matrix ~ 3D gap velocities
lam0 = vertcat(lam0{:}); % nb x 1 vector ~ 3D preload
area = kron(area,ones(3,1)); % nb x 1 vector ~ 3D quadrature weights

% Check for function handle to external forcing (explicit time dependence)
if nargin>5
    func_fex = varargin{2};
    % Test for correct dimension and reasonable output at start/end times
    fex_S = func_fex(t(1)); fex_E = func_fex(t(end));
    if size(fex_S,1)~=n || size(fex_S,2)~=1 || ...
            size(fex_E,1)~=n || size(fex_E,2)~=1
        error(['Function handle to external forcing must return ' ...
            'column vector of length ' num2str(n) ' for scalar input.']);
    elseif any(any(isnan([fex_S fex_E]))) || any(any(isinf([fex_S fex_E])))
        error(['Function handle to external forcing must not return ' ...
            'NaN or INF.']);
    end
else
    func_fex = @(t) zeros(n,size(t,2));
end

%% PREPARE OUTPUT

% Allocate function output variables
ntout = floor(nt/Sopt.dataStorageFrequency);
Q = zeros(n,ntout);
U = zeros(n,ntout);
LAM = zeros(nb,ntout);
CSTATE = zeros(nc,ntout);
Solinfo.NIT = zeros(1,ntout);
Solinfo.IEx = ones(1,ntout);
Solinfo.RESIDUUM = zeros(1,ntout);

% Prepare command window output
nPrintStats = Sopt.cmdWindowOutputFrequency+1;
kPrintStats = floor(linspace(1,nt,nPrintStats));
simProgress = linspace(0,1,nPrintStats);
simStopWatch = tic;

% Allocate solver statistics (results are averaged between outputs)
navg = ceil(nt/(nPrintStats-1));
numActiveContacts = zeros(1,nt);
numIterations = zeros(1,navg);
residuum = zeros(1,navg);
isConverged = ones(1,navg);
istep = 1; % counter for data storage
jstep = 1; % counter for command window output

%% START TIME STEP INTEGRATION

% Inform user about start of integration
if Sopt.verbosity>0
    fprintf(['Begin simulation with leapfrog scheme ' ...
        'for %.6g equidistant time steps (dt=%.5g).\n'],nt,dt)
    fprintf('--------------------\n');
end

% Auxiliary matrices for time integration scheme
MdtMinus = Mii/dt - Dii/2;
MdtPlusinv = (Mii/dt+Dii/2)\eye(size(Mii,1));
Kbbdtinv = 1/dt*((Kbb)\eye(size(Kbb,1)));

% Set tolerance for slip velocity
tolGamT = Sopt.tol*(norm(Kbbdtinv));

for k=1:nt
    %% SET POINT OF DEPARTURE FOR CURRENT ONE-STEP PROBLEM

    if k==1
        % Adopt initial values specified by user.
        % NOTE: Due to the shifted grid, initial coordinate and velocity
        % needed for this function actually correspond to different time
        % levels (coordinate at t=0 vs. velocity at t=-dt/2). This could be
        % fixed by setting up an equation system, from which the required
        % quantities are obtained for initial values given at the same
        % time instant. However, considering that typically we do not
        % know the precise initial values anyway and are interested in the
        % steady state, we prefer to simply accept inconsistent initial
        % values here.
        % ADDITIONAL NOTE: The qb_km1 here is only relevant for calculating
        % the boundary velocity. The boundary coordinate is not dependent
        % on this initial value.
        qi_k = q0(II);
        ui_km12 = u0(II);
        qb_km1 = q0(IB);
        ub_km12 = u0(IB);

        % Use zero initial guess for contact forces
        lam = zeros(nb,1);
    else
        % Adopt values from previous time level
        qi_k = qi_kp1;
        qb_km1 = qb_k;
        ui_km12 = ui_kp12;
        % ub_km12 = ub_km12; % COMMENTED OUT BECAUSE OF TRIVIAL IDENTITY
        % Adopt lam from previous time level, too
        % lam = lam; % COMMENTED OUT BECAUSE OF TRIVIAL IDENTITY
    end

    % Evaluate external force
    fex_k = real(feval(func_fex,t(k)));

    %% ESTIMATE ACTIVE SET

    % Predict gap assuming uniform motion
    gpre = qb_km1 + dt*ub_km12 + gimp(:,k);

    % Contacts are estimated as open depending on predicted positive gap
    % and previous non-positive normal contact force.
    lamNpre = lam(iN)+lam0(iN);
    isOpen = gpre(iN) > 0 & lamNpre <= 0;

    % The elimination of presumed sticking contacts is switched
    % on by specifying Sopt.eliminatePresumedStickingContacts = true.
    % Only then 'isSticking' may assume true values.
    if Sopt.eliminatePresumedStickingContacts
        lamTpre = sqrt(lam(iT1).^2 + lam(iT2).^2);
        gamTpre = sqrt( (ub_km12(iT1) + gamimp(iT1,k)).^2 + ...
            (ub_km12(iT2) + gamimp(iT2,k)).^2 );
        isSticking = lamTpre<mu.*lamNpre & gamTpre<tolGamT;
        isSticking(isOpen) = false;
    else
        isSticking = false(nc,1);
    end

    % The rest is treated as active.
    isActive = ~isOpen & ~isSticking;

    % Expand Booleans from scalar to 3D per contact
    ISSTICKING = logical(kron(isSticking,ones(3,1)));
    ISACTIVE = logical(kron(isActive,ones(3,1)));
    ISOPEN = logical(kron(isOpen,ones(3,1)));
    % Previous value of contact force is adopted for active contacts, while
    % zero (dynamic 'lam' plus preload) is used for the open ones.
    lam(ISOPEN) = -lam0(ISOPEN);

    % Update global 'c' vector
    c = Kbbdtinv*(fex_k(IB)-Kbi*qi_k) - qb_km1/dt + gamimp(:,k);
    %% COMPUTE ACTIVE CONTACT FORCES
    if any(~isOpen)
        if ~any(isSticking)
            % This is the simpler case where the set of presumed sticking
            % contacts is empty

            % Apply projected Jacobi over-relaxation solver
            [lamActive,output,isClosed,isSliding] = ...
                JORproj(Kbbdtinv(ISACTIVE,ISACTIVE),c(ISACTIVE)+...
                Kbbdtinv(ISACTIVE,ISOPEN)*lam(ISOPEN),...
                lam(ISACTIVE),mu(isActive),area(ISACTIVE),lam0(ISACTIVE),Sopt);

            % Store active contact forces
            lam(ISACTIVE) = lamActive;
        else
            % Here the elimination of presumed sticking contacts is applied
            if exist('isStickingOLD','var') && ...
                    exist('isActiveOLD','var') ...
                    && all(isStickingOLD==isSticking) && ...
                    all(isActiveOLD==isActive)
                % The appropriate Gred and epsAL are available
            else
                % The appropriate Gred and epsAL need to be calculated
                % lamSticking = -FART*lamActive - smell
                % gamActive = Gred*lamActive + cred
                FART = Kbbdtinv(ISSTICKING,ISSTICKING)\...
                    Kbbdtinv(ISSTICKING,ISACTIVE);
                Gred = Kbbdtinv(ISACTIVE,ISACTIVE) - ...
                    Kbbdtinv(ISACTIVE,ISSTICKING)*FART;
                epsAL = [];
            end
            smell = Kbbdtinv(ISSTICKING,ISSTICKING)\(c(ISSTICKING)+...
                Kbbdtinv(ISSTICKING,ISOPEN)*lam(ISOPEN));
            cred = c(ISACTIVE) - Kbbdtinv(ISACTIVE,ISSTICKING)*smell + ...
                Kbbdtinv(ISACTIVE,ISOPEN)*lam(ISOPEN);

            % Handle active set
            if any(ISACTIVE)
                % Apply projected Jacobi over-relaxation solver
                [lamActive,output,isClosed,isSliding,epsAL] = ...
                    JORproj(Gred,cred,lam(ISACTIVE),...
                    mu(isActive),area(ISACTIVE),lam0(ISACTIVE),Sopt,epsAL);
            else
                % All contacts presumed sticking; trivial solution
                lamActive = zeros(0,1);
                isClosed = false(0,1);
                isSliding = false(0,1);
                output.iter = 0;
                output.residuum = 0;
                output.isConverged = true;
            end

            % Store active and sticking contact forces
            lam(ISACTIVE) = lamActive;
            lam(ISSTICKING) = -FART*lamActive-smell;
        end

        % Store solver output for statistics
        numActiveContacts(jstep) = sum(isActive);
        numIterations(jstep) = output.iter;
        residuum(jstep) = output.residuum;
        isConverged(jstep) = output.isConverged;
    else
        % Mimic solver output in the trivial case of no active contacts
        numActiveContacts(jstep) = 0;
        numIterations(jstep) = 0;
        residuum(jstep) = 0;
        isConverged(jstep) = 1;
    end

    % These must be set in any case:
    isStickingOLD = isSticking;
    isActiveOLD = isActive;

    %% COMPUTE REMAINING UNKNOWNS OF ONE-STEP PROBLEM  

    % Evaluate gap and boundary velocity compatible with computed contact
    % forces
    gam = Kbbdtinv*lam + c;
    ub_km12 = gam - gamimp(:,k);

    % Evaluate boundary coordinates at t(k)
    qb_k = qb_km1 + ub_km12*dt;

    % Solve dynamic force balance for internal velocities at t(k)+dt/2
    ui_kp12  = MdtPlusinv * ...
        ( fex_k(II) + MdtMinus*ui_km12 - Kii*qi_k - Kib*qb_k );

    % Integrate inner coordinates to t(k+1)
    qi_kp1 = qi_k + ui_kp12*dt;

    %% OUTPUT

    % Store results at given data storage frequency
    if mod(k,Sopt.dataStorageFrequency) == 0
        % Coordinate, velocity, contact force and state output
        Q(:,istep) = [qb_k;qi_k];
        U(:,istep) = [ub_km12;ui_km12];
        LAM(:,istep) = lam;
        if any(isActive)
            idxActive = find(isActive);
            CSTATE(idxActive(isClosed&isSliding),istep) = 2;
            CSTATE(idxActive(isClosed&~isSliding),istep) = 1;
            CSTATE(isSticking,istep) = 1;
            % CSTATE(idxActive(~isClosed),istep) = 0; % covered by allocation
            % CSTATE(IIsep,istep) = 0; % covered by allocation
        end

        % Solver statistics output
        Solinfo.NIT(istep) = numIterations(jstep);
        Solinfo.RESIDUUM(istep) = residuum(jstep);
        Solinfo.IEx(istep) = isConverged(jstep);

        % Increment data storage counter
        istep = istep+1;
    end

    % Command window output
    if k == kPrintStats(1)
        if k==1
            if Sopt.verbosity>1
                fprintf(['Solver stats are printed for first step and' ...
                    ' then as average over intervals according to ' ...
                    'command window output frequency.']);
                fprintf('\n');
                if Sopt.verbosity>2
                    fprintf(['Residuum is measured' ...
                        ' in terms of L-infinity norm of vector of contact' ...
                        ' force difference between last two iterations.']);
                    fprintf('\n');
                end
            end
        end
        if Sopt.verbosity>0
            fprintf('Progress = %.2f %%.', simProgress(1)*100 );
        end
        if Sopt.verbosity>1
            fprintf(' Number of JORproj iterations = %2.1f.', ...
                mean(numIterations(1:jstep)));
            if Sopt.verbosity>2
                fprintf(' residuum = %.2g', mean(residuum(1:jstep)));
                fprintf(', number of active contacts = %2.1f.',...
                    mean(numActiveContacts(1:jstep)));
            end
        end
        if Sopt.verbosity>0
            fprintf('\n');
        end
        % Reset step counter and set trigger to next output
        jstep = 1;
        kPrintStats(1) = [];
        simProgress(1) = [];
    else
        % Increment step counter for printing solver statistics
        jstep = jstep+1;
    end
end
disp('--------------------');
disp('End of simulation.');
disp('COMPUTATIONAL EFFORT:');
disp(['Total elapsed time is ', ...
    num2str(toc(simStopWatch),'%16.1f'),' s']);
disp(['Average number of iterations is ', ...
    num2str(mean(Solinfo.NIT),'%2.1f'),'.']);
%______________________________________________________
end

function [lam,output,isClosed,isSliding,epsAL] = ...
    JORproj(G,c,lam,mu,area,lam0,Sopt,varargin)
%========================================================================
% This function implements the projected Jacobi over-relaxation method to
% iteratively solve the non-smooth equation system
%   lam - projC(lam - epsAL.*(G*lam+c)) = 0,
% where G is denoted Delassus matrix, projC denotes the projection onto the 
% set C, and epsAL>0 can be interpreted as Augmented Lagrangian parameter.
% The options structure Sopt must contain the relative error tolerance tol, 
% the scaling parameter epsALscl, and the maximum number of iterations 
% itermax. The residuum is measured in terms of the L-infinity norm of 
% the difference of lam between the last two iterations, and it is 
% the tolerance is relative to the L-infinity norm of lam. The scaling
% parameter increases (epsALscl>1) or decreases (0<epsALscl<1) the
% parameter epsAL with respect to the one recommended in [2].
%========================================================================

% Here we take care of the conversion between consistent nodal forces and
% stress. The initial guess of the contact force lam and preload lam0 are 
% divided by the respective area to go to stress level. Subsequently, the
% variable lam is in stress units. To this end, the coefficient matrix G
% must be right-multiplied by the diagonal matrix containing the respective
% areas, so that the new G*lam yields the same value as before. The 
% problem then is that G is no longer symmetric, which may lead to poor 
% convergence behavior of the projected Jacobi over-relaxation method. We 
% fix this by left-multiplying G by the same diagonal matrix, and also the
% c vector. This is allowed because effectively we change epsAL with this,
% which is arbitrary except for the requirement that epsAL>0, and this
% remains satisfied since the areas are strictly positive.
AREA = diag(area);
lam = AREA\lam;
lam0 = AREA\lam0;
G = AREA'*G*AREA;
c = AREA'*c;

% Select epsAL parameter
if nargin>=8 && isscalar(varargin{1}) && varargin{1}>0
    % Use specified parameter
    epsAL = varargin{1};
else
    % Set epsAL parameter as recommended in [2], Eqs. 4.59-4.60.
    G_diagonal_dominant = all((2*abs(diag(G))) >= sum(abs(G),2));
    if G_diagonal_dominant
        epsAL = Sopt.epsALscl./abs(diag(G));
    else
        epsAL = Sopt.epsALscl./sum(abs(G),2);
    end
end

% Carry out fixed point iterations in accordance with projected Jacobi
% over-relaxation method
converged = false; iter = 0;
while (~converged && iter <= Sopt.itermax)
    lamp = projC( lam -  epsAL.*(G*lam+c), mu, lam0);

    % Check for convergence
    residuum = max(abs( lamp - lam ));
    converged = ( residuum <= max(abs(lam))*Sopt.tol );

    % Increment loop counter
    iter = iter + 1;
    lam = lamp;
end

% Warn user if maximum number of iterations was reached before convergence.
if iter > Sopt.itermax && ~converged
    disp(['WARNING: ITERMAX=' num2str(iter) ' reached with residuum: ' ...
        num2str(residuum) '.']);
end

% Here we convert back from stress to force via the respective area
lam = AREA*lam;
lam0 = AREA*lam0;

% Store output
isClosed = lam(1:3:end)+lam0(1:3:end)>0;
isSliding = sqrt(lam(2:3:end).^2+lam(3:3:end).^2) >= ...
    mu.*(lam(1:3:end)+lam0(1:3:end));
isSliding(~isClosed) = false; % only closed contacts can slide
output.iter = iter-1;
output.residuum = residuum;
output.isConverged = converged;

end

function x = projC(x,mu,lam0)
% In the case of 3D contact, dry Coulomb friction in the tangential
% contact plane, and unilateral interaction in the normal contact
% direction are enforced.
xN = max(x(1:3:end)+lam0(1:3:end),0);
xT1 = x(2:3:end);
xT2 = x(3:3:end);
normxT = sqrt(xT1.^2+xT2.^2);
isSliding = normxT >= mu.*xN;
muxN_normxT = mu(isSliding).*xN(isSliding)./normxT(isSliding);
muxN_normxT(isnan(muxN_normxT)) = 0;
xT1(isSliding) = muxN_normxT .*xT1(isSliding);
xT2(isSliding) = muxN_normxT .*xT2(isSliding);
x(1:3:end) = xN-lam0(1:3:end);
x(2:3:end) = xT1;
x(3:3:end) = xT2;
end