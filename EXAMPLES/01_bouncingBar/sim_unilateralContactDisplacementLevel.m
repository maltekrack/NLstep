%========================================================================
% DESCRIPTION:
% This function is a variant of the sim_contact3D function in the SRC
% folder. The key difference is that this function is restricted to
% 1D unilateral contact, and the contact is enforced on displacement level
% (whereas NLstep otherwise enforces contact on velocity level in order to
% treat dry Coulomb friction). Further, this function permits only
% initially open contacts (initially preloaded contacts are not
% implemented), and no quadrature of contact stress is done (i.e. the area
% field of the nonlinear elements is ignored in this function, leading to 
% a unity weighting).
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
%
%   Q, U, LAM               matrices whose colums are the coordinates,
%                           velocities, and contact forces at the time 
%                           instants defined in t
%   CSTATE                  matrix containing status of each contact at the 
%                           time instants defined in t
%                               0: open; 1: closed
%   Solinfo                 Structure of solver info for each time instant
%       NIT                     number of iterations 
%       RESIDUUM                residuum norm
%       IEx                     1: converged; 0: not converged
%========================================================================
function [Q,U,LAM,CSTATE,Solinfo] = ...
    sim_unilateralContactDisplacementLevel(t,q0,u0,MBmodel,varargin)
disp('=================================================================');
disp('NLstep, Copyright (C) 2025 Malte Krack, Johann Gross');
disp('This program comes with ABSOLUTELY NO WARRANTY.');
disp('This is free software, and you are welcome to redistribute');
disp('it under certain conditions, see gpl-3.0.txt.');
disp('=================================================================');

%% SET AND CHECK OPTIONS

% Set default values
Sopt = struct( ...
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

% Check for 1D contact and evaluate imposed gaps
nl = MBmodel.nonlinear_elements;
if nb~=length(nl)
    error(['This function is for 1D contact only. ' ...
        'Every nonlinear element must be assigned to exactly one '...
        'contact boundary coordinate.']);
else
    % In the 1D case, number of contacts equals number of boundary
    % coordinates
    nc = nb;
end
imposedGap = cell(1,nc);
for inl = 1:length(nl)
    if ~strcmpi(nl{inl}.type,'unilateral')
        error('This function is for unilateral contact only.');
    end
    if isfield(nl{inl},'imposedGap')
        if isa(nl{inl}.imposedGap,'function_handle')
            imposedGap{inl} = feval(nl{inl}.imposedGap,t);
            if size(imposedGap{inl},1)~=1 || size(imposedGap{inl},2)~=nt
                error(['Imposed gap function handle must return row vector ' ...
                    'of the same length as the input time vector.']);
            end
        else
            if numel(nl{inl}.imposedGap)~=1
                error(['Imposed gap specified per unilateral contact element must' ...
                    ' be scalar.']);
            end
            imposedGap{inl} = repmat(nl{inl}.imposedGap,1,nt);
        end
    else
        error(['An imposed gap must be specified. Initially preloaded ' ...
            'contacts are not implemented in this particular function.']);
    end
end
imposedGap = vertcat(imposedGap{:}); % nc x nt matrix

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
Kbbinv = (Kbb)\eye(size(Kbb,1));

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
        qb_km1 = q0(IB); % ONLY NEEDED FOR BOUNDARY VELOCITY RECONSTRUCTION
        % ub_km12 = u0(IB); % NOT NEEDED IN THIS FUNCTION

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

    % Predict gap assuming zero contact force (Eq. 11 in [1]).
    gpre = Kbbinv*(fex_k(IB)-Kbi*qi_k) + imposedGap(k);

    % Contacts are estimated as open depending on predicted gap.
    isOpen = gpre > 0;

    % The rest is treated as active.
    isActive = ~isOpen;

    % Previous value of contact force is adopted for active contacts, while
    % zero is used for the separated ones.
    lam = lam.*isActive;

    %% COMPUTE ACTIVE CONTACT FORCES
    if any(isActive)
        % The complementary inequality contact constraints are formulated
        % as projective non-smooth equation system 
        %       lambda - projC(lambda - epsAL.*g) = 0,
        % where lambda is the contact force vector, g is the gap
        % vector, projC denotes the projection onto the admissible
        % set of contact forces, and epsAL>0 can be interpreted as
        % Augmented Lagrangian parameter. Contact kinematics, integration
        % scheme, and force balance at the boundary yield
        %       g = G*lambda + c,
        % where G is denoted Delassus matrix. Substituting this into the
        % above non-smooth equation system yields an implicit equation
        % system for lambda. This is iteratively solved using projected 
        % Jacobi over-relaxation.
        % 
        % In this function, the above equations are formulated on
        % displacement level. This is only possible in the frictionless
        % case. The rest of NLstep is formulated on velocity level to
        % treat dry Coulomb friction.

        % Evaluate G and c
        G = Kbbinv(isActive,isActive);
        c = Kbbinv(isActive,:)*(fex_k(IB) - Kbi*qi_k) + ...
            imposedGap(isActive,k);    

        % Apply projected Jacobi over-relaxation solver
        [lamActive,output,isClosed] = ...
            JORproj(G,c,lam(isActive),Sopt);

        % Store active contact forces
        lam(isActive) = lamActive;

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

    %% COMPUTE REMAINING UNKNOWNS OF ONE-STEP PROBLEM
    
    % Evaluate boundary coordinates compatible with computed contact forces
    % and static force balance at the boundary at t(k)
    qb_k  = Kbbinv*(lam - Kbi*qi_k + fex_k(IB));

    % Evaluate boundary velocity at t(k)-dt/2
    ub_km12  = 1/dt*(qb_k-qb_km1);

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
            CSTATE(idxActive(isClosed),istep) = 1;
            % CSTATE(idxActive(~isClosed),istep) = 0; % covered by allocation
            % CSTATE(isOpen,istep) = 0; % covered by allocation
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

function [lam,output,isClosed] = JORproj(G,c,lam,Sopt)
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

% Set epsAL parameter as recommended in [2], Eqs. 4.59-4.60.
G_diagonal_dominant = all((2*abs(diag(G))) >= sum(abs(G),2));
if G_diagonal_dominant
    epsAL = Sopt.epsALscl./abs(diag(G));
else
    epsAL = Sopt.epsALscl./sum(abs(G),2);
end

% Carry out fixed point iterations in accordance with projected Jacobi 
% over-relaxation method
converged = false; iter = 0;
while (~converged && iter <= Sopt.itermax)
    % In the case of unilateral contact, the set C is simply that of
    % non-negative real numbers. The projection can then simply be
    % implemented with the max function.
    lamp = max( lam -  epsAL.*(G*lam+c), 0*lam );

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

% Set output
isClosed = lam>0;
output.iter = iter-1;
output.residuum = residuum;
output.isConverged = converged;

end