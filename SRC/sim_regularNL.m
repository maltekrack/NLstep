%========================================================================
% Matlab function for the time step integration of a linear time-invariant
% mechanical system with local regular nonlinear elements and possible
% forcing with explicit time dependence, from a given set of initial values
% of displacement and velocities, using the Newmark's quadrature rules in
% conjunction with equilibrium averaging. Popular particular forms include
% the (conventional) Newmark method, the HHT (alpha-)method, the 
% Generalized-alpha method, and the Czekanski method (widely-used for
% contact).
%
% DYNAMIC FORCE BALANCE
% The differential equation system governing the motion of the mechanical
% system reads:
%
%       M * \ddot q + D * \dot q + K * q + fnl = fex(t),            (1)
%       fnl = \sum_j w_j * f_j ( w_j^T * q, q_j^T * \dot q ).       (2)
%
% Herein, \dot denotes derivative with respect to time t. q is the
% vector of n generalized coordinates. M, D, K are real n x n coefficient 
% matrices. Symmetry and/or positive definiteness is not presumed, which 
% can be useful for simulating massless-boundary models (singular M), 
% and accounting for gyroscopic or circulatoric forces (associated with 
% skew-symmetric part of D, K, respectively. The matrices are input as
% fields of the input variable model (model.M, model.D, model.K), which may
% be a structure, or an object of either the MechanicalSystem or the
% CMS_ROM class of NLvib (https://www.ila.uni-stuttgart.de/nlvib/).
% 
% LOCAL NONLINEAR ELEMENTS, 3D CONTACT CAPABILITIES
% The notion of local nonlinear elements is in accordance with that of
% NLvib. 
% NLvib is restricted to scalar functions f_j so that w_j are n x 1 
% vectors. The implementation here is more general: 3D contact models are 
% available, leading to a 3 x 1 vector function f_j and the corresponding 
% w_j is a n x 3 matrix. As in NLvib, frictional hysteresis is implemented 
% using an incremental evolution formulation which is resolved using a
% predictor-corrector scheme. For this, the coordinate at the previous time
% level is used (which is not available/used in most common integrators).
% In this case, hence, f_j is not an explicit function of displacement
% and velocity at the current time, in contrast to what Eq. (2) suggests. 
% Apart from this, only continuous functions f_j are allowed. If you 
% intend to model contact via set-valued laws, such as the (non-regularized)
% Coulomb-/Signorini-conditions, this function is not appropriate, and
% 'sim_contact3D.m' or 1D variants of it should be used. In the case of 3D
% contact, a scalar friction coefficient, and a scalar (penalty) stiffness 
% must be specified. Also, a scalar imposed normal gap or normal preload 
% can be specified. The contact is well-defined only if one of those is 
% specified. The nonlinear elements are defined in
% model.nonlinear_elements, which is a cell array of structures. For 
% element nl, the parameters are specified as fields of 
% model.nonlinear_elements{nl}, with the name type, and possible more
% parameters which are specific to the given type. The w_j is given by the
% field name force_direction. In the case of 3D contact only, it is
% optional (if you do not specify it, unit directions will be assumed, as
% in sim_contact3D.m), otherwise it is mandatory. For contact, important
% field names are frictionCoefficient, imposedGap, preload, and
% stiffness. Most of the names should be self-explanatory. In this 
% function, all of those must have fixed values (function handles are not
% allowed). The stiffness refers to the penalty parameter. The preload is 
% to be specified as consistent nodal force, and, similarly, the stiffness 
% is to be specified in the unit force divided by displacement. To 
% implement a uniform stiffness per area and/or an imposed normal pressure, 
% using a node-based quadrature, you have to multiply by the respective 
% area before. Note that this is in contrast to sim_contact3D.m, where the
% area field is treated within the function.
%
% TIME STEP INTEGRATION
% The equilibrium averaging concept is applied to (1)-(2), where the 
% inertia forces are replaced by
%           (1-alpha_m) * M * accE + alpha_m * M * accS.            (3)
% Herein, accE, accS are the acceleration vectors at current and previous
% time level, respectively, and alpha_m is the corresponding averaging
% parameter. Analogously, the remaining forces are averaged with a
% parameter alpha_f. The notation here follows closely that in the
% textbook [1]. For specific values of alpha_m, alpha_f, different
% popular methods arise such as the (conventional) Newmark method for
% alpha_f=0=alpha_m, or the HHT method.
% The integration method further relies on Newmark's quadrature rules,
%                 accE = a1*(qE-qS)-a2*uS-a3*accS,                  (4)
%                 velE = a4*(qE-qS)+a5*uS+a6+accS,                  (5)
% where velE, velS denote the velocity vector at current and previous
% time level, and qE, qS analogously for the displacement vector. The
% coefficients a1 through a6 are expressed using the parameters beta
% and gamma as commonly done in Newmark's method, as well as the time
% step dt,
%       a1 = 1/beta/dt^2,
%       a2 = 1/beta/dt,
%       a3 = (1-2*beta)/2/beta,
%       a4 = gamma/beta/dt,
%       a5 = 1-gamma/beta,
%       a6 = (1-gamma/2/beta)*dt.
% Different popular integration schemes of the Newmark family arise for
% specific settings of beta and gamma, such as the average constant
% acceleration variant (gamma = 1/2, beta = 1/4). Classical stability
% and accuracy results are available, see e.g. [1], FOR THE CASE
% of a linear time-invariant mechanical system and some additional
% assumptions on symmetry and positive definiteness of the coefficient
% matrices. The generalization to frictional-unilateral contact can be
% found e.g. in [2-3]. Those conditions DO NOT GENERALLY HOLD for the
% problem class of primary interest here, where general nonlinear behavior
% is allowed and the system may have a singular mass matrix.
% 
% SOLUTION OF THE ONE-STEP PROBLEM
% Inserting (3)-(5) into (1)-(2), one obtains an algebraic equation system
% for qE. This one-step problem is solved using the Newton method with
% a Cholesky factorization of the analytical Jacobian.
%
% REFERENCES
%   [1] M. Geradin, D. Rixen: Mechanical Vibrations. ISBN 978-1-118-90020-8
%   [2] Czekanski et al. https://doi.org/10.1002/cnm.411
%   [3] Czekanski et al. https://doi.org/10.1016/S0168-874X(01)00072-5
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
%   model                   model defining M,D,K and nonlinear_elements
%
% Optional input variables:
%   varargin{1} = Sopt   Options structure; only specified options 
%                           are changed from their default values
%       scheme                  This field specifies the scheme. Default is
%                               Newmark's average constant acceleration
%                               scheme. The user can set scheme.name to
%                                   'Newmark' and specify the parameter 
%                                       scheme.alpha (leading to 
%                                       alpha_f=0,
%                                       alpha_m=0, 
%                                       gamma=1/2+alpha, 
%                                       beta=1/4*(1+alpha)^2), or
%                                   'HHT' and specify the parameter 
%                                       scheme.alpha (leading to 
%                                       alpha_f=alpha, 
%                                       alpha_m,gamma,beta as in (b)), or
%                                   'GeneralizedAlpha' and specify the
%                                       parameter rho=scheme.rhoinfty
%                                       (leading to
%                                       alpha_m=(2*rho-1)/(rho+1), 
%                                       alpha_f=rho/(rho+1),
%                                       gamma=1/2+alpha_f-alpha_m,
%                                       beta=1/4*(1+alpha_f-alpha_m)^2), or
%                                   'Czekanski' to 
%                                       specify the parameter 
%                                       scheme.alpha_m only (leading to
%                                       alpha_f,beta,gamma set in
%                                       accordance with what is viewed 
%                                       optimal for contact in [2-3]).
%                               As alternative, the user does not set
%                               scheme.name and instead specifies 
%                                   scheme.alpha_m, scheme.alpha_f,
%                                   scheme.gamma, scheme.beta directly.
%       qprevious               coordinate vector at time level prior to q0
%                               (presuming equidistant time levels); needed 
%                               for consistent initialization of hysteretic
%                               nonlinear elements (default: q0)
%       fnlprevious             nonlinear force vector at time level
%                               prior to q0 (presuming equidistant time
%                               levels); needed for consistent
%                               initialization of hysteretic nonlinear
%                               elements (default: 0*q0)
%       tol                     scalar convergence tolerance; residuum 
%                               of the above mentioned algebraic equation
%                               system is measured in terms of the 
%                               L-infinity norm
%       itermax                 maximum number of Newton iterations per
%                               time step
%       dataStorageFrequency    set this to n>1 to store data only every n
%                               time levels
%       cmdWindowOutputFrequency    ... (as data storage but for command
%                                   window output)
%       verbosity               0: none; 1: basic; 2: iter; 3: detailed
%   varargin{2} = func_fex  function handle to external forcing with known
%                           explicit time dependence, fex(t)
%
% Output variables:
%   Q, U, FNL               matrices whose colums are the coordinates,
%                           velocities, and nonlinear forces at the time
%                           instants defined in t
%   Solinfo                 Structure of solver info for each time instant
%       NIT                     number of iterations
%       RESIDUUM                residuum norm
%       IEx                     1: converged; 0: not converged
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
function [Q,U,FNL,Solinfo] = sim_regularNL(t,q0,u0,model,varargin)

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
    'itermax',1e2,...
    'dataStorageFrequency',1,...
    'cmdWindowOutputFrequency',10,...
    'verbosity',1,...
    'qprevious',[],...
    'fnlprevious',[]);
% Default integration scheme is Newmark's average constant acceleration
% method (equivalent to trapezoidal quadrature rule)
alpha_m = 0;
alpha_f = 0;
beta = 1/4;
gamma = 1/2;
scheme = 'ConstantAverageAcceleration';

% Overwrite by user-specified values
if nargin>=5 && isstruct(varargin{1})
    tmp = varargin{1};
    new_fieldnames = fieldnames(tmp);
    for ij=1:length(new_fieldnames)
        if strcmpi(new_fieldnames{ij},'scheme')
            if isfield(tmp.scheme,'name') && ...
                    strcmpi(tmp.scheme.name,'Newmark') && ...
                    isfield(tmp.scheme,'alpha')
                scheme = 'Newmark';
                gamma = 1/2+tmp.scheme.alpha;
                beta = 1/4*(1+tmp.scheme.alpha)^2;
            elseif isfield(tmp.scheme,'name') && ...
                    strcmpi(tmp.scheme.name,'HHT') && ...
                    isfield(tmp.scheme,'alpha')
                scheme = 'HHT';
                alpha_f = tmp.scheme.alpha;
                gamma = 1/2+tmp.scheme.alpha;
                beta = 1/4*(1+tmp.scheme.alpha)^2;
            elseif isfield(tmp.scheme,'name') && ...
                    strcmpi(tmp.scheme.name,'GeneralizedAlpha') && ...
                    isfield(tmp.scheme,'rhoinfty')
                scheme = 'GeneralizedAlpha';
                rho = tmp.scheme.rhoinfty;
                alpha_m = (2*rho-1)/(rho+1);
                alpha_f = rho/(rho+1);
                gamma = 1/2+alpha_f-alpha_m;
                beta = 1/4*(1+alpha_f-alpha_m)^2;
           elseif isfield(tmp.scheme,'name') && ...
                    strcmpi(tmp.scheme.name,'czekanski') && ...
                    isfield(tmp.scheme,'alpha_m')
                scheme = 'Czekanski';
                alpha_m = tmp.scheme.alpha_m;
                if numel(alpha_m)~=1 || ~isreal(alpha_m) ...
                        || alpha_m>0.5 || alpha_m<-0.5
                    error(['alpha_m must be real scalar between -0.5 ' ...
                        'and +0.5.']);
                end
                alpha_f = (-2*alpha_m^2+alpha_m-1+...
                    sqrt(2*alpha_m^2-3*alpha_m+2)/(1-alpha_m));
                gamma = 0.5+alpha_f-alpha_m;
                beta = 0.25*(-2*alpha_m^2+...
                    alpha_m*(3+2*alpha_f)-2)/(alpha_m-1);
            elseif isfield(tmp.scheme,'alpha_m') && ...
                    isfield(tmp.scheme,'alpha_f') && ...
                    isfield(tmp.scheme,'beta') && ...
                    isfield(tmp.scheme,'gamma')
                scheme = 'EquilibriumAveragingNewmarkQuadrature';
                alpha_m = tmp.scheme.alpha_m;
                alpha_f = tmp.scheme.alpha_f;
                beta = tmp.scheme.beta;
                gamma = tmp.scheme.gamma;
            else
                error(['Please see header of this file for allowed ' ...
                    'settings of scheme field of options structure.']);
            end
        else
            Sopt.(new_fieldnames{ij}) = ...
                tmp.(new_fieldnames{ij});
        end
    end
end

% Check for reasonable values
if numel(alpha_m)~=1 || ...
        ~isreal(alpha_m) || isinf(alpha_m) || isnan(alpha_m)
    error('alpha_m must be real scalar.');
end
if numel(alpha_f)~=1 || ...
        ~isreal(alpha_f) || isinf(alpha_f) || isnan(alpha_f)
    error('alpha_f must be real scalar.');
end
if numel(beta)~=1 || ...
        ~isreal(beta) || isinf(beta) || isnan(beta)
    error('beta must be real scalar.');
end
if numel(gamma)~=1 || ...
        ~isreal(gamma) || isinf(gamma) || isnan(gamma)
    error('gamma must be real scalar.');
end
if numel(Sopt.tol)~=1 || ~isreal(Sopt.tol) || Sopt.tol<=0
    error('Sopt.tol must be positive real scalar.');
end
if numel(Sopt.itermax)~=1 || ~isreal(Sopt.itermax) || Sopt.itermax<1
    error('Sopt.itermax must be real scalar => 1.');
end
if numel(Sopt.dataStorageFrequency)~=1 || ...
        ~isreal(Sopt.dataStorageFrequency) || Sopt.dataStorageFrequency<1
    error('Sopt.dataStorageFrequency must be real scalar => 1.');
end
if numel(Sopt.cmdWindowOutputFrequency)~=1 || ...
        ~isreal(Sopt.cmdWindowOutputFrequency) || ...
        Sopt.cmdWindowOutputFrequency<1
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
    error('Time vector must contain at least two elements.')
end
nt = length(t);

% Check M,D,K matrices
if ~isa(model,'MechanicalSystem') && ~isa(model,'CMS_ROM') && ...
        ~(isfield(model,'M') && isfield(model,'D') && ...
        isfield(model,'K') && isfield(model,'nonlinear_elements'))
    error('Model must contain matrices M,D,K as fields.');
end
M = model.M;
D = model.D;
K = model.K;
n = size(M,1); % total number of generalized coordinates

% Check dimensions and reasonable values of coefficient matrices
if ~size(M,2)==n || ~size(D,1)==n || ~size(D,2)==n ...
        || ~size(K,1)==n || ~size(K,2)==n || ~isreal([M D K]) ...
        || any(any(isnan([M D K]))) || any(any(isinf([M D K])))
    error(['Coefficients M,D,K must be n x n matrices with neither ' ...
        'NaN nor INF, real entries.']);
end

% Check initial value dimensions and reasonable numbers
q0 = q0(:); u0 = u0(:); % enforce column vectors
if size(q0,1)~=n || size(u0,1)~=n
    error(['Initial coordinate and velocity must be provided as ' ...
        'vector of length matching the provided model.']);
elseif any(any(isnan([q0 u0]))) || any(any(isinf([q0 u0])))
    error(['Initial coordinate and velocity values must not be ' ...
        'NaN or INF.']);
end
if isempty(Sopt.qprevious)
    % No previous coordinate vector specified; set qprevious = q0
    qprevious = q0;
elseif numel(Sopt.qprevious)==n && ...
        ~(any(isnan(Sopt.qprevious(:))) || any(isinf(Sopt.qprevious(:))))
    qprevious = Sopt.qprevious(:); % enforce row vector
else
    error(['Specified previous coordinate vector must have appropriate ' ...
        'number of elements and reasonable values.']);
end
if isempty(Sopt.fnlprevious)
    % No previous coordinate vector specified; set qprevious = q0
    fnlprevious = 0*q0;
elseif numel(Sopt.fnlprevious)==n && ...
        ~(any(isnan(Sopt.fnlprevious(:))) || any(isinf(Sopt.fnlprevious(:))))
    fnlprevious = Sopt.fnlprevious(:); % enforce row vector
else
    error(['Specified previous coordinate vector must have appropriate ' ...
        'number of elements and reasonable values.']);
end

% Ensure force directions are provided for all nonlinear elements
NL = model.nonlinear_elements;
is3DcontactOnly = true;
isForceDirectionSpecified = true;
for inl=1:length(NL)
    if ~isfield(NL{inl},'force_direction')
        isForceDirectionSpecified = false;
    end
    if ~strcmpi(NL{inl}.type,'3D')
        is3DcontactOnly = false;
    end
end
if is3DcontactOnly && ~isForceDirectionSpecified
    % 3D contact only, presume unity force directions
    W = eye(n,3*length(NL));
    for inl=1:length(NL)
        NL{inl}.force_direction = W(:,(1:3)+3*(inl-1));
    end
elseif isForceDirectionSpecified
    % Check dimensions of specified force directions
    for inl=1:length(NL)
        if strcmpi(NL{inl}.type,'3D') && ...
                (size(NL{inl}.force_direction,1)~=n||...
                size(NL{inl}.force_direction,2)~=3)
            error(['In case of 3D contact, specified force direction ' ...
                'must be n x 3 matrix.']);
        elseif ~strcmpi(NL{inl}.type,'3D') && ...
                (size(NL{inl}.force_direction,1)~=n ||...
                size(NL{inl}.force_direction,2)~=1)
            error(['Force direction must be n x 1 vector ' ...
                '(except for 3D contact).']);
        end
    end
else
    error(['Force direction must be specified for all nonlinear ' ...
        'elements (except if you have 3D contact only).']);
end

% Check specified imposed gap / preload, friction coefficient, stiffness / 
% set default values
for inl = 1:length(NL)
    if strcmpi(NL{inl}.type,'3D')
        % Handle specified imposed gap
        if isfield(NL{inl},'imposedGap')
            if isa(NL{inl}.imposedGap,'function_handle')
                error(['Function handle not permitted for imposed gap ' ...
                    'in this function.']);
            elseif ~isscalar(NL{inl}.imposedGap)
                error(['Imposed gap specified per 3D contact element must' ...
                    ' be scalar.']);
            end
        else
            NL{inl}.imposedGap = [];
        end
        % Handle specified imposed gap velocity
        if isfield(NL{inl},'imposedGapVelocity')
            error('Imposed gap velocity not allowed in this function.');
        end
        % Handle specified preload
        if isfield(NL{inl},'preload')
            if isa(NL{inl}.preload,'function_handle')
                error('Preload must not be function handle.');
            elseif ~isscalar(NL{inl}.preload)
                error(['Preload specified per 3D contact element must' ...
                    ' be scalar.']);
            end
        else
            NL{inl}.preload = [];
        end
        % Handle specified friction coefficient
        if isfield(NL{inl},'frictionCoefficient')
            if isa(NL{inl}.frictionCoefficient,'function_handle')
                error('Friction coefficient must not be function handle.');
            elseif ~isscalar(NL{inl}.frictionCoefficient)
                error(['Friction coefficient specified per 3D contact element must' ...
                    ' be scalar.']);
            end
        else
            error(['A friction coefficient must be specified for 3D contact. ' ...
                'Set it to zero if you want frictionless contact.']);
        end
        % Handle specified (penalty) stiffness
        if isfield(NL{inl},'stiffness')
            if isa(NL{inl}.stiffness,'function_handle')
                error('Stiffness must not be function handle.');
            elseif ~isscalar(NL{inl}.stiffness)
                error(['Stiffness specified per 3D contact element must' ...
                    ' be scalar.']);
            end
        else
            error('A (penalty) stiffness must be specified for 3D contact.');
        end
    end
end

% Check for function handle to external forcing (explicit time dependence)
if nargin>=6
    func_fex = varargin{2};
    % Test for correct dimension and reasonable output at start/end times
    fex_S = func_fex(t(1)); fex_E = func_fex(t(end));
    if size(fex_S,1)~=n || size(fex_S,2)~=1 || ...
            size(fex_E,1)~=n || size(fex_E,2)~=1
        error(['Function handle to external forcing must return ' ...
            'column vector of length ' num2str(n) ' for scalar input' ...
            ' for specified model.']);
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
FNL = zeros(n,ntout);
Solinfo.NIT = zeros(1,ntout);
Solinfo.IEx = ones(1,ntout);
Solinfo.RESIDUUM = zeros(1,ntout);

% Prepare command window output
nPrintStats = Sopt.cmdWindowOutputFrequency+1;
kPrintStats = floor(linspace(2,nt,nPrintStats));
simProgress = linspace(0,1,nPrintStats);
simStopWatch = tic;

% Allocate solver statistics (results are averaged between outputs)
navg = ceil(nt/(nPrintStats-1));
numIterations = zeros(1,navg);
residuum = zeros(1,navg);
isConverged = ones(1,navg);
istep = 1; % counter for data storage
jstep = 1; % counter for command window output

%% START TIME STEP INTEGRATION

% Inform user about start of integration
if Sopt.verbosity>0
    fprintf(['Begin simulation with ' scheme ' scheme ' ...
        'for %.6g equidistant time steps (dt=%.5g).\n'],nt,dt)
    fprintf('--------------------\n');
end

% Auxiliary variables for time integration scheme: coefficients a1,...,a6 
% of Newmark quadrature rules as function of beta, gamma, dt
a1 = 1/beta/dt^2;
a2 = 1/beta/dt;
a3 = (1-2*beta)/2/beta;
a4 = gamma/beta/dt;
a5 = 1-gamma/beta;
a6 = (1-gamma/2/beta)*dt;
%   linear time-invaraint part of coefficient matrix of coordinate vector 
%   at next time level within algebraic equation system of one-step problem
S = a1*(1-alpha_m)*M+a4*(1-alpha_f)*D+(1-alpha_f)*K;

%% Loop over time steps
for k=1:nt
    %% SET POINT OF DEPARTURE FOR CURRENT ONE-STEP PROBLEM

    if k==1
        % Adopt initial coordinate and velocity vector specified by user
        qE = q0;
        uE = u0;

        % Evaluate nonlinear forces
        [fnlE,~] = nonlinear_forces(qE,uE,qprevious,fnlprevious,NL);

        % Determine initial forcing and accelerations
        fexE = real(feval(func_fex,t(k)));
        if rank(M)~=n
            % Mass matrix singular; presume zero initial accelerations
            accE = zeros(n,1);
        else
            % Mass matrix regular; initial accelerations follow from force
            % balance at time t=0
            accE = M\(fexE-D*uE-K*qE-fnlE);
        end

        % Store initial values, if applicable
        % NOTE: mod(k,Sopt.dataStorageFrequency) == 0 is the same as
        % Sopt.dataStorageFrequency==1, and we have istep=k=1;
        if mod(k,Sopt.dataStorageFrequency) == 0
            % Coordinate, velocity, nonlinear force
            Q(:,istep) = qE;
            U(:,istep) = uE;
            FNL(:,istep) = fnlE;

            % Increment storage counter
            istep = istep+1;
        end

        % Continue with next time level
        continue;
    else
        % Adopt values from previous time level
        qS = qE;
        uS = uE;
        accS = accE;
        fnlS = fnlE;
        fexS = fexE;
    end

    % Evaluate external force at time tE
    fexE = real(feval(func_fex,t(k)));

    %% SOLVE ONE-STEP PROBLEM

    % Predict displacement assuming uniform motion
    qpre = qS + uS*dt + 0.5*dt^2*accS;
    % qpre = qS + uS*dt + (0.5-beta)*dt^2*accS;

    % Evaluate right-hand side
    b = (1-alpha_f)*fexE+alpha_f*fexS...
        -(1-alpha_m)*M*(-a1*qS-a2*uS-a3*accS)...
        -(1-alpha_f)*D*(-a4*qS+a5*uS+a6*accS)...
        -alpha_m*M*accS...
        -alpha_f*D*uS...
        -alpha_f*K*qS;

    % Solve for displacement vector and nonlinear force vector at end of 
    % time step
    [qE,fnlE,output] = ...
        solveOneStepProblem(qpre,qS,uS,accS,a4,a5,a6,fnlS,S,b,...
        NL,alpha_f,Sopt);

    % Evaluate velocity and acceleration vectors according to Newmark
    % quadrature rules
    uE = a4*(qE-qS) + a5*uS + a6*accS;
    accE = a1*(qE-qS) - a2*uS - a3*accS;

    %% OUTPUT

    % Store solver output for statistics
    numIterations(jstep) = output.iter;
    residuum(jstep) = output.residuum;
    isConverged(jstep) = output.isConverged;

    % Store results at given data storage frequency
    if mod(k,Sopt.dataStorageFrequency) == 0
        % Coordinate, velocity, nonlinear force
        Q(:,istep) = qE;
        U(:,istep) = uE;
        FNL(:,istep) = fnlE;

        % Solver statistics output
        istep = istep+1;
        Solinfo.NIT(istep) = numIterations(jstep);
        Solinfo.RESIDUUM(istep) = residuum(jstep);
        Solinfo.IEx(istep) = isConverged(jstep);
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
                        ' in terms of L-infinity norm of residuum vector.']);
                    fprintf('\n');
                end
            end
        end
        if Sopt.verbosity>0
            fprintf('Progress = %.2f %%.', simProgress(1)*100 );
        end
        if Sopt.verbosity>1
            fprintf(' Number of Newton iterations = %2.1f.', ...
                mean(numIterations(1:jstep)));
            if Sopt.verbosity>2
                fprintf(' residuum = %.2g', mean(residuum(1:jstep)));
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
end

function [qE,fnlE,output] = ...
    solveOneStepProblem(qE,qS,uS,accS,a4,a5,a6,fnlS,S,b,NL,alpha_f,Sopt)
%========================================================================
% This function solves the algebraic equation system of the one-step
% problem, formulated as
%       S*qE + (1-alpha_f)*fnlE + alpha_f*fnlS = b,
% using the Netwon method with Cholesky factorization of the analytical
% Jacobian. Note that fnlE depends on qE explicitly and implicitly due to
% its dependence on uE (uE is related to qE via the quadrature rule).
%========================================================================
for ik=1:Sopt.itermax
    % Evaluate current iterate of velocity uE
    uE = a4*(qE-qS) + a5*uS + a6*accS;
    duE_dqE = a4;
    % Evaluate nonlinear forces and their gradients
    [fnlE,dfnlE_dqE,dfnlE_duE] = ...
        nonlinear_forces(qE,uE,qS,fnlS,NL);
    % Evaluate residual, set up Jacobian
    r = S*qE+(1-alpha_f)*fnlE+alpha_f*fnlS-b;
    dr_dqE = S+(1-alpha_f)*(dfnlE_dqE+dfnlE_duE*duE_dqE);
    % Newton update using Cholesky factorization
    Chol = chol(dr_dqE);
    qE = qE-Chol\(Chol'\r);
    % Evaluate new residuum
    r = (S*qE+(1-alpha_f)*fnlE+alpha_f*fnlS-b);
    % If L-infinity norm is below tolerance, terminate Newton iterations
    if max(abs(r))<Sopt.tol
        break;
    end
end
output.iter = ik;
output.residuum = max(abs(r));

% Check whether convergence criterion is met
if max(abs(r))>=Sopt.tol
    fprintf('noconv , residual %.2e\n', max(abs(r)));
    output.isConverged = 0;
else
    output.isConverged = 1;
end

end

function [fnlE,dfnlE_dqE,dfnlE_duE] = nonlinear_forces(qE,uE,qS,fnlS,NL)
%% Initialize nonlinear force vector and its derivatives
n = size(qE,1);
fnlE = zeros(n,1);
dfnlE_dqE = zeros(n,n);
dfnlE_duE = zeros(n,n);

% Evaluate all nonlinear elements
for inl=1:length(NL)

    % Determine force direction and calculate displacement and velocity of
    % nonlinear element
    w = NL{inl}.force_direction;
    qnl = w'*qE;
    unl = w'*uE;

    % Initialize derivativives
    % dfnldqnl = zeros(size(qnl,1)); % not needed as specified in all types
    dfnldunl = zeros(size(qnl,1));

    % Determine nonlinear force and its derivatives depending on type
    switch lower(NL{inl}.type)
        case 'unilateralspring'
            fnl = NL{inl}.stiffness*...
                (qnl+NL{inl}.gap).*...
                double(qnl+NL{inl}.gap<=0);
            dfnldqnl = NL{inl}.stiffness*...
                double(qnl+NL{inl}.gap<=0);
        case 'quadraticdamper'
            fnl = NL{inl}.damping*qnl.^2.*unl;
            dfnldqnl = NL{inl}.damping*2*qnl.*unl;
            dfnldunl = NL{inl}.damping*qnl.^2;
        case 'elasticdryfriction'
            % Stiffness and friction limit force parameters
            k = NL{inl}.stiffness;
            fc = NL{inl}.friction_limit_force;

            % Predict force assuming sticking contact
            fnlPre = w' * fnlS + k * (qnl-w'*qS);

            % Correct force distinguishing stick and stlip
            if abs(fnlPre)>fc
                % slip
                fnl = sign(fnlPre)*fc;
                dfnldqnl = 0;
            else
                % stick
                fnl = fnlPre;
                dfnldqnl = k;
            end
        case {'3d'}
            % Friction coefficient, (penalty) stiffness
            % NOTE: Isotropic behavior is presumed here by using scalar
            % values.
            mu = NL{inl}.frictionCoefficient;
            kc = NL{inl}.stiffness;

            % Handle specified imposed gap / preload
            if ~isfield(NL{inl},'preload')
                NL{inl}.preload = [];
            end
            if ~isfield(NL{inl},'imposedGap')
                NL{inl}.imposedGap = [];
            end
            if (~isempty(NL{inl}.preload)&&...
                    ~isempty(NL{inl}.imposedGap)) || ...
                    (isempty(NL{inl}.preload)&&isempty(NL{inl}.imposedGap))
                error('Either preload or imposed gap must be specified.');
            elseif ~isempty(NL{inl}.preload)
                fN0 = NL{inl}.preload;
            else
                fN0 = - kc*NL{inl}.imposedGap;
            end

            % Evaluate normal contact force
            qN = qnl(1); dqN = [1 0 0];
            fN = (kc*qN-fN0)*double(kc*qN-fN0<0);
            dfN = kc*double((kc*qN-fN0)<0)*dqN;

            % Determine tangential contact force
            if fN==0
                % separation
                fT = [0;0];
                dfT = zeros(2,3);
            else
                % contact; predict force assuming stick
                qT = qnl(2:3); dqT = [0 1 0; 0 0 1];
                fTpre = w(:,2:3)'*fnlS + kc*(qT-w(:,2:3)'*qS);
                dfTpre = kc*dqT;
                norm_fTpre = sqrt(fTpre'*fTpre);

                % Correct force distinguishing stick and slip
                if norm_fTpre<=-mu*fN
                    % stick
                    fT = fTpre;
                    dfT = dfTpre;
                else
                    % slip
                    fT = -mu*fN*fTpre/norm_fTpre;
                    dfT = - mu*fTpre/norm_fTpre*dfN ...
                        - mu*fN/norm_fTpre*dfTpre ...
                        + mu*fN*fTpre/norm_fTpre^3*(fTpre'*dfTpre);
                end
            end

            % Account for preload and assemble 3D contact force
            fnl = [fN + max(fN0,0);fT];
            dfnldqnl = [dfN;dfT];

        otherwise
            error(['Unknown nonlinear element ' ...
                NL{inl}.type '.']);
    end

    % Add contribution to nonlinear force and partial derivatives
    fnlE = fnlE + w*fnl;
    dfnlE_dqE = dfnlE_dqE + w*dfnldqnl*w';
    dfnlE_duE = dfnlE_duE + w*dfnldunl*w';
end

end
