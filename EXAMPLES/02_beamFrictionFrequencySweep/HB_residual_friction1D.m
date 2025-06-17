%========================================================================
% DESCRIPTION:
% This function is a variant of NLvib's HB_residual.m function, 
% implementing 1D dry frictional contact (Coulomb conditions) via the 
% Dynamic Lagrangian approach, formulated on velocity level, and
% restricting the harmonic set to the odd harmonics.
% 
% SETTING THE (CONTACT) PARAMETERS
% The structures describing the nonlinear elements, i.e., the entries of
% the cell array system.nonlinear_elements, must all have field 'type'
% equal to 'friction', and specify two mandatory parameters as fields, 
% namely the 'friction_limit_force' and the Dynamic Lagrangian parameter
% 'epsDL'. Trivial nonlinear force directions are presumed in the sense
% that friction conditions are imposed on the first nc coordinates of the
% model, where nc=length(system.nonlinear_elements).
% 
% ELIMINATION OF EVEN HARMONICS
% The even harmonics, including the zeroth, are removed from the set of
% unknowns and HB equations. This is justified as long as these are not 
% directly driven by external forcing, nor a corresponding symmetry 
% breaking bifurcation occurs. The elimination of the even harmonics 
% simplifies the problem. This also removes the ambiguity of the zeroth 
% harmonic in the case of permanently sticking contacts (where the slip
% displacement is arbitrary within the Coulomb limits). For each nonlinear 
% element, a value for the limit friction force and for the Dynamic 
% Lagrangian parameter 'epsDL' must be specified. 
% 
% NOTE: An alternative to the velocity-level formulation here is a
% displacement-level implementation by marching the incremental law
% governing the evolution of the hysteresis. This displacement-level
% marching procedure can be expected to yield faster harmonic convergence,
% at the cost that the contact force at different time levels cannot be 
% evaluated independently, resulting in a non-vectorized, sequential
% implementation.
%========================================================================
% This file is part of NLvib.
%
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
%
% COPYRIGHT AND LICENSING:
% NLvib Copyright (C) 2025  
%   Malte Krack  (malte.krack@ila.uni-stuttgart.de) 
%   Johann Gross (johann.gross@ila.uni-stuttgart.de)
%   University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY.
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
function [R,dR,Q,Fnl] = ...
    HB_residual_friction1D(X,system,H,N,analysis_type,varargin)
%% Handle input variables

% Check for 1D frictional contact and determine parameters
NL = system.nonlinear_elements;
fc = zeros(length(NL),1);
epsDL = zeros(length(NL),1);
for inl=1:length(NL)
    if ~strcmpi(NL{inl}.type,'friction')
        error('This function is for friction contact only.');
    end
    fc(inl) = NL{inl}.friction_limit_force;
    epsDL(inl) = NL{inl}.epsDL;
end

% System matrices
M = system.M;
D = system.D;
K = system.K;

% number of degrees of freedom
n = size(M,1);

% Conversion from sine-cosine to complex-exponential representation,
% accounting for non-trivial harmonic set
Hset = 1:2:H;
numH = numel(Hset);
ID = repmat(1:n,1,numH)+n*kron(Hset,ones(1,n));
IC = repmat(1:n,1,numH)+n*kron(0:2:2*(numH-1),ones(1,n)); IS = IC+n;
dX = eye(length(X));
Q = zeros(n*(H+1),1);   dQ = zeros(size(Q,1),size(dX,2));
Q(ID) = X(IC)-1i*X(IS); dQ(ID,:) = dX(IC,:)-1i*dX(IS,:);

% Handle analysis type
if nargin<=4 || isempty(analysis_type)
    % Default analysis: frequency response
    analysis_type = 'frf';
end
switch lower(analysis_type)
    case {'frf','frequency response'}
        % Frequency response analysis: X = [Qtrunc;Om]

        % Excitation 'excitation' is the fundamental harmonic of the
        % external forcing
        Fex1 = system.Fex1;
        % Setup excitation vector
        Fex = zeros(n*(H+1),1);
        h = 1;% Forcing of h-th harmonic
        Fex(h*n+(1:n)) = Fex1;
        dFex = zeros(size(Fex,1),length(X));

        % Excitation frequency
        Om  =  X(end);
        dOm = dX(end,:);

        % Scaling of dynamic force equilibrium
        if length(varargin)<2 || isempty(varargin{2})
            fscl = 1;
        else
            fscl = varargin{2};
        end
    otherwise
        error('So far only frequency response analysis is implemented.');
end

%% Evaluate linear part of dynamic force balance
Rlinc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*Q - Fex;
dRlinc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*dQ - dFex + ...
    ( -2*Om*kron(diag((0:H).^2),M) + 1i*kron(diag(0:H),D) )*Q*dOm;

%% Compute the Fourier coefficients of the nonlinear forces and the
% Jacobian using AFT in conjunction with the Dynamic Lagrangian approach
[Fnl,dFnl] = HB_nonlinear_forces_DLFT(Q,dQ,Om,dOm,Rlinc,dRlinc,H,N,...
    fc,epsDL);

%% Assemble the residual and the Jacobian
Rc = Rlinc + Fnl;
dRc = dRlinc + dFnl;

% Scale dynamic force equilibrium (useful for numerical reasons)
Rc = 1/fscl*(Rc);
dRc = 1/fscl*(dRc);

% Conversion from complex-exponential to sine-cosine representation, accounting for
% truncated harmonic set
R = zeros(size(X,1)-1,1); dR = zeros(size(X,1)-1,length(X));
R(IC) = real(Rc(ID)); dR(IC,:) = real(dRc(ID,:));
R(IS) = -imag(Rc(ID)); dR(IS,:) = -imag(dRc(ID,:));

end
%% Computation of the Fourier coefficients of the nonlinear forces and the 
% Jacobian using AFT
function [F,dF] = ...
    HB_nonlinear_forces_DLFT(Q,dQ,Om,dOm,Rlin,dRlin,H,N,fc,epsDL)

% Initialize output
F = zeros(size(Q));
dF = zeros(size(F,1),size(dQ,2));

% Specify equidistant time samples along one period
tau = (0:2*pi/N:2*pi-2*pi/N)';

% Loop over nonlinear contact elements
% NOTE: This loop could be vectorized exploiting independence of operations
% on individual nonlinear elements.
n = size(Q,1)/(H+1);
nc = length(fc);
for inl=1:nc
    % Index to Fourier coefficients associated with 'inl'
    INL = inl:n:size(Rlin,1);

    % Set up Fourier coefficients of gap velocity and argument of projection
    Gam = 1i*Om*(0:H)'.*Q(INL);
    dGam = repmat(1i*Om*(0:H)',1,size(dQ,2)).*dQ(INL,:) + ...
        repmat(1i*(0:H)'.*Q(INL),1,size(dQ,2)).*repmat(dOm,H+1,1);
    ARG = Rlin(INL) - epsDL*Gam;
    dARG = dRlin(INL,:) - epsDL*dGam;

    % Apply inverse discrete Fourier transform
    H_iDFT = exp(1i*tau*(0:H));
    arg  = real(H_iDFT*ARG);
    darg = real(H_iDFT*dARG);

    % Calculate nonlinear contact force as projection onto the admissible
    % set 'C', of contact forces, which is the interval [-fc,fc] in the 
    % case of Coulomb dry friction. Note the similarity of this with the
    % local function 'projC' in NLstep's 'sim_friction1D.m'.
    isSliding = abs(arg)>=fc(inl);
    lam = arg.*double(~isSliding) + fc(inl)*sign(arg).*double(isSliding);
    dlam = darg.*repmat(double(~isSliding),1,size(dQ,2));

    % Account for sign convention
    fnl = -lam;
    dfnl = -dlam;

    % Apply forward discrete Fourier transform using FFT and truncation to
    % half-spectrum notation
    Fnlc = fft(fnl(end-N+1:end))/N;
    dFnlc = fft(dfnl(end-N+1:end,:))/N;
    Fnl = [Fnlc(1);2*Fnlc(2:H+1)];
    dFnl = [dFnlc(1,:);2*dFnlc(2:H+1,:)];

    % Insert Fourier coefficients associated wth 'inl' into global vector
    F(INL) = Fnl;
    dF(INL,:) = dFnl;
end
end
