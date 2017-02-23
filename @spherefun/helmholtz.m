function u = helmholtz(f, K, m, n)
%HELMHOLTZ   Fast Helmholtz solver for the sphere.
%    U = HELMHOLTZ(F, K, N) solves U_xx + U_yy + U_zz + K^2U = F on the sphere
%    for U with a discretization of size N x N. F should be a SPHEREFUN and the
%    solution is returned as a SPHEREFUN.
%
%    HELMHOLTZ(F, K, M, N) same as HELMHOLTZ(F, K, N), but with a
%    discretization of size M x N.
%
%  Example:
%    K = 100; m = 1000; n = m;
%    f = spherefun( @(x,y,z) cos(x.*y.*z) );
%    u = spherefun.helmholtz(f, K, m, n);
%    plot( u )

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPERS NOTE:
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in
% the theta- and lambda-direction.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation decouples into n
% linear systems with pentadiagonal banded matrices, by exploiting
% symmetries these can reduced to tridiagonal systems.
%
% SOLVE COMPLEXITY:    O(M*N)  with M*N = total degrees of freedom

% TODO: 
% Make this code adaptive! 

% Solve standard Helmholtz equation. This parameter is kept for developers.
c = 1;

% If the call is helmholtz(f, K, m), then set n
if ( nargin < 4 )
    n = m;
end

if ( K == 0 )
    u = spherefun.poisson(f, 0, m, n);
    return
end

% Check for eigenvalues:
e = eig( [ -1 K^2 ; 1 0 ] ); % Is K = sqrt(l*(l+1)), where l is an integer?
e = e(e>0); e = e(abs( e - round(e) ) < 1e-13 );
if ( ~isempty( e ) )
    error('SPHEREFUN:HELMHOLTZ:EIGENVALUE',...
            'There are infinitely many solutions since K is an eigenvalue of the Helmholtz operator.')
end

% If m or n are non-positive then throw an error
if ( ( m <= 0 ) || ( n <= 0 ) )
    error('CHEBFUN:SPHEREFUN:HELMHOLTZ:badInput',...
        'Discretization sizes should be positve numbers');
end

% If m and n are 1, the solution is easy.  It's just a constant given by
% 1/K^2*mean(F).
if ( ( m == 1 ) && (n == 1 ) )
    f = spherefun(f);
    u = mean2(f)/K^2 + 0*f;
    return;
end

% Make m even so that the pole at theta=0 is always sampled.
m = m + mod(m,2);

% Construct useful spectral matrices:
Im = speye(m);

% Please note that DF1m here is different than trigspec.diff(m,1) because we
% take the coefficient space point-of-view and set the (1,1) entry to be
% nonzero.
DF1m = trigspec.diffmat(m, 1, 1);
DF2m = trigspec.diffmat(m, 2);
DF2n = trigspec.diffmat(n, 2);

% Multiplication for sin(theta).*cos(theta):
% Below is equivalent to
% Mcossin = spdiags(.25i*[-ones(m, 1) ones(m, 1)], [-2 2], m, m);
cfs = trigtech(@(theta) sin(pi*theta).*cos(pi*theta));
Mcossin = trigspec.multmat(m, cfs.coeffs);

% Multiplication for sin(theta)^2:
% Below is equivalent to
% Msin2 = spdiags(.5*[-.5*ones(m, 1) ones(m, 1) -.5*ones(m, 1)], [-2 0 2], m, m);
cfs = trigtech(@(theta) sin(pi*theta).^2);
Msin2 = trigspec.multmat(m, cfs.coeffs);

% Forcing term:
if ( isa(f, 'function_handle') )
    % Underlying discretization grid:
    lam0 = trigpts(n,[-pi, pi]); 
    th0 = trigpts(m,[-pi, pi]); 
    [rhs_lam, rhs_th] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_th);      % Get (trigvals,trigvals) of rhs
    F = trigtech.vals2coeffs(F);        % Get in Fourier basis
    F = trigtech.vals2coeffs(F.').';
elseif ( isa(f, 'spherefun') )
    F = coeffs2(f, n, m);
end

floorm = floor(m/2);
floorn = floor(n/2);

% Multiple rhs by sin(th)^2 and divide by K^2:
F = Msin2 * F / K^2;

% Want to solve
%    X L^T + X DF^T = F
% subject to zero integral constraints, i.e.,  w^T X_0  = 0.
% The matrix equation decouples because DF is diagonal.  Also
% We take advantage of even/odd symmetry to get a factor of 2 speed up.  
% TODO: if f is real there is another factor of 2 speed up possible.

% Matrix for solution's coefficients:
CFS = zeros(m, n);

% Form discretization of the theta-dependent operator:
L = c*(Msin2*DF2m + Mcossin*DF1m)/K^2 + Msin2;
scl = c*diag(DF2n)/K^2;

% We take advantage of even/odd symmetry to get a factor of 2 speed up.  
% TODO: if f is real there is another factor of 2 speed up possible.

%
% Solve for the negative even and negative odd modes:
%

% Matrices for the odd expansions in theta.
ii_o = 2-mod(floorm,2):2:m;
L_o = L(ii_o,ii_o);
Im_o = Im(ii_o,ii_o);
% Matrices for the even expansions in theta.
ii_e = 1+mod(floorm,2):2:m;
L_e = L(ii_e,ii_e);
Im_e = Im(ii_e,ii_e);
k_neg = floorn+1:-1:1;
for k = k_neg;
    CFS(ii_o,k) = (L_o + scl(k)*Im_o) \ F(ii_o,k);
    CFS(ii_e,k) = (L_e + scl(k)*Im_e) \ F(ii_e,k);
end

% Fill in the positive odd modes from the negative odd modes.
ii = floorn+2:1:n;
CFS(ii_o, ii) = (-1)^floorn*bsxfun(@times,(-1).^ii,-conj(CFS(ii_o, k_neg(2:numel(ii)+1))));

% Fill in the negative even modes from the negative even modes.
CFS(ii_e, ii) = (-1)^floorn*bsxfun(@times,(-1).^ii,-conj(CFS(ii_e, k_neg(2:numel(ii)+1))));

% Now, convert to a spherefun object:
u = spherefun.coeffs2spherefun(CFS);

end