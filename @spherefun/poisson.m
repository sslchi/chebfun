function u = poisson(f, const, m, n)
%POISSON   Fast Poisson solver for the sphere.
%   POISSON(F, C, N) solves laplacian(U) = F on the unit sphere, which in 
%   spherical coordinates (lam, th) is
%
%     sin(th)^2U_{th,th} + sin(th)cos(th)U_th + U_{lam,lam} = sin(th)^2*F
%
%   The equation is discretized on an N x N grid in spherical coordinates.
%   The integral of F is assumed to be zero, which is the compatibility
%   constraint for there to exist a solution to the Poisson problem on the
%   sphere. The mean value of the solution U is set to C.  This
%   function returns a SPHEREFUN representing the solution.
%
%   POISSON(F, C, M, N) is the same as POISSON(F, C, N), but with a
%   discretization of size M x N.
%
% EXAMPLE:
%   f = @(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th);
%   exact = @(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
%             sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th);
%   u = spherefun.poisson(f, 0, 100);
%   norm(spherefun(exact) - u)
%   mean2(u)

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

if ( nargin < 4 )
    n = m;
end

% Construct useful spectral matrices:
% Please note that DF1m is different than trigspec.diff(m,1) because we 
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
Im = speye(m);
scl = diag(DF2n); 

% There is a factor of 4 speed up here, by taking account of real 
% solution, and even/odd symmetry.

% Forcing term:
if ( isa(f, 'function_handle') )
    % Underlying discretization grid:
    lam0 = trigpts(n,[-pi, pi]); 
    th0 = trigpts(m,[-pi, pi]); 
    [rhs_lam, rhs_theta] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_theta);
    tol = 1e5*max(abs(F(:)))*chebfunpref().cheb2Prefs.chebfun2eps;
    F = trigtech.vals2coeffs(F);
    F = trigtech.vals2coeffs(F.').';
elseif ( isa(f, 'spherefun') )
    tol = 1e5*vscale(f)*chebfunpref().cheb2Prefs.chebfun2eps;
    F = coeffs2(f, n, m);
elseif ( isa( f, 'double' ) )
    tol = 1e5*chebfunpref().cheb2Prefs.chebfun2eps;
    F = f;       % Get trigcoeffs2 of rhs.
end

% First, let's project the rhs to have mean zero:
floorm = floor(m/2);
floorn = floor(n/2);
k = floorn + 1;
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
ii = 1:m;
meanF = en(ii)*F(ii, k)/en(floorn+1);

% Check that the mean of F is zero (or close enough).  If it is not then
% issue a warning
if ( abs(meanF) > tol )
    warning('CHEBFUN:SPHEREFUN:POISSON:meanRHS',...
       ['The integral of the right hand side may not be zero, which is '...
        'required for there to exist a solution to the Poisson '...
        'equation. Subtracting the mean off the right hand side now.']);
end        
F(floorm+1,k) = F(floorm+1,k)-meanF;

% Multiply the right hand side by (sin(theta)).^2
F = Msin2*F;

% Matrix for solution's coefficients:
CFS = zeros(m, n);

% Form discretization of the theta-dependent operator:
L = Msin2*DF2m + Mcossin*DF1m;
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
k_neg = floorn:-1:1;
for k = k_neg;
    CFS(ii_o,k) = (L_o + scl(k)*Im_o) \ F(ii_o,k);
    CFS(ii_e,k) = (L_e + scl(k)*Im_e) \ F(ii_e,k);
end
% Solve for the zero mode. Odd terms are easy
k = floorn+1;
CFS(ii_o,k) = (L_o + scl(k)*Im_o) \ F(ii_o,k);

% Solve for the zero mode. Even terms require some care since the operator L
% is singular
ii = [1:floor(m/4) floor(m/4)+2:numel(ii_e)];
CFS(ii_e, k) = [ en(ii_e) ; L_e( ii, :) ] \ [ 0 ; F(ii_e(ii), k) ];

% Fill in the positive odd modes from the negative odd modes.
ii = floorn+2:1:n;
CFS(ii_o, ii) = (-1)^floorm*bsxfun(@times,(-1).^ii,-conj(CFS(ii_o, k_neg(1:numel(ii)))));

% Fill in the negative even modes from the negative even modes.
CFS(ii_e, ii) = (-1)^floorm*bsxfun(@times,(-1).^ii,-conj(CFS(ii_e, k_neg(1:numel(ii)))));

u = spherefun.coeffs2spherefun( CFS ) + const; 

end