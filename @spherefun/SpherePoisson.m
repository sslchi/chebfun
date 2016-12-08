function SpherePoisson

f = @(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th);
exact = @(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
    sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th);
u = spherefun.poisson(f, 0, 100);

norm(spherefun(exact) - u)
mean2(u)

m = 101;
n = 101;
% n = 2*n+1;
const = 0;

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

% Underlying discretization grid:
lam0 = trigpts(n,[-pi, pi]);
th0 = trigpts(m,[-pi, pi]);

% Forcing term:
if ( isa(f, 'function_handle') )
    [rhs_lam, rhs_theta] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_theta);
    %tol = 1e5*max(abs(F(:)))*chebfunpref().cheb2Prefs.chebfun2eps;
    tol = chebfunpref().cheb2Prefs.chebfun2eps;
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
k = floor(n/2) + 1;
floorm = floor(m/2);
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
ii = 1:m;
meanF = en(ii)*F(ii, k)/en(floor(m/2)+1);

% Check that the mean of F is zero (or close enough).  If it is not then
% issue a warning
if ( abs(meanF) > tol )
    warning('CHEBFUN:SPHEREFUN:POISSON:meanRHS',...
        ['The integral of the right hand side may not be zero, which is '...
        'required for there to exist a solution to the Poisson '...
        'equation. Subtracting the mean off the right hand side now.']);
end
F(floor(m/2)+1,k) = F(floor(m/2)+1,k)-meanF;

% Multiply the right hand side by (sin(theta)).^2
F = Msin2*F;
[U, S, V] = svd( F, 0);
r = find(diag(S)/S(1,1)>1e-14,1,'last');
U = U(:,1:r); S = S(1:r,1:r); V = V(:,1:r); 
C = U*diag(sqrt(diag(S))); 
R = V*diag(sqrt(diag(S)));
R = R.';

% Form discretization of the theta-dependent operator:
L = Msin2*DF2m + Mcossin*DF1m;

% matrix equation is
% L*X + X*DF2n = F

% Solve
% L*X(:,[1:floor(n/2), floor(n/2)+2:n]) + X(:,[1:floor(n/2),
% floor(n/2)+2:n])*DF2n = F(:,[1:floor(n/2), floor(n/2)+2:n])
idx = [1:floor(n/2), floor(n/2)+2:n];
Dt = DF2n(idx,idx);
CFS = zeros(m,n);

% We want this, but we need to compute it using ADI:
% CFS(:,idx) =  lyap( L, Dt, -F(:,idx) );

% We now want to solve the above using factor independent ADI:
eA = eig(full(L));  % replace later
eB = eig(full(Dt)); % replace later
Int = [min(eA), max(eA), -max(eB), -min(eB)];

[~, Tinv, gam, ~] = mobiusTdisk(Int); %gam and cross ratio

% we estimate small 1/kp for use in ellipk
if abs(1/gam) < 1e-7
    kp = small_gam(gam);
else
    kp = sqrt(1-sqrt(1-(1/gam)^2));
end
[~, K] = ellipk( kp);
kp = sqrt(1-(1/gam)^2);%need this for ellipj function in fADIrank1

%approx max number ADI iterations
Nmax = ceil(1/pi^2*log(4/tol)*log(4*gam));

%--------Main loop: performs ADI ------------------------
YY =[];
ZZ =[];
DD = [];
A = L;
B = Dt;
for k = 1:r;
    [Z, D, Y] = rank1_fADI(C(:,k), R(k,idx), A, B, gam, Tinv, K, kp, Nmax);
    ZZ = [ZZ Z];
    YY = [YY Y];
    DD = [DD; D];
end

%share diag
ds = sign(DD);
d = sqrt(abs(DD));
ZZ = ZZ*diag(d); 
YY = conj(YY);
YY = (YY*diag(d.*ds));

CFS(:,idx) = ZZ*YY.';

% Now do the equation where we need the integral constraint:
% We will take X_{n/2+1,:} en = 0.

% Second, solve:
k = floor(n/2) + 1;
ii = [1:floorm floorm+2:m];
CFS(:, k) = [ en ; L( ii, :) ] \ [ 0 ; F(ii, k) ];
u = spherefun.coeffs2spherefun( CFS ) + const;

norm(spherefun(exact) - u)

end

%factored ADI on a rank 1 piece
function [Z, D, Y] = rank1_fADI(C, R, A, B, gam,Tinv, K,kp, N)

%compute ADI shifts
idx1 = 1:N;
[~, ~, dn] = ellipj((2*idx1-1)*K/(2*N), kp);
p1 = gam*dn; %optimal shift parameters
q =  -Tinv(-p1) ;
p = Tinv(p1) ;

%perform ADI
I = speye(size(A));
II = speye(size(B));
Z(:,1) = (A-p(1)*I)\C;
Y(1,:) = R/(B-q(1)*II);
for j = 1:N-1
    Z(:, j+1) = Z(:, j)+(A-p(j+1)*I)\((p(j+1)+q(j))*Z(:,j));
    Y(j+1,:) =  Y(j,:) + ((q(j+1)+p(j))*Y(j,:)/(B-q(j+1)*II));
end

%return pivots and avoid column-stacking Y
D = -(p+q).';
Y = Y.';
end

function g = small_gam(gamma)
% Taylor series estimate 1-sqrt(1-(1/gam)^2))
% when 1/gam is small.

x = 1/gamma;
cfs = 1./(sqrt(2)*[1  8 128 1024 ]);
x = [ x; x^3; x^5; x^7 ];
g = cfs*x;
end

function [T, Tinv, gam, M] = mobiusTdisk(I)
%given I = [a b c d] where [a b] and [c d] are two disjoint intervals
% on the real line, V are the four points [-gamma, -1 1 gamma]

a = I(1);
b = I(2);
c1 = I(3);
d1 = I(4);

%parameters
M = abs(c1-a)*abs(d1-b)/(abs(c1-b)*abs(d1-a));
gam = -1+2*M+2*sqrt(M^2-M);
A1 = -gam*a*(1-gam)+gam*(c1-gam*d1)+c1*gam-gam*d1;
B1 = -gam*a*(c1*gam-d1)-a*(c1*gam-gam*d1)-gam*(c1*d1-gam*d1*c1);
C1 = a*(1-gam)+gam*(c1-d1)+c1*gam-d1;
D11 = -gam*a*(c1-d1)-a*(c1-gam*d1)+c1*d1-gam*d1*c1;


T = @(z) (A1*z+B1)./(C1*z+D11);
Tinv = @(z) (D11*z-B1)./(-C1*z+A1);
end