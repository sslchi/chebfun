function [f, p] = nudct(c, x)
% NUDCT 

n = size(c, 1);
coeffs = c;
magic_number = 0.995; 
idx = abs(x)>magic_number; 

th = real(acos(x)/2/pi);

c(2:end) = c(2:end)/2;
c = [c(end:-1:2);c];
tol = eps; 
N = size(c,1);
gam = 1/2; 
K = 16; 
s = round(N*th);
t = mod(s, N) + 1;
[u, v] = constructAK( th, (n-1:-1:-(n-1))', K );
In = speye(size(c,1));
p = @(c) sum(u.*(In(t,:)*fft( ifftshift(flipud(repmat(c,1,K).*v),1),[],1 )),2);
f = p(c);

% Replace values near the endpoints: 
f(idx) = chebtech.clenshaw(x(idx),coeffs); 

end

function [s, t, gam] = FindAlgorithmicParameters( x, N )
% Find algorithmic parameters.
%
%  s/N = closest equispaced gridpoint to x
%  t/N = closest equispaced FFT sample to x
%  gam = perturbation parameter

s = round(N*x);
t = mod(s, N) + 1;
gam = norm( N*x - s, inf);
end

function [u, v] = constructAK( x, omega, K )
% Construct a low rank approximation to
%
%     A_{jk} = exp(-2*pi*1i*(x-s_j/N)*k), 0<=j,k<=N-1,
%
% where |x_j-j/N|<= gam <=1/2.  See [1].

N = size(omega,1);
[s, ~, gam] = FindAlgorithmicParameters( x, N );
er = N*x - s;
% scl = exp(-1i*pi*er);
u = (ChebP(K-1,er/gam)*Bessel_cfs(K, gam));
v = ChebP(K-1, 2*omega/N);
end

function cfs = Bessel_cfs(K, gam)
% The bivarate Chebyshev coefficients for the function f(x,y) = exp(-i*x.*y)
% on the domain [-gam, gam]x[0,2*pi] are given by Lemma A.2 of Townsend's
% DPhil thesis.
arg = -gam*pi/2;
[pp,qq] = meshgrid(0:K-1);
cfs = 4*(1i).^qq.*besselj((pp+qq)/2,arg).*besselj((qq-pp)/2, arg);
cfs(2:2:end,1:2:end) = 0;
cfs(1:2:end,2:2:end) = 0;
cfs(1,:) = cfs(1,:)/2;
cfs(:,1) = cfs(:,1)/2;
end

function T = ChebP( n, x )
% Evaluate Chebyshev polynomials of degree 0,...,n at points in x. Use the
% three-term recurrence relation:
N = size(x, 1);
T = zeros(N, n+1);
T(:,1) = 1;
if ( n == 0 )
    return
end
T(:,2) = x;
twoX = 2*x;
for k = 2:n
    T(:,k+1) = twoX.*T(:,k) - T(:,k-1);
end
end