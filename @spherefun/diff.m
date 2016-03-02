function f = diff( f, varargin )
% DIFF  Tangential derivative of a spherefun in Cartesian coordinates.
%
%  F = DIFF( F ) computes the first tangential derivative of F with respect
%  to x.  This is the projection of the surface gradient of f in the
%  x-direction
%  
%  F = DIFF( F, DIM )  computes the first tangential derivative of F. If
%  DIM = 1, the tangential derivative is taken in the x-direction. If DIM =
%  2, the tangential derivative is taken in the y-direction and if DIM = 3,
%  the tangential derivative is taken in the z-direction.
%
%  F = DIFF( F, DIM, K) computes the kth tangential derivatives of F in the
%  variable given by DIM.
%
%  See also GRADIENT, LAPLACIAN

% Check for empty:
if ( isempty( f ) )
    return
end 

% Parse user inputs:
if ( nargin == 1 )
    dim = 1;
    K = 1;
elseif ( nargin == 2 )
    dim = varargin{1};
    K = 1;
else
    dim = varargin{1};
    K = varargin{2};
end

if ( dim ~= 1 && dim ~= 2 && dim ~= 3 )
    error('SPHEREFUN:DIFF:DIM', 'Unrecognized coordinate dimension');
end

if ( abs( K - round(K) ) > eps )
    error('SPHEREFUN:DIFF:DIFFORDER', 'Fractional derivatives not allowed')
end
K = round( K );

% We are going to work at the tech level to make things faster.
[C, D, R] = cdr( f );

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C)+mod(length(C),2);

% The variable coefficients in the definitions of the derivatives means
% that the length of the columns and rows will increase by one wave number
% after taking the derivatives with respect to x and y. The z derivative
% only increases the columns wave number by 1. The means we need to pad the
% coefficients with on extra zero negative and positive coefficient before
% doing the computations.
if dim ~= 3
    m = length(R)+2;  % Pad rows
else
    m = length(R);
end
n = n + 2;  % Pad columns

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Alias will do the padding of the coefficients.
Ccfs = ctechs.alias(ctechs.coeffs,n);
Rcfs = rtechs.alias(rtechs.coeffs,m);

% Implement higher derivatives as repeated (iterated) differentiation
for j=1:K
    [Ccfs,D,Rcfs] = onediff(Ccfs, D, Rcfs, dim);
end
f.
% f = spherefun(real(F([n/2+1:end 1],:)));

end

% Computes one derivative of f for the given dimension
function F = onediff(C, D, R, dim)
% TODO: This code will not work for complex valued spherefuns, if we ever
% allow them.
% realf = isreal(f);

% Simplify f to avoid any extra work
% f = simplify(f);

n = length(C);
m = length(R);

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1],[-1 1],n,n);
Msinm = .5i*spdiags(ones(m,1)*[-1,1],[-1 1],m,m);
Mcosn = .5*spdiags(ones(n,1)*[1,1],[-1 1],n,n);
Mcosm = .5*spdiags(ones(m,1)*[1,1],[-1 1],m,m);

% Differentiation matrix
DFn = trigspec.diffmat(n,1); 
DFm = trigspec.diffmat(m,1); 

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Alias will do the padding of the coefficients.
ctechs.coeffs = ctechs.alias(ctechs.coeffs,n);
rtechs.coeffs = rtechs.alias(rtechs.coeffs,m);

Ccfs = ctechs.coeffs;
Rcfs = rtechs.coeffs;

% Compute the derivatives
dCdth = DFn*Ccfs;
dRdlam = DFm*Rcfs;

% dx and dy involve two terms while dz is only one so we handle the cases
% separately
if ( dim == 1 ) || ( dim == 2)
    if ( dim == 1 )            % x
        % Calculate -sin(lam)./sin(th) dfdlam
        C1 = Msinn \ Ccfs;
        R1 = -Msinm*dRdlam;

        % Calculate cos(lam)cos(th) dfdth
        C2 = Mcosn*dCdth;
        R2 = Mcosm*Rcfs;
    elseif ( dim == 2 )         % y
        % Calculate the C * D * R.' decomposition of cos(lam)./sin(th) dfdlam
        Ccfs = ctechs.coeffs;
        C1 = Msinn \ Ccfs;
        R1 = Mcosm*dRdlam;

        % Calculate the C * D * R.' decomposition of sin(lam)cos(th) dfdth
        C2 = Mcosn*dCdth;
        R2 = Msinm*Rcfs;
    end

    X = C1*D*R1.' + C2*D*R2.';

    % Compression plus may not preserve the expansion properties we want.
    % So we sample each piece add them together and construct a spherefun.
    % TODO: Fix this so everything is done in coefficient space, like this
    % f = f1 + f2;        
    % When constructing from samples, m must be even.
%     m = m + mod(m,2);
%     f = spherefun(sample(f1,m,n/2+1)+sample(f2,m,n/2+1));    
else
    % Calculate sin(th) dfdth
    C1 = -Msinn*dCdth;
    R1 = Rcfs;
    X = C1*D*R1.';
end    

F = trigtech.coeffs2vals(trigtech.coeffs2vals(X).').';
%     F = trigtech.coeffs2vals(trigtech.coeffs2vals(SphereFourierFilter(X)).').';
%     f = spherefun(real(F([n/2+1:end 1],:)));
%     f = spherefun(F);    

end

function X = SphereFourierFilter( X ) 
% Fourier filter on the sphere: 

[m, n] = size( X ); 

% Go to VALUES-FOURIER space:
X = trigtech.coeffs2vals( X );

% Construct mask: 
mask = ( ones(m,1)*abs(-floor(n/2):floor(n/2)-1) < n/2*abs(sin(pi*trigpts(m)))*ones(1,n)+1 );
mask(1, [floor(n/2) floor(n/2)+2]) = 0;

% Apply filter: 
X = mask.*X; 

% Convert back to FOURIER-FOURIER: 
X = trigtech.vals2coeffs( X ); 

end