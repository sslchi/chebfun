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
% if dim ~= 3
%     m = length(R)+2;  % Pad rows
% else
%     m = length(R);
% end
% n = n + 2;  % Pad columns
m = length(R);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% % Alias will do the padding of the coefficients.
% Ccfs = ctechs.alias(ctechs.coeffs,n);
% Rcfs = rtechs.alias(rtechs.coeffs,m);

Ccfs = ctechs.coeffs;
Rcfs = rtechs.coeffs;

idxPlus = f.idxPlus;
idxMinus = f.idxMinus;

% Implement higher derivatives as repeated (iterated) differentiation
for j=1:K
    % Alias will do the padding of the coefficients.
    n = n + 2;
    m = m + 2;
    Ccfs = ctechs.alias(Ccfs,n);
    Rcfs = rtechs.alias(Rcfs,m);
    [Ccfs,D,Rcfs,idxPlus,idxMinus] = onediff(Ccfs, D, Rcfs, idxPlus, idxMinus, dim);
end

ctechs = real(trigtech({'',Ccfs}));
f.cols.funs{1}.onefun = ctechs;
rtechs = real(trigtech({'',Rcfs}));
f.rows.funs{1}.onefun = rtechs;
f.idxPlus = idxPlus;
f.idxMinus = idxMinus;
f.pivotValues = 1./diag(D);

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs,[-1;1]);
f.rows.pointValues = feval(rtechs,[-1;1]); 

end

% Computes one derivative of f for the given dimension
function [C,D,R,idxPlus,idxMinus] = onediff(Ccfs, D, Rcfs, idxPlus, idxMinus, dim)
% TODO: This code will not work for complex valued spherefuns, if we ever
% allow them.
% realf = isreal(f);

% Simplify f to avoid any extra work
% f = simplify(f);

n = length(Ccfs);
m = length(Rcfs);

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1],[-1 1],n,n);
Msinm = .5i*spdiags(ones(m,1)*[-1,1],[-1 1],m,m);
Mcosn = .5*spdiags(ones(n,1)*[1,1],[-1 1],n,n);
Mcosm = .5*spdiags(ones(m,1)*[1,1],[-1 1],m,m);

% Differentiation matrix
DFn = trigspec.diffmat(n,1); 
DFm = trigspec.diffmat(m,1); 

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
    
    % Parity changes
    temp = idxPlus;
    idxPlus = idxMinus;
    idxMinus = temp;
    
    % idxPlus corresponds to the columns being even functions. We know that
    % the Fourier coefficients c_k, k=-floor(n/2):ceil(n/2)-1, of these 
    % columns then satisfy the property:
    %       c_{k} - c_{-k} = 0 and Im(c_k) = 0
    % Additionally, idxPlus corresponds to the rows only containing even
    % Fourier modes.  We will exploit these properties in compression plus
    % to make it faster and so that the results are guaranteed to also
    % satisfy these properties.
    
%     [Cp,dp,Rp] = compression_plus(C1(:,idxPlus), R1(:,idxPlus),...
%         C2(:,idxPlus), R2(:,idxPlus), D(idxPlus,idxPlus));

    if ~isempty( idxPlus )
        evenStart = 1 + mod(floor(m/2),2);
        
        B1c = real([C2(1,idxPlus);C1(2:n/2+1,idxPlus)+C1(n:-1:n/2+1,idxPlus)]);
        E1r = R1(evenStart:2:m,idxPlus);

        B2c = real([C2(1,idxPlus);C2(2:n/2+1,idxPlus)+C2(n:-1:n/2+1,idxPlus)]);
        E2r = R2(evenStart:2:m,idxPlus);

        [Bc,dp,Er] = compression_plus(B1c, E1r, B2c, E2r, D(idxPlus,idxPlus));

        Cp = real(0.5*[Bc(1,:);Bc(2:n/2,:);flip(Bc(2:n/2-1,:))]);
        Rp = zeros(m,size(Er,2));
        Rp(evenStart:2:m,:) = Er;
    else
        Cp = []; dp = []; Rp = [];
    end        
        
%     [Cm,dm,Rm] = compression_plus(C1(:,idxMinus), R1(:,idxMinus),...
%         C2(:,idxMinus), R2(:,idxMinus), D(idxMinus,idxMinus));
    if ~isempty( idxMinus )
        oddStart = 2 - mod(floor(m/2),2);
        
        B1c = imag(C1(2:n/2+1,idxMinus)-C1(n:-1:n/2+1,idxMinus));
        E1r = R1(oddStart:2:m,idxMinus);        

        B2c = imag(C2(2:n/2+1,idxMinus)-C2(n:-1:n/2+1,idxMinus));
        E2r = R2(oddStart:2:m,idxMinus);

        [Bc,dm,Er] = compression_plus(B1c, E1r, B2c, E2r, D(idxMinus,idxMinus));

        Cm = (0.5*1i)*real([zeros(1,size(Er,2));Bc(1:n/2,:);-flip(Bc(1:n/2-1,:))]);
        Rm = zeros(m,size(Er,2));
        Rm(oddStart:2:m,:) = Er;
    else
        Cm = []; dm = []; Rm = [];
    end

    C = [Cp Cm];
    R = [Rp Rm];
    
    idxPlus = 1:size(Cp,2);
    idxMinus = size(Cp,2)+1:size(Cp,2)+size(Cm,2);
    
    D = diag([dp;dm]);
    
    % Compression plus may not preserve the expansion properties we want.
    % So we sample each piece add them together and construct a spherefun.
    % TODO: Fix this so everything is done in coefficient space, like this
    % f = f1 + f2;        
    % When constructing from samples, m must be even.
%     m = m + mod(m,2);
%     f = spherefun(sample(f1,m,n/2+1)+sample(f2,m,n/2+1));    
else
    % Calculate sin(th) dfdth
    C = -Msinn*dCdth;
    R = Rcfs;
end    

%  F = trigtech.coeffs2vals(trigtech.coeffs2vals(X).').';
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

function [U,s,V] = compression_plus(Cf,Rf,Cg,Rg,D)

% Check for special cases
if isempty(Cf) && isempty(Cg)
    U = []; s = []; V = [];
    return
elseif isempty(Cf)
    U = Cg; s = diag(D); V = Rg;
    return
elseif isempty(Cg)
    U = Cf; s = diag(D); V = Rf;
    return
end

scl = diag(1./diag(D));
cols = [Cf Cg];
rows = [Rf Rg];

[Qcols, Rcols] = qr(cols);
[Qrows, Rrows] = qr(rows);

Z = zeros(length(scl));
D = [ scl, Z ; Z.', scl ];
[U, S, V] = svd(Rcols * D * Rrows.');
% If V is complex-valued, then conjugate: 
V = conj( V ); 
% Take diagonal from SIGMA:
s = diag(S);

% Compress the format if possible.
% [TODO]: What should EPS be in the tolerance check below?

% Remove singular values that fall below eps*vscale: 
idx = find( s > 10*eps, 1, 'last');

U = Qcols*U(:,1:idx);
V = Qrows*V(:,1:idx);
s = s(1:idx);

end

function f = projectOntoBMCI( f )
% PROJECTONTOBMCI  Projection onto BMC-I symmetry.
%
% g = projectOntoBMCI(f) is the orthogonal projection of f onto BMC-I 
% symmetry, i.e., a function that is
% 1. even in theta for every even wave number in lambda;
% 2. odd in theta for every odd wave number in lambda;
% Additionally, for all but k=0 wavenumber lambda the resulting projection
% enforces the spherefun is zero at the poles. 
%
% The projection is orthogonal, i.e., the correction matrix to fix up the
% structure has the smallest possible Frobenius norm.

% Even part
feven = f;
feven.cols = feven.cols(:,feven.idxPlus);
feven.rows = feven.rows(:,feven.idxPlus);
feven = projectOntoEvenBMCI( feven );

% Odd part
fodd = f;
fodd.cols = fodd.cols(:,fodd.idxMinus);
fodd.rows = fodd.rows(:,fodd.idxMinus);
fodd = projectOntoOddBMCI( fodd );

% Put pieces back together.
f.cols(:,f.idxPlus) = feven.cols;
f.rows(:,f.idxPlus) = feven.rows;

f.cols(:,f.idxMinus) = fodd.cols;
f.rows(:,f.idxMinus) = fodd.rows;

end

function f = projectOntoEvenBMCI( f )
% Project a spherefun to have even BMC-I symmetry, i.e., a spherefun that
% is pi-periodic in lambda and even in theta. The projection is orthogonal,
% i.e., the correction matrix to fix up the structure has the smallest
% possible Frobenius norm.

% Nothing to project
if isempty( f )
    return;
end

% Operate on the column coefficients first to project them onto even
% functions.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
[m, n] = size(X); 

isevenM = false;
if mod(m,2) == 0
    X(1,:) = 0.5*X(1,:);
    X = [X;X(1,:)];
    m = m+1;
    isevenM = true;
end

% Goal is to enforce the expansion is even in theta

evenModes = 1:n;
if f.nonZeroPoles
    zeroMode = 1;
    % Need to handle the zero mode in lambda separately as the expansion
    % only has to be forced to be even in theta.

    % Let 
    % I = eye(m); A = I - fliplr(I); A = A(1:(m-1)/2,:);
    % then we want to compute to find C with smallest two-norm such that
    % A*(X + C) = 0
    % The solution is 
    % C = A'*((A*A')\(A*X))
    % However, because of the form of the matrix equation the solution can
    % be worked out to just be
    C = 0.5*(X(1:m,zeroMode)-X(m:-1:1,zeroMode));

    % Update coeff matrix: 
    X(:,zeroMode) = X(:,zeroMode) - C; 

    % The result of the code now needs to operate on the remaining even,
    % non-zero modes.
    evenModes = 2:n;
end

if ~isempty( evenModes )
    % Second do the even, non-zero modes in lambda

    % First enforce these are even as above.
    C = 0.5*(X(1:m,evenModes)-X(m:-1:1,evenModes));
    X(:,evenModes) = X(:,evenModes) - C;

    % Now enforce these are zero at the poles (evenness is preserved)
    % Letting 
    % A = [[ones(1,m); (-1).^waveNumbers];
    % We want to find the C with smallest two-norm such that 
    % A*(X + C) = 0
    % The solution is 
    % C = A'*((A*A')\(A*X))
    % However, we can again work out the solution in close form because of
    % the special structure of the matrix equations.
    X(1:2:m,evenModes) = bsxfun(@minus,...
        X(1:2:m,evenModes),(2/(m+1))*sum(X(1:2:m,evenModes),1));
    X(2:2:m,evenModes) = bsxfun(@minus,...
        X(2:2:m,evenModes),(2/(m-1))*sum(X(2:2:m,evenModes),1));    
end

% If m is even we need to remove the mode that was appended 
if ( isevenM )
    X(1,:) = (X(1,:)+X(end,:));
    X = X(1:m-1,:);
end

ctechs = real(trigtech({'',X}));
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an even BMCI
% function should only contain even wave numbers. The projection is to
% simply zero out the odd wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X,1); 
zeroMode = floor(n/2)+1;
oddModes = [fliplr(zeroMode-1:-2:1) zeroMode+1:2:n];
X(oddModes,:) = 0;
rtechs = real(trigtech({'',X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs,[-1;1]);
f.rows.pointValues = feval(rtechs,[-1;1]); 

end

function [C, R] = projectOntoOddBMCI( C, R )
% Project a spherefun to have odd BMC-I symmetry, i.e., a spherefun that is
% pi-anti-periodic in lambda and even in theta. The projection is
% orthogonal, i.e., the correction matrix to fix up the structure has the
% smallest possible Frobenius norm.

% Nothing to project
if isempty( C ) || isempty( R )
    return;
end

% Operate on the column coefficients first to project them onto odd
% functions.
X = C;

% Get size: 
m = size(X,1);

isevenM = false;
if mod(m,2) == 0
    X(1,:) = 0.5*X(1,:);
    X = [X;X(1,:)];
    m = m+1;
    isevenM = true;
end

% I = eye(m); A = I + fliplr(I); 
% A = A(1:(m-1)/2+1,:); A((m-1)/2+1,(m-1)/2+1) = 1;
% 
% % Solution to underdetermined system A*(X + Y) = 0 with smallest Frobenius
% % norm: 
% C = A\(A*X);
% % Update coeff matrix: 
% X = X - C; 

% Update coeff matrix: 
X = X - 0.5*(X(1:m,:)+X(m:-1:1,:)); 

% If m is even we need to remove the mode that was appended 
if ( isevenM )
    X(1,:) = (X(1,:)+X(end,:));
    X = X(1:m-1,:);
end
C = X;

% Now operate on the rows. The coefficients for the rows of an odd BMCI
% function should only contain odd wave numbers. The projection is to
% simply zero out the even wave numbers.
X = R;
n = size(X,1); 
zeroMode = floor(n/2)+1;
evenModes = [fliplr(zeroMode-2:-2:1) zeroMode:2:n];
X(evenModes,:) = 0;

R = X;

end