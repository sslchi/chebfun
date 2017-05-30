function [u, ZZ, DD, YY] = PoissonDisk_parity(f, n)
% uses ADI to find soln to poisson's eqn. This useful for finding
% a low rank factorization without calling constructor/forming an
% explicit 2D coeff matrix. However, 
% it can be less accurate due to some estimations that are req'd. 

    tol = 1e-15*f.vscale;  %need a better way to set tol
    
%%
% 1) work on even-periodic part of the problem: 

    %even-periodic (incl. pole)
    id = f.idxPlus;
    if isempty(id)
        up = diskfun(); 
        ZZ1 = []; 
        YY1 = []; 
        DD1 = []; 
    else
    fp = f;
    fp.cols = fp.cols(:, id);
    fp.rows = fp.rows(:, id);
    fp.pivotValues = fp.pivotValues(id, :);
    fp.pivotLocations = fp.pivotLocations(id, :); 
    fp.idxPlus = 1:length(id);
    fp.idxMinus = [];
    
    %call ADI: 
    [ZZ, DD, YY] = poisson_diskADI(fp, n);
    
    %build a diskfun 
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    up = fp;
    up.cols = c;
    up.rows = r; 
    up.pivotValues = 1./diag(DD);
    up.idxPlus = 1:length(up.pivotValues); 
    up.idxMinus = []; 
    up.pivotLocations = nan(length(up.pivotValues), 2); 
    up.nonZeroPoles = 0;                       
    
    %check for a nonzero pole 
    %the ADI algorithm results in information about the pole
    % being distributed across all rows and cols. 
    % We assume in all other diskfun routines that all
    % pole information is contained in first row/col only. We must correct
    % this structure. 
    
%     uv = sample(up); 
%    [vp, poleconst] = checkPole(uv(1,:), tol);
%     if poleconst==0
%         error('pole is not constant-valued')
%     end
    vp = feval(up,0,0); 
    if abs(vp) < tol;                      
       up.nonZeroPoles = 0; 
       %enforce that each rank-1 term is zero at pole
       up = projectOntoBMCII_zp(up);
    else
       up.nonZeroPoles = 1; 
                             
    %build a diskfun upzp = (up-vp) from coeffs. This will be zero at origin.  
    mm = length(up.rows); 
    nn = length(up.cols);
    uvcol = zeros(nn,1); 
    uvrow = zeros(mm, 1); 
    uvcol(1) =-vp ; 
    uvrow((mm-mod(mm,2))/2+1) = vp;
    c = real(chebfun(chebtech2({'',uvcol})));
    r = real(chebfun(uvrow,[-pi,pi],'coeffs','periodic'));
    
    %note: this method potentially overestimates the rank by 2. 
    % use compression plus instead
%     upzp = up; 
%     upzp.cols = [c up.cols]; 
%     upzp.rows = [r up.rows];
%     upzp.pivotValues = [vp; up.pivotValues];
    
% (This is generic compression plus. TODO: since this is 
% always a low_rank(X)+rank_1 , highly structured,
% plus can probably be done faster than this.)

 fScl = diag(1./up.pivotValues);
 gScl = 1./vp;
 cols = [c up.cols];
 rows = [r up.rows];

[Qcols, Rcols] = qr(cols);
[Qrows, Rrows] = qr(rows);

Zt = zeros(length(gScl), length(fScl));
Dt = [gScl, Zt ; Zt.', fScl];
[U, S, V] = svd(Rcols * Dt * Rrows.');
s = diag(S);

vup = vscale(up); 
vvp = abs(vp);

vscl = 2*max(vup, vvp); 
% Remove singular values that fall below eps*vscale: 
idx = find( s > eps * vscl, 1, 'last');

    U = U(:,1:idx);
    V = V(:,1:idx);
    s = s(1:idx);
    upzp = up;
    upzp.cols = Qcols * U;
    upzp.rows = Qrows * conj(V); 
    upzp.pivotValues = 1./s;                               
 
    %now  make each rank-1 term possess BMC-II symmetry
    upzp = projectOntoBMCII_zp(upzp);
    %now build even part of u with correct structure
    up = upzp; 
    up.cols = [-c upzp.cols]; 
    up.rows = [r up.rows];
    up.pivotValues = [vp; up.pivotValues];
    end  
    
                                                              
    ZZ1 = ZZ; 
    DD1 = DD; 
    YY1 = YY; 
    end
    
    %%
    % odd-antiperiodic part  
    id = f.idxMinus;
    if isempty(id)
       um = diskfun(); 
    else
    fm = f;
    fm.cols = fm.cols(:, id);
    fm.rows = fm.rows(:, id);
    fm.pivotValues = fm.pivotValues(id);
    fm.pivotLocations = fm.pivotLocations(id, :);
    fm.idxMinus = 1:length(id);
    fm.idxPlus = [];
    fm.nonZeroPoles = 0;
    
    
    [ZZ, DD, YY] = poisson_diskADI(fm, n);
    
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    um = fm;
    um.cols = c;
    um.rows = r; 
    um.pivotValues = 1./diag(DD);
    um.idxMinus = 1:length(um.pivotValues); 
    um.idxPlus = []; 
    um.pivotLocations = nan(length(up.pivotValues), 2); 
    um = projectOntoBMCII(um); 
    ZZ = [ZZ1 ZZ]; 
    DD = diag([diag(DD1);diag(DD)]); 
    YY = [YY1 YY]; 
    end

       
    %%
    %combine results and build diskfun for solution
    
    pivots = [up.pivotValues;um.pivotValues];
    cols = [up.cols um.cols ];
    rows = [up.rows um.rows];

    numPlus = length(up.pivotValues);
    idxPlus = 1:numPlus;
    numMinus = length(um.pivotValues);
    idxMinus = (numPlus+1):(numPlus+numMinus);
    
    u = up;
    u.cols = cols;
    u.rows = rows;
    u.pivotValues = pivots;
    u.pivotLocations = nan(length(u.pivotValues), 2); 
    u.idxPlus = idxPlus;
    u.idxMinus = idxMinus;
    u.nonZeroPoles = up.nonZeroPoles;
   
   
end



function f = projectOntoBMCII_zp(f)
% enforces that f, which is zero-valued at the origin, has BMC-II symmetry.
% this is needed because ADI enforces BMC symmetry but spreads information
% about the origin across all rank 1 terms. 
% we only need the correction the even part of the solution u. 
% this is a slightly specialized version of @diskfun\projectOntoBMCII. 


% Operate on the column coefficients first.
X = f.cols.funs{1}.onefun.coeffs;

% Get size: 
[m, n] = size(X); 

 % First we enforce that the every function in r is even: 
 X(2:2:end, :) = 0; 
 
% now we enforce that every function in r is zero at origin:
evenModes = 1:n;
% to enforce that f(r) = 0 for each column function, 
% we want C with smallest 2-norm such that  A*(X(1:2:end, 2:end)+C) = 0, 
% where A =[1 -1 1 -1..] . 
% The solution is 
% C = A'*((A*A')\(A*X)). 
% odd modes in r are zero, they won't contribute.
even = 1:2:m;
Xe = X(even, evenModes); 
factor = 1/length(even)*(sum(Xe(1:2:end, :),1)-sum(Xe(2:2:end, :),1));
C = ((-1*ones(length(even), 1)).^((2:length(even)+1)'))*factor; 
%now add C to X
X(even, evenModes) = Xe-C; 
     
ctechs = chebtech2({'',X}); 
f.cols.funs{1}.onefun = ctechs;

% Now operate on the rows. The coefficients for the rows of an even BMC-II
% function should only contain even wave numbers. The projection is to
% simply zero out the odd wave numbers.
X = f.rows.funs{1}.onefun.coeffs;
n = size(X, 1);
zeroMode = floor(n/2) + 1;
oddModes = [fliplr(zeroMode-1:-2:1) zeroMode+1:2:n];
X(oddModes, :) = 0;
rtechs = real(trigtech({'', X}));
f.rows.funs{1}.onefun = rtechs;

% Weird feval behavior in chebfun requires this
f.cols.pointValues = feval(ctechs, [-1; 1]);
f.rows.pointValues = feval(rtechs, [-1; 1]); 

end

function [pole, constVal] = checkPole(val, tol)
% Check that the values at the pole are constant.

% Take the mean of the values that are at the pole.
pole = mean(val);

% Compute their standard deviation
stddev = std(val);

% If the standard deviation does not exceed the 1e8*tolearnce then the pole
% is "constant".
% TODO: Get a better feel for the tolerance check.
if ( (stddev <= 1e8*tol) || (stddev < eps) )
    constVal = 1;
else
    constVal = 0;
end

end

%ADI function
function [ZZ, DD, YY] = poisson_diskADI(f, n)
%performs ADI to find the CDR factorization of a diskfun u that satisfies 
% poisson's equation with homogeneous b.c. and rhs f, with n*(2n+1) deg of freed.

%%
tol = 1e-14; %this should be user-selected: ||X||_2 < tol*sigma_1(X); 
m = n;
n = 2*n+1;

%get coeffs and set up RHS 
[C, D, R] = cdr(f); 
C = chebtech2.alias(C.coeffs, n); 
R = trigtech.alias(R.coeffs, m); 
[~, r] = size(C);  
ds = sign(diag(D));
d = sqrt(abs(diag(D)));
C = C*diag(d); 
R = (R*diag(d.*ds)).';

%set up operators
Im = speye(m, m);

% S converts T coeffs to U coeffs
dg = .5*ones(n - 2, 1);
S = spdiags([1 0 ; .5 0 ; dg -dg], [0 2], n, n);

%M multiplies by (1-r.^2) in T
M = spdiags ([-1*ones(n, 1)/4 ones(n,1)/2 -1*ones(n, 1)/4]...
    , [-2 0 2], n, n);
M(2,2) = 1/4; 
M(3,1) = -1/2; 


%Mr2 multiplies by r.^2 in T
Mr2 = spdiags ([ones(n, 1)/4 ones(n,1)/2 ones(n, 1)/4]...
    , [-2 0 2], n, n);
Mr2(2,2) = 3/4; 
Mr2(3,1) = 1/2; 

%Mr2c multiplies by r.^2 in C1
Mr2c = Mr2; 
Mr2c(1,1) = 1/4; 
Mr2c(2,2) = 1/2; 
Mr2c(3,1) = 1/4; 

%Mpc multiplies by a reqd poly in C1 (see paper)
Mpc = mypolyM(n); 

%D1 is first derivative in C1
D1 = (1:n)'-1; 
D1 = spdiags(D1, 1, n, n); 


%D2 is [(1-r.^2).*second derivative] in C1
idx = (1:n).';
idx = idx-1; 
D2 = -1*spdiags([(-.5*idx+.5*idx.^2), [0;0;-.5*(idx(1:end-2)+2)-.5*(idx(1:end-2)+2).^2]], [0;2], n, n); 
D2(1:2, :) = [0 0 3 zeros(1, n-3); 0 0 0 6 zeros(1, n-4)]; 

%left multiply operator
A = Mr2c*(D2)+Mpc*D1-4*S*Mr2;

%right multiply operator: B is D2 for trig
B = -1*((1:m/2).^2).';
B = [B(end:-1:1); 0; B(1:end-1)] ; 
B = spdiags(B, 0, m, m); 

SM = S*M; 

%adjust columns of RHS to account for basis shift: 
C =  M\(Mr2*C);

%----constants associated with shift parameters----------------
%approximate min and max eigenvalues of A: 
[Amin, Amax] = estimate_Aeigs(n); 
Int = [Amin, Amax,0, -full(B(1,1))];

[T, Tinv, gam, ~] = mobiusTdisk(Int); %gam and cross ratio


%we estimate small 1/kp for use in ellipk
if abs(1/gam) < 1e-7
    kp = small_gam(gam); 
else
kp = sqrt(1-sqrt(1-(1/gam)^2));
end
[~, K] = ellipk( kp);
kp = sqrt(1-(1/gam)^2);%need this for ellipj function in fADIrank1

%approx max number ADI iterations

%A is mildy nonnormal, so we need to add the effect of the condition
% number for V that diagonalizes A. it grows less than quadratically with
% n, so we esitmate the log(cond(V)) with a linear fit

KV = ceil(1.274376695304042*log(n) -0.990672481311612);

%get shift parameters: 

N = ceil(1/pi^2*(log(4/tol)+KV)*log(4*gam));  % approx number ADI iterations/degree of rat function
idx = 1:N;
[~, ~, dn] = ellipj((2*idx-1)*K/(2*N), kp);
p = gam*dn; %optimal shift parameters
q = -Tinv(-p);
p = Tinv(p);
%%
% Now do factored ADI
%--------Main loop: performs factored ADI ------------------------
YY =[];
ZZ =[];



ZZ(:,1:r) = (A-SM*p(1))\(SM*C);  
YY(1:r,:) = R/(B-q(1)*Im);
        for j = 1:N-1
            ZZ(:, j*r+1:(j+1)*r) = ZZ(:, (j-1)*r+1:j*r)+...
                (A-SM*p(j+1))\((p(j+1)+q(j))*SM*ZZ(:,(j-1)*r+1:j*r));
            YY(j*r+1:(j+1)*r,:) =  YY((j-1)*r+1:j*r,:) +...
                ((q(j+1)+p(j))*YY((j-1)*r+1:j*r,:)/(B-q(j+1)*Im)); 
        end
        
%return pivots 
d = diag(-(p+q)); 
DD = kron(d,speye(r)); 
YY = YY.';

%-----compress results ----------


%share diag
ds = sign(diag(DD));
d = sqrt(abs(diag(DD)));
ZZ = ZZ*diag(d); 
YY = YY*diag(d.*ds);


%compression step
[QZ, RZ] = qr( ZZ, 0); 
[QY, RY] = qr( YY, 0);

% XX * YY^T = (QX*RX) * (QY*RY)^T = QX*(RX*RY^T)*QY^T
inner_piece = RZ*RY'; 
[U, S, V ] = svd( inner_piece ); 
tol = eps; 
idx = find(diag(S)>tol*S(1,1), 1, 'last'); 
U = U(:,1:idx); 
V = V(:,1:idx); 
S = S(1:idx,1:idx); 


ZZ = QZ * U; 
DD = S; 
YY = QY * V; 
%ZZ is in Chebyshev-Dirichlet: 
%multiply by M to get back into Chebyshev basis.
ZZ = M*ZZ;



function M = mypolyM(m)
%multiplication by -5x^3+ x in C1 basis
    
a1 = [0;-2.75;0; -1.25];
a1 = [a1 ; zeros(m - 4, 1)];

M = spdiags([-1.25*ones(m,1) -2.75*ones(m, 1) -2.75*ones(m, 1) -1.25*ones(m,1)]/2, [-3 -1 1 3], m, m); 
sz = length(a1)-2; 
sub = 1:sz;
H = 0*speye(sz); 
H(2,1) = -1.25/2;
H(1,2) = -1.25/2;
    
M(sub, sub)= M(sub, sub) - H;
end

function g = small_gam(gamma) 
% Taylor series estimate 1-sqrt(1-(1/gam)^2))
% when 1/gam is small. 

x = 1/gamma;
cfs = 1./(sqrt(2)*[1  8 128 1024 ]);
x = [ x; x^3; x^5; x^7 ];
g = cfs*x; 
end

function [Aleft, Aright] = estimate_Aeigs(n)
%cheaply approximates the eigenvalues for A

%estimate leftmost bound on interval containing eigenvalues
grd = 3.992788438847006;
int = -1.135681960892715;
es = -1e-10;
lin = @(x) -exp(log(x)*grd + int-es); 
Aleft = lin(n); 

%estimate rightmost bound on interval containing eigenvalues

grd = -0.224921441849691;
int = -1.592726925430925;
es = 1e-2;
lin = @(x) -exp(log(x)*grd + int-es);
Aright = lin(n); 

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


end


