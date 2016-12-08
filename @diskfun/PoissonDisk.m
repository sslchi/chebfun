function [u, ZZ, DD, YY] = poisson_disk_ADI(f, n)
%performs ADI to find the CDR factorization of a diskfun u that satisfies 
% poisson's equation with homogeneous b.c. and rhs f, with n deg of freed.
%
% Input: f is a diskfun, DOF for BMC function = (2*n+1)*n; choose n even.
%
%TO DO: 
%
% 1. separate and solve even-pi-periodic and odd-pi-antiperiodic pieces
% individually. Need to deal with pole correctly.
%
% 2. create a CUR2diskfun command; to do it right requires the above
%
%
% 3. add in code for non-homogeneous bcs. 

%%
tol = 1e-15;
m = n;
n = 2*n+1;

%get coeffs and set up RHS 
[C, D, R] = cdr(f); 
C = chebtech2.alias(C.coeffs, n); 
R = trigtech.alias(R.coeffs, m); 
[~, r] = size(C); 
%dv = diag(D); 
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

[~, Tinv, gam, ~] = mobiusTdisk(Int); %gam and cross ratio

%mu = exp(pi^2/2/log(16*cr)); %used to bound # rank 1 ADI steps in
%factor-independent ADI; not used in this code. 
 
%we estimate small 1/kp for use in ellipk
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
for k = 1:r; 
        [ Z,D, Y] = rank1_fADI(C(:,k), R(k,:), A, B, Im, gam,...
            Tinv, K,kp, SM, Nmax); 
        ZZ = [ZZ Z]; 
        YY = [YY Y];        
        DD = [DD; D];
end


%-----compress results and build a diskfun----------

%share diag
ds = sign(DD);
d = sqrt(abs(DD));
ZZ = ZZ*diag(d); 
YY = conj(YY);
YY = (YY*diag(d.*ds));


%compression step
[QZ, RZ] = qr( ZZ, 0); 
[QY, RY] = qr( YY, 0);

% % XX * YY^T = (QX*RX) * (QY*RY)^T = QX*(RX*RY^T)*QY^T
inner_piece = RZ*RY'; 
[U, S, V ] = svd( inner_piece ); 
tol = 1e-15; 
idx = find(diag(S)>tol*S(1,1), 1, 'last');
U = U(:,1:idx); 
V = V(:,1:idx); 
S = S(1:idx,1:idx); 


ZZ = QZ * U; 
DD = S; 
YY = QY * V; 
% ZZ is in Chebyshev-Dirichlet: 
% multiply by M to get back into Chebyshev basis.
ZZ = M*ZZ;

%build a diskfun (currently this is not quite correct)
c = real(chebfun(chebtech2({'',ZZ})));
r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
u = f;
u.cols = c;
u.rows = r; 
u.pivotValues = 1./diag(DD);
u.idxPlus = 1:length(u.pivotValues); 
u.idxMinus = []; 
%u.pivotLocations = not sure what goes here
u.nonZeroPoles =1; 

%----------------------------------------

%factored ADI on a rank 1 piece
function [Z,D, Y] = rank1_fADI(C, R,A, B, Im, gam,Tinv, K,kp,  SM, N)

%compute ADI shifts
idx1 = 1:N;
[~, ~, dn] = ellipj((2*idx1-1)*K/(2*N), kp);
p1 = gam*dn; %optimal shift parameters
q =  -Tinv(-p1) ; 
p = Tinv(p1) ; 

%perform ADI
Z(:,1) = (A-SM*p(1))\(SM*C);  
Y(1,:) = R/(B-q(1)*Im);
        for j = 1:N-1
            Z(:, j+1) = Z(:, j)+(A-SM*p(j+1))\((p(j+1)+q(j))*SM*Z(:,j));
            Y(j+1,:) =  Y(j,:) + ((q(j+1)+p(j))*Y(j,:)/(B-q(j+1)*Im)); 
        end
        
%return pivots and avoid column-stacking Y
D = -(p+q).';
Y = Y.';
end

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
%NOTE: this needs improvement.

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






