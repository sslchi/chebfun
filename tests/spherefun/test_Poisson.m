function pass = test_Poisson( ) 
% Check correctness of Poisson solver on the sphere: 

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

for k = [1 2 3 4]
    
    if k == 1
        m = 60; n = 40;
    elseif k == 2
        m = 61; n = 41;
    elseif k == 3
        m = 62; n = 42;
    else
        m = 63; n = 43;
    end

    nxt = 1;  
    for L = 1:3
        for M=0:L 

        f = spherefun.sphharm(L, M); 

        u = spherefun.poisson(-L*(L+1)*f, 0, m, n);

        pass(nxt,k) = ( norm( u - f, 2 ) < tol );

        nxt = nxt + 1; 
        end
    end

end

% Flatten pass array to so we can add a few more tests:
pass = pass(:);
nxt = numel(pass) + 1;

% Discretization sizes:
m = 40; 
n = 40;
 
% Example 1: 
f = spherefun(@(lam,th) -6*cos(lam).*cos(th).*sin(th)); 
exact = spherefun(@(lam,th) sin(th).*cos(th).*cos(lam));
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol ); 
nxt = nxt + 1;

% Example 2: 
f = spherefun(@(lam,th) -4*(3*cos(th)+5*cos(3*th)).*sin(lam).*sin(th)); 
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2);
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol ); 
nxt = nxt + 1;

% Inline example: 
f = spherefun(@(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th));
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
            sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th));
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol );
nxt = nxt + 1;

% Check that the mean is zero.
pass(nxt) = ( abs(mean2(u)) < tol );
nxt = nxt + 1;

% Combination of odd/even spherical harmonics
f = -5*6*(spherefun.sphharm(5,0) + spherefun.sphharm(5,3) + spherefun.sphharm(5,-4)) +...
    -6*7*(spherefun.sphharm(6,0) + spherefun.sphharm(6,-5) + spherefun.sphharm(6,4)) + ...
    -7*8*(spherefun.sphharm(7,2) + spherefun.sphharm(7,-4));
exact = (spherefun.sphharm(5,0) + spherefun.sphharm(5,3) + spherefun.sphharm(5,-4)) +...
        (spherefun.sphharm(6,0) + spherefun.sphharm(6,-5) + spherefun.sphharm(6,4)) + ...
        (spherefun.sphharm(7,2) + spherefun.sphharm(7,-4));
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol );
nxt = nxt + 1;

% Test m and n odd
m = 21; n = 21;
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol );
nxt = nxt + 1;

% Test m even and n odd
m = 22; n = 23;
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol );
nxt = nxt + 1;

% Test m odd and n even
m = 19; n = 22;
u = spherefun.poisson(f, 0, m, n);
pass(nxt) = ( norm(u - exact, inf) < tol );
nxt = nxt + 1;

% Check that the code properly deals with a right hand side without a 
% mean of zero.
f = spherefun(@(x,y,z) 1 + x);
warning('off','CHEBFUN:SPHEREFUN:POISSON:meanRHS');
u = spherefun.poisson(f, 0, 10);
warning('on','CHEBFUN:SPHEREFUN:POISSON:meanRHS');
% This f - mean2(f)
g = spherefun(@(x,y,z) x);
% Solution with the mean of f set to zero
v = spherefun.poisson(g, 0, 10);
pass(nxt) = ( norm(u - v) < tol );
nxt = nxt + 1;

% Check that the code allows the mean of the solution to be set to
% something other than zero.
f = spherefun(@(x,y,z) x.*y.*z );
u = spherefun.poisson(f, 1, 10);
pass(nxt) = ( abs(mean2(u)-1) < tol );


end