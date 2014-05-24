% Test file for fourtech/cumsum.m
function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = fourtech();

%%
% Spot-check antiderivatives for a couple of functions.  We verify that the
% fourtech antiderivatives match the true ones up to a constant by checking 
% that the standard deviation of the difference between the two on a large 
% random grid is small. We also check that feval(cumsum(f), -1) == 0 each 
% time.

% Note that fourtech's cumsum only works when mean of the function is zero.
% Thus in all these tests, the functions have this property.  There is one
% test at the end that verifies an error message is given when cumsum is
% applied to a function without zero mean.
 
k = 2; a = 0;
f = testclass.make(@(x) sin(k*pi*(x-a)).*cos((k+1)*pi*(x-a)), [], [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k+1)*cos(pi*(x-a))-cos((2*k+1)*pi*(x-a)))/(2*(pi+2*k*pi));
err = feval(F, x) - F_ex(x);
tol = 10*F.vscale.*F.epslevel;
pass(1) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

k = 200; a = 0.17;
f = testclass.make(@(x) sin(k*pi*(x-a)).*cos((k+1)*pi*(x-a)), [], [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k+1)*cos(pi*(x-a))-cos((2*k+1)*pi*(x-a)))/(2*(pi+2*k*pi));
err = feval(F, x) - F_ex(x);
tol = 10*F.vscale.*F.epslevel;
pass(2) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

k1 = 5; a1 = -0.33; k2 = 40; a2 = 0.17;
f = testclass.make(@(x) sin(k1*pi*(x-a1)).*cos((k1+1)*pi*(x-a1)) + 1i*sin(k2*pi*(x-a2)).*cos((k2+1)*pi*(x-a2)), [], [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k1+1)*cos(pi*(x-a1))-cos((2*k1+1)*pi*(x-a1)))/(2*(pi+2*k1*pi)) + 1i*((2*k2+1)*cos(pi*(x-a2))-cos((2*k2+1)*pi*(x-a2)))/(2*(pi+2*k2*pi));
err = feval(F, x) - F_ex(x);
tol = 10*F.vscale.*F.epslevel;
pass(3) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = testclass.make(@(x) sin(4*x).^2, [], [], pref);
F = testclass.make(@(x) 0.5*x - 0.0625*sin(8*x), [], [], pref);
G = cumsum(f);
err = G - F;
tol = 10*G.vscale.*G.epslevel;
values = err.coeffs2vals(err.coeffs); 
pass(5) = (std(values) < tol) && (abs(feval(G, -1)) < tol);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a 
% constant.

f = testclass.make(@(x) x.*(x - 1).*sin(x) + 1, [], [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
tol = 10*g.vscale.*g.epslevel;
pass(6) = (norm(err, inf) < 100*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
tol = 10*h.vscale.*h.epslevel;
pass(7) = (std(err) < tol)  && (abs(feval(h, -1)) < tol);

%%
% Check operation for array-valued chebtech objects.

f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], [], pref);
F_exact = testclass.make(@(x) [(-cos(x)) (x.^3/3) (exp(1i*x)/1i)], [], [], pref);
F = cumsum(f);
err = std(feval(F, x) - feval(F_exact, x));
tol = 10*max(F.vscale.*F.epslevel);
pass(8) = (norm(err, inf) < tol)  && all(abs(feval(F, -1)) < tol);
  

end
