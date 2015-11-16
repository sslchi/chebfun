%% chebpde.m -- an executable m-file for solving a partial differential equation..
% Automatically created in CHEBGUI by user nhale.
% Created on November 13, 2015 at 14:29.

%% Problem description.
% Solving
%   u_t = 0.1*u" + u',
% for x in [-3,3] and t in [0,6], subject to
%   neumann at x = -3
% and
%   dirichlet at x = 3

%% Problem set-up
% Create an interval of the space domain...
dom = [-3 3];
%...and specify a sampling of the time domain:
t = 0:.1:6;

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) 0.1.*diff(u,2)+diff(u) - u.^2;

% Assign boundary conditions.
bc.left = 'dirichlet';
bc.right = 'dirichlet';

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial condition.
u0 = sin(pi.*x).*(x./6+1./2).^2;

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'Ylim', [-1,1], 'plot', 'on');

%% Call pde15s to solve the problem.
[t, u] = pde15s(pdefun, t, u0, bc, opts);