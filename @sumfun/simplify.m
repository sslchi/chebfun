function f = simplify(f, tol)
%SIMPLIFY  Simplifies the smoothPart of a SINGFUN object.
%  F = SIMPLIFY(F) returns a SINGFUN which has a simplified smoothPart by
%  calling SIMPLIFY in SMOOTHFUN and has exponents less than 1 by absorbing 
%  boundary roots with integer multiplicity into the smoothPart.
%
%  F = SIMPLIFY(F, TOL) does the same thing but uses TOL supplied by the user
%  and pass it to SIMPLIFY in SMOOTHFUN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.funs = cellfun(@(f) simplify(f, tol), f.funs, 'UniformOutput', false);

end
