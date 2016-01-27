function val = feval(f, x)
%FEVAL   Evaluate a SINGFUN.
%   FEVAL(F, X) evaluates the SINGFUN F at the given points X.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% For evaluation, the underlying SMOOTHFUN, i.e., f.smoothPart, is first
% evaluated at X and then the values computed are scaled by the singular
% factors.

% Evaluate each of the funs:
val = cellfun(@(f) feval(f, x), f.funs, 'UniformOutput', false);
val = sum(cell2mat(val), 2);

end
