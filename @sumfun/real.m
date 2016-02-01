function f = real(f)
%REAL   Real part of a SINGFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL, CONJ.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty arguments:
if ( isempty(f) )
    return;
end

% Compute the real part of the smooth part of F.
f.funs = cellfun(@real, f.funs, 'UniformOutput', false);

end
