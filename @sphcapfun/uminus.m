function g = uminus(f)
%UMINUS   Unary minus for a DISKFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = f; 
g.diskFunction = -f.diskFunction; 

end