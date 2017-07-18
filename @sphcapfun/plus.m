function h = plus(f, g)
%+   Plus for DISKFUN objects.
% F + G adds F and G. F and G can be scalars or DISKFUN objects.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

domf = f.domain; 
domg = g.domain; 

if ( norm( domf-domg ) > 10*eps )
    error('Domains of sphcapfun objects are inconsistent.')
end

hf = plus( f.diskFunction, g.diskFunction ); 

h = sphcapfun;
h.domain = domf; 
h.diskFunction = hf; 

end
