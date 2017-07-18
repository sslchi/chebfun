function h = times(f, g)
%.*   Pointwise multiplication for DISKFUN objects.
%   F.*G multiplies DISKFUN objects F and G. Alternatively F or G could be 
%   a double.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

domf = f.domain; 
domg = g.domain; 
if ( norm( domf-domg ) > 10*eps )
    error('Domains of sphcapfun objects are inconsistent.')
end

if ( isa(f, 'sphcapfun') && isa(g, 'sphcapfun') ) 
    h = f; 
    h.diskFunction = times( f.diskFunction, g.diskFunction ); 
elseif ( isa(f, 'sphcapfun') )
    h = f; 
    h.diskFunction = times( f.diskFunction, g ); 
elseif ( isa(g, 'sphcapfun') )
    h = g; 
    h.diskFunction = times( f, g.diskFunction ); 
end


end