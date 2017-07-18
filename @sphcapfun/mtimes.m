function h = mtimes(f, g)
%*   Pointwise multiplication for DISKFUN objects.
%   c*F or F*c multiplies a DISKFUN F by a scalar c.
%
%   F*G computes the integral of F(s,t)G(l,s) over s, and this is the 
%   continuous analogue of matrix-matrix multiplication.
%
% See also DISKFUN/TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'sphcapfun') && isa(g, 'sphcapfun') ) 
    h = f; 
    h.diskFunction = mtimes( f.diskFunction, g.diskFunction ); 
elseif ( isa(f, 'sphcapfun') )
    h = f; 
    h.diskFunction = mtimes( f.diskFunction, g ); 
elseif ( isa(g, 'sphcapfun') )
    h = g; 
    h.diskFunction = mtimes( f, g.diskFunction ); 
end

end