function f = flipud(g)
%FLIPUD   Flip/reverse a SPHCAPFUN over the x-axis.
%
%   G = FLIPUD(F) returns a SPHCAPFUN G that is flipped over the 
%   x-axis, that is G(x,y) = F(x, -y).
%
% See also SPHCAPFUN/FLIPLR, SPHCAPFUN/FLIPDIM, SPHCAPFUN/ROTATE. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( g ) ) 
    f = diskfun;
    return
end 
f = g; 
f.rows = flipud(g.rows); 

end
