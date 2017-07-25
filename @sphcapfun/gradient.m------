function G = gradient( f ) 
%GRADIENT   Numerical gradient of a SPHCAPFUN. 
%   G = GRADIENT(F) returns the numerical gradient of the
%   SPHCAPFUN F as a SPHCAPFUNV G.
%
% See also SPHCAPFUNV/DIV, SPHCAPFUNV/CURL, SPHCAPFUN/CURL, SPHCAPFUNV/VORTICITY

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if isempty( f )
    G = diskfunv;
    return
end

fx = diff(f, 1);   % diff in x-variable
fy = diff(f, 2);   % diff in y-variable 

G = diskfunv(fx,fy);

end

