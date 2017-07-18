function g = coeffs2sphcapfun(X, dom)
%COEFFS2DISKFUN   Convert a matrix of Chebyshev-Fourier coefficients to a 
%                 diskfun. 
% 
%   F = coeffs2diskfun(X) returns a diskfun object F that has a
%   Chebyshev-Fourier matrix of coefficients X.  This is useful for
%   computing quantities on the disk with the function F.
% 
% See also DISKFUN/COEFFS2

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 ) 
    dom = [-pi pi 0 2*pi/3];
end

f = diskfun.coeffs2diskfun( X ); 
g = sphcapfun;
g.domain = dom;
g.diskFunction = f; 

end