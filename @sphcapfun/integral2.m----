function varargout = integral2(varargin)
%INTEGRAL2  Double integral of a SPHCAPFUN over its domain.
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   SPHCAPFUN.
%
%
%   I = INTEGRAL2(F, [a b c d]) integrates F in polar coordinates
%   over the region [a b] x [c d], where a and b are angular values 
%   in the interval [-pi  pi], and  c and d are radial values in the 
%   interval [0 1].
%
% See also SPHCAPFUN/INTEGRAL, SPHCAPFUN/SUM2, SPHCAPFUN/QUAD2D.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = quad2d(varargin{:});

end
