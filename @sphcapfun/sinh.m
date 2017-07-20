function varargout = sinh(varargin)
%SINH   Hyperbolic sine of a SPHCAPFUN.
%
% See also SPHCAPFUN/SIN and SPHCAPFUN/COSH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sinh@separableApprox(varargin{:});

end
