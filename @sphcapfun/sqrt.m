function varargout = sqrt(varargin)
%SQRT   Square root of a SPHCAPFUN object.
%   SQRT(F) returns the square root of a positive SPHCAPFUN F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sqrt@separableApprox(varargin{:});

end
