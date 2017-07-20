function varargout = exp(varargin)
%EXP  Exponential of a SPHCAPFUN
%   EXP(F) returns the exponential of a SPHCAPFUN. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = exp@separableApprox(varargin{:});

end
