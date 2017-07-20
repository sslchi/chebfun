function varargout = power(varargin)
%.^	     Pointwise power of a SPHCAPFUN. 
%
% F.^G returns a SPHCAPFUN F to the scalar power G, a scalar F to the
% SPHCAPFUN power G, or a SPHCAPFUN F to the SPHCAPFUN power G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = power@separableApprox(varargin{:});

end
