function varargout = transpose(varargin)
%.'   SPHCAPFUN transpose.
%   F.' is the non-conjugate transpose of the underlying quasimatrix for 
%   the SPHCAPFUN object F. Given F in polar coordinates, F(theta, r), 
%   F.' = F(r, theta). To transpose a SPHCAPFUN with respect to Cartesian
%   coordinates, use the command FLIPXY.
% 
% See also SPHCAPFUN/CTRANSPOSE and SPHCAPFUN/FLIPXY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = transpose@separableApprox(varargin{:});

end
