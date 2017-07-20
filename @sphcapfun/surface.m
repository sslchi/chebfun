function varargout = surface(varargin)
%SURFACE  Plot surface of a SPHCAPFUN.
%   SURFACE(F) plots the SPHCAPFUN object F on the surface of the unit disk.
%
%   SURFACE(X, Y, F, ...) calls separableApprox/SURF.  See this function for
%   details.
%
%   SURFACE(..., 'PropertyName', PropertyValue,...) sets the value of the specified
%   surface property. Multiple property values can be set with a single
%   statement.
%
%   H = SURFACE(...) returns a handle to the figure.
%
% See also SPHCAPFUN/PLOT and SPHCAPFUN/SURF.


% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surface@separableApprox(varargin{:});

end
