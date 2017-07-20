function varargout = lu(varargin)
%LU   LU factorization of a SPHCAPFUN.
%
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a SPHCAPFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHCAPFUN:LU:notSupported',...
        'LU factorization of a SPHCAPFUN is not supported.');

end
