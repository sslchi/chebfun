function varargout = cumprod(varargin)
%CUMPROD  Indefinite product integral of a SPHCAPFUN. 
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a SPHCAPFUN user.
%
% See also SPHCAPFUN/CUMSUM, SPHCAPFUN/SUM, SPHCAPFUN/PROD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHCAPFUN:CUMPROD:notSupported',...
        'Indefinite product integral of a SPHCAPFUN not supported.');

end
