function F = vertcat(varargin)
%VERTCAT   Vertical concatenation of SPHCAPFUN objects.
%   K = VERTCAT(F, G, H) is the vertical concatenation of SPHCAPFUN objects F, 
%   G, and H. The output K is a SPHCAPFUNV object.
% 
%   [F ; G] is equivalent to VERTCAT(F, G).
%
%   VERTCAT(F) returns the SPHCAPFUN F. 
% 
% See also SPHCAPFUNV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1};
    
elseif ( nargin == 2 )
    if ( all(cellfun(@(F) isa(F,'diskfun'), varargin)) )
        F = diskfunv(varargin{:});
        
    else
        error('SPHCAPFUN:vertcat:tooManyComponents', ...
            'Only SPHCAPFUN objects are valid to concatenate.');
    end
    
else
    error('SPHCAPFUN:vertcat:tooManyInputs', ...
        'Can only vertically concatenate two SPHCAPFUN objects.');
end
    
end
