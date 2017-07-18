function vscl = vscale(f)
%VSCALE   Vertical scale of a DISKFUN.
%   VSCL = VSCALE(F) returns the vertical scale of a DISKFUN as determined
%   by evaluating on a coarse tensor-product grid. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

vscl = vscale( f.diskFunction ); 

end