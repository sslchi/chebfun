function G = curl( f ) 
%CURL   Numerical curl of a scalr SPHCAPFUN times the normal to the disk.
%   G = CURL(F) returns a SPHCAPFUNV G representing the numerical
%   curl of the scalar SPHCAPFUN F times the normal to the disk, i.e.
%   curl([0 0 F]). This is equivalent to taking the gradient of F and 
%   crossing it with the normal component to the disk. Note that this is
%   different from the curl of a diskfunv (i.e., curl(V), where V is a
%   vector-valued function): see SPHCAPFUNV/curl.
%
% See also SPHCAPFUN/GRADIENT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if isempty( f )
    G = diskfunv;
    return;
end

% Curl of F is curl(0, 0,F) = (F_y, -F_x)
G = diskfunv(diff(f, 2), -diff(f, 1));

end

