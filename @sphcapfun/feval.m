function y = feval(varargin)
%FEVAL   Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, X, Y) evaluates a diskfun F at a point (X,Y) in Cartesian
%   cooridnates, where X and Y are doubles.
%
%   Y = FEVAL( F, THETA, R, 'polar') evaluates a diskfun F in polar
%   coordinates (THETA,R).  Here THETA and R are doubles representing the
%   central angle (in radians) and radius in polar coordinates and must be
%   points in the unit disk.
%
%   Y = FEVAL(F, c), where c is a complex-valued chebfun representing a
%   contour, evaluates F along the contour.
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

varargin{1} = varargin{1}.diskFunction;
y = feval( varargin{:} ); 

end
