function varargout = surf(f, varargin)
%SURF   Surface plot of a SPHCAPFUN.
%   SURF(F) plots the SPHCAPFUN object F on the surface of the unit disk.
%
%   SURF(X, Y, F, ...) calls separableApprox/SURF.  See this function for
%   details.
%
%   SURF(..., 'PropertyName', PropertyValue,...) sets the value of the specified
%   surface property. Multiple property values can be set with a single
%   statement.
%
%   H = SURF(...) returns a handle to the figure.
%
% See also SPHCAPFUN/PLOT and SEPARABLEAPPROX/SURF.

% Empty check:
if ( isempty(f) )
    h = surf([]);
    if ( nargout == 1 )
        varargout = { h };
    end
    return
end

% How dense to make the samples?
minPlotNum = 200;
defaultOpts = {'facecolor', 'interp', 'edgealpha', .5, 'edgecolor', 'none'};

% Number of points to plot
j = 1; 
argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'numpts')
        minPlotNum = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isempty(argin) )
    argin = {};
end

if ( isa(f,'sphcapfun') )
    
    if ( ((nargin == 1) || ...
            ( (nargin > 1) && ~isempty(argin) && ...
            ~isa(argin{1}, 'separableApprox')) || ...
            ((nargin == 3) && isempty(argin))) || ...
            ((nargin == 2) && plotEarth == true) )
        % surf(f,...)
        
        dom = f.domain;
        l = linspace(dom(1), dom(2), minPlotNum);
        t = linspace(dom(3), dom(4), minPlotNum);
%         t = linspace(0, 1, minPlotNum );
        C = fevalm(f, l, t); 
%         t = linspace(dom(3), dom(4), minPlotNum);
        [ll, tt] = meshgrid(l, t);
        
        % Adjust tt if the domain is colatitude
        if (1) % iscolat( f )
            tt = pi/2 - tt;
%             tt = .5 - tt;
        end
        [xx,yy,zz] = sph2cart(ll,tt,ones(size(ll)));
        
        % Make some corrections to C for prettier plotting.
        if ( norm(C - C(1,1),inf) < 1e-10 )
            % If vals are very close up to round off then the color scale is
            % hugely distorted. This fixes that.
            [n, m] = size(C);
            C = C(1,1)*ones(n, m);
        end
        h = surf(xx, yy, zz, C, defaultOpts{:}, argin{:});
        xlim([-1 1])
        ylim([-1 1]) 
        zlim([-1 1])
        xlabel('x'), ylabel('y'), zlabel('z')
        % Make the aspect ratio equal if the plot is not currently being
        % held.
        if ( ~ishold )
            daspect([1 1 1]);
        end        
    else
        % Pass this along to the surf function in separableApprox.
        h = surf@separableApprox(f, varargin{:});
    end
end

if ( nargout > 0 )
    varargout = { h };
end

end
