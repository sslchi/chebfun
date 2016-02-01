function data = plotData(f, g)
%PLOTDATA   Useful data values for plotting a SINGFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the smooth part of F and then scales
%   it by the singular factors given in the EXPONENTS of F.
%
%   DATA = PLOTDATA(F, G) is similar but for plot calls of the form PLOT(F, G),
%   where both F and G are SINGFUN objects.
%
%   DATA = PLOTDATA(F, G, H) is for plots of the form PLOT3(F, G, H). In this
%   instance, DATA also contains fields zLine and zPoints for the data
%   corresponding to H.
%
% See also PLOT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    tmp = cellfun(@plotData, f.funs, 'UniformOutput', false);
    data = tmp{1};

    xLine = cell2mat(cellfun(@(f) f.xLine, tmp, 'UniformOutput', false).');
    data.xLine = unique(xLine);
    data.yLine = feval(f, data.xLine);

    xPoints = cell2mat(cellfun(@(f) f.xPoints, tmp, 'UniformOutput', false).');
    data.xPoints = unique(xPoints);
    data.yPoints = feval(f, data.xPoints);
    
elseif ( nargin == 2 )
        % PLOT(F, G)
       
    tmpF = cellfun(@plotData, f.funs, 'UniformOutput', false);
    data = tmpF{1};
    tmpG = cellfun(@plotData, g.funs, 'UniformOutput', false);

    xF = cell2mat(cellfun(@(f) f.xLine, tmpF, 'UniformOutput', false).');
    xF = unique(xF);
    xG = cell2mat(cellfun(@(f) f.xLine, tmpG, 'UniformOutput', false).');
    xG = unique(xG);
    x = unique([xF, xG]);
    
    data.xLine = feval(f, x);
    data.yLine = feval(g, x);

    xPoints = cell2mat(cellfun(@(f) f.xPoints, tmpF, 'UniformOutput', false).');
    data.xPoints = unique(xPoints);
    data.yPoints = feval(f, data.xPoints);
end

end
