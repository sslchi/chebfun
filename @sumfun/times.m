function h = times(f, g)
%.*   Multiply SINGFUNS with SINGFUNS and SMOOTHFUNS
%   F.*G multiplies SINGFUN objects F and G or a SINGFUN by a scalar/SMOOTHFUN
%   if either F or G is a scalar/SMOOTHFUN.
%
% See also LDIVIDE, RDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Empty arguments:
if ( isempty(f) || isempty(g) )
    h = [];
    return
end

%%

if ( ~isa(g, 'sumfun') )
    h = f;
    h.funs = cellfun(@(f) times(f, g), f.funs, 'UniformOutput', false);
elseif ( ~isa(f, 'sumfun') )
    h = g;
    h.funs = cellfun(@(g) times(f, g), g.funs, 'UniformOutput', false);
else 
    h = f;
    funs = {};
    for k = 1:size(f.funs)
        funs = [funs cellfun(@(f) times(f, g), f{k}.funs, 'UniformOutput', false)];
    end
    h.funs = funs;
end

end
