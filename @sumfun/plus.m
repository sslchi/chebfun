function f = plus(f, g)
%+   Addition of SINGFUN objects with SINGFUNs and SMOOTHFUNs.
%   F + G adds F and G, where F and G may be SINGFUN objects or scalars.
%
% See also MINUS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    f = [];
    return
end

if ( ~isa(f, 'sumfun') )
    f = sumfun(g.funs{1}.make(f));
elseif ( ~isa(g, 'sumfun') )
    g = sumfun(f.funs{1}.make(g));
end

f.funs = [f.funs, g.funs];
f = merge(f);

end
