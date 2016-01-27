function f = merge(f)

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = numel(f.funs);

exps = zeros(n, 2);
for k = 1:n
    if ( isa(f.funs{k}, 'singfun') ) 
        exps(k,:) = get(f.funs{k}, 'exponents');
    end
end

% Find compatable exponents
d1 = bsxfun(@minus, exps(:,1), exps(:,1)');
mask1 = triu(d1 - round(d1) == 0 - eye(n));
d2 = bsxfun(@minus, exps(:,2), exps(:,2)');
mask2 = triu(d2 - round(d2) == 0 - eye(n));
[i, j] = find( mask1 & mask2 );

% Add compatable funs
remove = [];
while ( ~isempty(i) )
    i1 = i(1); j1 = j(1);
    i(1) = []; j(1) = [];
    i(i == j1) = [];
    
    f.funs{i1} = f.funs{i1} + f.funs{j1};
    remove = [remove ; j1];     %#ok<AGROW>
end
f.funs(remove) = [];

% Remove resulting zero funs:
if ( numel(f.funs) > 1 )
    isz = cellfun(@iszero, f.funs);
    if ( all(isz) )
        f.funs = {f.funs{1}.make(0)};
    else
        f.funs(isz) = [];
    end
end

% 
if ( numel(f.funs) == 1 )
    f = f.funs{1};
end

end
