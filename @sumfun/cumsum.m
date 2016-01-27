function f = cumsum(f, dim)


f.funs{1}
f.funs{2}
if ( nargin == 1 || dim == 1 )
    cumsum(f.funs{1})
    cumsum(f.funs{2})
    
    f.funs = cellfun(@cumsum, f.funs, 'UniformOutput', false);
else
    % TODO:
    error
end

end
