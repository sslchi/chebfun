function deltas = getDeltas(f)
% TODO: Document.

deltas = [];
for k = 1:numel(f)
    fk = f.blocks{k};
    if ( isa(fk, 'chebfun') && isdelta(fk) )
         deltas = [deltas, get(fk, 'deltas')];
    end
end

if ( ~isempty(deltas) )
    deltas = struct('loc', deltas(1,:), 'mag', deltas(2,:));
else
    deltas = struct('loc',[], 'mag', []);
end

end