function d = mergeDomains(varargin)
%MERGEDOMAINS  Merge domains of CHEBMATRIX blocks.
%   D = MERGEDOMAINS(D1,D2,...) merges domains (union of breakpoints, while
%   checking endpoints). Each input is inspected.
%       If numeric scalar, it is ignored.
%       If a numeric vector, it is taken as a domain.
%       If a chebfun or linBlock, its domain is extracted.
%
% See also DOMAIN.MERGE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Start by getting all the domains from VARARGIN:
doms = cell(size(varargin));
for argCounter = 1:nargin
    item = varargin{argCounter};
    if ( isnumeric(item) )
        if ( length(item) > 1)
            doms{argCounter} = item;
        end
        % else do nothing
    elseif ( isa(item, 'chebfun') )
        doms{argCounter} = domain(item);
    else % We must have a linBlock, whose domain we can access directly
        doms{argCounter} = item.domain;
    end
end

% Throw away empties
doms(cellfun(@isempty, doms)) = [];

% If we're left with one domain, return it!
if ( length(doms) == 1 )
    d = doms{1};
    return
    
% If we're left with an empty domain, return it as well.
elseif ( isempty(doms) )
    d = [];
    return
end

%%
% Often in the CHEBMATRIX case we're dealing with objects that live on the same
% domains. The DOMAIN.MERGE() method is pretty expensive, so we check whether we
% can get away with a cheaper version of merging the domains -- namely, if
% they're all of the same size, we check their diffs, and if they agree, we're
% happy!

% Get the lengths of all the domains:
domLengths = cellfun(@length, doms);

% Do the length of the domains agree?
lengthMatch = ~any(domLengths - domLengths(1));
if ( lengthMatch )
    % Compare the other entries of DOMS with the first entry.
    
    % A row vector containing the first domain:
    dom1 = doms{1};
    
    % Function for comparing the domains
    diffFun = @(dom) norm(dom - dom1);
    
    % Subtract the first domain from all the domain cell, and check whether we
    % get any non-zero entries. If so, the domains do not match. Otherwise they
    % match.    
    domainMatch = ~any(cellfun(diffFun, doms));
    
    if ( domainMatch )
        % The domains match! Return the first entry of DOMS:
        d = dom1;
    else
        % Domains do not match, so we need to merge all the domains:
        d = domain.merge(doms{:});
    end
    
else
    % Lengths do not match, so we need to merge all the domains:
    d = domain.merge(doms{:});
end

end
