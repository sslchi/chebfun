function val = get(f, propName)
%GET   GET method for the ballfun class.
%   P = GET(F, PROP) returns the property P of the ballfun object F 
%   specified in the string PROP. Valid entries for the string PROP are:
%   'COLS'
%   'ROWS' 
%   'TUBES'
%   'CORE'
%   'DOMAIN'

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop through an array of ballfunS objects.
if ( numel(f) > 1 )
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

% Get the properties.
switch ( propName )
    case 'cols'
        val = f.cols;
    case 'rows'
        val = f.rows;
    case 'tubes'
        val = f.tubes;        
    case 'core'
        val = f.core;
    case 'domain'
        val = f.domain;
    otherwise
        error('CHEBFUN:ballfun:get:propName', ...
            [propName,' is not a valid ballfun property.'])
end