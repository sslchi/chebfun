function f = compose(f, op, varargin)
%COMPOSE   compose command for SPHCAPFUN objects. 
% 
%   F = COMPOSE(F, OP) returns the SPHCAPFUN that approximates OP(F).
% 
%   F = COMPOSE(F, OP, G) returns the SPHCAPFUN that approximates OP(F, G).
%
%   F = COMPOSE(F, G) with a CHEBFUN G with one column returns a SPHCAPFUN that
%   approximates G(F).  If G has two columns, the result is a SPHCAPFUNV. For a
%   CHEBFUN2 or CHEBFUN2V, the composition is interpreted as G(real(F),
%   imag(F)).
%
% This command is a wrapper for the SPHCAPFUN constructor. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(op, 'chebfun') )
    % Composition OP(f) of SPHCAPFUN f and CHEBFUN OP.
    
    if ( length(op.domain) > 2 )
        % If OP has several pieces, OP(SPHCAPFUN) might be inaccurate.
        warning('CHEBFUN:SPHCAPFUN:compose:pieces', ...
            ['The composition of a CHEBFUN with several pieces and a SPHCAPFUN\n', ...
            'might be inaccurate.']);
    end
    
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op.domain, tol) )
        error('CHEBFUN:SPHCAPFUN:COMPOSE:DomainMismatch', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    nColumns = size(op, 2);
    if ( nColumns == 1 )
        % Call constructor:
        f = diskfun(@(x,y) op(feval(f, x, y)));
        
    elseif ( nColumns == 2 )
        % Extract columns of the CHEBFUN OP:
        op1 = op(:,1);
        op2 = op(:,2);
        
        % Call constructor:
        f = diskfunv(@(x,y) op1(feval(f, x, y)), @(x,y) op2(feval(f, x, y)));
        
    else
        % The CHEBFUN object OP has a wrong number of columns.
        error('CHEBFUN:SPHCAPFUN:COMPOSE:Columns', ...
            'The CHEBFUN object must have 1 or 2 columns.')
        
    end
    
elseif ( isa(op, 'chebfun2') )
    % Composition OP(f) of the SPHCAPFUN f and CHEBFUN2 OP, interpreted as
    % OP(real(f), imag(f)).  For now SPHCAPFUNs are real, but this might change
    % in the future.
        
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op.domain(1:2), tol) )
        error('CHEBFUN:SPHCAPFUN:COMPOSE:DomainMismatch2', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = diskfun(@(x,y) op(feval(real(f), x, y), feval(imag(f), x, y)));
    
elseif ( isa(op, 'chebfun2v') )
    % Composition OP(f) of the SPHCAPFUN object f and the CHEBFUN2V OP with two
    % components, interpreted as OP(real(f), imag(f)).
    % For now SPHCAPFUNs are real, but this might change in the future.
        
    % Check that OP has two components:
    if ( op.nComponents ~= 2 )
        error('CHEBFUN:SPHCAPFUN:COMPOSE:C2VnComponents', ...
            'The Chebfun2v objects must have two components.')
    end
    
    % Get the components:
    op1 = op(1);
    op2 = op(2);
    
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op1.domain(1:2), tol) )
        error('CHEBFUN:SPHCAPFUN:COMPOSE:DomainMismatch2v', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = diskfunv(@(x,y) op1(feval(real(f), x, y), feval(imag(f), x, y)), ...
        @(x,y) op2(feval(real(f), x, y), feval(imag(f), x, y)));
    
elseif ( nargin == 2 && nargin(op) == 1)
    % OP has one input variable.
    f = diskfun(@(x,y) op( feval(f, x, y, 'polar') ), 'polar');
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     % promote
        g = diskfun(@(x,y) g + 0*x); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = diskfun(@(x,y) f + 0*x); 
    end
    
    f = diskfun(@(x,y) op( feval(f, x, y, 'polar'), feval(g, x, y, 'polar') ), 'polar'); 
else
    % Not sure what to do, error: 
    error('CHEBFUN:SPHCAPFUN:COMPOSE:OP', 'NARGIN(OP) not correct.')
end

end 

