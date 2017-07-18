function g = constructor(g, op, varargin)
%CONSTRUCTOR   The main SPHCAPFUN constructor.
% This code is for generating a SPHCAPFUN object that represents a function
% on the unit disk. A SPHCAPFUN object is a real-valued function that is
% represented as a sum of rank 1 outerproducts of univariate functions in
% `doubled-up' polar coordinates.
%
% The algorithm for constructing a SPHCAPFUN comes in two phases:
%
% PHASE 1: The first phase attempts to determine the numerical rank of the
% function by performing Gaussian elimination with special 2x2 pivoting
% matrices on a tensor grid of sample values. GE is performed until the
% sample matrix is approximated. At the end of this stage we have candidate
% pivot locations and pivot elements.
%
% PHASE 2: The second phase attempts to resolve the corresponding column
% and row slices by sampling along the slices and performing GE (pivoting
% at 2x2 matrices) on the skeleton. Sampling along each slice is increased
% until the Fourier/Chebyshev coefficients of the slice fall below machine
% precision.
%
% See also SPHCAPFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )          % SPHCAPFUN( )
    return
end

% Parse the inputs:
[op, dom, pref, fixedRank, vectorize] = parseInputs(op, varargin{:});

op = @(lam, th) op(lam, (th-dom(3))/(dom(4)-dom(3)));

g = sphcapfun; 
g.diskFunction = diskfun(op, [dom(1:2) 0 1], pref);
g.domain = dom; 

end

% 
% % Return op if construction is from coefficients which is handled in
% % parseInputs.
% if ( isa(op, 'sphcapfun') )  
%     g = op;
%     % Fix the rank:
%     g = fixTheRank(g, fixedRank);
%     return
% end
% 
% % Preferences:
% tech        = chebtech2;  
% tpref       = tech.techPref;
% minSample   = 4;
% maxSample   = tpref.maxLength;
% cheb2Prefs  = pref.cheb2Prefs;
% sampleTest  = cheb2Prefs.sampleTest;
% maxRank     = cheb2Prefs.maxRank;
% pseudoLevel = cheb2Prefs.chebfun2eps;
% 
% alpha = 100; % Default value for coupling parameter
% 
% if ( isa(op, 'char') )     % SPHCAPFUN( CHAR )
%     op = str2op( op );
% end
% 
% % TODO: 
% % 1. Need to allow for different domains.
% % 2. Add non-adaptive construction
% % 3. Add tensor-product.
% % 4. Construction from samples/values: need a check about orientation of
% % the sample.
% 
% % Deal with constructions from numeric data:
% if ( isa(op, 'double') )    % SPHCAPFUN( DOUBLE )
%     g = constructFromDouble(op, dom, alpha, pref);
%     % Fix the rank:
%     g = fixTheRank(g, fixedRank);
%     return
% end
% 
% %
% % Construction is from a function handle.
% %
% 
% % Check for op = @(lam,th) constant function
% [ll, tt] = meshgrid(dom(1:2), dom(3:4));
% if ( ~vectorize && numel(op(ll,tt)) == 1 )
%     op1 = op;
%     op = @(ll, tt) op1(ll, tt) + 0*ll;
% end
% 
% factor  = 8; % Ratio between size of matrix and no. pivots.
% isHappy = 0; % If we are currently unresolved.
% failure = 0; % Reached max discretization size without being happy.
% 
% while ( ~isHappy && ~failure )
%     %
%     % Setup Phase I: GE with block 2-by-2 pivoting to determine the
%     % numerical rank and pivot locations.  Sampling is done on
%     % Fourier-Chebyshev grid.
%     
%     grid = minSample;          
%     happyRank = 0;             % Happy with phase one? 
%     strike = 1;
%     while ( ~happyRank && ~failure && strike < 3)
%         grid = 2*grid;
% 
%         % Sample function on a tensor product grid.
%         [x, y] = getPoints(grid, grid);
%         [xx, yy] = meshgrid(x, y);
%         F = evaluate(op, xx, yy, vectorize);
%         
%         if ( ~isreal( F ) ) 
%             warning('SPHCAPFUN:CONSTRUCTOR:COMPLEX', ...
%                     ['Only real-valued sphcapfuns are currently supported. The '...
%                      'imaginary part is being set to zero now.'])
%              F = real( F );   
%         end
% 
%         [tol, vscale] = getTol(F, 2*pi/grid, 1/grid, dom, pseudoLevel);
%         pref.chebfuneps = tol;
%         
%         % Does the function blow up or evaluate to NaN?:
%         if ( isinf(vscale) )
%             error('CHEBFUN:SPHCAPFUN:constructor:inf', ...
%                 'Function returned INF when evaluated');
%         elseif ( any(isnan(F(:)) ) )
%             error('CHEBFUN:SPHCAPFUN:constructor:nan', ...
%                 'Function returned NaN when evaluated');
%         end
%         
%         % Do GE
%         [pivotIndices, pivotArray, removePoles, happyRank] = ...
%             PhaseOne(F, tol, alpha, factor);
% 
%         if ( grid > factor*(maxRank-1) )
%             warning('SPHCAPFUN:CONSTRUCTOR:MAXRANK', ... 
%                                     'Unresolved with maximum rank.');
%             failure = 1;
%         end
%         
%         % If the function is 0+noise then stop after three strikes.
%         if ( max(abs(pivotArray(1,:))) < 1e4*tol )
%             strike = strike + 1;
%         end
%     end
% 
%     % Do Phase 2: resolve along the column and row slices.
%     [cols, pivots, rows, pivotLocations, idxPlus, idxMinus, isHappy, failure] = ...
%         PhaseTwo(op, pivotIndices, pivotArray, grid, dom, vscale, ...
%         maxSample, removePoles, vectorize, pref);
%     
%     g.cols = chebfun(cols, dom(3:4)-[1 0], pref);
%     g.rows = chebfun(rows, dom(1:2), 'trig', pref);
%     if ( all(pivots) == 0 )
%         pivots = inf;
%     end
%     g.pivotValues = pivots;
%     g.domain = dom;
%     g.idxPlus = idxPlus;
%     g.idxMinus = idxMinus;
%     g.nonZeroPoles = removePoles;
%     g.pivotLocations = adjustPivotLocations(pivotLocations, pivotArray); 
%   
%     % Sample Test:
%     if ( sampleTest )
%         % Wrap the op with evaluate in case the 'vectorize' flag is on:
%         sampleOP = @(th,r) evaluate(op, th, r, vectorize);
%         % Evaluate at points in the domain:
%         pass = g.sampleTest( sampleOP, tol, vectorize);
%         if ( ~pass )
%             % Increase minSamples and try again.
%             minSample = 2*minSample;
%             isHappy = 0;
%         end
%     end
% end
% 
% % Simplifying rows and columns after they are happy.
% g = simplify( g, pseudoLevel );
% 
% % Fix the rank, if in nonadaptive mode.
% g = fixTheRank( g , fixedRank );
% 
% % Project onto BMC-II symmetry so the function is smooth on the disk.
% g = projectOntoBMCII( g );  
% 
% end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [op, dom, pref, fixedRank, vectorize] = parseInputs(op, varargin)

if ( isa(op, 'char') )     % SPHCAPFUN( CHAR )
    op = str2op(op);
end

% If the operator has one argument, then throw an error
if ( isa(op, 'function_handle') )
    % If the operator has one argument, then throw an error
    if ( nargin(op) <= 1 )
        error('CHEBFUN:SPHCAPFUN:CONSTRUCTOR:toFewInputArgs',...
            'The function %s must accept 2 or 3 input arguments.',op);
    % If f is defined in terms of x,y,z; then convert it to
    % (longitude,latitude).
    elseif ( nargin(op) == 3 )
        % Wrap op so it can be evaluated in spherical coordinates
        op = @(lam, th) spherefun.sphf2cartf(op, lam, th, 0);
    end
end

% Get the domain: (Always first if given)
% If domain is empty take it to be co-latitude.
dom = [-pi, pi, 0, 2*pi/3]; 
fixedRank = NaN;
fixedLength = 0;

while ( numel(varargin) > 0 && isnumeric(varargin{1}) )
    d = varargin{1};
    varargin(1) = [];
    
    if ( numel(d) == 4 )                 % SPHCAPFUN(OP, [A B C D])
        dom = d;
%         if ( norm(dom(:)' - [-pi pi 0 pi]) > 0 )
%             error('CHEBFUN:SPHCAPFUN:CONSTRUCTOR:domain',...
%                 ['The only domain presently supported in sphcapfun is [-pi pi]x[0 pi] in '...
%                 'intrinstic (spherical) coordinates, which corresponds to colatitude.']);
%         end
    elseif ( numel(d) == 2 )             % SPHCAPFUN(OP, [M N])
        % Interpret this as the user wants a degree (M,N)
        % sphcapfun
        fixedLength = 1;
        m = d(1); 
        n = d(2);        
    elseif ( numel(d) == 1 )             % SPHCAPFUN(OP, K)
        fixedRank = d;
    else
        error('CHEBFUN:SPHCAPFUN:CONSTRUCTOR:domain',... 
              ['A domain is rarely given for sphcapfun, ', ... 
              'but it needs to be given by four corner values',... 
              'in intrinstic coordinates.'])
    end
end

if ( fixedLength )  % Check that m and n are positive integers
    if ( ( m <= 0 ) || ( n <= 0 ) || ( abs(round(m)-m)  > eps ) || ...
            ( abs(round(n)-n) > eps ) )
        error('CHEBFUN:SPHCAPFUN:constructor:parseInputs:domain2', ...
            ['When constructing with fixed lengths, the values '...
             'for the lengths must be positive integers.']);
    end
end

if ( ( fixedRank < 0 ) || ( abs(round(fixedRank)-fixedRank) > eps ) )
        error('CHEBFUN:SPHCAPFUN:constructor:parseInputs:domain3', ...
            ['When constructing with a fixed rank, the value must '...
             'be a positive integer.']);
end    

% Preferences structure given?
isPref = find(cellfun(@(p) isa(p, 'chebfunpref'), varargin));
if ( any(isPref) )
    pref = varargin{isPref};
    varargin(isPref) = [];
else
    pref = chebfunpref();
end

isEpsGiven = find(cellfun(@(p) strcmpi(p, 'eps'), varargin));
if ( isEpsGiven )
    pseudoLevel = varargin{isEpsGiven+1};
    varargin(isEpsGiven+(0:1)) = [];
else
    pseudoLevel = 0;
end
pref.cheb2Prefs.chebfun2eps = max(pref.cheb2Prefs.chebfun2eps, pseudoLevel);

% Look for vectorize flag:
vectorize = find(cellfun(@(p) strncmpi(p, 'vectori', 7), varargin));
if ( vectorize )
    varargin(vectorize) = [];
    vectorize = true;
else
    vectorize = false;
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~vectorize && ~isnumeric(op) ) % another check
    [vectorize, op] = vectorCheck(op, dom, pref.chebfuneps);
end

isCoeffs = find(cellfun(@(p) strcmpi(p, 'coeffs'), varargin));
if ( isCoeffs )
    varargin(isCoeffs) = [];
    op = sphcapfun.coeffs2sphcapfun(op);
end

% Deal with SPHCAPFUN(OP, [M N]) now that all the other things are set.
if ( fixedLength )
    [x, y] = getPoints(m, n, dom);
    [xx, yy] = meshgrid(x, y);
    % Handle the special case of the input being a sphcapfun.  We can't call
    % evaluate here because, we have to use feval(op,xx,yy).
    if ( isa(op,'sphcapfun') )
        op = feval(op, xx, yy);
    else
        op = evaluate(op, xx, yy, vectorize);
    end    
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = fixTheRank( g , fixedRank )
% Fix the rank of a SPHCAPFUN. Used for nonadaptive calls to the constructor.

if ( fixedRank < 0 )
    error('CHEBFUN:SPHCAPFUN:constructor:fixTheRank:negative', ...
        'Nonadaptive rank should be positive.')
elseif ( fixedRank > 0 )
    if ( length(g.pivotValues) > fixedRank )
        % Truncate things:
        g.cols = g.cols(:,1:fixedRank);
        g.rows = g.rows(:,1:fixedRank);
        g.pivotValues = g.pivotValues(1:fixedRank);
        g.idxPlus = g.idxPlus( g.idxPlus <= fixedRank );
        g.idxMinus = g.idxMinus( g.idxMinus <= fixedRank );
    elseif ( length(g.pivotValues) < fixedRank )
        % Pad things with zero columns:
        zcols = chebfun(0, g.cols.domain);
        zrows = chebfun(0, g.rows.domain, 'trig');
        for jj = length(g.pivotValues) : fixedRank - 1
            g.cols = [g.cols zcols];
            g.rows = [g.rows zrows];
            g.pivotValues = [g.pivotValues ; 0];
        end
    end
elseif ( fixedRank == 0 )
    g.cols = chebfun(0, g.cols.domain);
    g.rows = chebfun(0, g.rows.domain, 'trig'); 
    g.pivotValues = Inf;
    g.idxPlus = [];
    g.idxMinus = 1;
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vectorize, op] = vectorCheck(op, dom, pseudoLevel)
% Check for cases: @(x,y) x*y, and @(x,y) x*y'

vectorize = false;

if isa(op,'sphcapfun')
    return;
end

% Evaluate at a 2-by-2 grid on the interior of the domain.
[xx, yy] = meshgrid( dom(1:2)/3 + diff(dom(1:2))/3,...
                     dom(3:4)/2 + diff(dom(3:4))/3);
try
    A = op(xx, yy);
catch
    throwVectorWarning();
    vectorize = true;
    return
end
if ( isscalar(A) )
    op = @(x,y) op(x,y) + 0*x + 0*y;
    A = op(xx, yy);
end
B = zeros(2);
for j = 1:2
    for k = 1:2
        B(j,k) = op(xx(j,k), yy(j,k));
    end
end
if ( any(any( abs(A - B) > min( 1000*pseudoLevel, 1e-4 ) ) ) )
    % Function handle probably needs vectorizing.
    % Give user a warning and then vectorize.
    throwVectorWarning();
    vectorize = true;
end
    function throwVectorWarning()
        warning('CHEBFUN:SPHCAPFUN:constructor:vectorize',...
            ['Function did not correctly evaluate on an array.\n', ...
            'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
            'Use the ''vectorize'' flag in the SPHCAPFUN constructor\n', ...
            'call to avoid this warning message.']);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = str2op( op )
% OP = STR2OP(OP), finds the dependent variables in a string and returns an
% op handle than can be evaluated.

depvar = symvar( op );
if ( numel(depvar) > 2)
    error('CHEBFUN:SPHCAPFUN:constructor:str2op:depvars', ...
        'Too many dependent variables in string input.');
elseif ( numel(depvar) == 1 )
    % Treat as a complex variable:
    op = eval(['@(' real(depvar{1}) + 1i*imag(depvar{1}) ')' op]);
elseif ( numel(depvar) == 2 )
    op = eval(['@(' depvar{1} ',' depvar{2} ')' op]);
else
    op = eval(['@(' depvar{1} ',' depvar{2} ',' depvar{3} ')' op]);
end

end
