function out = sum(f)
%SUM   Definite integral of a SINGFUN on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.
%
% See also CUMSUM, DIFF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers:
%
% The main algorithm:
%
% When the smoothPart of F is a CHEBTECH, that is, it can be written as
% Chebyshev sum, the integral
%
% I = \int_{-1}^{1} F dx = \sum_{0}^{n-1} c_r M_r,
%
% where M_r = \int_{-1}^{1} (1+x)^a(1-x)^b T_r(x) dx is the rth Jacobi moment.
%
% The computation of M_r is treated differently for different a and b:
%
% (I) when a == b, M_r are the Gegenbauer moments, which bear a closed-form
% solution as indicated in [2] and [4].
%
% (II) when a ~= b, M_r are the general Jacobi moments, which can be obtained
% using a three-term recursive relation discussed in [3].
%
% This way, all quadratures in SINGFUN, along with those in CHEBTECH are now
% entirely of Clenshaw-Curtis style.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful References:
%
% [1]. K. Xu and M. Javed, Singfun Working Note, August 2013
%
% [2]. Hunter, D., and Nikolov, G., Gaussian Quadrature of Chebyshev
% Polynomials, J. Comput. Appl. Math. 94 (1998), 123-131.
%
% [3]. Piessens, R., and Branders, M., The Evaluation and Application of Some
% Modified Moments, BIT 13 (1973), 443-450.
%
% [4]. Sommariva, A., Fast construction of Fejer and Clenshaw–Curtis rules for
% general weight functions, Computers & Mathematics with Applications 65
% (2012), 682-693.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = sum(cellfun(@sum, f.funs));


end
