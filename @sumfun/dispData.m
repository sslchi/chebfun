function data = dispData(f)
%DISPDATA   Useful information for DISPLAY at higher levels.
%   DATA = DISPDATA(F) extracts useful information from the given SINGFUN F and
%   the information DATA will be used by DISPLAY at higher levels. Currently, 
%   the only information it extracts is exponents.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

nFuns = numel(f.funs);
data{1}.name = '  sumfun';
data{1}.data = ['  ' int2str(nFuns) ' pieces'];

% More information can be appended to DATA:

end
