function f = flipud(g)
%FLIPUD   Flip/reverse a DISKFUN over the x-axis.
%
%   G = FLIPUD(F) returns a DISKFUN G that is flipped over the 
%   x-axis, that is G(x,y) = F(x, -y).
%
% See also DISKFUN/FLIPLR, DISKFUN/FLIPDIM, DISKFUN/ROTATE. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( g ) ) 
    f = diskfun;
    return
end 
f = g; 
% flipud command for trigfuns will not work correctly if # modes is even.
 m = length(g.rows); 
 r = rank(g); 
 if (mod(m,2)==0)
       R = g.rows; 
       rtechs = R.funs{1}.onefun;
       % Alias to get the extra zero in place.
       rtechs.coeffs = rtechs.alias(rtechs.coeffs, m+1); 
       g.rows.funs{1}.onefun = rtechs;
       f.rows = flipud(g.rows); 
       R = f.rows; 
       rtechs = R.funs{1}.onefun;
       % Alias to remove zero.
       rtechs.coeffs = rtechs.alias(rtechs.coeffs, m); 
       f.rows.funs{1}.onefun = rtechs;
       
 else
       f.rows = flipud(g.rows);
 end   
        
        
end
