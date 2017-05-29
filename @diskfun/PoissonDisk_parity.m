function [u, ZZ, DD, YY] = poisson_disk_ADIparity(f, n)
%does ADI to find soln to Poisson's eqn. 
%uses parity prop and returns a diskfun

%this is in progress; need to figure out what to do about pivot locations 
% to make it work right. 


    
    %even-periodic (incl. pole)
    id = f.idxPlus;
    fp = f;
    fp.cols = fp.cols(:, id);
    fp.rows = fp.rows(:, id);
    fp.pivotValues = fp.pivotValues(id, :);
    fp.pivotLocations = fp.pivotLocations(id, :); 
    fp.idxPlus = 1:length(id);
    fp.idxMinus = [];
    
    %call ADI: 
    [ZZ, DD, YY] = poisson_disk_ADI_II(fp, n);
    
    %build a diskfun (currently this is not quite correct)
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    up = fp;
    up.cols = c;
    up.rows = r; 
    up.pivotValues = 1./diag(DD);
    up.idxPlus = 1:length(up.pivotValues); 
    up.idxMinus = []; 
    %up.pivotLocations = nan(length(up.pivotValues), 2); %was not constructed with pivots; need a nonempty
                            %value here
    up.nonZeroPoles = 0; %need to deal with this!                        
   %check for a nonzero pole
 
    vp = feval(up,0,0); 
    if abs(vp) > 1e-14
       up.nonZeroPoles = 1; 
    else
       up.nonZeroPoles = 0; 
    end                        
                            
    % in other routines,
    % it is assumed that all pole information is stored in first row and 
    % col, so we need to adjust our representation.
    % Unfortunately, I'm not sure how this can be done in coeff space using
    % the current compression plus algorithm (relies on fact that pole is
    % already stored correctly).
    
   
    
    ZZ1 = ZZ; 
    DD1 = DD; 
    YY1 = YY; 
    %odd-antiperiodic
    id = f.idxMinus;
    fm = f;
    fm.cols = fm.cols(:, id);
    fm.rows = fm.rows(:, id);
    fm.pivotValues = fm.pivotValues(id);
    fm.pivotLocations = fm.pivotLocations(id, :);
    fm.idxMinus = 1:length(id);
    fm.idxPlus = [];
    fm.nonZeroPoles = 0;
    
    [ZZ, DD, YY] = poisson_disk_ADI_II(fm, n);
    
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    um = fm;
    um.cols = c;
    um.rows = r; 
    um.pivotValues = 1./diag(DD);
    um.idxMinus = 1:length(up.pivotValues); 
    um.idxPlus = []; 
    %um.pivotLocations = nan(length(up.pivotValues), 2); 
    ZZ = [ZZ1 ZZ]; 
    DD = diag([diag(DD1);diag(DD)]); 
    YY = [YY1 YY]; 
    
    %
    %combine results
    
    pivots = [up.pivotValues;um.pivotValues];
    cols = [up.cols um.cols ];
    rows = [up.rows um.rows];
    %locations = [up.pivotLocations;um.pivotLocations];

    numPlus = length(up.pivotValues);
    idxPlus = 1:numPlus;
    numMinus = length(um.pivotValues);
    idxMinus = (numPlus+1):(numPlus+numMinus);
    
    u = up;
    u.cols = cols;
    u.rows = rows;
    u.pivotValues = pivots;
    %u.pivotLocations = locations;
    u.pivotLocations = nan(length(u.pivotValues), 2); 
    u.idxPlus = idxPlus;
    u.idxMinus = idxMinus;
    u.nonZeroPoles = up.nonZeroPoles;
   
    %rotate required due to SVD
    u = flipud(u); 
end