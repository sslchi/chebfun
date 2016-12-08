function u = poisson_disk_ADIparity(f, n)
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
    [~,ZZ, DD, YY] = poisson_disk_ADI(fp, n);
    
    %build a diskfun (currently this is not quite correct)
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    up = fp;
    up.cols = c;
    up.rows = r; 
    up.pivotValues = 1./diag(DD);
    up.idxPlus = 1:length(up.pivotValues); 
    up.idxMinus = []; 
    %u.pivotLocations = not sure what goes here
    %up.nonZeroPoles =; %think more about this
    
   
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
    
    [~, ZZ, DD, YY] = poisson_disk_ADI(fp, n);
    
    c = real(chebfun(chebtech2({'',ZZ})));
    r = real(chebfun(YY,[-pi,pi],'coeffs','periodic'));
    um = fm;
    um.cols = c;
    um.rows = r; 
    um.pivotValues = 1./diag(DD);
    um.idxMinus = 1:length(up.pivotValues); 
    um.idxPlus = []; 
    %u.pivotLocations = not sure what goes here

    
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
    u.idxPlus = idxPlus;
    u.idxMinus = idxMinus;
    u.nonZeroPoles = up.nonZeroPoles;
    
end