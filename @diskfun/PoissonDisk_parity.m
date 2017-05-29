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
    up.pivotLocations = nan(length(up.pivotValues), 2); 
    up.nonZeroPoles = 0;                       
   %check for a nonzero pole
 
    vp = feval(up,0,0); 
    if abs(vp) > 1e-14
       up.nonZeroPoles = 1; 
    else
       up.nonZeroPoles = 0; 
    end                        
    
    %the ADI algorithm results in information about the pole
    % being distributed across all rows and cols. 
    % We assume in all other diskfun routines that all
    % pole information is contained in first row/col only.
    %  Not sure how to fix this without 
    % using value space. 
%     m = length(f.rows); 
%     nn = length(f.cols); 
%     uv = sample(up, m +mod(m,2), nn);
%     uv = sample(up); 
%     uv = diskfun(uv); 
    vp = diskfun(@(x,y) 0*x-vp); %build const. diskfun value at pole. 
    %upnew = uv-vp; 
    %up = uv; %test whether correct pole info fixes the accuracy!
    %upzp = @separableApprox.plus(up,vp, 1); %separate into zero-pole part
                                    % using compression plus.
   
    
fScl = diag(1./up.pivotValues);
gScl = diag(1./vp.pivotValues);
cols = [up.cols, vp.cols];
rows = [up.rows, vp.rows];

[Qcols, Rcols] = qr(cols);
[Qrows, Rrows] = qr(rows);

Zt = zeros(length(fScl), length(gScl));
Dt = [ fScl, Zt ; Zt.', gScl ];
[U, S, V] = svd(Rcols * Dt * Rrows.');
% Take diagonal from SIGMA:
s = diag(S);

% Compress the format if possible.
% [TODO]: What should EPS be in the tolerance check below?
vup = vscale(up); 
vvp = vscale(vp);

vscl = 2*max(vup, vvp); 
% Remove singular values that fall below eps*vscale: 
idx = find( s > 10*eps * vscl, 1, 'last');

if ( isempty(idx) )
    % Return 0 separableApprox
    h = 0*f;
else
    U = U(:,1:idx);
    V = V(:,1:idx);
    s = s(1:idx);
    upzp = up;
    upzp.cols = Qcols * U;
    upzp.rows = Qrows * conj(V);
    % [TODO]: PivotValues have very little meaning after this compression step.
    % For now we assign the singular values as the pivot values. 
    upzp.pivotValues = 1./s;
end                               
                                    
                                                             
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
    um.idxMinus = 1:length(um.pivotValues); 
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
   
    u = flipud(u); 
end