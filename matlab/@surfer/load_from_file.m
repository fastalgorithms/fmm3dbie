function obj = load_from_file(ff)
    dat = load(ff);
    norder  = dat(1);
    npatches = dat(2);
    npols = (norder+1)*(norder+2)/2;
    npts = npatches*npols;
    srcvals = dat(3:end);
    srcvals = reshape(srcvals,[npts,12]);
    srcvals = srcvals';
    obj = surfer(npatches,norder,srcvals);
    
end
