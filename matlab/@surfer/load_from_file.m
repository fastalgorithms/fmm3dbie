function obj = load_from_file(ff)
% LOAD_FROM_FILE  Create surfer object by loading a .go3 file
%
% S = load_from_file('file.go3')
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
