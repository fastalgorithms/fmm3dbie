function [ffrm] = get_second_fundamental_form(S,p,varargin)
%
%
%
%
%  Input arguments:
%    * S: surfer object
%    * opts: options struct (optional)
%
%  Output arguments:
%

    opts = [];
    if(nargin == 2)
       opts = varargin{1};
    end
     
% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;


    ffrm = zeros(4,npts);
    mex_id_ = 'get_second_fundamental_form(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], io double[xx])';
[ffrm] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ffrm, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 4, npts);

    ffrm = reshape(ffrm,[2,2,npts]);    

end

%-------------------------------------------------
%
%%
%%   Laplace dirichlet routines
%
%
%-------------------------------------------------

