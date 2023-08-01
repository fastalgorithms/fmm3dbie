function [surf_grad] = get_surface_grad(S,p,varargin)
%
%
%
%  surf_grad = get_surfgrad(S,p)
%     This subroutine evaluates the surface gradient of a given function
%     p and returns the gradient in cartesian coordinates (default) or 
%     in the dxyz/du, and dxyz/dv basis 
%
%  Input arguments:
%    * S: surfer object
%    * p: input function on surface
%    * opts: options struct (optional)
%        opts.iscartesian (true): return surface gradient in
%           cartesian coordiantes if true, else return
%           and dxyz/du and dxyz/dv components
%
%  Output arguments:
%    * surf_grad: double (3,S.npts) or double (2,npts)
%        depending on the flag opts.iscartesian
%

    opts = [];
    if(nargin == 3)
       opts = varargin{1};
    end

    iscartesian = true;
    if(isfield(opts,'iscartesian'))
        iscartesian = opts.iscartesian;
    end
     
% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    if(isreal(p))
        nd0 = 1;
        p = reshape(p,[1,npts]);
    else
        nd0 = 2;
        p = [real(p(:)) imag(p(:))].';
    end
    n2 = 2;


    nduse = nd0*n2;



    dp = zeros(nduse,npts);
    mex_id_ = 'get_surf_grad(i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[xx], io double[xx])';
[dp] = fmm3dbie_routs(mex_id_, nd0, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, p, dp, 1, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, nd0, npts, nduse, npts);
    
    if(nd0 ~= 1)
        dp = dp(1:2:4,1:npts) + 1j*dp(2:2:4,1:npts); 
    end
    if(~iscartesian)
        surf_grad = dp;
    else
        surf_grad = repmat(dp(1,1:npts),[3,1]).*srcvals(4:6,1:npts) + repmat(dp(2,1:npts),[3,1]).*srcvals(7:9,1:npts);
    end
end

%
%
%
%
%


%-------------------------------------------------
%
%%
%%   Helmholtz dirichlet routines
%
%
%-------------------------------------------------

