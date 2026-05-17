function [tritree] = gettritree(srccoefs, targinfo, rfac, opts)
%
%  tritree = gettritree(srccoefs, targinfo, rfac, opts)
%
%  Syntax: 
%   [tritree] = gettritree(srccoefs, targinfo) 
%   [tritree] = gettritree(srccoefs, targinfo, rfac) 
%   [tritree] = gettritree(srccoefs, targinfo, rfac, opts) 
%
%  This subroutine constructs a tree on the standard simplex
%  T_{0} = {(u,v): 0<u<1, 0<v<1, u+v<1}
%
%  such that for the all targets, there exists a sub-partition of
%  the triangle satisfying T_{0} = \cup_{j=1}^{N_{i}} T_{j,i}
%
%  with each d(t_{i},T_{j,i})>=rfac*R_{j,i} where R_{j,i} is the radius
%  of the bounding sphere centered at the centroid of T_{j,i}
%
%  Input arguments:
%    srccoefs: (9,npols) koornwinder polynomial coefficients of
%                        [r; du; dv]
%    * targinfo: target info  
%       targinfo.r = (3,nt) target locations
%    * rfac: radius defining the target specific partition (optional)
%       default, rfac = 3
%    * opts: options struct
%        opts.ntrimax - maximum number of triangles in the partition (3000)
%        opts.nlevmax - maximum depth of tree (6)
%    
%
   
    if (nargin < 3), rfac = 3; end
    if (nargin < 4), opts = []; end

    ntrimax = 3000;    if isfield(opts, 'ntrimax'), ntrimax = opts.ntrimax; end
    nlevmax = 8;       if isfield(opts, 'nlevmax'), nlevmax = opts.nlevmax; end

    npols = size(srccoefs, 2);
    npatches = 1;
    norder = (-3 + sqrt(8*npols + 1))/2;
    
    if isa(targinfo, 'struct')
        targs = targinfo.r;
    else
        targs = targinfo;
    end

    ntarg = size(targs, 2);
    
    itargptr = 1;
    ntargptr = ntarg;

    ichild_start = zeros(ntrimax,1);
    da = zeros(ntrimax,1);
    tricm = zeros(3,ntrimax);
    trirad = zeros(ntrimax,1);
    tverts = zeros(6,ntrimax);
    itrirel = zeros(ntarg,ntrimax);

    ier = 0;
    ntri = 0;
    nlev = 0;
    mex_id_ = 'gettritree(i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], io int64_t[x], io int64_t[x], io int64_t[x], io double[x], io double[xx], io double[x], io double[xx], io int64_t[xx], io int64_t[x])';
[ntri, nlev, ichild_start, da, tricm, trirad, tverts, itrirel, ier] = kern_routs(mex_id_, npatches, norder, npols, srccoefs, ntarg, targs, itargptr, ntargptr, ntrimax, nlevmax, rfac, ntri, nlev, ichild_start, da, tricm, trirad, tverts, itrirel, ier, 1, 1, 1, 9, npols, 1, 3, ntarg, 1, 1, 1, 1, 1, 1, 1, ntrimax, ntrimax, 3, ntrimax, ntrimax, 6, ntrimax, ntarg, ntrimax, 1);
   
    tritree = [];
    if ier > 0
       warning('FMM3DBIE.GETTRITREE: too few triangles allocated, tritree not formed');
       return;
    else
      tritree.ichild_start = ichild_start(1:ntri);
      tritree.da = da(1:ntri);
      tritree.tricm = tricm(:,1:ntri);
      tritree.trirad = trirad(1:ntri);
      tritree.tverts = reshape(tverts(:,1:ntri), [2,3,ntri]);
      tritree.itrirel = itrirel(:,1:ntri);
    end
end


%
%
%
