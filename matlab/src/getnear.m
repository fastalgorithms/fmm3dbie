function [rsc] = getnear(S, targinfo, rfac)
%
%  getnear(S, varargin)
%
%  Syntax: 
%   [rsc] = getnear(S) 
%   [rsc] = getnear(S, targinfo)
%   [rsc] = getnear(S, targinfo, rfac)
%
%  This subroutine returns the row sparse compressed struct
%  for near interactions between the surface S, and targets.
%
%  The row sparse compressed representation is a sparse format
%  between target indices and patch indices of the surface. 
%  It additionally includes an indexing array which points to where 
%  quadrature corrections for interaction between target and patch
%  are stored.
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * targinfo: target info (optional) 
%       targinfo.r = (3,nt) target locations
%      targinfo = S if unspecified
%    * rfac: radius defining near neighborhood (optional)
%        target t_{i} is close to patch \Gamma_{j} if 
%        d(t_{i}, c_{j}) <= rfac*R_{j}
%        where c_{j} is the centroid of \Gamma_{j}, and 
%        R_{j} is the radius of the bounding sphere
%        centered at c_{j}
%
%  Output arguments:
%    * rsc: row sparse compressed struct
%       rsc.row_ptr: (nt+1,1)
%         pointer to col_ind array where list of relevant source
%         patches for target i starts
%       rsc.col_ind: (nnz,1)
%         list of source patches relevant for all targets,
%         sorted by target number
%       rsc.iquad: (nnz+1,1)
%         iquad(i) is the location in quadrature correction array 
%         where quadrature for interaction corresponding to
%         col_ind(i) starts
%      
%    

    [~,~,norders,ixyzs,iptype,~] = extract_arrays(S);
    npatches = S.npatches;
    npts = S.npts;

    if nargin < 2
       targs = S.r;
    else
       if isa(S, 'struct') || isa(S,'surfer')
           targs = targinfo.r;
        else
            targs = targinfo;
        end
    end
    ndtarg = size(targs,1); 
    ntarg = size(targs,2); 

    if nargin < 3
        [iptype_c, iptype_n] = groupcounts(iptype);
        [~, ind] = max(iptype_c);
        iptype_avg = iptype_n(ind);

        [norder_c, norder_n] = groupcounts(norders);
        [~, ind] = max(norder_c);
        norder_avg = norder_n(ind);

        % get nearfield definition
        rfac = 0;
        rfac0 = 0;
        mex_id_ = 'get_rfacs(i int[x], i int[x], io double[x], io double[x])';
[rfac, rfac0] = fmm3dbie_routs(mex_id_, norder_avg, iptype_avg, rfac, rfac0, 1, 1, 1, 1);
    end
    
    cms = S.cms;
    rads = S.rads;

    rad_near = rads*rfac;

    nnz = 0;
    n3 = 3;
    % find the number of near field points
    mex_id_ = 'findnearmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, nnz, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, 1);
    row_ptr = zeros(ntarg+1,1);
    col_ind = zeros(nnz,1);

    ntp1 = ntarg+1;
    nnzp1 = nnz+1;
    % find indices of near field points
    mex_id_ = 'findnear(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x], io int[x])';
[row_ptr, col_ind] = fmm3dbie_routs(mex_id_, cms, npatches, rad_near, ndtarg, targs, ntarg, row_ptr, col_ind, n3, npatches, 1, npatches, 1, ndtarg, ntarg, 1, ntp1, nnz);

    % get location in wnear_ij array where quadrature for col_ind(i) starts for a single kernel. 
    iquad = zeros(nnz+1,1);
    npp1 = npatches+1;
    mex_id_ = 'get_iquad_rsc(i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x])';
[iquad] = fmm3dbie_routs(mex_id_, npatches, ixyzs, npts, nnz, row_ptr, col_ind, iquad, 1, npp1, 1, 1, ntp1, nnz, nnzp1);
    
    rsc = [];
    rsc.row_ptr = row_ptr;
    rsc.col_ind = col_ind;
    rsc.iquad = iquad;

end
