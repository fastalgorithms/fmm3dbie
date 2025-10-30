function [row_ptr, col_ind, iquad] = getnear(S, targinfo, rfac)
%
%  getnear
%    This subroutine returns the indices of the near quadrature corrections
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
        iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
        norder_avg = floor(sum(norders)/(npatches+0.0d0));
        
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

end