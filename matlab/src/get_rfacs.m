function [rfac, rfac0] = get_rfacs(norder_avg, iptype_avg)
% get multiple of centroid radius that defines the nearfield
    rfac = 0;
    rfac0 = 0;
    mex_id_ = 'get_rfacs(i int[x], i int[x], io double[x], io double[x])';
[rfac, rfac0] = fmm3dbie_routs(mex_id_, norder_avg, iptype_avg, rfac, rfac0, 1, 1, 1, 1);
end

