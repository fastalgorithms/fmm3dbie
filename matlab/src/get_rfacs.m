function [rfac, rfac0] = get_rfacs(norder_avg, iptype_avg)
    rfac = 0;
    rfac0 = 0;
    mex_id_ = 'get_rfacs(i int64_t[x], i int64_t[x], io double[x], io double[x])';
[rfac, rfac0] = kern_routs(mex_id_, norder_avg, iptype_avg, rfac, rfac0, 1, 1, 1, 1);
end

