function novers = kernel3d_getnear_overs(S, t, eps, zk, kernel_order)
%KERNEL3D_GETNEAR_OVERS   Get per-patch oversampling orders for a kernel.
%
%   novers = kernel3d_getnear_overs(S, t, eps, zk, kernel_order)
%
%   Builds the near-interaction struct via getnear(S, t), populates the
%   required Q.wavenumber and Q.kernel_order fields, then calls
%   get_oversampling_parameters(S, Q, eps).
%
%   Input:
%     S            - surfer object (source surface)
%     t            - target info struct with field .r (3 x nt), or surfer
%     eps          - precision requested
%     zk           - wavenumber (0 for Laplace)
%     kernel_order - singularity order: -1 (weak), 0 (pv), 1 (hs)
%
%   Output:
%     novers - (npatches x 1) array of oversampling orders

rsc              = getnear(S, t);
rsc.targinfo     = t;
rsc.wavenumber   = zk;
rsc.kernel_order = kernel_order;
novers           = get_oversampling_parameters(S, rsc, eps);

end
