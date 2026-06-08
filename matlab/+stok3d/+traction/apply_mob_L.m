function pot = apply_mob_L(sigma, S, icomps, area, centroid, rmoi_inv)
% Apply the generalized ones-matrix L to sigma (matrix-free).
%
% L[sigma]|_{Gamma_i} = rint/area + cross(Ri * rint_torque, x-xc)
%   rint        = int sigma dS
%   rint_torque = int (x-xc) x sigma dS
%   Ri          = rmoi_inv(:,:,i)

npts  = size(sigma, 2);
ncomp = length(area);
pot   = zeros(3, npts);
pts   = @(ic) S.ixyzs(icomps(ic)) : S.ixyzs(icomps(ic+1))-1;

for ic = 1:ncomp
    dx  = S.r(:, pts(ic)) - centroid(:,ic);
    wi  = S.wts(pts(ic));
    sig = sigma(:, pts(ic));

    rint        = sig * wi;
    rint_torque = cross(dx, sig, 1) * wi;

    omega = rmoi_inv(:,:,ic) * rint_torque;
    pot(:, pts(ic)) = rint/area(ic) + cross(omega*ones(1,numel(pts(ic))), dx);
end
end
