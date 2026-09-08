function v = mixed_matvec(x, srfrs, alpha, eps, Qcor, Sovers, ...
    area1, centroid1, rmoi_inv1)
%STOK3D.MIXED_MATVEC  Matrix-free apply for the Stokes mixed BC system.
%
% System: [S'+L    (aS+D)' ] [sigma1]   S1 targets (traction BC, rigid body)
%         [S       aS+D    ] [sigma2]   S2 targets (velocity BC)
%
% One stfmm3d call (ifppregtarg=3), all sources to all targets.
% L stabilizes the rigid-body null space on S1.

S1 = srfrs(1);  S2 = srfrs(2);
n1 = 3*S1.npts;

sigma1 = reshape(x(1:n1),     3, S1.npts);
sigma2 = reshape(x(n1+1:end), 3, S2.npts);

S1over = Sovers{1}{1,1};  xi1 = Sovers{2}{1,1};
S2over = Sovers{1}{2,2};  xi2 = Sovers{2}{2,2};

sig1w = reshape(xi1*sigma1(:), 3, S1over.npts) .* S1over.wts';
sig2w = reshape(xi2*sigma2(:), 3, S2over.npts) .* S2over.wts';

srcinfo = struct('sources', [S1over.r,            S2over.r        ], ...
                 'stoklet',  [sig1w,               alpha*sig2w     ], ...
                 'strslet',  [zeros(3,S1over.npts), sig2w          ], ...
                 'strsvec',  [zeros(3,S1over.npts), S2over.n       ]);

out = stfmm3d(eps, srcinfo, 0, [S1.r, S2.r], 3);

row1 = trac_from_fmm(out, S1.n, 1, S1.npts);
row2 = reshape(out.pottarg(:, S1.npts+1:end), 3, S2.npts);

% Stabilize S1 block against rigid-body null space
icomps1 = [1; S1.npatches+1];
L1 = stok3d.traction.apply_mob_L(sigma1, S1, icomps1, area1, centroid1, rmoi_inv1);
row1 = row1 + L1;

v = [row1(:); row2(:)] + Qcor * x;
end

function t = trac_from_fmm(out, n, i1, nt)
idx = i1 : i1+nt-1;
pre = reshape(out.pretarg(idx), 1, nt);
g11 = reshape(out.gradtarg(1,1,idx), 1, nt);
g12 = reshape(out.gradtarg(1,2,idx), 1, nt);
g13 = reshape(out.gradtarg(1,3,idx), 1, nt);
g21 = reshape(out.gradtarg(2,1,idx), 1, nt);
g22 = reshape(out.gradtarg(2,2,idx), 1, nt);
g23 = reshape(out.gradtarg(2,3,idx), 1, nt);
g31 = reshape(out.gradtarg(3,1,idx), 1, nt);
g32 = reshape(out.gradtarg(3,2,idx), 1, nt);
g33 = reshape(out.gradtarg(3,3,idx), 1, nt);
t = zeros(3, nt);
t(1,:) = n(1,:).*(-pre+2*g11) + n(2,:).*(g12+g21) + n(3,:).*(g13+g31);
t(2,:) = n(1,:).*(g21+g12)    + n(2,:).*(-pre+2*g22) + n(3,:).*(g23+g32);
t(3,:) = n(1,:).*(g31+g13)    + n(2,:).*(g32+g23)    + n(3,:).*(-pre+2*g33);
end
