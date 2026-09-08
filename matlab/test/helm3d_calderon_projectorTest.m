%
% Tests the Helmholtz layer potential quadratures, and in particular the
% hypersingular (ipv = 2) quadrature for D', through the Calderon identity
% P x = x applied to the Cauchy data of a point source well separated from
% the surface. The Calderon projector is P = C + I/2 with
%
%     C = [ -D   S  ]
%         [ -D'  S' ]
%
% built as kernel3d('h','trans_sys',zk,[-1 1; -1 1]).
%
% The psi row contains D', whose error constant is much larger than that of
% S, D and S', hence the separate tolerances.
%

run ../startup.m

eps    = 1e-7;
norder = 8;
na     = 2;
zk     = 1.1;

src = [];  src.r = [3.5; 2.1; -1.7];

coefs = [-1, 1; -1, 1];

tol_phi = 1e-3;
tol_psi = 2e-1;

for iptype = [1, 11]

    S    = geometries.ellipsoid([1,1,1.2], na*[1,1,1], [], norder, iptype);
    npts = S.npts;

    assert(mean(sum(S.r .* S.n, 1)) > 0, 'normals are not outward');
    assert(min(vecnorm(S.r - src.r)) > 0.5, 'source is too close to the surface');

    fprintf('\n=== iptype = %d, norder = %d, npatches = %d, npts = %d ===\n', ...
        iptype, norder, S.npatches, npts);

    kern = kernel3d('h', 'trans_sys', zk, coefs);

    Cmat = surfermat(S, kern, eps);
    Pmat = Cmat + 0.5*eye(size(Cmat));

    phi = helm3d.kern(zk, src, S, 's');       phi = phi(:);
    psi = helm3d.kern(zk, src, S, 'sprime');  psi = psi(:);

    x  = reshape([phi.'; psi.'], [], 1);
    Px = Pmat*x;

    res = reshape(Px - x, 2, []).';
    err_phi = norm(res(:,1)) / norm(phi);
    err_psi = norm(res(:,2)) / norm(psi);
    fprintf('  P x = x   phi row : rel err = %.3e\n', err_phi);
    fprintf('            psi row : rel err = %.3e\n', err_psi);

    assert(err_phi < tol_phi, ...
        'iptype %d: P x = x failed in the phi row (%.3e >= %.3e)', ...
        iptype, err_phi, tol_phi);
    assert(err_psi < tol_psi, ...
        'iptype %d: P x = x failed in the psi row (%.3e >= %.3e)', ...
        iptype, err_psi, tol_psi);
end

fprintf('\nhelm3d_calderon_projectorTest: all checks passed\n');
