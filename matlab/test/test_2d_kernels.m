% Test the 2D kernels (lap2d, helm2d, flex2d) on a flat surfer.
%
% For f = 1 on the unit disk the volume potentials have closed forms:
%
%   Laplace:   N[1](x)  = (1 - r^2)/4
%              dN/dn    = -1/2                              on r = 1
%   Helmholtz: V[1](x)  = (i pi/(2k)) J0(kr) H1(k) - 1/k^2
%              dV/dn    = -(i pi/2) J1(k) H1(k)             on r = 1
%
% Flexural kernels are checked for
% self-convergence under refinement.
%
% Boundary targets are tagged with patch_id / uvs_targ, with the disk
% patches rotated so that their boundary edge is the v = 0 edge.

test_lap2d_volume();
test_helm2d_volume();
test_flex2d_convergence();
test_fmm();


function [S, tbd] = disk_and_boundary(norder)
S = aligned_disk([3,3,3], norder);
th = linspace(0, 2*pi, 41); th = th(1:end-1);
tbd = []; tbd.r = [cos(th); sin(th); zeros(size(th))];
tbd.n = [cos(th); sin(th); zeros(size(th))];
tbd = add_patch_ids(S, tbd);
end


function S = aligned_disk(npars, norder)
% unit disk with each boundary patch's v = 0 edge on the circle
S = geometries.disk([], [], npars, norder);
vs = cell2mat(S.end_pt_verts(:).');
S = align_patches(S, 1 - vecnorm(vs), []);
end


function t = add_patch_ids(S, t)
% on-surface patch ids and uvs for boundary targets
[~, pid, uvs, ~, flags] = get_closest_pts(S, t);
assert(all(flags == 0), 'closest point search failed');
assert(max(abs(uvs(2,:))) < 1e-10, 'boundary targets not on v = 0 edge');
t.patch_id = pid(:);
t.uvs_targ = uvs;
end


function test_lap2d_volume()

tol = 1e-7;
eps = 1e-10;
[S, tbd] = disk_and_boundary(8);
f = ones(S.npts, 1);

% on-surface (volume-to-volume)
A = surfermat(S, kernel3d('lap2d', 's'), eps);
r2 = S.r(1,:).^2 + S.r(2,:).^2;
uex = (1 - r2(:))/4;
err = norm(A*f - uex, inf) / norm(uex, inf);
fprintf('lap2d  s  (v2v): rel err %5.2e\n', err);
assert(err < tol);

% normal derivative on the boundary curve
B = surferkernevalmat(S, kernel3d('lap2d', 'sp'), tbd, eps);
err = norm(B*f + 0.5, inf) / 0.5;
fprintf('lap2d  sp (v2b): rel err %5.2e\n', err);
assert(err < tol);

end


function test_helm2d_volume()

tol = 1e-7;
eps = 1e-10;
zk  = 2.3;
[S, tbd] = disk_and_boundary(8);
f = ones(S.npts, 1);

A = surfermat(S, kernel3d('helm2d', 's', zk), eps);
r = sqrt(S.r(1,:).^2 + S.r(2,:).^2).';
uex = (1i*pi/(2*zk))*besselj(0, zk*r)*besselh(1, 1, zk) - 1/zk^2;
err = norm(A*f - uex, inf) / norm(uex, inf);
fprintf('helm2d s  (v2v): rel err %5.2e\n', err);
assert(err < tol);

B = surferkernevalmat(S, kernel3d('helm2d', 'sp', zk), tbd, eps);
dex = -(1i*pi/2)*besselj(1, zk)*besselh(1, 1, zk);
err = norm(B*f - dex, inf) / abs(dex);
fprintf('helm2d sp (v2b): rel err %5.2e\n', err);
assert(err < tol);

end


function test_flex2d_convergence()

tol = 1e-6;
eps = 1e-10;
zk  = 1.7;
nu  = 0.3;

% boundary targets need n, d, d2 for the plate kernels
th = linspace(0, 2*pi, 21); th = th(1:end-1);
tbd0 = []; tbd0.r = [cos(th); sin(th); zeros(size(th))];
tbd0.n  = [cos(th);  sin(th); zeros(size(th))];
tbd0.d  = [-sin(th); cos(th); zeros(size(th))];
tbd0.d2 = [-cos(th); -sin(th); zeros(size(th))];

types = {'s', 'clamped_plate_bcs', 'supported_plate_bcs', 'free_plate_bcs'};
vals = cell(2, numel(types));
for iref = 1:2
    S = aligned_disk([], 6 + 2*iref);
    f = ones(S.npts, 1);
    tbd = add_patch_ids(S, tbd0);
    for it = 1:numel(types)
        K = kernel3d('flex2d', types{it}, zk, nu);
        vals{iref, it} = surferkernevalmat(S, K, tbd, eps) * f;
    end
end

for it = 1:numel(types)
    err = norm(vals{1,it} - vals{2,it}, inf) / norm(vals{2,it}, inf);
    fprintf('flex2d %-20s: self-conv %5.2e\n', types{it}, err);
    assert(err < tol);
end

end


function test_fmm()
% fmm vs direct sum, skipped if fmm2d is not on the path

if ( exist(['fmm2d.' mexext], 'file') ~= 3 )
    fprintf('fmm2d not found, skipping fmm tests\n');
    return
end

tol = 1e-8;
S = geometries.disk([], [], [], 4);
sigma = randn(S.npts, 1) + 1i*randn(S.npts, 1);
th = linspace(0, 2*pi, 31); th = th(1:end-1);
t = []; t.r = 1.5*[cos(th); sin(th); zeros(size(th))];
t.n = [cos(th); sin(th); zeros(size(th))];
t.d = [-sin(th); cos(th); zeros(size(th))];

kerns = {kernel3d('lap2d', 's'), kernel3d('lap2d', 'sp'), ...
         kernel3d('helm2d', 's', 2.3), kernel3d('helm2d', 'sp', 2.3), ...
         kernel3d('helm2d', 's2trans', 2.3, [1; 2]), ...
         kernel3d('flex2d', 's', 1.7), ...
         kernel3d('flex2d', 'clamped_plate_bcs', 1.7), ...
         kernel3d('flex2d', 'supported_plate_bcs', [1.7, 0.5i], 0.3), ...
         kernel3d('flex2d', 'clamped_plate_bcs', [1.7, 0])};

for i = 1:numel(kerns)
    K = kerns{i};
    assert(~isempty(K.fmm), 'missing fmm for %s %s', K.name, K.type);
    pf = K.fmm(1e-12, S, t, sigma);
    pd = K.eval(S, t)*sigma;
    err = norm(pf - pd) / norm(pd);
    fprintf('%-7s %-20s fmm err %5.2e\n', K.name, K.type, err);
    assert(err < tol);
end

assert(isempty(kernel3d('flex2d', 'free_plate_bcs', 1.7, 0.3).fmm));
assert(isempty(kernel3d('flex2d', 's', 0).fmm));

end
