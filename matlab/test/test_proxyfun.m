% Test byindex.proxy_points, byindex.proxyfun and byindex.flamfuns. The
% FLAM tests only run if rskelf and ifmm are on the path.
run ../startup.m
rng(42);

%% Now run the tests

test_proxy_geometry(2);
test_proxy_geometry(3);
test_stretched_box();
test_tight_klb();
test_em3d_kern_out();
test_kern_in();
test_components();
test_offsurface_groups();
test_nan();
test_occ_split_planes();
test_flam_rskelf();
test_2d();


function test_proxy_geometry(dim)
% Rings (dim = 2) and spheres: frames, weights, far-field compression,
% rescaling and unsupported fields.

eps = 1e-4;
if dim == 2, K = kernel3d('lap2d', 's'); else, K = kernel3d('l', 's'); end
l = ones(dim, 1);  ctr = zeros(dim, 1);
p = byindex.proxy_points(K, l, eps, struct('nsrc', 24, 'nproxy_max', 768, 'rmax', 4));

assert(numel(p.width) == dim && all(p.w > 0), 'wrong proxy geometry');
assert(max(abs(vecnorm(p.n) - 1)) < 1e-12, 'proxy normals must be unit vectors');
assert(max(abs(sum(p.n.*p.du, 1))) < 1e-12, 'proxy tangents must be orthogonal');
assert(all(sum(cross(p.du, p.dv).*p.n, 1) > 0), 'du x dv must be along n');
if dim == 2
    assert(all(p.r(3,:) == 0) && isequal(p.d, p.du), 'wrong ring frame');
    assert(abs(sum(p.w) - 2*pi*p.radius) < 1e-8, 'wrong ring length');
else
    assert(abs(sum(p.w) - 4*pi*p.radius^2) < 1e-8, 'wrong shell area');
end

S = box_sources(0.6*l, 128);
check_far(S, K, p, l, eps);

% Reusing the surface for a box scaled by 1/2 and shifted scales the
% ring weights by length and the shell weights by area
sc = 0.5;  ctr = (1:dim).';  pad = [ctr; zeros(3-dim, 1)];
S.r = sc*S.r + pad;  S.wts = sc^2*S.wts;
pxy = struct('r', sc*p.r + pad);
x = S.r(1:dim,:);
C = byindex.proxyfun(S, K, [], p, 'c', x, x, 1:3, [], sc*l, ctr);
R = byindex.proxyfun(S, K, [], p, 'r', x, x, 1:3, [], sc*l, ctr);
s3 = struct('r', S.r(:,1:3));
ref = K.eval(s3, pxy).*S.wts(1:3).';
assert(norm(C - ref, 'fro') < 1e-12*norm(ref, 'fro'), 'outgoing scaling failed');
ref = K.eval(pxy, s3).*(sc^(dim-1)*p.w).';
assert(norm(R - ref, 'fro') < 1e-12*norm(ref, 'fro'), 'incoming weight scaling failed');

bad = K;  bad.targ_fields = {'kappa'};
caught = false;
try
    byindex.proxy_points(bad, l, eps);
catch ME
    caught = contains(ME.message, 'kappa');
end
assert(caught, 'unsupported proxy fields must be rejected');

end


function test_stretched_box()
% A box flattened along z: equal-area points on each shell, and the
% far field of a target-derivative kernel compressed both ways.

eps = 1e-4;
l = [1; 1; 0.5];
K = kernel3d('h', 'dp', 1.5);
p = byindex.proxy_points(K, l, eps, struct('rmax', 4));
assert(isequal(p.width, l), 'wrong proxy box');

i0 = [0 cumsum(p.nproxy)];
for k = 1:p.nshell
    w = p.w(i0(k)+1:i0(k+1));
    assert(max(w)/min(w) - 1 < 1e-6, 'shell %d: points not equal-area', k);
end

check_far(box_sources(l, 500), K, p, l, eps);

end


function test_tight_klb()
% klb to 1e-10, near the ID cutoff of eps/100.

l = ones(3, 1);
K = kernel3d('bel', 'klb');
p = byindex.proxy_points(K, l, 1e-10, struct('rmax', 4));
check_far(box_sources(l, 500), K, p, l, 1e-10);

end


function test_em3d_kern_out()
% em3d densities are in the (du, dv) frame, which is not smooth on the
% proxies, so the outgoing kernel takes Cartesian currents.

eps = 1e-4;
zk = 1.1;
K = kernel3d('em', 'nrccie-bc', zk, 0.5);
Kout = kernel3d('em', 'nrccie-eval', zk)*@cartesian;
l = ones(3, 1);
p = byindex.proxy_points(Kout, l, eps, struct('rmax', 4));
p.kern_out = Kout;

S = box_sources(l, 64);
far = far_targets(l, 32);
x = S.r;
C = byindex.proxyfun(S, K, [], p, 'c', x, x, 1:3*S.npts, [], l, zeros(3, 1));
check_id(K.eval(S, far).*repelem(S.wts, 3).', C, 10*eps);

end


function test_kern_in()
% An incoming kernel that oscillates more than the outgoing one sets the
% resolution, and reproduces its far field on the targets.

eps = 1e-4;
l = ones(3, 1);
Kout = kernel3d('l', 's');
Kin = kernel3d('h', 's', 4);
opts = struct('rmax', 4);
p0 = byindex.proxy_points(Kout, l, eps, opts);
opts.kern_in = Kin;
p = byindex.proxy_points(Kout, l, eps, opts);
assert(sum(p.nproxy) > sum(p0.nproxy), 'kern_in must refine the proxy surface');
assert(isequal(p.kern_in, Kin), 'proxy.kern_in not kept');

S = box_sources(l, 300);
far = far_targets(l, 48);
R = byindex.proxyfun(S, Kin, [], p, 'r', S.r, S.r, 1:S.npts, [], l, zeros(3, 1));
check_id(Kin.eval(far, S).', R.', 10*eps);

end


function test_offsurface_groups()
% Two off-surface target groups, the second needing normals, with a 2x1
% kernel array: one proxy copy per group, and rows from both in a box.

eps = 1e-4;  tol = 1e-3;

S    = geometries.sphere(2, 8, [0;0;0], 4, 1);
kern = [kernel3d('h', 's', 1.2); kernel3d('h', 'sp', 1.2)];
l = 1;
proxy = byindex.proxy_points(kern, l, eps);

ntarg = 100;
ctr = [2; 0; 0];
tg = struct('r', {[], []}, 'n', {[], []});
for k = 1:2
    th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
    tg(k).n = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
    tg(k).r = ctr + (2 + rand(1,ntarg)) .* tg(k).n;
end
rx = [tg.r];

% columns: sources on S near ctr
srcids   = find(vecnorm(S.r - ctr) < 0.5);
src2targ = byindex.kernbyindexeval((1:2*ntarg).', srcids, S, kern, tg);
[Kpxy, nbr] = byindex.proxyfun(S, kern, tg, proxy, 'c', rx, S.r, ...
    srcids, 1:2*ntarg, l, ctr);
assert(isempty(nbr), 'far targets should be removed from nbr');
assert(size(Kpxy, 1) == 2*size(proxy.r, 2), 'one proxy copy per target group');
check_id(src2targ, Kpxy, tol);

% rows: targets of both groups in a box at the origin
ctr = [0; 0; 0];
for k = 1:2
    tg(k).r = 0.5*rand(1,ntarg) .* tg(k).n;
end
rx = [tg.r];
targ2src = byindex.kernbyindexeval((1:2*ntarg).', (1:S.npts).', S, kern, tg);
Kpxy = byindex.proxyfun(S, kern, tg, proxy, 'r', rx, S.r, ...
    1:2*ntarg, 1:S.npts, l, ctr);
assert(isequal(size(Kpxy), [2*ntarg, size(proxy.r, 2)]), 'wrong rc = r size');
check_id(targ2src.', Kpxy.', tol);

end


function test_components()
% kern_out/kern_in with one-component blocks applied to each component,
% blocks shared across surfers, and the stacked rskelf block.

K = kernel3d('l', 's');
p = struct('r', 2*eye(3), 'n', eye(3), 'du', eye(3), 'dv', eye(3), ...
    'w', [2; 3; 4], 'width', ones(3, 1), 'radius', 1.5, 'ifaddtrans', true, ...
    'kern_out', [K; 2*K], 'kern_in', [K, 2*K, 3*K]);
S = struct('r', rand(3, 4)/4, 'npts', 4, 'wts', (1:4).'/4);
T = S;  T.r = rand(3, 4)/4;
srfs = [S, T];
allsrc = struct('r', [S.r, T.r]);
targ = struct('r', rand(3, 3)/4);
main = kernel3d(@(s, t) kron(K.eval(s, t), [1 2 3; 4 5 6]));
rx = repelem(targ.r, 1, 2);  cx = repelem(allsrc.r, 1, 3);
jj = [24 1 14 8];  ii = [6 1 4 2];

C = byindex.proxyfun(srfs, main, targ, p, 'c', rx, cx, jj, [], 1, zeros(3, 1));
B = kron(K.eval(allsrc, p), eye(3)).*repelem([S.wts; T.wts].', 3);
ref = [B; 2*B];
assert(norm(C - ref(:,jj), 'fro') < 1e-12, 'outgoing components/blocks failed');

R = byindex.proxyfun(srfs, main, targ, p, 'r', rx, cx, ii, [], 1, zeros(3, 1));
B = kron(K.eval(p, targ), eye(2)).*repelem(p.w.', 2);
ref = [B, 2*B, 3*B];
assert(norm(R - ref(ii,:), 'fro') < 1e-12, 'incoming components/blocks failed');

main = kernel3d(@(s, t) kron(K.eval(s, t), [1 2; 3 4]));
x = repelem(allsrc.r, 1, 2);  slf = [16 1 10 4];
C = byindex.proxyfun(srfs, main, [], p, 'c', x, x, slf, [], 1, zeros(3, 1));
R = byindex.proxyfun(srfs, main, [], p, 'r', x, x, slf, [], 1, zeros(3, 1));
B = byindex.proxyfun(srfs, main, [], p, x, slf, [], 1, zeros(3, 1));
assert(norm(B - [C; R.'], 'fro') < 1e-12, 'nonsymmetric rskelf stacking failed');

p.kern_out = kernel3d(@(s, t) kron(K.eval(s, t), [1 2 3]));
caught = false;
try
    byindex.proxyfun(srfs, main, [], p, 'c', x, x, slf, [], 1, zeros(3, 1));
catch ME
    caught = strcmp(ME.identifier, 'BYINDEX:proxyfun:proxyKernel');
end
assert(caught, 'incompatible physical components must be rejected');

end


function test_nan()
% NaN kernels error at once, and radpp kernels are NaN beyond their fit.

K = kernel3d(@(s, t) nan(size(t.r, 2), size(s.r, 2)));
caught = false;
try
    byindex.proxy_points(K, 1, 1e-4);
catch ME
    caught = strcmp(ME.identifier, 'BYINDEX:proxy_points:nan');
end
assert(caught, 'NaN kernels must be rejected');

K = kernel3d.radpp(@lap2d.green, 6, 'sg');
v = K.eval(struct('r', zeros(3, 1)), struct('r', [3 7; 0 0; 0 0]));
assert(all(isfinite(v(1:3))) && all(isnan(v(4:6))), 'radpp must be NaN beyond its fit');

end


function test_occ_split_planes()
[x, y] = meshgrid(linspace(-1, 1, 9));
S = struct('r', [x(:).'; y(:).'; zeros(1, 81)], ...
    'npts', 81, 'norders', 1, 'wts', ones(81, 1)/81);
Q = sparse(1, 3, 1, 81, 81);
K = kernel3d.ones();
K.getquad = @(varargin) Q;
opts = struct('type', 'rskelf', 'occ', 8, 'occ_max', 9, ...
    'proxyopts', struct('nsrc', 8, 'nproxy_max', 32));
funs = byindex.flamfuns(S, K, 1e-4, [], opts);
assert(funs.occ == 9, 'split-plane points changed the occupancy estimate');

end


function test_flam_rskelf()
% rskelf from byindex.flamfuns on a flattened ellipsoid with the Debye
% pair: I + R + W, klb outgoing and the full kernel incoming.

if ~has_flam()
    fprintf('FLAM not found, skipping test_flam_rskelf\n');
    return
end

eps = 1e-4;  tol = 1e-3;
S = geometries.ellipsoid([1 1 0.4], [], [], 3);
KF = kernel3d('bel', 'rlb') + kernel3d.ones()/S.area;
opts = struct('type', 'rskelf', 'occ', 64, ...
    'kern_out', kernel3d('bel', 'klb'), 'kern_in', KF);
funs = byindex.flamfuns(S, KF, eps, [], opts);
F = rskelf(funs.Afun, funs.x, funs.occ, eps, funs.pxyfun, flam_opts('rskelf'));

A = funs.Afun(1:S.npts, 1:S.npts);
sig = randn(S.npts, 1);
u = A*sig;
err = norm(rskelf_mv(F, sig) - u)/norm(u);
assert(err < tol, 'rskelf: err %.2e > tol %.2e', err, tol);

end


function test_2d()
% A flat 2 x 1 rectangle: rings with points evenly spaced in length, and
% rskelf and ifmm (with targets of shape 2 x k x nch) from
% byindex.flamfuns, which build 2D trees with unequal box sides.

eps = 1e-4;  tol = 1e-3;
K = kernel3d('lap2d', 's');
l = [2; 1];
p = byindex.proxy_points(K, l, eps, struct('rmax', 6));
i0 = [0 cumsum(p.nproxy)];
for k = 1:p.nshell
    w = p.w(i0(k)+1:i0(k+1));
    assert(max(w)/min(w) - 1 < 1e-6, 'ring %d: points not evenly spaced', k);
end
check_far(box_sources(l, 300), K, p, l, eps);

if ~has_flam()
    fprintf('FLAM not found, skipping FLAM in test_2d\n');
    return
end
S = geometries.square([8 4], 4).scale([2 1 1]);
sig = randn(S.npts, 1);

funs = byindex.flamfuns(S, K, eps, [], struct('type', 'rskelf', 'occ', 64));
assert(size(funs.x, 1) == 2, 'expected a 2D tree');
F = rskelf(funs.Afun, funs.x, funs.occ, eps, funs.pxyfun, flam_opts('rskelf'));
u = funs.Afun(1:S.npts, 1:S.npts)*sig;
err = norm(rskelf_mv(F, sig) - u)/norm(u);
assert(err < tol, 'rskelf 2D: err %.2e > tol %.2e', err, tol);

targs = [2.5 + 2*rand(1, 16, 25); 2*rand(1, 16, 25) - 1];
funs = byindex.flamfuns(S, K, eps, targs, struct('occ', 64));
assert(size(funs.rx, 1) == 2 && size(funs.cx, 1) == 2, 'expected a 2D tree');
F = ifmm(funs.Afun, funs.rx, funs.cx, funs.occ, eps, funs.pxyfun, flam_opts('ifmm'));
u = surferkernevalmat(S, K, struct('r', [targs(:,:); zeros(1, 400)]), eps)*sig;
err = norm(ifmm_mv(F, funs.xinterp*sig, funs.Afun) + funs.Q*sig - u)/norm(u);
assert(err < tol, 'ifmm 2D: err %.2e > tol %.2e', err, tol);

end


function check_far(S, K, p, l, eps)
% Proxy blocks from both sides of a box reproduce the far field of K.
dim = numel(l);
far = far_targets(l, 48);
x = S.r(1:dim,:);
C = byindex.proxyfun(S, K, [], p, 'c', x, x, 1:S.npts, [], l, zeros(dim, 1));
R = byindex.proxyfun(S, K, [], p, 'r', x, x, 1:S.npts, [], l, zeros(dim, 1));
check_id(K.eval(S, far).*S.wts.', C, 10*eps);
check_id(K.eval(far, S).', R.', 10*eps);
end


function check_id(A, B, tol)
% An ID of the columns of B reproduces the columns of A.
[~, R, p] = qr(B, 0);
d = abs(diag(R));
k = sum(d > tol/10*d(1));
assert(k < size(B, 2), 'test must compress some columns');
T = R(1:k,1:k)\R(1:k,k+1:end);
E = A(:,p(k+1:end)) - A(:,p(1:k))*T;
assert(norm(E, 'fro') < tol*norm(A, 'fro'), 'proxy ID does not reproduce the far field');
end


function S = box_sources(l, n)
% n sources in the box with sides l, with frames and weights.
dim = numel(l);
S = random_frames([l(:).*(rand(dim, n) - 0.5); zeros(3-dim, n)], dim);
S.npts = n;
S.wts = (0.5 + rand(n, 1))/n;
end


function T = far_targets(l, n)
% n targets between 2 and 3 box sizes from the center, with frames.
dim = numel(l);
v = randn(dim, n);  v = v./vecnorm(v);
T = random_frames([l(:).*(2 + rand(1, n)).*v; zeros(3-dim, n)], dim);
end


function S = random_frames(r, dim)
n = size(r, 2);
nv = [randn(dim, n); zeros(3-dim, n)];  nv = nv./vecnorm(nv);
du = cross(nv, randn(3, n));  du = du./vecnorm(du);
if dim == 2, du = cross(nv, repmat([0; 0; 1], 1, n)); end
S = struct('r', r, 'n', nv, 'du', du, 'dv', cross(nv, du), ...
    'd', du, 'mean_curv', randn(1, n));
end


function M = cartesian(s)
% Tangential (du, dv) currents and charge to Cartesian currents and charge
n = size(s.r(:,:), 2);
ru = s.du(:,:)./vecnorm(s.du(:,:));
rv = cross(s.n(:,:), ru);
M = zeros(4, 3, n);
M(1:3,1,:) = reshape(ru, 3, 1, n);
M(1:3,2,:) = reshape(rv, 3, 1, n);
M(4,3,:) = 1;
end


function opts = flam_opts(routine)
opts = struct('proxybylevel', contains(fileread(which(routine)), 'proxybylevel'));
end


function tf = has_flam()
tf = exist('rskelf', 'file') == 2 && exist('ifmm', 'file') == 2;
end
