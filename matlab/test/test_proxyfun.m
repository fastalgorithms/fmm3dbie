% Test byindex.proxy_points and byindex.proxyfun.
%
% The proxy functions are called with the arguments FLAM would pass, and
% the skeletons they give are checked against the interactions they
% stand in for. FLAM is not needed, but if rskelf and ifmm are on the
% path, small rskelf and ifmm examples are also run.

run ../startup.m

rng(42);

%% Now run the tests

test_proxy_points_calling_sequences();
test_proxy_points_rejects();
test_rskelf_addtrans();
test_ifmm_mixed_opdims();
test_ifmm_offsurface_groups();
test_proxybylevel();
test_flam_ifmm_mixed_opdims();


function test_proxy_points_calling_sequences()
% All calling sequences, and the shape of the output.

eps  = 1e-6;
kern = kernel3d('l', 's');

p1 = byindex.proxy_points(kern, 1, eps);
p2 = byindex.proxy_points(kern, 1, eps, struct('nshell_max', 1, 'ifaddtrans', true));

assert(p1.nshell == 1 && p2.nshell == 1, 'Laplace should need one shell');
assert(p2.ifaddtrans && ~p1.ifaddtrans, 'ifaddtrans not stored');
npxy = p1.nproxy*p1.nshell;
assert(isequal(size(p1.r), [3 npxy]) && isequal(size(p1.n), [3 npxy]) && ...
    numel(p1.w) == npxy, 'wrong proxy sizes');
assert(abs(p1.radius - 1.5) < 1e-12, 'default radius should be 1.5 box sizes');
assert(abs(sum(p1.w) - 4*pi*1.5^2) < 1e-10, 'weights should sum to the sphere area');

% Tighter tolerances need more proxy points
p5 = byindex.proxy_points(kern, 1, 1e-3);
assert(p5.nproxy < p1.nproxy, 'nproxy should grow as eps decreases');

end


function test_proxy_points_rejects()
% Kernels needing fields other than n at the proxies, and proxy surfaces
% that would contain non-neighbor boxes, are rejected.

caught = false;
try
    byindex.proxy_points(kernel3d('bel', 'rlb'), 1, 1e-6);
catch
    caught = true;
end
assert(caught, 'mean_curv at proxy targets should raise an error');

caught = false;
try
    byindex.proxy_points(kernel3d('l', 's'), 1, 1e-6, struct('crad', 2));
catch
    caught = true;
end
assert(caught, 'crad > sqrt(3) should raise an error');

% Missing eps: a warning within a factor of ten, an error otherwise.
% Too few proxy points first.
kern = kernel3d('l', 's');
wstate = warning('off', 'BYINDEX:proxy_points:nproxy');
cleanup = onCleanup(@() warning(wstate));
lastwarn('');
p = byindex.proxy_points(kern, 1, 1e-6, struct('nproxy_max', 500));
[~, id] = lastwarn;
assert(strcmp(id, 'BYINDEX:proxy_points:nproxy') && p.nproxy <= 500, ...
    'expected a warning and nproxy <= nproxy_max');

caught = false;
try
    byindex.proxy_points(kern, 1, 1e-6, struct('nproxy_max', 400));
catch
    caught = true;
end
assert(caught, 'missing eps by more than a factor of ten should raise an error');

% Too few shells: log|x-y| needs more than one shell
kern = kernel3d(@logg, 's');
p = byindex.proxy_points(kern, 1, 1e-5);
assert(p.nshell > 1, 'log kernel should need more than one shell');

warning('off', 'BYINDEX:proxy_points:nshell');
lastwarn('');
p = byindex.proxy_points(kern, 1, 1e-5, struct('nshell_max', 2));
[~, id] = lastwarn;
assert(strcmp(id, 'BYINDEX:proxy_points:nshell') && p.nshell == 2, ...
    'expected a warning and nshell = nshell_max');

caught = false;
try
    byindex.proxy_points(kern, 1, 1e-5, struct('nshell_max', 1));
catch
    caught = true;
end
assert(caught, 'one shell for the log kernel should raise an error');

end


function test_rskelf_addtrans()
% rskelf calling sequence with a non-symmetric kernel and
% proxy.ifaddtrans: the skeleton compresses both the columns and the rows
% of slf, and nbr keeps only the points inside the proxy surface.

eps = 1e-6;  tol = 1e-5;

S    = geometries.sphere(2, 4, [0;0;0], 4, 1);
kern = kernel3d('l', 'd');
ctr  = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kern, l, eps, struct('ifaddtrans', true));

dist    = vecnorm(S.r - ctr);
srcids  = find(dist < 0.5);
targids = find(dist > 2.0);
nbrin   = find(dist > 0.5 & dist < 1.4);

src2targ = byindex.kernbyindex(targids, srcids, S, kern);
targ2src = byindex.kernbyindex(srcids, targids, S, kern);

[Kpxy, nbr] = byindex.proxyfun(S, kern, [], proxy, S.r, srcids, ...
    [nbrin, targids], l, ctr);
assert(isequal(sort(nbr), sort(nbrin)), 'nbr filtered incorrectly');
npxy = size(proxy.r, 2);
assert(size(Kpxy, 1) == 2*npxy, 'ifaddtrans should give %d rows, got %d', ...
    2*npxy, size(Kpxy, 1));

check_lr('rskelf addtrans, columns', src2targ,    Kpxy, tol);
check_lr('rskelf addtrans, rows',    targ2src.', Kpxy, tol);

end


function test_ifmm_mixed_opdims()
% ifmm calling sequence with mixed opdims, so rows and columns have
% different points: rc = 'c' and rc = 'r' on surfer 1.

eps = 1e-6;  tol = 1e-5;

S1 = geometries.sphere(2, 4, [0;0;0], 4, 1);
S2 = geometries.sphere(2, 4, [8;0;0], 4, 1);
srfrs = [S1, S2];
kerns = mixed_kerns();

ctr = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kerns, l, eps);
rx = [repelem(S1.r, 1, 3), S2.r];
cx = [S1.r, S2.r];

near = find(vecnorm(S1.r - ctr) < 0.5);
far  = find(vecnorm(S1.r - ctr) > 2.0);

% columns: sources near ctr, targets far from it
cnear = near;
rfar  = [reshape(3*far + (-2:0).', 1, []), 3*S1.npts + (1:S2.npts)];
src2targ = byindex.kernbyindex(rfar, cnear, srfrs, kerns);
[Kpxy, nbr] = byindex.proxyfun(srfrs, kerns, [], proxy, 'c', rx, cx, cnear, rfar, l, ctr);
assert(isempty(nbr), 'far points should be removed from nbr');
assert(size(Kpxy, 1) == 4*size(proxy.r, 2), 'wrong number of proxy rows');
check_lr('ifmm mixed opdims, rc = c', src2targ, Kpxy, tol);

% rows: targets near ctr, sources far from it
rnear = reshape(3*near + (-2:0).', 1, []);
cfar  = [far, S1.npts + (1:S2.npts)];
targ2src = byindex.kernbyindex(rnear, cfar, srfrs, kerns);
Kpxy = byindex.proxyfun(srfrs, kerns, [], proxy, 'r', rx, cx, rnear, cfar, l, ctr);
assert(size(Kpxy, 2) == 2*size(proxy.r, 2), 'wrong number of proxy columns');
check_lr('ifmm mixed opdims, rc = r', targ2src.', Kpxy.', tol);

end


function test_ifmm_offsurface_groups()
% Two off-surface target groups, the second needing normals, with a 2x1
% kernel array.

eps = 1e-6;  tol = 1e-5;

S    = geometries.sphere(2, 4, [0;0;0], 4, 1);
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
check_lr('offsurface groups, rc = c', src2targ, Kpxy, tol);

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
check_lr('offsurface groups, rc = r', targ2src.', Kpxy.', tol);

end


function test_proxybylevel()
% With one argument, proxyfun returns the proxy function for that level,
% rebuilt for a box of size width/2^(lvl-1).

eps = 1e-6;  tol = 1e-5;

S    = geometries.sphere(2, 4, [0;0;0], 4, 1);
kern = kernel3d('h', 's', 4);
width = 4;
proxy = byindex.proxy_points(kern, width, eps);
pxyfun = @(varargin) byindex.proxyfun(S, kern, [], proxy, varargin{:});

lvl = 3;
pxyfunlvl = pxyfun(lvl);
assert(isa(pxyfunlvl, 'function_handle'), 'expected a function handle for lvl');

plvl = byindex.proxy_points(kern, width/2^(lvl-1), eps);
assert(plvl.nproxy < proxy.nproxy, 'smaller boxes should need fewer proxies');

ctr = [2; 0; 0];  l = width/2^(lvl-1);
dist    = vecnorm(S.r - ctr);
srcids  = find(dist < l/2);
targids = find(dist > 1.6*l);

src2targ = byindex.kernbyindex(targids, srcids, S, kern);
Kpxy = pxyfunlvl(S.r, srcids, targids, l, ctr);
assert(size(Kpxy, 1) == size(plvl.r, 2), 'level proxy function should use plvl');
check_lr('proxybylevel', src2targ, Kpxy, tol);

% Levels up to proxy.lvl keep the surface instead of rebuilding it
proxy.lvl = 2;
pxyfun = @(varargin) byindex.proxyfun(S, kern, [], proxy, varargin{:});
Kpxy = feval(pxyfun(2), S.r, srcids, targids, l, ctr);
assert(size(Kpxy, 1) == size(proxy.r, 2), 'level 2 should keep proxy');

end


function test_flam_ifmm_mixed_opdims()
% ifmm from byindex.flamfuns with mixed opdims and oversampled sources:
% each point appears once per component and the matrix is not square.

if ~has_flam(), return; end

eps = 1e-4;  tol = 1e-3;

S1 = geometries.sphere(1, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(1, 2, [3;0;0], 3, 1);
srfrs = [S1, S2];
kerns = mixed_kerns();

funs = byindex.flamfuns(srfrs, kerns, eps);
assert(size(funs.xinterp, 1) > size(funs.xinterp, 2), 'expected oversampling');
F = ifmm(funs.Afun, funs.rx, funs.cx, 32, eps, funs.pxyfun);

Afull = surfermat(srfrs, kerns, eps, struct('unif_nover', true));
sig = randn(size(funs.Q, 2), 1);
u   = Afull*sig;
err = norm(ifmm_mv(F, funs.xinterp*sig, funs.Afun) + funs.Q*sig - u)/norm(u);
assert(err < tol, 'ifmm mixed opdims: err %.2e > tol %.2e', err, tol);

end


function tf = has_flam()
% Helper: true if FLAM's rskelf and ifmm are on the path.
tf = exist('rskelf', 'file') == 2 && exist('ifmm', 'file') == 2;
if ~tf
    st = dbstack(1);
    fprintf('FLAM not found, skipping %s\n', st(1).name);
end
end


function kerns = mixed_kerns()
% Helper: 2x2 kernel array with opdims [3 1], [3 1]; [1 1], [1 1]
kerns(2,2) = kernel3d();
kerns(1,1) = grad(kernel3d('l', 's'));
kerns(1,2) = grad(kernel3d('l', 'd'));
kerns(2,1) = kernel3d('l', 'sp');
kerns(2,2) = kernel3d('l', 's');
end


function varargout = logg(s, t)
% Helper: G = log|x-y| with target derivatives, helm3d.green format
x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);  z = t(3,:).' - s(3,:);
r2 = x.^2 + y.^2 + z.^2;
X = cat(3, x, y, z);
varargout{1} = log(r2)/2;
varargout{2} = X./r2;
H = zeros([size(x) 3 3]);
for i = 1:3, for j = 1:3
    H(:,:,i,j) = (i==j)./r2 - 2*X(:,:,i).*X(:,:,j)./r2.^2;
end, end
varargout{3} = H;
end


function check_lr(label, src2targ, Kpxy, tol)
% Helper: verify the interpolative decomposition of Kpxy (column-pivoted
% QR, as in FLAM) also compresses src2targ.
[~, R, jpvt] = qr(Kpxy, 0);
d  = abs(diag(R));
k  = max(sum(d > tol/10 * d(1)), 1);
T  = R(1:k, 1:k) \ R(1:k, k+1:end);
sk = jpvt(1:k);
rd = jpvt(k+1:end);
err = norm(src2targ(:, rd) - src2targ(:, sk)*T, 'fro') / norm(src2targ, 'fro');
assert(err < tol, '%s: err %.2e > tol %.2e', label, err, tol);
end
