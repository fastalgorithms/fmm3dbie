function obj = radgreen(green, type, varargin)
%KERNEL3D.RADGREEN   Kernel from a radial Green's function, with exact
%   eval and near quadrature from a fit built at getquad time.
%
%   K = kernel3d(@green, type, ...)
%
%   GREEN has the helm3d.green format, or returns an isotropic tensor
%   (val (nt,ns,3,3), grad (nt,ns,3,3,3)) for the types 's', 'efie',
%   'curl' and 'mfie'. It is checked to be radial (or isotropic) and its
%   derivatives against finite differences.
%
%   TYPE and options are as in KERNEL3D.RADPP, including 'fmm' for an FMM
%   of GREEN (not for tensor greens). 'ipv' overrides the
%   automatic choice; kernel_order is max(ipv)-1. Kernels on the same
%   GREEN (up to a factor) combine under +, -, scalar * and interleave into
%   one kernel, and kernels of one built-in type add as that type on the
%   sum of their greens.
%
%   Example:
%      g = @(s,t) helm3d.green(zk, s, t);
%      K = kernel3d(g, 'cp', [1 1i]);
%
% See also KERNEL3D.RADPP

if ischar(green)
    obj = combine(green, type, varargin{:});
    return
end
if nargin < 2 || isempty(type), type = 's'; end

% a positional coefs or nu goes to radpp as is
npos = find(cellfun(@ischar, varargin), 1) - 1;
if isempty(npos), npos = numel(varargin); end
pos = varargin(1:npos);  varargin(1:npos) = [];

% the singularity of green, detected or given
sing = [];
k = find(strcmpi(varargin(1:2:end), 'singularity'));
if ~isempty(k)
    sing = varargin{2*k};
    varargin(2*k-1:2*k) = [];
    if ischar(sing) && strcmpi(sing, 'log'), sing = 0; end
    if ~isnumeric(sing) || ~isscalar(sing)
        error('KERNEL3D.RADGREEN: ''singularity'' must be a power or ''log''.');
    end
end
istensor = numel(green(zeros(3,1), [1;0;0])) == 9;
if istensor
    prof  = @(r, nd) tprofiles(green, r, nd);
    check = @(nd, name) tcheck(green, nd, name);
    k = find(strcmpi(varargin(1:2:end), 'dim'));
    if ~isempty(k) && varargin{2*k} ~= 3
        error('KERNEL3D.RADGREEN: a tensor green needs dim = 3.');
    end
    if ~isempty(k), varargin(2*k-1:2*k) = []; end
    varargin = [varargin, {'dim', 13}];
else
    prof  = @(r, nd) gprofiles(green, r, nd);
    check = @(nd, name) gcheck(green, nd, name);
end
if isempty(sing), sing = gsing(@(r) sum(abs(prof(r, 0)), 2)); end

part = struct('profiles', prof, 'check', check, 'handle', green, ...
              'sing', sing, 'tensor', istensor);
obj = kernel3d.radpp(mkgreen({part}, 1), [], type, pos{:}, varargin{:});
% built-in types: kernels of one type on different greens add as that
% type on the sum of the greens
obj.params.key = [];
if ischar(type)
    p = obj.params;
    obj.params.key = {p.spec.name, p.kcoefs, p.nu, p.dim};
end

end


function G = mkgreen(parts, c)
% the green sum_k c(k)*parts{k}
keep = c ~= 0;
if any(keep), parts = parts(keep);  c = c(keep); end
q = parts{1};
G = struct('profiles', q.profiles, 'check', q.check, 'sing', q.sing, ...
           'tensor', q.tensor, 'parts', {parts}, 'c', c);
if numel(parts) == 1
    if c ~= 1, G.profiles = @(r, nd) c*q.profiles(r, nd); end
    return
end
G.profiles = @(r, nd) sumprof(parts, c, r, nd);
G.check = @(nd, name) checkall(parts, nd, name);
try
    G.sing = gsing(@(r) sum(abs(G.profiles(r, 0)), 2));
catch
    G.sing = max(cellfun(@(q) q.sing, parts));
end
end


function f = sumprof(parts, c, r, nd)
f = 0;
for k = 1:numel(parts), f = f + c(k)*parts{k}.profiles(r, nd); end
end


function checkall(parts, nd, name)
for k = 1:numel(parts), parts{k}.check(nd, name); end
end


function f = gprofiles(green, r, nd)
% on-axis profiles [G | G_x | G_xx G_yy | G_xxx G_xyy]
r = r(:);
out = cell(1, nd+1);
[out{:}] = green(zeros(3,1), [r.'; zeros(2, numel(r))]);
f = out{1}(:);
if nd >= 1, f = [f, out{2}(:,1,1)]; end
if nd >= 2, f = [f, out{3}(:,1,1,1), out{3}(:,1,2,2)]; end
if nd >= 3, f = [f, out{4}(:,1,1,1,1), out{4}(:,1,1,2,2)]; end
end


function p = gsing(v)
% power p of the singularity G ~ r^-p at 0, from the slopes of log|G| at
% small r. log r, bounded and smooth greens give p = 0 (weakly singular):
% the slope of log|log r| tends to 0 only like 1/log r, hence the 0.25
% margin, which still counts r^-1/2 as p = 1
r = 10.^-(4:0.5:10).';
v = v(r);
sl = diff(log(abs(v) + realmin))./diff(log(r));
sl = sl(end-3:end);
if any(~isfinite(sl)) || max(sl) - min(sl) > 0.15
    error(['KERNEL3D.RADGREEN: cannot determine the singularity of ' ...
           'green; pass ''singularity'', p.']);
end
p = max(0, ceil(-median(sl) - 0.25));
end


function f = tprofiles(green, r, nd)
% on-axis profiles [A B | A' B'] of G_ij = A d_ij + B u_i u_j
r = r(:);
out = cell(1, nd+1);
[out{:}] = green(zeros(3,1), [r.'; zeros(2, numel(r))]);
A = out{1}(:,1,2,2);
f = [A, out{1}(:,1,1,1) - A];
if nd >= 1
    Ap = out{2}(:,1,2,2,1);
    f = [f, Ap, out{2}(:,1,1,1,1) - Ap];
end
end


function tcheck(green, nd, name)
% outputs and sizes, tensor transformation under rotation, and grad
% against finite differences of val
names = {'val', 'grad'};
need  = strjoin(names(1:nd+1), ', ');
err0  = 'KERNEL3D.RADGREEN: ';

s = randn(3, 2);  t = randn(3, 4) + [2; 0; 0];
out = cell(1, nd+1);
try
    [out{:}] = green(s, t);
catch ME
    error([err0 'type ''%s'' needs [%s] from green: %s'], name, need, ...
          ME.message);
end
for k = 1:nd+1
    sz = [4 2 3*ones(1, k+1)];
    if ~isequal(size(out{k}), sz)
        error([err0 '%s of green should be %s, not %s.'], names{k}, ...
              mat2str(sz), mat2str(size(out{k})));
    end
end

% G(Q x, Q y) = Q G(x, y) Q^T, and invariance under a shift
y = [0.2; -0.1; 0.4];  x = y + [0.5; 0.6; -0.3];
[Q, ~] = qr(randn(3));
G1 = reshape(green(y, x), 3, 3);
G2 = reshape(green(Q*y + 0.7, Q*x + 0.7), 3, 3);
if norm(Q*G1*Q.' - G2, 'fro') > 1e-10*norm(G1, 'fro')
    error([err0 'green is not an isotropic tensor.']);
end

% grad against finite differences of val
if nd >= 1
    h = 1e-5;
    o = cell(1, 2);
    [o{:}] = green(y, x);
    Ga = reshape(o{2}, 9, 3);
    F = zeros(9, 3);
    for k = 1:3
        e = zeros(3, 1);  e(k) = h;
        F(:,k) = reshape(green(y, x + e) - green(y, x - e), 9, 1)/(2*h);
    end
    rel = norm(Ga - F, 'fro')/norm(Ga, 'fro');
    if rel > 1e-6
        error([err0 'grad of green does not match finite differences ' ...
               'of val (%.1e).'], rel);
    end
end
end


function gcheck(green, nd, name)
% enough outputs, of the right sizes, radial, and consistent derivatives
names = {'val', 'grad', 'hess', 'third'};
need  = strjoin(names(1:nd+1), ', ');
err0  = 'KERNEL3D.RADGREEN: ';

s = randn(3, 2);  t = randn(3, 4) + [2; 0; 0];
out = cell(1, nd+1);
try
    [out{:}] = green(s, t);
catch ME
    error([err0 'type ''%s'' needs [%s] from green: %s'], name, need, ...
          ME.message);
end
for k = 1:nd+1
    sz = [4 2 3*ones(1, k-1)];
    if ~isequal(size(out{k}), sz)
        error([err0 '%s of green should be %s, not %s.'], names{k}, ...
              mat2str(sz), mat2str(size(out{k})));
    end
end

% radial: same r in several directions and after a shift
y = [0.2; -0.1; 0.4];  r0 = 0.8;
U = randn(3, 8);  U = U./vecnorm(U);
o1 = cell(1, min(nd, 1) + 1);  o2 = o1;
[o1{:}] = green(y, y + r0*U);
[o2{:}] = green(y + 1.3, y + 1.3 + r0*U);
v = o1{1};
sc = max(abs(v));
if (max(abs(v - v(1))) > 1e-10*sc) || (max(abs(o2{1} - v)) > 1e-10*sc)
    error([err0 'green is not radial.']);
end
if nd >= 1
    g = reshape(o1{2}, 8, 3);
    gr = sum(g.*U.', 2);
    if norm(g - gr.*U.', 'fro') > 1e-10*norm(g, 'fro')
        error([err0 'grad of green is not radial.']);
    end
end

% each derivative against finite differences of the one below
x = y + 0.9*U(:,1);  h = 1e-5;
for k = 1:nd
    ok = cell(1, k+1);
    [ok{:}] = green(y, x);
    A = reshape(ok{k+1}, [], 3);
    F = zeros(size(A));
    for j = 1:3
        e = zeros(3, 1);  e(j) = h;
        op = cell(1, k);  om = op;
        [op{:}] = green(y, x + e);
        [om{:}] = green(y, x - e);
        F(:, j) = (op{k}(:) - om{k}(:))/(2*h);
    end
    rel = norm(A - F, 'fro')/max(norm(A, 'fro'), realmin);
    if rel > 1e-6
        error([err0 '%s of green does not match finite differences ' ...
               'of %s (%.1e).'], names{k+1}, names{k}, rel);
    end
end
end


%%% arithmetic on green kernels

function out = combine(op, varargin)
% kernels may be given as a kernel3d array or a cell array
switch op
    case 'compatible'
        out = compatible(varargin{1});
    case 'summable'
        Ks = tocell(varargin{1});
        out = compatible(Ks) || samekey(Ks);
    case 'sum'      % sum(K1, a, K2, b) = a*K1 + b*K2
        if samekey(varargin([1 3]))
            out = keysum(varargin{:});
        else
            out = lincomb(varargin{:});
        end
    case 'scale'    % scale(K, a)
        if ~isempty(varargin{1}.params.key)
            out = keysum(varargin{1}, varargin{2});
        else
            out = lincomb(varargin{1}, varargin{2});
        end
    case 'interleave'
        out = blocks(varargin{1});
    otherwise
        error('KERNEL3D.RADGREEN: unknown operation ''%s''.', op);
end
end


function tf = isgreen(Ks)
tf = all(cellfun(@(K) strcmp(K.name, 'radgreen') && isempty(K.diag), Ks(:)));
end


function tf = compatible(Ks)
% green kernels with no diag term, the same dim and proportional greens
Ks = tocell(Ks);
tf = isgreen(Ks);
for k = 2:numel(Ks)
    if ~tf, return; end
    tf = Ks{k}.params.dim == Ks{1}.params.dim && ...
         ~isempty(ratio(Ks{1}.params.green, Ks{k}.params.green));
end
end


function tf = samekey(Ks)
% green kernels of one built-in type
tf = isgreen(Ks) && ~isempty(Ks{1}.params.key) && ...
     all(cellfun(@(K) isequal(K.params.key, Ks{1}.params.key), Ks(:)));
end


function tf = samegreen(g1, g2)
% equal handles, or anonymous functions with the same text and the same
% captured values
tf = isequal(g1, g2);
if tf, return; end
try
    f1 = functions(g1);  f2 = functions(g2);
    tf = strcmp(f1.type, 'anonymous') && strcmp(f2.type, 'anonymous') && ...
         strcmp(f1.function, f2.function) && isequal(f1.workspace, f2.workspace);
catch
    tf = false;
end
end


function lam = ratio(G1, G2)
% lam with G2 = lam*G1, or [] if the greens are not proportional
lam = [];
if numel(G1.parts) ~= numel(G2.parts), return; end
c1 = G1.c(:);  c2 = zeros(size(c1));
for k = 1:numel(G2.parts)
    i = find(cellfun(@(q) samegreen(q.handle, G2.parts{k}.handle), G1.parts), 1);
    if isempty(i), return; end
    c2(i) = G2.c(k);
end
l = c2(1)/c1(1);
if norm(c2 - l*c1) <= 1e-14*norm(c2), lam = l; end
end


function K = keysum(K1, a, K2, b)
% a*K1 (+ b*K2) of one built-in type: that type on a*G1 (+ b*G2)
Ks = {K1};  cs = a;
if nargin > 2, Ks = {K1, K2};  cs = [a b]; end
parts = {};  c = [];
for k = 1:numel(Ks)
    G = Ks{k}.params.green;
    for j = 1:numel(G.parts)
        i = find(cellfun(@(q) samegreen(q.handle, G.parts{j}.handle), parts), 1);
        if isempty(i)
            parts{end+1} = G.parts{j};  c(end+1) = 0;  
            i = numel(parts);
        end
        c(i) = c(i) + cs(k)*G.c(j);
    end
end
p = K1.params;
ipv = p.ipv;                   % the worse of the two, as for other sums
for k = 2:numel(Ks), ipv = max(ipv, Ks{k}.params.ipv); end
ipv = {'ipv', ipv};
% the fmm of the sum of the greens
fmm = {};
if all(cellfun(@(Kk) ~isempty(Kk.params.fmm), Ks))
    fs = cellfun(@(Kk) Kk.params.fmm, Ks, 'UniformOutput', false);
    fmm = {'fmm', radpp_sum(cs, fs{:})};
end
K = kernel3d.radpp(mkgreen(parts, c), [], p.spec, 'dim', p.dim, ipv{:}, ...
        'coefs', p.kcoefs, 'nu', p.nu, 'opts', p.fitopts, ...
        'zk', max(cellfun(@(K) K.zk, Ks)), fmm{:}, ...
        'ifcomplex', any(cellfun(@(K) logical(K.ifcomplex), Ks)) || ~isreal(cs));
K.params.key = p.key;
% otherwise the sum of the fmms of the kernels, if they have them
if isempty(fmm) && all(cellfun(@(Kk) isa(Kk.fmm, 'function_handle'), Ks))
    fs = cellfun(@(Kk) Kk.fmm, Ks, 'UniformOutput', false);
    K.fmm = radpp_sum(cs, fs{:});
end
end


function K = lincomb(K1, a, K2, b)
% a*K1 (+ b*K2), weights padded to the higher derivative order
Ks = {K1};  cs = a;
if nargin > 2, Ks = {K1, K2};  cs = [a b]; end
if numel(Ks) == 2 && ~isequal(K1.opdims, K2.opdims)
    error('KERNEL3D.RADGREEN: kernel dimensions must agree.');
end
p  = Ks{1}.params;
nd = max(cellfun(@(K) K.params.nderiv, Ks));
nc = ncompof(nd, p.dim);
sp = specs(Ks);
cw = cs;                       % proportional greens: G_k = lam*G_1
for k = 2:numel(Ks)
    cw(k) = cs(k)*ratio(p.green, Ks{k}.params.green);
end

spec = p.spec;
spec.name = 'custom';
spec.nderiv = nd;
spec.wfun = @(t) wsum(t, sp, cw, nc, 'wfun');
spec.vfun = [];
if any(arrayfun(@(q) ~isempty(q.vfun), sp))
    spec.vfun = @(t) wsum(t, sp, cw, nc, 'vfun');
end
spec.targ_fields = unique([sp.targ_fields]);

ipv = Ks{1}.params.ipv;
for k = 2:numel(Ks), ipv = max(ipv, Ks{k}.params.ipv); end
K = rebuild(Ks, spec, ipv, cw);

% fmms are combined as in the generic arithmetic
if all(cellfun(@(Kk) isa(Kk.fmm, 'function_handle'), Ks))
    fs = cellfun(@(Kk) Kk.fmm, Ks, 'UniformOutput', false);
    K.fmm = @(eps, src, targ, sigma) fmmsum(fs, cs, eps, src, targ, sigma);
end
end


function pot = fmmsum(fs, cs, eps, src, targ, sigma)
pot = 0;
for k = 1:numel(fs)
    pot = pot + cs(k)*fs{k}(eps, src, targ, sigma);
end
end


function W = wsum(t, sp, cs, nc, f)
W = 0;
for k = 1:numel(sp)
    if ~isempty(sp(k).(f))
        W = W + cs(k)*pad(sp(k).(f)(t), nc);
    end
end
if isequal(W, 0)
    W = zeros(size(t.r,2), nc, size(sp(1).ent, 1));
end
end


function K = blocks(Ks)
% interleave a block array of green kernels into one kernel
Ks = tocell(Ks);
[m, n] = size(Ks);
rd = cellfun(@(K) K.opdims(1), Ks(:,1)).';
cd = cellfun(@(K) K.opdims(2), Ks(1,:));
r0 = [0 cumsum(rd)];  c0 = [0 cumsum(cd)];
opdims = [r0(end) c0(end)];

p  = Ks{1}.params;
nd = max(cellfun(@(K) K.params.nderiv, Ks(:)));
nc = ncompof(nd, p.dim);

% entries of the result in column-major order: block, local entry
[I, J] = ndgrid(1:opdims(1), 1:opdims(2));
ne = numel(I);
map = zeros(ne, 2);
ipv = zeros(1, ne);
for e = 1:ne
    k = find(I(e) > r0, 1, 'last');  l = find(J(e) > c0, 1, 'last');
    q = Ks{k,l}.params;
    le = find(q.ent(:,1) == I(e) - r0(k) & q.ent(:,2) == J(e) - c0(l));
    map(e,:) = [sub2ind([m n], k, l), le];
    ipv(e) = q.ipv(le);
end

sp = specs(Ks(:));
lam = cellfun(@(K) ratio(p.green, K.params.green), Ks(:));
spec = p.spec;
spec.name = 'custom';
spec.nderiv = nd;
spec.opdims = opdims;
spec.ent = [I(:) J(:)];
spec.wfun = @(t) wblocks(t, sp, lam, map, nc, 'wfun');
spec.vfun = [];
if any(arrayfun(@(q) ~isempty(q.vfun), sp))
    spec.vfun = @(t) wblocks(t, sp, lam, map, nc, 'vfun');
end
spec.targ_fields = unique([sp.targ_fields]);
K = rebuild(Ks(:).', spec, ipv, lam);

% block fmm, as in the generic interleave
if all(cellfun(@(Kk) isa(Kk.fmm, 'function_handle'), Ks(:)))
    K.fmm = @(eps, src, targ, sigma) fmmblocks(Ks, r0, c0, opdims, eps, ...
                                               src, targ, sigma);
end
end


function pot = fmmblocks(Ks, r0, c0, opdims, eps, src, targ, sigma)
% sigma is (opdims(2)*ns, ndens), one column per density
[m, n] = size(Ks);
ns = size(src.r(:,:), 2);
ndens = numel(sigma)/(opdims(2)*ns);
nt = size(targ.r(:,:), 2);
sig = reshape(sigma, opdims(2), ns, ndens);
pot = zeros(opdims(1), nt, ndens);
for k = 1:m
    for l = 1:n
        sl = sig(c0(l)+1:c0(l+1), :, :);
        pkl = Ks{k,l}.fmm(eps, src, targ, reshape(sl, [], ndens));
        pot(r0(k)+1:r0(k+1), :, :) = pot(r0(k)+1:r0(k+1), :, :) + ...
            reshape(pkl, r0(k+1) - r0(k), nt, ndens);
    end
end
pot = reshape(pot, [], ndens);
end


function W = wblocks(t, sp, lam, map, nc, f)
nt = size(t.r, 2);
W = zeros(nt, nc, size(map, 1));
Wb = cell(numel(sp), 1);
for b = unique(map(:,1)).'
    if ~isempty(sp(b).(f)), Wb{b} = lam(b)*pad(sp(b).(f)(t), nc); end
end
for e = 1:size(map, 1)
    if ~isempty(Wb{map(e,1)})
        W(:,:,e) = Wb{map(e,1)}(:,:,map(e,2));
    end
end
end


function K = rebuild(Ks, spec, ipv, cs)
p = Ks{1}.params;
K = kernel3d.radpp(p.green, [], spec, 'ipv', ipv, 'dim', p.dim, ...
        'opts', p.fitopts, ...
        'zk', max(cellfun(@(K) K.zk, Ks)), ...
        'ifcomplex', any(cellfun(@(K) logical(K.ifcomplex), Ks)) || ~isreal(cs));
K.params.key = [];
end


function c = tocell(Ks)
if iscell(Ks), c = Ks; return; end
c = cell(size(Ks));
for k = 1:numel(Ks), c{k} = Ks(k); end
end


function sp = specs(Ks)
sp = Ks{1}.params.spec;
for k = 2:numel(Ks), sp(k) = Ks{k}.params.spec; end
end


function W = pad(W, nc)
W = cat(2, W, zeros(size(W,1), nc - size(W,2), size(W,3)));
end


function nc = ncompof(nd, dim)
if dim == 3, nc = [1 4 10 20]; else, nc = [1 3 6 10]; end
nc = nc(nd+1);
end
