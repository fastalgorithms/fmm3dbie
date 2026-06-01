% Test surfermat, surfermatapply, surferkernevalmat, and surferkerneval.

%% Now run the tests

test_apply_paths();
test_eval_paths();
test_objover_modes();
test_block_kernels();


function test_apply_paths()
% surfermat/surfermatapply: dense, corrections, nonsmoothonly paths and rfac passthrough.

tol_apply = 1e-6;

S     = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
zk    = 1;
kerns = kernel3d('h','c',zk,[1,1]);
srfrs = S;
eps   = 1e-10;

tic
Smat = surfermat(srfrs, kerns, eps);
Smat = Smat - 0.5*eye(size(Smat));

[syscors, novers] = surfermat(srfrs, kerns, eps, struct('corrections',  1));
syscors  = syscors  - 0.5*speye(size(syscors));
toc

skern = kernel3d('h','s',zk);
src   = []; src.r = [2;0;1];
rhs   = -skern.eval(src, merge(srfrs));

sysapply = @(d) surfermatapply(srfrs, kerns, d, eps, novers, syscors);

assert(norm(sysapply(rhs) - Smat*rhs) / norm(Smat*rhs) < tol_apply, 'corrections apply error');

end


function test_eval_paths()
% surferkernevalmat/surferkerneval: smooth, corrections, nonsmoothonly,
% dense paths and rfac passthrough.

tol_eval = 1e-6;

S     = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
zk    = 1;
kerns = kernel3d('h','c',zk,[1,1]);
srfrs = S;
eps   = 1e-10;

kernseval = kerns(1,:);

skern = kernel3d('h','s',zk);
src   = []; src.r = [2;0;1];
rhs   = -skern.eval(src, merge(srfrs));

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; zeros(1,nplot^2)];

tic
uscat_smth  = surferkerneval(srfrs, kernseval, rhs, targs, eps);

opts_cor = []; opts_cor.corrections = 1;
[evalcors, novers_eval] = surferkernevalmat(srfrs, kernseval, targs, eps, opts_cor);
opts_cor.corrections = evalcors;  opts_cor.objover = novers_eval;
uscat_cor   = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_cor);

evalmat     = surferkernevalmat(srfrs, kernseval, targs, eps);
uscat_dense = evalmat * rhs;
toc

uref = uscat_smth;
assert(norm(uref - uscat_cor)   / norm(uref) < tol_eval, 'corrections eval error');
assert(norm(uref - uscat_dense) / norm(uref) < tol_eval, 'dense evalmat error');

end


function test_objover_modes()
% objover input modes (plain vector, ifreturnovers, cell-vector broadcast)
% for surfermat/apply and surferkernevalmat/eval with multiple surfers.
% Also checks on-surface eval against dense evalmat.

tol_apply = 1e-6;
tol_eval  = 1e-6;

zk  = 1;
eps = 1e-10;

S1 = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), [1]);
S2 = slicesurfer(geometries.ellipsoid([3,1,1],  [3,3,3],[],4), [1]);
srfrs    = [S1, S2];
kerns    = kernel3d('h','c',zk,[1,1]);
kernseval = kerns(1,:);

skern  = kernel3d('h','s',zk);
src_ms = []; src_ms.r = [5;0;0];
rhs    = -skern.eval(src_ms, merge(srfrs));

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; zeros(1,nplot^2)];

tic
[syscors,  novers]         = surfermat(srfrs, kerns, eps, struct('corrections',1));
[syscors2, objover_precomp] = surfermat(srfrs, kerns, eps, struct('corrections',1,'ifreturnovers',1));
syscors  = syscors  - 0.5*speye(size(syscors));
syscors2 = syscors2 - 0.5*speye(size(syscors2));

[evalcors,  novers_eval]         = surferkernevalmat(srfrs, kernseval, targs, eps, struct('corrections',1));
[evalcors2, objover_eval_precomp] = surferkernevalmat(srfrs, kernseval, targs, eps, ...
    struct('corrections',1,'ifreturnovers',1));
toc

% surfermatapply
ref    = surfermatapply(srfrs, kerns, rhs, eps, novers, syscors);

novers_11     = novers{1,1};
novers_11cell = cell(size(novers));
for ii = 1:numel(novers), novers_11cell{ii} = novers_11; end
ref_11 = surfermatapply(srfrs, kerns, rhs, eps, novers_11cell, syscors);

so_vec = objover_precomp{1}(1,:);  xi_vec = objover_precomp{2}(1,:);

assert(iscell(objover_precomp) && numel(objover_precomp)==2 && iscell(objover_precomp{1}), ...
    'surfermat ifreturnovers=1: second output must be a 2-element cell of cells');
assert(norm(surfermatapply(srfrs, kerns, rhs, eps, novers_11,     syscors)  - ref_11) == 0, ...
    'surfermatapply plain-vector objover gave different result');
assert(norm(surfermatapply(srfrs, kerns, rhs, eps, objover_precomp,  syscors2) - ref)    / norm(ref)    < tol_apply, ...
    'surfermatapply precomputed objover error');
assert(norm(surfermatapply(srfrs, kerns, rhs, eps, {so_vec,xi_vec},  syscors2) - ref)    / norm(ref)    < tol_apply, ...
    'surfermatapply cell-vector broadcast objover error');

% surferkerneval
opts_ref = []; opts_ref.corrections = evalcors; opts_ref.objover = novers_eval;
ref_eval = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ref);

novers_eval_1     = novers_eval{1};
novers_eval_1cell = cellfun(@(~) novers_eval_1, novers_eval, 'UniformOutput', false);
opts_ref_1 = []; opts_ref_1.corrections = evalcors; opts_ref_1.objover = novers_eval_1cell;
ref_eval_1 = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ref_1);

so_eval = objover_eval_precomp{1}(:);  xi_eval = objover_eval_precomp{2}(:);
opts_ep = []; opts_ep.corrections = evalcors2; opts_ep.objover = objover_eval_precomp;
opts_eb = []; opts_eb.corrections = evalcors2; opts_eb.objover = {so_eval, xi_eval};

assert(iscell(objover_eval_precomp) && numel(objover_eval_precomp)==2 && iscell(objover_eval_precomp{1}), ...
    'surferkernevalmat ifreturnovers=1: second output must be a 2-element cell of cells');
assert(norm(surferkerneval(srfrs, kernseval, rhs, targs, eps, ...
    struct('corrections', evalcors, 'objover', novers_eval_1)) - ref_eval_1) == 0, ...
    'surferkerneval plain-vector objover gave different result');
assert(norm(surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ep) - ref_eval) / norm(ref_eval) < tol_eval, ...
    'surferkerneval precomputed objover error');
assert(norm(surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_eb) - ref_eval) / norm(ref_eval) < tol_eval, ...
    'surferkerneval cell-vector broadcast objover error');

% on-surface eval: single surfer, corrections path vs dense evalmat
kern_on      = kernel3d('h', 'c', zk, [1,1]);
kernseval_on = kern_on(1,:);
rhs_on       = -skern.eval(src_ms, S1);

evalmat_on = surferkernevalmat(S1, kernseval_on, S1, eps);
[evalcors_on, objover_on] = surferkernevalmat(S1, kernseval_on, S1, eps, ...
    struct('corrections', 1, 'ifreturnovers', 1));

opts_on = []; opts_on.corrections = evalcors_on; opts_on.objover = objover_on;
assert(norm(surferkerneval(S1, kernseval_on, rhs_on, S1, eps, opts_on) - evalmat_on*rhs_on) / norm(evalmat_on*rhs_on) < tol_eval, ...
    'surferkerneval on-surface vs dense evalmat mismatch');

end


function test_block_kernels()
% Interleaved, block (2x2), and stacked kernel assembly via surfermat/apply/eval.

tol_blk = 1e-13;

S     = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
srfrs = S;
eps   = 1e-10;

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; zeros(1,nplot^2)];

opts_ns = []; opts_ns.nonsmoothonly = 1;
zks = 0.9;

tr_kern = kernel3d('h','trans_sys',zks);
s_kern  = kernel3d('h','s',zks);
d_kern  = kernel3d('h','d',zks);
sp_kern = kernel3d('h','sp',zks);
dp_kern = kernel3d('h','dp',zks);

trmat = surfermat(srfrs, tr_kern, eps, opts_ns);
smat  = surfermat(srfrs, s_kern,  eps, opts_ns);
dmat  = surfermat(srfrs, d_kern,  eps, opts_ns);
spmat = surfermat(srfrs, sp_kern, eps, opts_ns);
dpmat = surfermat(srfrs, dp_kern, eps, opts_ns);

trmat2 = zeros(size(trmat));
trmat2(1:2:end, 1:2:end) = dmat;
trmat2(1:2:end, 2:2:end) = smat;
trmat2(2:2:end, 1:2:end) = dpmat;
trmat2(2:2:end, 2:2:end) = spmat;

assert(norm(trmat-trmat2,1) < tol_blk*norm(trmat,1), 'interleaved surfermat incorrect');

s_kern0 = kernel3d('l','s');
d_kern0 = kernel3d('l','d');

s0mat = surfermat(srfrs, s_kern0, eps, opts_ns);
d0mat = surfermat(srfrs, d_kern0, eps, opts_ns);
tmat  = surfermat(srfrs, kernel3d([s_kern,s_kern0;d_kern,d_kern0]), eps, opts_ns);

assert(norm(tmat(1:2:end,1:2:end) - smat,  1) < tol_blk*norm(smat,  1), 'block surfermat: (1,1) mismatch');
assert(norm(tmat(2:2:end,1:2:end) - dmat,  1) < tol_blk*norm(dmat,  1), 'block surfermat: (2,1) mismatch');
assert(norm(tmat(1:2:end,2:2:end) - s0mat, 1) < tol_blk*norm(s0mat, 1), 'block surfermat: (1,2) mismatch');
assert(norm(tmat(2:2:end,2:2:end) - d0mat, 1) < tol_blk*norm(d0mat, 1), 'block surfermat: (2,2) mismatch');

rhs0    = randn(S.npts,1);
rhs1    = randn(S.npts,1);
rhs_blk = [rhs0(:).'; rhs1(:).']; rhs_blk = rhs_blk(:);

a0 = surfermatapply(srfrs, s_kern0, rhs0, eps);
a1 = surfermatapply(srfrs, s_kern,  rhs1, eps);
b0 = surfermatapply(srfrs, d_kern0, rhs0, eps);
b1 = surfermatapply(srfrs, d_kern,  rhs1, eps);
ab = surfermatapply(srfrs, kernel3d([s_kern0,s_kern;d_kern0,d_kern]), rhs_blk, eps);

assert(norm(ab(1:2:end) - (a0+a1)) < tol_blk*norm(a0+a1), 'block surfermatapply: row 1 mismatch');
assert(norm(ab(2:2:end) - (b0+b1)) < tol_blk*norm(b0+b1), 'block surfermatapply: row 2 mismatch');

c0 = surferkerneval(srfrs, s_kern0, rhs0, targs, eps);
c1 = surferkerneval(srfrs, s_kern,  rhs1, targs, eps);
d0 = surferkerneval(srfrs, d_kern0, rhs0, targs, eps);
d1 = surferkerneval(srfrs, d_kern,  rhs1, targs, eps);
cd = surferkerneval(srfrs, kernel3d([s_kern0,s_kern;d_kern0,d_kern]), rhs_blk, targs, eps);

assert(norm(cd(1:2:end) - (c0+c1)) < tol_blk*norm(c0+c1), 'block surferkerneval: row 1 mismatch');
assert(norm(cd(2:2:end) - (d0+d1)) < tol_blk*norm(d0+d1), 'block surferkerneval: row 2 mismatch');

k1    = 2*(kernel3d('stok','s') + kernel3d('stok','d'));
k2    = kernel3d([s_kern, d_kern, s_kern]);
smat1 = surfermat(srfrs, k1, eps);
smat2 = surfermat(srfrs, k2, eps);
smat  = surfermat(srfrs, kernel3d([k1;k2]), eps);

for i = 1:3
    assert(norm(smat(i:4:end,:) - smat1(i:3:end,:), 1) < tol_blk*norm(smat1,1), ...
        'stacked surfermat: row block %d mismatch', i);
end
assert(norm(smat(4:4:end,:) - smat2, 1) < tol_blk*norm(smat2,1), 'stacked surfermat: bottom block mismatch');

end
