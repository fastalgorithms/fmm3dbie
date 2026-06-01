% Test surfermat, surfermatapply, surferkernevalmat, and surferkerneval.
% Covers: apply paths, eval paths, rfac passthrough, objover input modes,
% on-surface eval, and block/stacked kernel assembly.

tol_apply = 1e-6;
tol_eval  = 1e-6;
tol_blk   = 1e-13;

S     = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
zk    = 1;
kerns = kernel3d('h','c',zk,[1,1]);
srfrs = S;
eps   = 1e-10;

skern = kernel3d('h','s',zk);
src   = []; src.r = [2;0;1];

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; zeros(1,nplot^2)];

%% surfermat / surfermatapply paths

tic
Smat = surfermat(srfrs, kerns, eps);
Smat = Smat - 0.5*eye(size(Smat));

[syscors, novers, rfac_sys] = surfermat(srfrs, kerns, eps, struct('corrections',  1));
syscors  = syscors  - 0.5*speye(size(syscors));

[sysnsmth, novers2]         = surfermat(srfrs, kerns, eps, struct('nonsmoothonly', 1));
sysnsmth = sysnsmth - 0.5*speye(size(sysnsmth));
toc

rhs = -skern.eval(src, merge(srfrs));

assert(iscell(rfac_sys) && isequal(size(rfac_sys), [numel(srfrs), numel(srfrs)]), ...
    'surfermat rfac must be a cell of size (nsurfers x nsurfers)');
assert(~isnan(rfac_sys{1,1}), 'surfermat rfac{1,1} should be non-NaN');

sysapply  = @(d) surfermatapply(srfrs, kerns, d, eps, novers,  syscors);
sysapply2 = @(d) surfermatapply(srfrs, kerns, d, eps, novers2, sysnsmth, struct('usematlab',0));

assert(norm(sysapply(rhs)  - Smat*rhs) / norm(Smat*rhs) < tol_apply, 'corrections apply error');
assert(norm(sysapply2(rhs) - Smat*rhs) / norm(Smat*rhs) < tol_apply, 'nonsmoothonly apply error');

% rfac passthrough (layer_eval path only): must be bit-identical
ref_apply2 = sysapply2(rhs);

opts_rfac = struct('usematlab', 0, 'rfac', rfac_sys);
assert(norm(surfermatapply(srfrs, kerns, rhs, eps, novers2, sysnsmth, opts_rfac) - ref_apply2) == 0, ...
    'surfermatapply with cell opts.rfac gave different result');

opts_rfac_scalar = struct('usematlab', 0, 'rfac', rfac_sys{1,1});
assert(norm(surfermatapply(srfrs, kerns, rhs, eps, novers2, sysnsmth, opts_rfac_scalar) - ref_apply2) == 0, ...
    'surfermatapply with scalar opts.rfac gave different result');

%% surferkerneval / surferkernevalmat paths

kernseval = kerns(1,:);

tic
uscat_smth  = surferkerneval(srfrs, kernseval, rhs, targs, eps);

opts_cor = []; opts_cor.corrections = 1;
[evalcors, novers_eval, rfac_eval] = surferkernevalmat(srfrs, kernseval, targs, eps, opts_cor);
opts_cor.corrections = evalcors;  opts_cor.objover = novers_eval;
uscat_cor   = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_cor);

evalmat     = surferkernevalmat(srfrs, kernseval, targs, eps);
uscat_dense = evalmat * rhs;

uscat_fmm   = surferkerneval(srfrs, kernseval, rhs, targs, eps, struct('usematlab',0));

opts_ns = []; opts_ns.nonsmoothonly = 1;
[nsmth, novers_ns] = surferkernevalmat(srfrs, kernseval, targs, eps, opts_ns);
opts_ns2 = []; opts_ns2.usematlab = 0;
opts_ns2.corrections = nsmth;  opts_ns2.objover = novers_ns;
uscat_nsmth = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ns2);
toc

uref = uscat_fmm;
assert(norm(uref - uscat_smth)  / norm(uref) < tol_eval, 'smooth-only eval error');
assert(norm(uref - uscat_cor)   / norm(uref) < tol_eval, 'corrections eval error');
assert(norm(uref - uscat_dense) / norm(uref) < tol_eval, 'dense evalmat error');
assert(norm(uref - uscat_nsmth) / norm(uref) < tol_eval, 'nonsmoothonly eval error');

assert(iscell(rfac_eval) && numel(rfac_eval) == numel(srfrs), ...
    'surferkernevalmat rfac must be a cell of length nsurfers');
assert(~isnan(rfac_eval{1}), 'surferkernevalmat rfac{1} should be non-NaN');

% rfac passthrough
opts_ns_rfac = struct('usematlab', 0, 'corrections', nsmth, 'objover', novers_ns, 'rfac', rfac_eval);
assert(norm(surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ns_rfac) - uscat_nsmth) == 0, ...
    'surferkerneval with cell opts.rfac gave different result');

opts_ns_rfac_scalar = struct('usematlab', 0, 'corrections', nsmth, 'objover', novers_ns, 'rfac', rfac_eval{1});
assert(norm(surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ns_rfac_scalar) - uscat_nsmth) == 0, ...
    'surferkerneval with scalar opts.rfac gave different result');

%% objover input modes: surfermat/apply and surferkernevalmat/eval, multiple surfers

S1 = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), [1]);
S2 = slicesurfer(geometries.ellipsoid([3,1,1],  [3,3,3],[],4), [1]);
srfrs_ms  = [S1, S2];
kerns_ms  = kernel3d('h','c',zk,[1,1]);
kernseval_ms = kerns_ms(1,:);

src_ms = []; src_ms.r = [5;0;0];
rhs_ms = -skern.eval(src_ms, merge(srfrs_ms));

tic
[syscors_ms,  novers_ms]       = surfermat(srfrs_ms, kerns_ms, eps, struct('corrections',1));
[syscors_ms2, objover_precomp] = surfermat(srfrs_ms, kerns_ms, eps, struct('corrections',1,'ifreturnovers',1));
syscors_ms  = syscors_ms  - 0.5*speye(size(syscors_ms));
syscors_ms2 = syscors_ms2 - 0.5*speye(size(syscors_ms2));

[evalcors_ms,  novers_eval_ms]       = surferkernevalmat(srfrs_ms, kernseval_ms, targs, eps, struct('corrections',1));
[evalcors_ms2, objover_eval_precomp] = surferkernevalmat(srfrs_ms, kernseval_ms, targs, eps, ...
    struct('corrections',1,'ifreturnovers',1));
toc

% surfermatapply: plain-vector, precomputed, cell-vector broadcast
ref_ms = surfermatapply(srfrs_ms, kerns_ms, rhs_ms, eps, novers_ms, syscors_ms);

novers_ms_11     = novers_ms{1,1};
novers_ms_11cell = cell(size(novers_ms));
for ii = 1:numel(novers_ms), novers_ms_11cell{ii} = novers_ms_11; end
ref_ms_11 = surfermatapply(srfrs_ms, kerns_ms, rhs_ms, eps, novers_ms_11cell, syscors_ms);

so_vec = objover_precomp{1}(1,:);  xi_vec = objover_precomp{2}(1,:);

assert(iscell(objover_precomp) && numel(objover_precomp)==2 && iscell(objover_precomp{1}), ...
    'surfermat ifreturnovers=1: second output must be a 2-element cell of cells');
assert(norm(surfermatapply(srfrs_ms, kerns_ms, rhs_ms, eps, novers_ms_11,     syscors_ms)  - ref_ms_11) == 0, ...
    'surfermatapply plain-vector objover gave different result');
assert(norm(surfermatapply(srfrs_ms, kerns_ms, rhs_ms, eps, objover_precomp,  syscors_ms2) - ref_ms) / norm(ref_ms) < tol_apply, ...
    'surfermatapply precomputed objover error');
assert(norm(surfermatapply(srfrs_ms, kerns_ms, rhs_ms, eps, {so_vec,xi_vec},  syscors_ms2) - ref_ms) / norm(ref_ms) < tol_apply, ...
    'surfermatapply cell-vector broadcast objover error');

% surferkerneval: plain-vector, precomputed, cell-vector broadcast
opts_ref_ms = []; opts_ref_ms.corrections = evalcors_ms; opts_ref_ms.objover = novers_eval_ms;
ref_eval_ms = surferkerneval(srfrs_ms, kernseval_ms, rhs_ms, targs, eps, opts_ref_ms);

novers_eval_ms_1     = novers_eval_ms{1};
novers_eval_ms_1cell = cellfun(@(~) novers_eval_ms_1, novers_eval_ms, 'UniformOutput', false);
opts_ref_ms_1 = []; opts_ref_ms_1.corrections = evalcors_ms; opts_ref_ms_1.objover = novers_eval_ms_1cell;
ref_eval_ms_1 = surferkerneval(srfrs_ms, kernseval_ms, rhs_ms, targs, eps, opts_ref_ms_1);

so_eval = objover_eval_precomp{1}(:);  xi_eval = objover_eval_precomp{2}(:);
opts_ep = []; opts_ep.corrections = evalcors_ms2; opts_ep.objover = objover_eval_precomp;
opts_eb = []; opts_eb.corrections = evalcors_ms2; opts_eb.objover = {so_eval, xi_eval};

assert(iscell(objover_eval_precomp) && numel(objover_eval_precomp)==2 && iscell(objover_eval_precomp{1}), ...
    'surferkernevalmat ifreturnovers=1: second output must be a 2-element cell of cells');
assert(norm(surferkerneval(srfrs_ms, kernseval_ms, rhs_ms, targs, eps, ...
    struct('corrections', evalcors_ms, 'objover', novers_eval_ms_1)) - ref_eval_ms_1) == 0, ...
    'surferkerneval plain-vector objover gave different result');
assert(norm(surferkerneval(srfrs_ms, kernseval_ms, rhs_ms, targs, eps, opts_ep) - ref_eval_ms) / norm(ref_eval_ms) < tol_eval, ...
    'surferkerneval precomputed objover error');
assert(norm(surferkerneval(srfrs_ms, kernseval_ms, rhs_ms, targs, eps, opts_eb) - ref_eval_ms) / norm(ref_eval_ms) < tol_eval, ...
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

%% Block and stacked kernel assembly

opts_ns_blk = []; opts_ns_blk.nonsmoothonly = 1;
zks = 0.9;

tr_kern = kernel3d('h','trans_sys',zks);
s_kern  = kernel3d('h','s',zks);
d_kern  = kernel3d('h','d',zks);
sp_kern = kernel3d('h','sp',zks);
dp_kern = kernel3d('h','dp',zks);

trmat = surfermat(srfrs, tr_kern, eps, opts_ns_blk);
smat  = surfermat(srfrs, s_kern,  eps, opts_ns_blk);
dmat  = surfermat(srfrs, d_kern,  eps, opts_ns_blk);
spmat = surfermat(srfrs, sp_kern, eps, opts_ns_blk);
dpmat = surfermat(srfrs, dp_kern, eps, opts_ns_blk);

trmat2 = zeros(size(trmat));
trmat2(1:2:end, 1:2:end) = dmat;
trmat2(1:2:end, 2:2:end) = smat;
trmat2(2:2:end, 1:2:end) = dpmat;
trmat2(2:2:end, 2:2:end) = spmat;

assert(norm(trmat-trmat2,1) < tol_blk*norm(trmat,1), 'interleaved surfermat incorrect');

s_kern0 = kernel3d('l','s');
d_kern0 = kernel3d('l','d');

s0mat = surfermat(srfrs, s_kern0, eps, opts_ns_blk);
d0mat = surfermat(srfrs, d_kern0, eps, opts_ns_blk);

tmat = surfermat(srfrs, kernel3d([s_kern,s_kern0;d_kern,d_kern0]), eps, opts_ns_blk);

assert(norm(tmat(1:2:end,1:2:end) - smat,  1) < tol_blk*norm(smat,  1), 'block surfermat: (1,1) mismatch');
assert(norm(tmat(2:2:end,1:2:end) - dmat,  1) < tol_blk*norm(dmat,  1), 'block surfermat: (2,1) mismatch');
assert(norm(tmat(1:2:end,2:2:end) - s0mat, 1) < tol_blk*norm(s0mat, 1), 'block surfermat: (1,2) mismatch');
assert(norm(tmat(2:2:end,2:2:end) - d0mat, 1) < tol_blk*norm(d0mat, 1), 'block surfermat: (2,2) mismatch');

rhs0 = randn(S.npts,1);
rhs1 = randn(S.npts,1);
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

k1   = 2*(kernel3d('stok','s') + kernel3d('stok','d'));
k2   = kernel3d([s_kern, d_kern, s_kern]);
smat1 = surfermat(srfrs, k1, eps);
smat2 = surfermat(srfrs, k2, eps);
smat  = surfermat(srfrs, kernel3d([k1;k2]), eps);

for i = 1:3
    assert(norm(smat(i:4:end,:) - smat1(i:3:end,:), 1) < tol_blk*norm(smat1,1), ...
        'stacked surfermat: row block %d mismatch', i);
end
assert(norm(smat(4:4:end,:) - smat2, 1) < tol_blk*norm(smat2,1), 'stacked surfermat: bottom block mismatch');
