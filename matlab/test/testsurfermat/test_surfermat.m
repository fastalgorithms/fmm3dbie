% Test that surfermat, surfermatapply, surferkernevalmat, and surferkerneval
% all agree to within discretization tolerance.
%
% Three paths through the system matrix are exercised:
%   Smat      - dense matrix (reference)
%   sysapply  - matlab smooth apply + corrections sparse mat
%   sysapply2 - layer_eval apply + nonsmoothonly sparse mat
%
% Three paths through the evaluation matrix are exercised:
%   evalmat  - dense matrix (reference)
%   uscat2   - corrections-path apply
%   uscat4   - layer_eval apply  (also used as reference for the rest)
%   uscat5   - nonsmoothonly-path apply

tol_apply = 1e-6;   % relative tolerance for apply consistency checks
tol_eval  = 1e-6;   % relative tolerance for evaluation consistency checks

%% Geometry and kernel

S    = geometries.ellipsoid([1,1,1.1],[3,3,3],[],6);
S = slicesurfer(S,[1]);
zk   = 1;
kerns = kernel3d('h','c',zk,[1,1]);
srfrs = S;
eps   = 1e-10;

%% Build system matrices: dense, corrections-sparse, nonsmoothonly-sparse

tic
Smat = surfermat(srfrs, kerns, eps);
Smat = Smat - 0.5*eye(size(Smat));

[syscors, novers]  = surfermat(srfrs, kerns, eps, struct('corrections',   1));
syscors = syscors - 0.5*speye(size(syscors));

[sysnsmth, novers2] = surfermat(srfrs, kerns, eps, struct('nonsmoothonly', 1));
sysnsmth = sysnsmth - 0.5*speye(size(sysnsmth));
toc

% RHS: exterior Helmholtz point-source

skern  = kernel3d('h','s',zk);
smerge = merge(srfrs);
src    = []; src.r = [2;0;1];
rhs    = -skern.eval(src, smerge);

assert(numel(rhs) == size(Smat,1), ...
    'rhs length (%d) must match system matrix rows (%d)', numel(rhs), size(Smat,1));

%% Verify that both apply paths match the dense matrix-vector product

sysapply  = @(dens) surfermatapply(srfrs, kerns, dens, eps, novers,  syscors);
sysapply2 = @(dens) surfermatapply(srfrs, kerns, dens, eps, novers2, sysnsmth, struct('usematlab',0));

matlabapply_err = norm(sysapply(rhs)  - Smat*rhs) / norm(Smat*rhs);
layerapply_err  = norm(sysapply2(rhs) - Smat*rhs) / norm(Smat*rhs);

assert(matlabapply_err < tol_apply, ...
    'corrections apply error %.2e exceeds tolerance %.2e', matlabapply_err, tol_apply);
assert(layerapply_err < tol_apply, ...
    'nonsmoothonly apply error %.2e exceeds tolerance %.2e', layerapply_err, tol_apply);

%% Build off-surface evaluation matrices and compare all eval paths

kernseval = kerns(1,:);

nplot = 30;
xx = linspace(-2,6,nplot);  yy = 2*linspace(-2,2,nplot);
[XX,YY] = meshgrid(xx,yy);
targs = []; targs.r = [XX(:).'; YY(:).'; 0*YY(:).'];

tic

% Smooth quadrature only, no correction
uscat_smth = surferkerneval(srfrs, kernseval, rhs, targs, eps);

% Corrections sparse mat path
opts_cor = []; opts_cor.corrections = 1;
[evalcors, novers_eval] = surferkernevalmat(srfrs, kernseval, targs, eps, opts_cor);
opts_cor.corrections = evalcors;  opts_cor.novers = novers_eval;
uscat_cor = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_cor);

% Dense evaluation matrix
evalmat  = surferkernevalmat(srfrs, kernseval, targs, eps);
uscat_dense = evalmat * rhs;

% layer_eval path (reference for eval comparisons)
uscat_fmm = surferkerneval(srfrs, kernseval, rhs, targs, eps, struct('usematlab',0));

% Nonsmoothonly sparse mat + layer_eval path
opts_ns = []; opts_ns.nonsmoothonly = 1;
[nsmth, novers_ns] = surferkernevalmat(srfrs, kernseval, targs, eps, opts_ns);
opts_ns2 = []; opts_ns2.usematlab = 0;
opts_ns2.corrections = nsmth;  opts_ns2.novers = novers_ns;
uscat_nsmth = surferkerneval(srfrs, kernseval, rhs, targs, eps, opts_ns2);

toc

uref = uscat_fmm;
eval_err_smth  = norm(uref - uscat_smth)  / norm(uref);
eval_err_cor   = norm(uref - uscat_cor)   / norm(uref);
eval_err_dense = norm(uref - uscat_dense) / norm(uref);
eval_err_nsmth = norm(uref - uscat_nsmth) / norm(uref);

assert(eval_err_smth < tol_eval, ...
    'smooth-only eval error %.2e exceeds tolerance %.2e', eval_err_smth, tol_eval);
assert(eval_err_cor < tol_eval, ...
    'corrections eval error %.2e exceeds tolerance %.2e', eval_err_cor, tol_eval);
assert(eval_err_dense < tol_eval, ...
    'dense evalmat error %.2e exceeds tolerance %.2e', eval_err_dense, tol_eval);
assert(eval_err_nsmth < tol_eval, ...
    'nonsmoothonly eval error %.2e exceeds tolerance %.2e', eval_err_nsmth, tol_eval);


%%

opts = [];
opts.nonsmoothonly = 1;
zks = [0.9];
% tr_kern = kernel3d('h','trans_sys_diff',zks);
% 
% s_kern  = kernel3d('h','s_diff',zks);
% d_kern  = kernel3d('h','d_diff',zks);
% sp_kern = kernel3d('h','sp_diff',zks);
% dp_kern = kernel3d('h','dp_diff',zks);

tr_kern = kernel3d('h','trans_sys',zks);

s_kern  = kernel3d('h','s',zks);
d_kern  = kernel3d('h','d',zks);
sp_kern = kernel3d('h','sp',zks);
dp_kern = kernel3d('h','dp',zks);

trmat = surfermat(srfrs, tr_kern, eps,opts);

smat  = surfermat(srfrs, s_kern, eps,opts);
dmat  = surfermat(srfrs, d_kern, eps,opts);
spmat = surfermat(srfrs, sp_kern, eps,opts);
dpmat = surfermat(srfrs, dp_kern, eps,opts);
trmat2 = zeros(size(trmat));
trmat2(1:2:end, 1:2:end) = dmat;
trmat2(1:2:end, 2:2:end) = smat;
trmat2(2:2:end, 1:2:end) = dpmat;
trmat2(2:2:end, 2:2:end) = spmat;

assert(norm(trmat-trmat2,1)<1e-14*norm(trmat,1), 'interleaved surfermat inccorrect')

%%

s_kern0 = kernel3d('l','s');
d_kern0 = kernel3d('l','d');

smat  = surfermat(srfrs, s_kern, eps,opts);
dmat  = surfermat(srfrs, d_kern, eps,opts);

s0mat  = surfermat(srfrs, s_kern0, eps,opts);
d0mat  = surfermat(srfrs, d_kern0, eps,opts);

tmat = surfermat(srfrs, kernel3d([s_kern,s_kern0;d_kern,d_kern0]), eps,opts);

tol_blk = 1e-13;
assert(norm(tmat(1:2:end,1:2:end) - smat,  1) < tol_blk*norm(smat,  1), ...
    'block surfermat: (1,1) block mismatch');
assert(norm(tmat(2:2:end,1:2:end) - dmat,  1) < tol_blk*norm(dmat,  1), ...
    'block surfermat: (2,1) block mismatch');
assert(norm(tmat(1:2:end,2:2:end) - s0mat, 1) < tol_blk*norm(s0mat, 1), ...
    'block surfermat: (1,2) block mismatch');
assert(norm(tmat(2:2:end,2:2:end) - d0mat, 1) < tol_blk*norm(d0mat, 1), ...
    'block surfermat: (2,2) block mismatch');

rhs0 = randn(S.npts,1);
rhs1 = randn(S.npts,1);
rhs = [rhs0(:).';rhs1(:).']; rhs = rhs(:);

a0 = surfermatapply(srfrs, s_kern0, rhs0,eps);
a1 = surfermatapply(srfrs, s_kern, rhs1,eps);
b0 = surfermatapply(srfrs, d_kern0, rhs0,eps);
b1 = surfermatapply(srfrs, d_kern, rhs1,eps);

ab = surfermatapply(srfrs, kernel3d([s_kern0,s_kern;d_kern0,d_kern]), rhs,eps);

assert(norm(ab(1:2:end) - (a0+a1)) < tol_blk*norm(a0+a1), ...
    'block surfermatapply: row 1 mismatch');
assert(norm(ab(2:2:end) - (b0+b1)) < tol_blk*norm(b0+b1), ...
    'block surfermatapply: row 2 mismatch');

c0 = surferkerneval(srfrs, s_kern0, rhs0,targs,eps);
c1 = surferkerneval(srfrs, s_kern, rhs1,targs,eps);
d0 = surferkerneval(srfrs, d_kern0, rhs0,targs,eps);
d1 = surferkerneval(srfrs, d_kern, rhs1,targs,eps);

cd = surferkerneval(srfrs, kernel3d([s_kern0,s_kern;d_kern0,d_kern]), rhs,targs,eps);

assert(norm(cd(1:2:end) - (c0+c1)) < tol_blk*norm(c0+c1), ...
    'block surferkerneval: row 1 mismatch');
assert(norm(cd(2:2:end) - (d0+d1)) < tol_blk*norm(d0+d1), ...
    'block surferkerneval: row 2 mismatch');

%%
k1 = 2*(kernel3d('stok','s')+kernel3d('stok','d'));
k2 = kernel3d([s_kern,d_kern, s_kern]);

smat1 = surfermat(srfrs,k1,eps);
smat2 = surfermat(srfrs,k2,eps);

smat = surfermat(srfrs, kernel3d([k1;k2]),eps);

for i = 1:3
    assert(norm(smat(i:4:end,:) - smat1(i:3:end,:),1) < tol_blk*norm(smat1,1), ...
        'stacked surfermat: row block %d mismatch', i);
end
assert(norm(smat(4:4:end,:) - smat2,1) < tol_blk*norm(smat2,1), ...
    'stacked surfermat: bottom block mismatch');
