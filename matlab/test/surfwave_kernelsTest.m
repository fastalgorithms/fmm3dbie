% SURFWAVE_KERNELSTEST  Compare Fortran and MATLAB surface wave kernels.
%
% Applies the Fortran near quadrature weights and the MATLAB smooth rule to
% smooth densities at far targets and checks that the integrals agree.
% Kernels tested:
%   flexural:  gs, gphi, gphi_bilap, s3d_gphi, s3d, s3d_sum
%   capillary: gs, gphi, lap_gphi, s3d_gphi, s3d, s3d_sum, gs_sprime,
%              gphi_sprime
%   gravity:   gs, gphi
% and that kernel3d capillary, iceflex and gravity types construct.

% run ../startup.m

rng(42);

norder = 8;
eps    = 1e-12;
tol    = 1e-9;
nfail  = 0;

S = flat_surfer(norder);
[srcvals, srccoefs, norders, ixyzs, iptype, wts] = extract_arrays(S);
npatches = S.npatches;
npts     = S.npts;

ntarg = 12;
th    = linspace(0, 2*pi, ntarg+1); th(end) = [];
rad   = 40;

tn = [cos(th); sin(th); ones(1,ntarg)];
tn = tn ./ vecnorm(tn);
targinfo          = [];
targinfo.r        = [rad*cos(th); rad*sin(th); zeros(1,ntarg)];
targinfo.du       = repmat([1;0;0], 1, ntarg);
targinfo.dv       = repmat([0;1;0], 1, ntarg);
targinfo.n        = tn;
targinfo.patch_id = -ones(ntarg,1);
targinfo.uvs_targ = zeros(2,ntarg);

targinfo_ev    = targinfo;
targinfo_ev.d  = repmat([1;0;0], 1, ntarg);
targinfo_ev.d2 = zeros(3, ntarg);

[rsc, nquad] = all_pairs_rsc(S, ntarg);

srcinfo = [];
srcinfo.r  = S.r;  srcinfo.n = S.n;
srcinfo.du = S.du; srcinfo.dv = S.dv;
srcinfo.d  = S.du; srcinfo.d2 = zeros(3,npts);

ref = @(K) K(srcinfo, targinfo) .* repmat(wts(:).', ntarg, 1);

xs = S.r(1,:).';  ys = S.r(2,:).';
sig = [ones(npts,1), ...
       xs + 2*ys, ...
       exp(-((xs-0.3).^2 + (ys-0.2).^2))];

lap3d  = @(s,t) surfwave.flex.lap3dkern(s.r, t.r);
algmat = ref(lap3d);

fprintf('\n=== flexural ===\n');
alpha_f = 23; gamma_f = 0.5; nu = 0.33;
[rts_f, ejs_f] = surfwave.flex.find_roots_flex(alpha_f, gamma_f);
zp_f = complex([rts_f(:); ejs_f(:); alpha_f; gamma_f; nu]);

flexgw = @(iker) surfwave.flex.getnearquad_flex(npatches, norders, ...
    ixyzs, iptype, npts, srccoefs, srcvals, targinfo, targinfo.patch_id, ...
    targinfo.uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, rsc.col_ind, ...
    rsc.iquad, rsc.rfac0, zp_f, nquad, iker);

flexk = @(ktype) @(s,t) surfwave.flex.kern(s, t, ktype, nu, rts_f, ejs_f);
flex_s3dg = flexk('s3d_gphi');

nfail = nfail + check('flex gs       ', flexgw(1), ref(flexk('gs_s')),      S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('flex gphi     ', flexgw(2), ref(flexk('gphi_s')),    S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('flex gphi bila', flexgw(4), ref(flexk('gphi_bilap')),S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('flex s3d gphi ', flexgw(5), ref(flexk('s3d_gphi')),  S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('flex s3d lap  ', flexgw(6), ref(lap3d), S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('flex s3d sum  ', flexgw(7), ref(@(s,t) lap3d(s,t) + flex_s3dg(s,t)), S, rsc, ntarg, sig, tol, algmat);

fprintf('\n=== capillary ===\n');
beta = 1.3; gamma_c = -0.7;
[rts_c, ejs_c] = surfwave.capillary.find_roots_capillary(beta, gamma_c);
zp_c = complex([rts_c(:); ejs_c(:)]);

capgw = @(iker) surfwave.capillary.getnearquad_capillary(npatches, ...
    norders, ixyzs, iptype, npts, srccoefs, srcvals, targinfo, ...
    targinfo.patch_id, targinfo.uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, ...
    rsc.col_ind, rsc.iquad, rsc.rfac0, zp_c, nquad, iker);

capk = @(ktype) @(s,t) surfwave.capillary.kern(rts_c, ejs_c, s, t, ktype);
cap_s3dg = capk('s3d_gphi');
lap3d    = @(s,t) surfwave.flex.lap3dkern(s.r, t.r);

nfail = nfail + check('cap gs        ', capgw(0), ref(capk('gs_s')),        S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap gphi      ', capgw(1), ref(capk('gphi_s')),      S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap lap gphi  ', capgw(3), ref(capk('lap_gphi')),    S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap s3d gphi  ', capgw(5), ref(capk('s3d_gphi')),    S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap s3d lap   ', capgw(6), ref(lap3d), S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap s3d sum   ', capgw(7), ref(@(s,t) lap3d(s,t) + cap_s3dg(s,t)), S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap sprime gs ', capgw(8), ref(capk('gs_sprime')),   S, rsc, ntarg, sig, tol, algmat);
nfail = nfail + check('cap sprime gp ', capgw(9), ref(capk('gphi_sprime')), S, rsc, ntarg, sig, tol, algmat);

fprintf('\n=== gravity ===\n');
g = 1.7;
zp_g = complex([g; 0; 0; 0; 0; 0]);

gravgw = @(iker) surfwave.gravity.getnearquad_gravity(npatches, ...
    norders, ixyzs, iptype, npts, srccoefs, srcvals, targinfo, ...
    targinfo.patch_id, targinfo.uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, ...
    rsc.col_ind, rsc.iquad, rsc.rfac0, zp_g, nquad, iker, ...
    (iker==1)*1 + (iker~=1)*g);

gravk = @(ktype) @(s,t) surfwave.gravity.kern(g, s, t, ktype);

nfail = nfail + check('grav gs       ', gravgw(0), ref(gravk('gs_s')),   S, rsc, ntarg, sig, tol);
nfail = nfail + check('grav gphi     ', gravgw(1), ref(gravk('gphi_s')), S, rsc, ntarg, sig, tol);

fprintf('\n=== kernel3d construction ===\n');
nfail = nfail + check_ctor('capillary', {{'s'},{'lap'},{'s3d'},{'s3d_lap'}, ...
    {'s3d_sum'},{'sp'},{'d'},{'dp'}}, ...
    @(ty) kernel3d('capillary', ty{1}, rts_c, ejs_c), srcinfo, targinfo_ev);
nfail = nfail + check_ctor('iceflex', {{'gs'},{'gphi'},{'gphi_bilap'},{'s3d'}, ...
    {'s3d_gphi'},{'s3d_sum'},{'gs_v2b'},{'gphi_v2b'}}, ...
    @(ty) kernel3d('iceflex', ty{1}, nu, rts_f, ejs_f), srcinfo, targinfo_ev);
nfail = nfail + check_ctor('gravity', {{'s'},{'grad'}}, ...
    @(ty) kernel3d('gravity', ty{1}, g), srcinfo, targinfo_ev);

fprintf('\n');
assert(nfail == 0, 'surfwave_kernelsTest: %d check(s) failed.', nfail);
fprintf('surfwave_kernelsTest: all checks passed.\n');

function S = flat_surfer(norder)
%FLAT_SURFER  Two flat triangular patches in the z = 0 plane.
uv    = koorn.rv_nodes(norder);
npols = size(uv, 2);
verts = {[0 0 0; 1 0 0; 0 1 0].', [1 1 0; 0 1 0; 1 0 0].'};
srcvals = zeros(12, 2*npols);
for ip = 1:2
    V  = verts{ip};
    du = V(:,2) - V(:,1);
    dv = V(:,3) - V(:,1);
    n  = cross(du, dv);  n = n/norm(n);
    idx = (ip-1)*npols + (1:npols);
    srcvals(1:3,  idx) = V(:,1) + du*uv(1,:) + dv*uv(2,:);
    srcvals(4:6,  idx) = repmat(du, 1, npols);
    srcvals(7:9,  idx) = repmat(dv, 1, npols);
    srcvals(10:12,idx) = repmat(n,  1, npols);
end
S = surfer(2, norder, srcvals, 1);
end

function [rsc, nquad] = all_pairs_rsc(S, ntarg)
%ALL_PAIRS_RSC  RSC pattern listing every (target, patch) pair as near.
npatches = S.npatches;
ixyzs    = S.ixyzs(:);
npols    = ixyzs(2:end) - ixyzs(1:end-1);

row_ptr = (0:npatches:npatches*ntarg).' + 1;
col_ind = repmat((1:npatches).', ntarg, 1);
nnz     = numel(col_ind);
iquad   = [1; 1 + cumsum(npols(col_ind))];

rsc         = [];
rsc.row_ptr = row_ptr;
rsc.col_ind = col_ind;
rsc.iquad   = iquad;
rsc.nnz     = nnz;
rsc0        = getnear(S, S.r(:,1));
rsc.rfac0   = rsc0.rfac0;
nquad       = iquad(end) - 1;
end

function nf = check(name, wnear, refmat, S, rsc, ntarg, sig, tol, algmat)
%CHECK  Compare the two rules as operators applied to smooth densities.

wnear = wnear(:).';
ixyzs = S.ixyzs(:);
got   = zeros(size(refmat));
for i = 1:ntarg
    for kk = rsc.row_ptr(i):rsc.row_ptr(i+1)-1
        ip   = rsc.col_ind(kk);
        cols = ixyzs(ip):ixyzs(ip+1)-1;
        got(i, cols) = wnear(rsc.iquad(kk):rsc.iquad(kk+1)-1);
    end
end

pf  = got    * sig;
pm  = refmat * sig;
den = max(max(abs(pm), [], 1), realmin);
err = max(max(abs(pf - pm), [], 1) ./ den);

pmc  = conj(refmat) * sig;
denc = max(max(abs(pmc), [], 1), realmin);
errc = max(max(abs(pf - pmc), [], 1) ./ denc);

ew = norm(got(:) - refmat(:), inf) / max(norm(refmat(:), inf), realmin);

if ~isfinite(err) || err > tol
    rat = got(1,1:min(4,size(got,2))) ./ refmat(1,1:min(4,size(refmat,2)));
    if nargin >= 9 && ~isempty(algmat)
        D = refmat - got;
        c = (algmat(:)'*D(:)) / (algmat(:)'*algmat(:));
        resid = norm(D(:) - c*algmat(:), inf) / max(norm(D(:), inf), realmin);
        fprintf(2, '        ref-got vs Laplace: c = %10.6g%+10.6gi   resid %.2e  %s\n', ...
            real(c), imag(c), resid, ...
            ternary(resid < 1e-6, '<-- MISSING ALGEBRAIC TERM', ''));
    end
    fprintf(2, '  FAIL  %s  integral %.3e  (tol %.1e)   [entrywise %.1e]\n', ...
        name, err, tol, ew);
    fprintf(2, '        vs conj(ref): %.3e   %s\n', errc, ...
        ternary(errc < tol, '<-- CONJUGATED', ''));
    fprintf(2, '        first ratios got/ref: %s\n', ...
        strjoin(arrayfun(@(z) sprintf('%.4g%+.4gi', real(z), imag(z)), ...
        rat, 'UniformOutput', false), '  '));
    nf = 1;
else
    fprintf('  ok    %s  integral %.3e   [entrywise %.1e]\n', name, err, ew);
    nf = 0;
end
end

function nf = check_ctor(family, types, ctor, srcinfo, targinfo)
%CHECK_CTOR  Every type constructs, and eval returns the declared opdims.
nf = 0;
ns = size(srcinfo.r, 2);
nt = size(targinfo.r, 2);
for i = 1:numel(types)
    ty = types{i};
    try
        K = ctor(ty);
        A = K.eval(srcinfo, targinfo);
        want = [K.opdims(1)*nt, K.opdims(2)*ns];
        if ~isequal(size(A), want)
            fprintf(2, '  FAIL  %s %-22s opdims %s, eval gave %s\n', family, ty{1}, ...
                mat2str(want), mat2str(size(A)));
            nf = nf + 1;
        else
            fprintf('  ok    %s %-22s [%d %d]\n', family, ty{1}, K.opdims(1), K.opdims(2));
        end
    catch ME
        fprintf(2, '  FAIL  %s %-22s %s\n', family, ty{1}, ME.message);
        nf = nf + 1;
    end
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
