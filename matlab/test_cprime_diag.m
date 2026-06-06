% Diagnostic: isolate which part of cprime kern disagrees with fmm.
% We test alpha=1,beta=0 (pure sprime) and alpha=0,beta=1 (pure dprime)
% separately, comparing kern matrix-vector product vs fmm.

run ../startup.m
rng(42);

ns = 30; nt = 20;

% Well-separated sources/targets
src.r = randn(3,ns); src.r = src.r./vecnorm(src.r)*0.4;
src.n = randn(3,ns); src.n = src.n./vecnorm(src.n);
targ.r = randn(3,nt); targ.r = targ.r./vecnorm(targ.r)*3.0;
targ.n = randn(3,nt); targ.n = targ.n./vecnorm(targ.n);

sigma = randn(3,ns);
wts   = ones(1,ns)*0.01;

function pot = kern_matvec(src, targ, sigma, wts, alpha, beta)
    ns = size(src.r,2); nt = size(targ.r,2);
    K = stok3d.kern(src, targ, 'cprime', alpha, beta);
    % K: (3,nt,3,ns); sigma: (3,ns); wts: (1,ns)
    sw = sigma.*wts;
    pot = squeeze(sum(K.*reshape(sw,[1,1,3,ns]),[3,4]));
end

function pot = fmm_eval(src, targ, sigma, wts, alpha, beta)
    sw = sigma.*wts;
    coefs = [alpha; beta];
    pot = stok3d.fmm(1e-10, src, targ, 'cprime', sw, coefs);
end

% Test 1: pure sprime (alpha=1, beta=0)
p_kern = kern_matvec(src,targ,sigma,wts,1,0);
p_fmm  = fmm_eval(src,targ,sigma,wts,1,0);
e1 = norm(p_kern(:)-p_fmm(:))/norm(p_fmm(:));
fprintf('alpha=1,beta=0 (sprime): kern vs fmm err = %.3e\n', e1);

% Test 2: pure dprime (alpha=0, beta=1)
p_kern = kern_matvec(src,targ,sigma,wts,0,1);
p_fmm  = fmm_eval(src,targ,sigma,wts,0,1);
e2 = norm(p_kern(:)-p_fmm(:))/norm(p_fmm(:));
fprintf('alpha=0,beta=1 (dprime): kern vs fmm err = %.3e\n', e2);

% Test 3: combined
p_kern = kern_matvec(src,targ,sigma,wts,0.7,1.3);
p_fmm  = fmm_eval(src,targ,sigma,wts,0.7,1.3);
e3 = norm(p_kern(:)-p_fmm(:))/norm(p_fmm(:));
fprintf('alpha=0.7,beta=1.3 (combined): kern vs fmm err = %.3e\n', e3);

% Test 4: check dprime vs FD of d kernel in target normal direction
h = 1e-5;
targ_p = targ; targ_p.r = targ.r + h*targ.n;
targ_m = targ; targ_m.r = targ.r - h*targ.n;
K_d_p = stok3d.kern(src,targ_p,'d');
K_d_m = stok3d.kern(src,targ_m,'d');
dprime_fd = (K_d_p - K_d_m)/(2*h);  % (3,nt,3,ns) FD approx of d/dn_x D

% dprime from kern (beta=1 part):
K_dp = stok3d.kern(src,targ,'cprime',0,1);
e4 = norm(K_dp(:)-dprime_fd(:))/norm(dprime_fd(:));
fprintf('dprime kern vs FD of D: err = %.3e\n', e4);

% Test 5: check sprime vs FD of s kernel
K_s_p = stok3d.kern(src,targ_p,'s');
K_s_m = stok3d.kern(src,targ_m,'s');
sprime_fd = (K_s_p - K_s_m)/(2*h);

K_sp = stok3d.kern(src,targ,'cprime',1,0);
e5 = norm(K_sp(:)-sprime_fd(:))/norm(sprime_fd(:));
fprintf('sprime kern vs FD of S: err = %.3e\n', e5);
