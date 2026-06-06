% Test: compare stok3d.kern 'cprime' against stok3d.fmm 'cprime'
% for a small random set of sources and targets.
%
% cprime kernel matrix acts on density sigma as:
%   pot_i(x) = sum_is sum_j kern(i,it,j,is) * sigma(j,is) * wts(is)
% whereas fmm evaluates the integral directly.
%
% We skip the diagonal (self-interaction) which has a jump.

rng(42);
ns = 20;
nt = 15;

% Random sources/targets well-separated from each other
src.r = randn(3, ns)*0.3 + [3;0;0];
src.n = randn(3, ns); src.n = src.n ./ vecnorm(src.n);
targ.r = randn(3, nt)*0.3 + [-3;0;0];
targ.n = randn(3, nt); targ.n = targ.n ./ vecnorm(targ.n);

sigma = randn(3, ns);
wts   = rand(1, ns)*0.1 + 0.05;

alpha = 1.7;
beta  = 0.9;
coefs = [alpha; beta];

% --- kernel matrix evaluation ---
K = stok3d.kern(src, targ, 'cprime', alpha, beta);
% K is (3, nt, 3, ns); contract with sigma.*wts
pot_kern = squeeze(sum(K .* reshape(sigma.*wts, [1,1,3,ns]), [3,4]));
% pot_kern: (3, nt)

% --- FMM evaluation ---
pot_fmm = stok3d.fmm(1e-10, src, targ, 'cprime', sigma.*wts, coefs);

err = norm(pot_kern(:) - pot_fmm(:)) / norm(pot_fmm(:));
fprintf('Relative error kern vs fmm: %.3e\n', err);

if err < 1e-4
    fprintf('AGREE: both implementations match.\n');
else
    fprintf('DISAGREE: there is a bug in one of them.\n');
    % Show per-target breakdown
    for it = 1:nt
        eit = norm(pot_kern(:,it) - pot_fmm(:,it));
        fprintf('  targ %2d: |kern-fmm| = %.3e\n', it, eit);
    end
end
