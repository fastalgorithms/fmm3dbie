% Plane wave scattering from a perfectly conducting torus, using the
% generalized Debye source representation [1]. Checked against NRCCIE.
% This demo takes a few minutes and needs FLAM.
%
% Representation:
%   E = ik S[J] - grad S[r] - curl S[M]
%   H = ik S[M] - grad S[q] + curl S[J]
%   J = ik (grad_G L^{-1} r - n x grad_G L^{-1} q) + a_1 h_1 + a_2 h_2
%   M = n x J
% L is the surface Laplacian and h_1, h_2 are harmonic fields.
%
% Equations for (r, q, a), are rescaled versions of:
%   n.H = -n.H_in,  S_0[div_G E_tan] = -S_0[div_G E_in,tan],
%   <E_tan, h_j> = -<E_in, h_j>
%
% L^{-1} uses the parametrix K, L^{-1} r = K sigma with
% (I + R + W) sigma = r, as in [2]. I + R + W is factored with rskelf
% and the gradients of K are applied with ifmm.
%
% [1] C. L. Epstein and L. Greengard, Debye sources and the numerical
%     solution of the time harmonic Maxwell equations, Comm. Pure Appl.
%     Math. 63 (2010), 413-463.
% [2] T. Goodwill, J. Hoskins, Z. Gimbutas and B. Wu, A parametrix for
%     the surface Stokes equation, arXiv:2602.20395 (2026).

assert(exist('rskelf', 'file') == 2, 'pec_debye_demo needs FLAM on the path');

%% Geometry and kernels
S = geometries.startorus([2,1,0.2], 3, [1,1,1], [12 8], 8);

n = S.npts;  w = S.wts(:).';  w3 = repelem(w, 3).';
zk = 1.2;  eps = 1e-8;
% Points are at most dmax apart. FLAM's proxy shells extend a little
% farther.
dmax = norm(max(S.r, [], 2) - min(S.r, [], 2));
rmax = 4*dmax;

N  = @(nn) reshape([0*nn(1,:); nn(3,:); -nn(2,:); -nn(3,:); 0*nn(1,:); ...
                    nn(1,:); nn(2,:); -nn(1,:); 0*nn(1,:)], 3, 3, []);   % n x
P  = @(nn) eye(3) - reshape(nn, 3, 1, []).*reshape(nn, 1, 3, []);
pg = @(v) reshape(v, size(v,1), 1, []);

tic;
% grad K, with lap2d.green = -K
wg = @(t) cat(2, zeros(size(t.r(:,:),2), 1, 3), ...
              repmat(reshape(-eye(3), 1, 3, 3), size(t.r(:,:),2), 1, 1));
Kg = kernel3d.radpp(@lap2d.green, rmax, wg, 'nderiv', 1, 'opdims', [3 1], ...
    'ifcomplex', false);
KF  = kernel3d('bel','rlb') + kernel3d.ones()/S.area;          % I + R + W

% [J;r;q] -> [E;H], plus the exterior jumps (r n - J)/2 and (q n - n x J)/2
nxJq = @(s) [cat(2, N(s.n(:,:)), zeros(3, 2, size(s.r(:,:),2))); ...
             repmat([0 0 0 0 1], 1, 1, size(s.r(:,:),2))];
jumps = @(t) 0.5*[cat(2, -eye(3) + 0*N(t.n(:,:)), pg(t.n(:,:)), 0*pg(t.n(:,:))); ...
                  cat(2, -N(t.n(:,:)), 0*pg(t.n(:,:)), pg(t.n(:,:)))];
Kev = kernel3d('em','nrccie-eval',zk);
Kfo = Kev*[eye(3) zeros(3,2); 0 0 0 1 0] + ...
      [zeros(3) -eye(3); eye(3) zeros(3)]*(Kev*nxJq);
Kfo.src_fields = union(Kfo.src_fields, {'n'});
Kfo.fmm = @(eps, s, t, sigma) debye_fmm(eps, zk, s, t, sigma);
Kfld = Kfo + kernel3d.eye(jumps);
Kfld.targ_fields = union(Kfld.targ_fields, {'n'});

% S_0[div e]
wd   = @(t) repmat(reshape([zeros(1,3); eye(3)], 1, 4, 3), size(t.r(:,:),2), 1, 1);
Kdiv = kernel3d.radpp(@lap3d.green, rmax, wd, 'nderiv', 1, 'opdims', [1 3], ...
    'ifcomplex', false);
Kdiv.fmm = @div_fmm;
fprintf('kernels: %.1fs\n', toc);

%% Surface Laplacian inverse and harmonic fields
tic;
% Outgoing, the proxies only need the parametrix, not the derivatives in
% R. Incoming, the derivatives act on the targets in the box, so the
% proxy surface is also built for KF.
Kp = kernel3d('bel','klb');
fF = byindex.flamfuns(S, KF, eps, [], struct('type', 'rskelf', ...
    'kern_out', Kp, 'kern_in', KF, 'occ_max', 16384));
FF = rskelf(fF.Afun, fF.x, fF.occ, eps, fF.pxyfun);
fprintf('rskelf of I+R+W: %.1fs\n', toc);

tic;
% Outgoing, the proxies use K itself. Incoming, the gradient is at the
% targets in the box.
fG = byindex.flamfuns(S, Kg, eps, [], struct('type', 'ifmm', ...
    'kern_out', Kp, 'kern_in', Kg));
FG = ifmm(fG.Afun, fG.rx, fG.cx, fG.occ, eps, fG.pxyfun, struct('store','a'));
fprintf('ifmm of grad K: %.1fs\n', toc);

% grad_G L^{-1} r - n x grad_G L^{-1} q
gradK    = @(y) ifmm_mv(FG, fG.xinterp*y, fG.Afun) + fG.Q*y;
PN       = cat(2, P(S.n), -N(S.n));
beltrami = @(r, q) reshape(pagemtimes(PN, reshape(permute(reshape( ...
    gradK(rskelf_sv(FF, [r(:) q(:)])), 3, n, 2), [1 3 2]), 6, 1, [])), 3, []);

% Harmonic fields, with phi = atan2(y,x):
%   h_1 = grad_G phi - grad_G L^{-1} (L phi)
%   h_2 = n x h_1
tic;
x1 = S.r(1,:);  y1 = S.r(2,:);  rho2 = x1.^2 + y1.^2;
g = [-y1; x1; 0*x1]./rho2;
nHn = (2*x1.*y1.*(S.n(1,:).^2 - S.n(2,:).^2) + ...
       2*(y1.^2 - x1.^2).*S.n(1,:).*S.n(2,:))./rho2.^2;
Lphi = -nHn - 2*S.mean_curv(:).'.*sum(S.n.*g);
h = g - S.n.*sum(S.n.*g) - beltrami(Lphi, 0*Lphi);
Hv = [h(:), reshape(cross(S.n, h), [], 1)];
Hv = Hv/chol(Hv.'*(w3.*Hv));                                    % orthonormal
ng = size(Hv, 2);
fprintf('harmonic fields: %.1fs\n', toc);

%% Near corrections
tic;
[Cf, ovf] = surfermat(S, Kfld, eps, struct('corrections', true));
[Cd, ovd] = surfermat(S, Kdiv, eps, struct('corrections', true));
fprintf('near corrections: %.1fs\n', toc);

%% Solve
% x = [r q] at each point, then a. bc gives n.H, S_0[div_G E_tan] and
% <E_tan, h_j> from [E;H], and sc scales them to I + compact.
tang = @(E) reshape(E - S.n.*sum(S.n.*E), [], 1);
bc   = @(EH) [reshape([sum(S.n.*EH(4:6,:)); ...
                       surfermatapply(S, Kdiv, tang(EH(1:3,:)), eps, ovd, Cd).'], [], 1); ...
              Hv.'*(w3.*tang(EH(1:3,:)))];
sc   = [repmat([2; 4], n, 1); -2*ones(ng, 1)];
e0   = [repmat([0; 1], n, 1); zeros(ng, 1)];
dens = @(x) reshape([1i*zk*beltrami(x(1:2:2*n), x(2:2:2*n)) + ...
                     reshape(Hv*x(2*n+1:end), 3, n); x(1:2:2*n).'; x(2:2:2*n).'], [], 1);
A = @(x) sc.*bc(reshape(surfermatapply(S, Kfld, dens(x), eps, ovf, Cf), 6, [])) + ...
         e0*(w*x(1:2:2*n))/S.area;

[Ein, Hin] = em3d.planewave(zk, [pi/3, 0], [1; 0.2], S);
rhs = -sc.*bc([Ein; Hin]);

tic;
x = gmres(A, rhs, [], 1e-5, 200);
tgmres = toc

%% Check against NRCCIE away from the surface
targ = randn(3, 40);  targ = 4*targ./vecnorm(targ);
tic;
EH = reshape(surferkerneval(S, Kfo, dens(x), targ, eps), 6, []);
fprintf('Debye reference field evaluation: %.1fs\n', toc);
tic;
dn = em3d.pec.solver(S, Ein, Hin, eps, zk, 1.0, struct('eps_gmres', 1e-10));
fprintf('NRCCIE reference solve: %.1fs\n', toc);
[En, Hn] = em3d.pec.eval(S, dn, struct('r', targ), eps, zk, 1.0);
fprintf('E vs NRCCIE %.2e, H vs NRCCIE %.2e\n', ...
    norm(EH(1:3,:)-En,'fro')/norm(En,'fro'), norm(EH(4:6,:)-Hn,'fro')/norm(Hn,'fro'));

%% Plot solutions and harmonic fields
r = x(1:2:2*n).';  q = x(2:2:2*n).';
J = reshape(dens(x), 5, []);  J = J(1:3,:);
funs = {r, 'r'; q, 'q'; J, 'J'; reshape(Hv(:,1), 3, n), 'h_1'};

figure(1); clf
subplot(1,2,1); plot(S, real(r)); colorbar; title('Re r')
subplot(1,2,2); plot(S, real(q)); colorbar; title('Re q')

figure(2); clf
for k = 1:size(funs, 1)
    f = funs{k,1};
    err = max(surf_fun_error(S, f), [], 1)/max(abs(f(:)));
    subplot(2,2,k); plot(S, log10(err)); colorbar
    title(['log_{10} relative tail, ' funs{k,2}])
end


figure(4); clf
for j = 1:ng
    hj = reshape(Hv(:,j), 3, n);
    subplot(1,ng,j); plot(S, vecnorm(hj)); hold on
    quiver(S, hj, 4, 'k'); hold off
    cb = colorbar; ylabel(cb, sprintf('|h_%d|', j))
    title(sprintf('h_%d', j))
end

%% Plot the electric field
nplt = 150;
xx = linspace(-5, 5, nplt);
[X, Y] = meshgrid(xx, xx);
tplt = [X(:).'; Y(:).'; 0*X(:).'];
ext = ~surferinterior(S, tplt);
sigma = dens(x);

tic;
Etot = NaN(3, nplt^2);
EHs = reshape(surferkerneval(S, Kfo, sigma, tplt(:,ext), eps), 6, []);
Etot(:,ext) = EHs(1:3,:) + em3d.planewave(zk, [pi/3, 0], [1; 0.2], tplt(:,ext));
% The PEC total electric field is normal on the surface. The Debye
% source r is not the physical surface charge, so evaluate the trace.
EHsurf = reshape(surfermatapply(S, Kfld, sigma, eps, ovf, Cf), 6, []);
rho = sum(S.n.*(EHsurf(1:3,:) + Ein));
fprintf('field evaluation: %.1fs\n', toc);

figure(5); clf
subplot(1,2,1)
pcolor(X, Y, reshape(real(Etot(1,:)), nplt, nplt)); shading interp; hold on
plot(S, real(S.n(1,:).*rho)); hold off
axis equal tight; colorbar
title('$\Re E_x$', 'Interpreter', 'latex')
set(gca, 'fontsize', 14)
subplot(1,2,2)
pcolor(X, Y, reshape(vecnorm(Etot), nplt, nplt)); shading interp; hold on
plot(S, abs(rho)); hold off
axis equal tight; colorbar
title('$|E|$', 'Interpreter', 'latex')
set(gca, 'fontsize', 14)


function pot = debye_fmm(eps, zk, s, t, sigma)
% [E;H] from [J;r;q] with one emfmm3d call, E and H as two density sets
ns = size(s.r(:,:), 2);  nt = size(t.r(:,:), 2);
sig = reshape(sigma, 5, ns, []);  nds = size(sig, 3);
J = sig(1:3,:,:);  M = cross(repmat(s.n(:,:), 1, 1, nds), J, 1);
src = [];
src.sources = s.r(:,:);
src.nd = 2*nds;
src.h_current = permute(cat(3, -M, J), [3 1 2]);
src.e_current = permute(cat(3, 1i*zk*J, 1i*zk*M), [3 1 2]);
src.e_charge  = permute(cat(3, -sig(4,:,:), -sig(5,:,:)), [3 2 1]);
U = emfmm3d(eps, zk, src, t.r(:,:), 1, 0, 0);
E = reshape(U.E, 2*nds, 3, nt);
pot = reshape(cat(1, permute(E(1:nds,:,:), [2 3 1]), ...
    permute(E(nds+1:end,:,:), [2 3 1])), 6*nt, nds);
end

function pot = div_fmm(eps, s, t, sigma)
% S_0[div e] as the potential of dipoles -e, with the real and imaginary
% parts as two density sets
ns = size(s.r(:,:), 2);  nt = size(t.r(:,:), 2);
e = reshape(sigma, 3, ns, []);  nds = size(e, 3);
src = [];
src.sources = s.r(:,:);
src.nd = 2*nds;
src.dipoles = -permute(cat(3, real(e), imag(e)), [3 1 2]);
U = lfmm3d(eps, src, 0, t.r(:,:), 1);
p = reshape(U.pottarg, 2*nds, nt);
pot = (p(1:nds,:) + 1i*p(nds+1:end,:)).';
end
