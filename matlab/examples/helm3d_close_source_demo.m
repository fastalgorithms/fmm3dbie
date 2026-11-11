% Exterior Helmholtz Dirichlet problem on a stellarator with a close source.
% Assumes pwd is the directory this script is in.

run('../startup.m');

%% Geometry and a source just inside the boundary

zk = 1;
norder = 20;
nordcheck = 6;
resolve_tol = 1e-6;
nrefmax = 12;
source_gap = 0.01;

S = geometries.stellarator([5,15],norder,11);
npatches0 = S.npatches;
[~,inode] = max(S.r(1,:));
anchor = S.r(:,inode);
normal = S.n(:,inode);
source = [];
source.r = anchor - source_gap*normal;
rhsfun = @(T) helm3d.kern(zk,source,T,'s').';
f = @(T) [T.r; rhsfun(T)];

%% Resolve the surface and Dirichlet data through order 6, then balance

tstart = tic;
S = resolve_fun(S,resolve_tol,f,nordcheck,nrefmax,Inf);
err = surf_fun_error(S,f(S),Inf,nordcheck);
rerr = max(err(1:3,:),[],'all');
rhserr = max(err(4,:));
assert(max(err(:)) < resolve_tol, 'Surface or Dirichlet data are not resolved');
fprintf('Order-6 tails: coordinates %.3e, rhs %.3e\n',rerr,rhserr);
S = oversample(S,nordcheck);
rhs = rhsfun(S).';
fprintf('Geometry: %d -> %d patches, %d solver nodes (%.1f s)\n', ...
    npatches0,S.npatches,S.npts,toc(tstart));

figure(1); clf
subplot(1,2,1)
plot(S,log10(abs(rhs))); colorbar
hold on
wireframe(S)
plot3(source.r(1),source.r(2),source.r(3),'r*');
hold off
axis equal tight
title('log_{10}|rhs|');
subplot(1,2,2)
plot(S,log10(surf_fun_error(S,rhs)))
colorbar
hold on
wireframe(S)
hold off
axis equal tight
title('log_{10} rhs tails');

%% Solve with a combined single and double layer representation

target_offsets = source_gap*[1;2;5;10;50];
targets = [];
targets.r = anchor + normal*target_offsets.';
inside = surferinterior(S,[source.r,targets.r]);
assert(inside(1) && ~any(inside(2:end)), ...
    'The source must be interior and the check targets exterior');

quad_tol = 1e-7;
rep_pars = [-1i*zk,1];
solve_opts = [];
solve_opts.ifout = 1;
solve_opts.eps_gmres = 1e-9;
tstart = tic;
[sigma,errs,rres] = helm3d.dirichlet.solver(S,rhs,quad_tol,zk,rep_pars,solve_opts);
fprintf('GMRES: %d iterations, relative residual %.3e (%.1f s)\n', ...
    numel(errs),rres,toc(tstart));
assert(rres < 10*solve_opts.eps_gmres, 'GMRES did not converge');

%% Compare with the exact point-source field, including close targets

u = helm3d.dirichlet.eval(S,sigma,targets,quad_tol,zk,rep_pars);
uex = helm3d.kern(zk,source,targets,'s');
relerr = abs(u-uex)./abs(uex);
fprintf('Worst relative field error: %.3e\n',max(relerr));

figure(2); clf
subplot(1,2,1)
plot(S,log10(abs(sigma))); colorbar
hold on
wireframe(S)
hold off
axis equal tight
title('log_{10}|density|');
subplot(1,2,2)
loglog(target_offsets,relerr,'o-'); grid on
xlabel('Distance to the boundary');
ylabel('Relative error');
