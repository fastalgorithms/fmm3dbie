fname = '../../geometries/meshes/simple_torus.msh';
opts = [];
opts.nrefine = 1;
Sall = multiscale_mesher(fname, 5, opts);
% Pick out most refined surface
S = Sall{opts.nrefine+1};

% set data to be 1 on a collection of patches defined by ipatch_ind
rhs = zeros(S.npts,1);
ipatch_ind = [1:15 100:115];
isrc_ind = [];
for i = ipatch_ind
    isrc_ind = [isrc_ind S.ixyzs(i):(S.ixyzs(i+1)-1)];
end
rhs(isrc_ind) = 1;

% solve
zk = 1.1;
eps = 1e-7;
[dens, errs, rres] = solver(S, 'helm', 'neu', rhs, eps, zk);

% post process fields at 100 random targets
targinfo = [];
ntarg = 100;
targinfo.r = rand(3,ntarg);
targinfo.r(3,:) = targinfo.r(3,:) + 10;
u = eval_fields(S, 'helm', 'neu', dens, targinfo, eps, zk);
