%
% Level restriction in resolve_fun on a rotated rounded cube.
%
% The surface is the superellipsoid |x|^q + |y|^q + |z|^q = 1, rotated so
% that its edges cut through the base patches of the cube-sphere instead of
% lying on their boundaries. resolve_fun then refines thin strips along the
% edges next to large, unsplit face patches. The faces are flat, so r, n and
% the (smooth, harmonic) Dirichlet data are resolved there at any tolerance,
% but the double-layer density varies near the edges and is not.
%
% For each case, resolve_fun refines on [r; n; data] with (lvlr = 'a') or
% without (lvlr = 'n') level restriction, the interior Dirichlet problem is
% solved with the Laplace double layer (surfermatapply + gmres), and the
% potential is compared with the exact solution at interior points near
% the surface, near the edges and away from them.
%
% Expect the unrestricted error to stall near 1e-4 as the tolerance
% tightens, while the level-restricted mesh reaches ~1e-6 at about the
% cost of the unrestricted 1e-5 run.
%
run ../startup.m

%% Parameters

q         = 200;            % edge rounding radius ~ 1.4/q
angles    = [0.3, 0.5, 0.7];% rotation about x, y, z
norder    = 20;             % order of the initial surface
na        = 2;              % base patches per cube-sphere edge
iptype    = 1;
geom_tol  = 1e-8;           % order-norder resolution of the exact surface
nordcheck = 8;              % solve order
nrefmax   = 15;
quad_tol  = 1e-10;
gmres_tol = 1e-12;
depths    = [0.02, 0.005];  % target depths below the surface
ntarg     = 1000;

% (lvlr, resolve tolerance) pairs
cases = {'n', 1e-4; 'a', 1e-4; 'n', 1e-5};

%% Geometry

Rx = [1 0 0; 0 cos(angles(1)) -sin(angles(1)); 0 sin(angles(1)) cos(angles(1))];
Ry = [cos(angles(2)) 0 sin(angles(2)); 0 1 0; -sin(angles(2)) 0 cos(angles(2))];
Rz = [cos(angles(3)) -sin(angles(3)) 0; sin(angles(3)) cos(angles(3)) 0; 0 0 1];
R = Rz*Ry*Rx;

tic
S0 = rounded_cube(R,q,na,norder,iptype,geom_tol);
fprintf('initial order-%d surface: %d patches (%.1f s)\n',norder,S0.npatches,toc);

% Exact interior solution (harmonic) and the data resolve_fun refines on.
uex = @(r) (r(1,:) + 0.5*r(2,:).*r(3,:) + 0.3*(r(1,:).^2 - r(3,:).^2)).';
f = @(T) [T.r; T.n; uex(T.r).'];

%% Targets: points on the exact surface pushed inward along the normal

rng(0);
p = randn(3,ntarg); p = p./vecnorm(p);
x = cube_map(R*p,0*p,0*p,q);
nrm = x.^(q-1); nrm = nrm./vecnorm(nrm);
s2 = sort(abs(x),1,'descend');
isedge = s2(2,:) > 0.9;           % within about 0.1 of an edge
targ = zeros(3,0); isedget = false(1,0);
for d = depths
    targ = [targ, R.'*(x - d*nrm)]; %#ok<AGROW>
    isedget = [isedget, isedge]; %#ok<AGROW>
end
targs = struct('r',targ);
utarg = uex(targs.r);

%% Solve

dkern = kernel3d('l','d');
ncase = size(cases,1);
res = struct('lvlr',cases(:,1),'tol',cases(:,2));

figure(1); clf
for k = 1:ncase
    lvlr = cases{k,1};
    tol = cases{k,2};

    tic
    S = resolve_fun(S0,tol,f,nordcheck,nrefmax,Inf,struct('lvlr',lvlr));
    res(k).tres = toc;
    res(k).npatches = S.npatches;

    S = oversample(S,nordcheck);
    res(k).npts = S.npts;

    subplot(1,ncase,k)
    wireframe(S,struct('wfill',1))
    axis equal off
    title(sprintf('lvlr=%s, tol %.0e: %d patches',lvlr,tol,S.npatches));
    drawnow

    tic
    [cors,objover] = surfermat(S,dkern,quad_tol,struct('corrections',true));
    res(k).tcor = toc;
    res(k).nnz = nnz(cors);

    % Interior limit is D_pv + jump*I, with D_pv[1] = jump = -/+ 1/2.
    Dapply = @(sig) surfermatapply(S,dkern,sig,quad_tol,objover,cors);
    jump = sign(median(Dapply(ones(S.npts,1))))/2;

    tic
    [sigma,~,~,iter] = gmres(@(sig) Dapply(sig) + jump*sig,uex(S.r),[], ...
        gmres_tol,200);
    res(k).tsolve = toc;
    res(k).iter = iter(end);

    u = surferkerneval(S,dkern,sigma,targs,quad_tol);
    err = abs(u(:) - utarg)/max(abs(utarg));
    res(k).err_face = max(err(~isedget));
    res(k).err_edge = max(err(isedget));

    fprintf(['lvlr=%s, tol %.0e: %d patches, %d pts, nnz %.3g, %d iters, ' ...
        'err faces %.2e edges %.2e (%.0f s)\n'],lvlr,tol,res(k).npatches, ...
        res(k).npts,res(k).nnz,res(k).iter,res(k).err_face,res(k).err_edge, ...
        res(k).tres + res(k).tcor + res(k).tsolve);
    clear S cors objover Dapply sigma
end

%% Summary

fprintf('\n%-5s %-6s %8s %9s %10s %8s | %9s %9s\n','lvlr','tol', ...
    'patches','npts','nnz','time s','err face','err edge');
for k = 1:ncase
    r = res(k);
    fprintf('%-5s %-6.0e %8d %9d %10.3g %8.0f | %9.2e %9.2e\n',r.lvlr,r.tol, ...
        r.npatches,r.npts,r.nnz,r.tres + r.tcor + r.tsolve,r.err_face,r.err_edge);
end


function S = rounded_cube(R,q,na,norder,iptype,geom_tol)
% Rotated superellipsoid as a radial map of the unit sphere. Sphere
% patches are split until the exact surface is resolved at order norder.
Ssph = geometries.sphere(1,na,[0;0;0],norder,iptype);
for iter = 1:40
    [x,du,dv] = cube_map(R*Ssph.r,R*Ssph.du,R*Ssph.dv,q);
    x = R.'*x; du = R.'*du; dv = R.'*dv;
    n = cross(du,dv,1);
    n = n./vecnorm(n);
    n = n.*sign(sum(n.*Ssph.n,1));
    S = surfer(Ssph.npatches,norder,[x;du;dv;n],iptype);
    err = surf_fun_error(S,[S.r; S.n],Inf);
    ibad = find(max(err,[],1) > geom_tol);
    if isempty(ibad), return, end
    Ssph = split_patches(Ssph,ibad);
end
error('rounded_cube: surface not resolved to %g at order %d',geom_tol,norder);
end


function [x,du,dv] = cube_map(p,dpu,dpv,q)
% Radial projection of the unit sphere onto sum |x_i|^q = 1 (q even), with
% the tangents carried along by the chain rule.
N = sum(p.^q,1).^(1/q);
gN = p.^(q-1)./N.^(q-1);
x = p./N;
du = dpu./N - p.*sum(gN.*dpu,1)./N.^2;
dv = dpv./N - p.*sum(gN.*dpv,1)./N.^2;
end
