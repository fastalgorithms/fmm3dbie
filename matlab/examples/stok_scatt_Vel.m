% example of viscous flow around an ellipsoid with no-slip boundary
% condition (called "velocity" boundary condition, where total 
% velocity = 0) in the exterior "scattering" setting 
% (flow prescribed at infinity, scattered flow decays, total flow satisfies
% the boundary condition)
%
% a good representation for this problem is a combined layer potential,
% which is what the built in solver will use 
%
% Askham 1/14/25, derived from Barnett's helm_scatt_Dir.m

% velocity at infinity 
velinf = [1,0.5,0].'; 

% geometry discretization

ord = 6;          % panel order
refine = 1;       % integer refinement factor (increase 1 to 2 to check conv)
iptype = 11;      % 11 selects quad patches with product GL nodes 
S = geometries.ellipsoid([1.5,0.5,1.0], ceil(refine*[5 2 3]), [], ord, ...
    iptype);

% for visualization
% note: we set up a background flow field that has zero z component. 
% because of domain symmetry, the flow field on the slice z=0 will have 
% zero z component. we can visualize the flow on this slice as streamlines
% that are essentially in 2D

dx = 0.03; L=3; g = -L:dx:L; n = numel(g); z0=0.0; [xx yy] = meshgrid(g,g);
sheet = []; sheet.r = [xx(:) yy(:) z0+0*xx(:)]';  % sheet of test points (3*n^2)

test_pt = []; test_pt.r = [2;2;2];

% find the points on sheet that are outside ellipsoid 
npts = S.npts;
tol0 = 1e-4;
Dlap = lap3d.eval(S,'double',ones(npts,1),sheet,tol0);
out = abs(Dlap) < abs(Dlap+1);
sheet.r = sheet.r(:,out);

% set up and solve BVP
tol = 1e-4;
rhs = repmat(-velinf,1,npts);

fprintf('solving...\n'); 
start = tic; [dens,gmres_rres] = stok3d.solver(S, 'vel', rhs, tol);
t1 = toc(start); fprintf('%5.2e s : solve time\n',t1);
fprintf('eval soln...\n'); 
uscat = nan(3,n,n);
start = tic; uscat(:,reshape(out,size(xx))) = stok3d.eval(S, 'vel', dens, sheet, tol);
t1 = toc(start); fprintf('%5.2e s : eval time\n',t1);
start = tic; uscat_test = stok3d.eval(S, 'vel', dens, test_pt, tol);
t1 = toc(start); fprintf('%5.2e s : eval time\n',t1);

%

figure(1); clf
utot = uscat + velinf;
[sx,sy] = meshgrid(g(2),g(2:5:end));
streamline(g,g,reshape(utot(1,:),size(xx)),reshape(utot(2,:),size(xx)),sx,sy);
axis equal
hold on
plot(S)
title('streamlines')
fprintf('%5.2e : max of z component of velocity on this slice\n', ...
    max(abs(utot(3,:))));

%

figure(2); clf;
errs = surf_fun_error(S,dens);
plot(S,errs(1,:))
colorbar
title('error estimate for density')

%

figure(3); clf 
semilogy(gmres_rres,'b-x')
title('gmres relative residual')
