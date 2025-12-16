% Two spheres test to get 4 digits of accuracy in the solution
a = 1;
d = 0.05;
c0 = [2*a + d; 0; 0];

norder = 6;
na = 3;
naref = na*2;
S1 = geometries.sphere(a, na, [], norder, 1);
S2 = geometries.sphere(a, na, c0, norder, 1);

S = merge([S1, S2]);
rhs = [ones(S1.npts,1); zeros(S2.npts,1)];

S1ref = geometries.sphere(a, naref, [], norder, 1);
S2ref = geometries.sphere(a, naref, c0, norder, 1);

Sref = merge([S1ref, S2ref]);
rhs_ref = [ones(S1ref.npts,1); zeros(S2ref.npts,1)];
plot(S);

% error estimate based on curvature
erra = norm((S.mean_curv - a).*sqrt(S.wts))/S.area;
fprintf('error in resolving curvature = %d\n', erra);

%%
rep_pars = [1.0; 1.0];
[densities, errs] = lap3d.solver(S, 'dir', rhs, 1e-8, rep_pars);

%%
figure(2)
clf
plot(S, densities)

figure(3)
clf
plot(S, rhs)

%%
dcoefs = vals2coefs(S, densities);
nppatch = S.npts/S.npatches;
dcoefs = reshape(dcoefs, [nppatch, S.npatches]);
dsub = sum(abs(dcoefs(end-norder+1:end,:)),1);

figure(4);
plot(S, dsub.')
%%
[densities_ref, errs] = lap3d.solver(Sref, 'dir', rhs_ref, 1e-8, rep_pars);


dcoefs_ref = vals2coefs(Sref, densities_ref);
dcoefs_ref = reshape(dcoefs_ref, [nppatch, Sref.npatches]);
dsub_ref = sum(abs(dcoefs_ref(end-norder+1:end,:)),1);

figure(5);
plot(Sref, dsub_ref.')

%%
x = S1.r(1,:).';
y = S1.r(2,:).';
z = S1.r(3,:).';
c1 = sum(densities(1:S1.npts).*S1.wts.*x);

xref = S1ref.r(1,:).';
yref = S1ref.r(2,:).';
c1_ref = sum(densities_ref(1:S1ref.npts).*S1ref.wts.*xref);

err1 = abs(c1 - c1_ref)./(abs(c1));
fprintf('error in polarization charge on obstacle 1=%d\n', err1);
