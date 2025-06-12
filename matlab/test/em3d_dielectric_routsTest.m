%
% This file tests the Maxwell pec 
%
%
run ../startup.m
S = geometries.sphere(1, 12, [0;0;0], 6, 1);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

om = 1.1;

ep0 = 1.4;
mu0 = 1;
ep1 = 1;
mu1 = 1;
rep_params = complex([ep0, mu0, ep1, mu1]);

zkin = om*sqrt(ep1*mu1);
zkout = om*sqrt(ep0*mu0);


eps = 1e-8;

xyz_in = [0.17; 0.23; -0.11];
xyz_out = [3.1; 1.0; 1.0];
src_info = [];
src_info.r = xyz_out;
src_info.edips = rand(3,1) + 1j*rand(3,1);
src_info.hdips = rand(3,1) + 1j*rand(3,1);


[Einc, Hinc] = em3d.incoming_sources(zkin, src_info, S, 'ehd');

%
opts = [];
opts.eps_gmres = 1e-8;
[densities, errs, rres, Q] = em3d.dielectric.solver(S, Einc, Hinc, eps, om, rep_params, opts);

%%
%
targ_info = [];
targ_info.r = xyz_in;

[E, H] = em3d.dielectric.eval(S, densities, targ_info, eps, om, rep_params);


[E_ex, H_ex] = em3d.incoming_sources(zkin, src_info, targ_info, 'ehd');

u_comp = [E; H];
u_ex = [E_ex; H_ex];
utest = u_comp + u_ex;


fprintf('Error in fields iterative solver=%d\n',norm(utest(:))/norm(u_ex(:)));

%% Now post process on surface

au = densities(:,1);
av = densities(:,2);
bu = densities(:,3);
bv = densities(:,4);

avec = au.'.*S.dru + av.'.*S.drv;
bvec = bu.'.*S.dru + bv.'.*S.drv;

[adiv, au2, av2] = get_surf_div_ortho(S, avec);
[bdiv, bu2, bv2] = get_surf_div_ortho(S, bvec);


%%

j1 = bvec*sqrt(ep1);
rho1 = bdiv*sqrt(ep1)/1j/zkin;

d1 = complex(zeros(4,S.npts));
d1(1:3,:) = j1;
d1(4,:) = rho1.';

j2 = avec*sqrt(mu1);
rho2 = adiv*sqrt(mu1)/1j/zkin;
d2 = complex(zeros(4,S.npts));
d2(1:3,:) = j2;
d2(4,:) = rho2.';
%%


[E11, H11] = em3d.pec.eval(S, d1, targ_info, eps, zkin, 1);
[E12, H12] = em3d.pec.eval(S, d2, targ_info, eps, zkin, 1);

E2 = sqrt(mu1)*(E11 + H12);
H2 = sqrt(ep1)*(H11 - E12);


[E11, H11] = em3d.pec.eval(S, d1, S, eps, zkin, 1);
[E12, H12] = em3d.pec.eval(S, d2, S, eps, zkin, 1);

H11 = H11 + cross(S.n, j1, 1)/2;
H12 = H12 + cross(S.n, j2, 1)/2;
E11 = E11 - rho1.'.*S.n/2;
E12 = E12 - rho2.'.*S.n/2;

E2surf = sqrt(mu1)*(E11 + H12);
H2surf = sqrt(ep1)*(H11 - E12);

[E_ex, H_ex] = em3d.incoming_sources(zkin, src_info, S, 'ehd');

function [jdiv, ju, jv] = get_surf_div_ortho(S, j)
    jproju = sum(j.*S.du,1).';
    jprojv = sum(j.*S.dv,1).';
    
    ju = zeros(size(jproju));
    jv = zeros(size(jprojv));
    g = zeros(size(jproju));
    for i=1:S.npatches
        a = S.ffforminv{i};
        a11 = squeeze(a(1,1,:));
        a12 = squeeze(a(1,2,:));
        a21 = squeeze(a(2,1,:));
        a22 = squeeze(a(2,2,:));
        b = S.ffform{i};
        b11 = squeeze(b(1,1,:));
        b12 = squeeze(b(1,2,:));
        b21 = squeeze(b(2,1,:));
        b22 = squeeze(b(2,2,:));

        iinds = S.ixyzs(i):S.ixyzs(i+1)-1;
        ju(iinds) = a11.*jproju(iinds) + a12.*jprojv(iinds);
        jv(iinds) = a21.*jproju(iinds) + a22.*jprojv(iinds);
        g(iinds) = sqrt(b11.*b22 - b12.*b21);
    end
    juuse = ju.*g;
    jvuse = jv.*g;

    juc = vals_to_coefs_surface_fun(S, juuse);
    jvc = vals_to_coefs_surface_fun(S, jvuse);
    uv = koorn.rv_nodes(S.norders(1));
    [~, dersu, dersv] = koorn.ders(S.norders(1), uv);

    [~, npols] = size(uv);
    juc = reshape(juc, [npols, S.npatches]);
    jvc = reshape(jvc, [npols, S.npatches]);

    ju_u = dersu.'*juc;
    jv_v = dersv.'*jvc;
    
    jdiv = ju_u + jv_v;
    jdiv = jdiv(:);
    jdiv = jdiv./g;

end