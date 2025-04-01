nch = 16;

S = geometries.disk([],[],[4 4 nch/4],8);

cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) circle(t),nch,cparams);
chnkr = sort(chnkr);

targinfo.r = [chnkr.r; 0*chnkr.r(1,:,:)];
[sxyz, uvsloc, dists, flags] = get_closest_pts(S,targinfo) ;

figure(1)
plot(S,rand(S.npatches,1))
hold on
plot3(targinfo.r(1,:),targinfo.r(2,:),targinfo.r(3,:),'x-')

figure(2)
plot(targinfo.r(1,:),targinfo.r(2,:),'x-')
hold on
plot(sxyz(1,:),sxyz(2,:),'x-')
 
alpha = 1;
beta = 600+0.2i;
gamma = 1;
[rts,ejs] = polynya3d.find_roots(alpha,beta,gamma);
zpars = [rts;ejs];

eps = 1e-6;
tic
A = polynya3d.matgen(S, zpars, eps);
toc

%%

f = exp(-10*((S.r(1,:) - 0.1).^2 + (S.r(2,:) - 0.1).^2 + S.r(3,:).^2 ) / 0.1).';

figure(3)
errs = surf_fun_error(S, f); plot(S, log(errs));
colorbar

figure(4)
plot(S,f); shading interp
colorbar

u = A.'*f;

figure(5)
plot(S,real(u)); shading interp
colorbar

%%

load('gong.mat')
sound(y)