% 

nch = 24;

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

max(dists) % max distance between pt on chnkr and its nearest pt on patch

us = uvsloc(1,:);
vs = uvsloc(2,:);
d1 = sqrt(us.^2+vs.^2);
d2 = sqrt((us-1).^2+vs.^2);
d3 = sqrt(us.^2+(vs-1).^2);
min([d1 d2 d3]) % closest distance to the vertex of a triangle

