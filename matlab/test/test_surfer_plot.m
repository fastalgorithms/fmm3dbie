% tester for surfer plotting.
% assumes pwd is the directory this script is in

run ../startup.m
close('all')
B = surfer.load_from_file('../../geometries/sphere_768_o03.go3');
figure
clf
tic, plot(B); toc;
colorbar

[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(B);

novers = randi([4,6],B.npatches,1);
[Bover,xmat] = oversample(B,novers);
figure
clf
plot(Bover);
colorbar

B2 = geometries.sphere(1, 4, [0;0;0], 6, 11);
figure
clf
plot(B2)
colorbar

B2 = geometries.sphere(1, 4, [0;0;0], 6, 12);
figure
clf
plot(B2, B2.r(1,:));
colorbar
