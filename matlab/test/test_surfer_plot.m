run ../startup.m
close('all')
B = surfer.load_from_file('~/git/fmm3dbie/geometries/sphere_768_o03.go3');
figure
clf
tic, plot(B); toc;

[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(B);

novers = randi([4,6],B.npatches,1);
[Bover,xmat] = oversample(B,novers);
figure
clf
plot(Bover);

B2 = geometries.sphere(1, 4, [0;0;0], 6, 11);
figure
clf
plot(B2)

B2 = geometries.sphere(1, 4, [0;0;0], 6, 12);
figure
clf
plot(B2, B2.r(1,:));