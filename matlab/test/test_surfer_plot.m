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

B2 = surfer.sphere(4, 1, 1, 11);
figure
clf
plot(B2)

B2 = surfer.sphere(4, 1, 1, 12);
figure
clf
plot(B2, B2.r(1,:));