addpath(genpath('~/git/fmm3dbie/matlab'))
close('all')
B = surfer.load_from_file('~/git/fmm3dbie/geometries/sphere_768_o03.go3');
figure
clf
tic, plot(B); toc;

return
[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(B);

novers = randi([4,6],B.npatches,1);
[Bover,xmat] = oversample(B,novers);
figure
clf
plot(Bover);
