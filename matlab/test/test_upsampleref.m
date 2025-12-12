close all
clear all
% Test interpolation from coarser to finer mesh
a = 1;
d = 0.05;

norder = 8;
na = 3;
naref = 2*na;

patch_ids = [1,7,30];
% patch_ids = [];
iptype = 1;
S1 = geometries.sphere(a, na, [], norder, iptype);
[srcvals1,srccoefs1,norders1,ixyzs1,iptype1,wts1] = ...
    extract_arrays(S1);

[Sover,ipatchid,uvstot,interp_mat,ivec,jvec,vvec] = split_patches(S1,patch_ids);

[srcvalsover,srccoefsover,nordersover,ixyzsover,iptypeover,wtsover] = ...
    extract_arrays(Sover);
%
idx = 5;
f_f = srcvalsover(idx,:)';
f_c = srcvals1(idx,:)';
A = get_interp_mat(S1,patch_ids);
error_interp = norm(A*f_c-f_f,inf);
