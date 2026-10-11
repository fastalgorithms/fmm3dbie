function p = probe_ptinfo()
%PROBE_PTINFO   Single random point used to probe function handles for
%   their output dimensions. Carries every per-point surfer field, so a
%   handle may use any of them.
p = [];
for f = {'r','n','du','dv','dru','drv','d','d2'}
    p.(f{1}) = randn(3,1);
end
p.uvs_targ = rand(2,1);
for f = {'wts','mean_curv'}
    p.(f{1}) = rand(1,1);
end
p.patch_id = 1;
% fundamental forms, stored on a surfer as one (2,2,npts) cell per patch
I  = [dot(p.du,p.du) dot(p.du,p.dv); dot(p.du,p.dv) dot(p.dv,p.dv)];
II = randn(2); II = (II + II.')/2;
p.ffform = {I}; p.sfform = {II}; p.ffforminv = {inv(I)};
end
