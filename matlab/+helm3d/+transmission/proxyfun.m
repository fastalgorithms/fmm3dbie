function [Kpxy,nbr] = proxyfun(x_all, slf_all, nbr_all, proxy_dict, l, ctr,zks,S)

% Convert system source indices back to original indices
jpts = idivide(int64(slf_all(:)-1), int64(2))+1;
[juni,~,ijuni] = unique(jpts);
% Convert system neighbour indices back to original indices
nbr_pts = idivide(int64(nbr_all(:)-1), int64(2))+1;
[nbr, ~, ~] = unique(nbr_pts);
slf = juni;

proxy = proxy_dict.proxy;
weigt = l*proxy_dict.weigt;
norms = proxy_dict.norms;
% Shift and scale precomputed proxy surface
proxy = bsxfun(@plus, proxy*l, ctr');

% Allocate a stacked proxy matrix of all kernels to compress
i = length(proxy);
j = length(slf);

nzks = numel(zks);
Kpxy = zeros(4*i*nzks, 2*j);

% Exterior kernels to compress
srcinfo = [];
srcinfo.r = S.r(:,slf);
srcinfo.n = S.n(:,slf);
targinfo = [];
targinfo.r = proxy;
targinfo.n = norms;

area = S.wts';

for iz=1:nzks
    z_k = zks(iz);

    Kmat = helm3d.kern(z_k,srcinfo,targinfo,'sprime');
    Kmat = bsxfun(@times, Kmat, sqrt(area(j)));
    Kmat = bsxfun(@times, weigt.', Kmat);


    D = helm3d.kern(z_k,srcinfo,targinfo,'d');

    D = bsxfun(@times, D, sqrt(area(j)));
    D = bsxfun(@times, weigt.', D);

    ioff = (iz-1)*4*i;
    iend = 4*i;
    Kpxy((2:4:iend)+ioff,1:2:end) = Kmat;
    Kpxy((1:4:iend)+ioff,2:2:end) = Kmat;
    Kpxy((4:4:iend)+ioff, 1:2:end) = D;
    Kpxy((3:4:iend)+ioff, 2:2:end) = D;
end
ctruse = ctr(:);
dxyz = abs(S.r(1:3,nbr)-ctruse(1:3))/l;
nbr = nbr(max(dxyz) < 2.5);



% select relevant columns 
ijuni2 = (ijuni-1)*2 + mod(slf_all(:)-1, 2)+1;
Kpxy = Kpxy(:, ijuni2);


dx = x_all(1, nbr_all) - ctr(1);
dy = x_all(2, nbr_all) - ctr(2);
dz = x_all(3, nbr_all) - ctr(3);
dist = max(abs([dx;dy;dz]));
nbr = nbr_all(dist/l < 2.5);


end
