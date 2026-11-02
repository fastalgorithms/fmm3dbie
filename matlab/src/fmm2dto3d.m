function U = fmm2dto3d(fmm, eps, srcuse, pg, targ, pgt)
%FMM2DTO3D converts a 2d fmm into a 3d format
%
% Syntax: U = fmm2dto3d(fmm, eps, srcuse, pg, targ, pgt)

if ( pg ~= 0 || isfield(srcuse, 'dipoles') )
    error('FMM2DTO3D:pg', 'Only charges and target outputs.');
end
nd = 1;
if ( isfield(srcuse, 'nd') ), nd = srcuse.nd; end
nt = size(targ, 2);
srcuse.sources = srcuse.sources(1:2,:);

V = fmm(eps, srcuse, 0, targ(1:2,:), pgt);

U = [];
U.pottarg = reshape(V.pottarg, nd, nt);
if ( pgt > 1 )
    U.gradtarg = zeros(nd, 3, nt, 'like', V.gradtarg);
    U.gradtarg(:,1:2,:) = reshape(V.gradtarg, nd, 2, nt);
end
if ( pgt > 2 )
    U.hesstarg = zeros(nd, 6, nt, 'like', V.hesstarg);
    U.hesstarg(:,[1 4 2],:) = reshape(V.hesstarg, nd, 3, nt);
end

end
