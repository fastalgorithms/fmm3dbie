function U = fmmgreen(eps, zk, srcuse, pg, targ, pgt)
%FLEX2D.FMMGREEN   FMM for flex2d.green with the interface of HFMM3D, for
% KERNEL3D.RADPP. Charges only, pg = 0; Hessians ordered xx yy zz xy xz yz.
%
% Syntax: U = flex2d.fmmgreen(eps, zk, srcuse, pg, targ, pgt)
%
% see also FLEX2D.FMM

if ( pg ~= 0 || isfield(srcuse, 'dipoles') )
    error('FLEX2D:fmmgreen:pg', 'Only charges and target outputs.');
end
nd = 1;
if ( isfield(srcuse, 'nd') ), nd = srcuse.nd; end
ns = size(srcuse.sources, 2);
nt = size(targ, 2);
q = reshape(srcuse.charges, nd, ns);
s = struct('r', srcuse.sources(1:2,:));
t = targ(1:2,:);

U = [];
U.pottarg = zeros(nd, nt);
if ( pgt > 1 ), U.gradtarg = zeros(nd, 3, nt); end
if ( pgt > 2 ), U.hesstarg = zeros(nd, 6, nt); end
for j = 1:nd
    out = cell(1, pgt);
    [out{:}] = flex2d.fmm(eps, zk, s, t, 'green', q(j,:).');
    U.pottarg(j,:) = out{1};
    if ( pgt > 1 ), U.gradtarg(j,1:2,:) = reshape(out{2}, 1, 2, nt); end
    if ( pgt > 2 ), U.hesstarg(j,[1 4 2],:) = reshape(out{3}, 1, 3, nt); end
end
if ( nd == 1 )
    U.pottarg = reshape(U.pottarg, 1, nt);
    if ( pgt > 1 ), U.gradtarg = reshape(U.gradtarg, 3, nt); end
    if ( pgt > 2 ), U.hesstarg = reshape(U.hesstarg, 6, nt); end
end

end
