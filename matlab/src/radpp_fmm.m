function pot = radpp_fmm(fmmfun,p,eps,srcinfo,targinfo,sigma)
%RADPP_FMM   Fast multipole method for a radpp or radgreen kernel, from
%   an FMM for its Green's function G.
%
% Syntax:
%   pot = radpp_fmm(fmmfun, p, eps, srcinfo, targinfo, sigma)
%
% Input:
%   fmmfun   - FMM for G with the interface of HFMM3D without the
%              wavenumber, U = fmmfun(eps, srcuse, pg, targ, pgt), with
%              U.pottarg(i) = sum_j c_j G(t_i,y_j) + v_j . grad_y G(t_i,y_j),
%              e.g. @(eps,s,pg,t,pgt) hfmm3d(eps, zk, s, pg, t, pgt).
%              Kernels with second derivatives need U.hesstarg.
%   p        - the kernel parameters, K.params
%   eps      - precision requested
%   srcinfo  - ptinfo struct: .r (3,:), .n (3,:)
%   targinfo - ptinfo struct: .r (3,:) and the fields the kernel needs
%   sigma    - density (opdims(2)*ns x ndens), already scaled by
%              quadrature weights
%
% Output:
%   pot  - potential (opdims(1)*nt x ndens)

m = p.opdims(1);  n = p.opdims(2);
ns = size(srcinfo.r(:,:), 2);
ndens = numel(sigma)/(n*ns);
sigma = reshape(sigma,n,ns,ndens);

dim = p.dim;
nt = size(targinfo.r(:,:),2);
targ = zeros(3,nt);
targ(1:size(targinfo.r,1),:) = targinfo.r(:,:);

% the weights of each entry on the derivatives of the charge and dipole
% fields, (nt,ncomp,nent)
W = p.wfun(targinfo);
V = [];
if p.isrc, V = p.vfun(targinfo); end
nent = size(p.ent,1);

% the density components with charges and with dipoles
isq = false(n,1);  isd = false(n,1);
for e = 1:nent
    l = p.ent(e,2);
    isq(l) = isq(l) || nnz(W(:,:,e)) > 0;
    if p.isrc, isd(l) = isd(l) || nnz(V(:,:,e)) > 0; end
end
lq = find(isq);  ld = find(isd);
nq = numel(lq);  ndip = numel(ld);
nd = (nq + ndip)*ndens;
if nd == 0
    pot = zeros(m*nt,ndens);
    return
end

% the fmm densities: charges sigma_l, then dipoles sigma_l n_y, for each
% density
srcuse = [];
srcuse.sources = zeros(3,ns);
srcuse.sources(1:size(srcinfo.r,1),:) = srcinfo.r(:,:);
srcuse.nd = nd;
srcuse.charges = zeros(nq + ndip,ndens,ns);
srcuse.charges(1:nq,:,:) = permute(sigma(lq,:,:),[1 3 2]);
srcuse.charges = reshape(srcuse.charges,nd,ns);
if ndip > 0
    sd = zeros(nq + ndip,ndens,ns);
    sd(nq+1:end,:,:) = permute(sigma(ld,:,:),[1 3 2]);
    ny = zeros(3,ns);
    ny(1:size(srcinfo.n,1),:) = srcinfo.n(:,:);
    srcuse.dipoles = reshape(sd,nd,1,ns).*reshape(ny,1,3,ns);
end

pgt = p.nderiv + 1;
U = fmmfun(eps,srcuse,0,targ,pgt);

% the field derivatives in the order of RADPP_ROTATE, (nt,ncomp,nd), from
% the hessian ordered xx, yy, zz, xy, xz, yz
ih = [1 4 5 2 6 3];
if dim == 2, ih = [1 4 2]; end
F = zeros(nt,p.ncomp,nd);
F(:,1,:) = reshape(reshape(U.pottarg,nd,nt).',nt,1,nd);
if pgt >= 2
    g = permute(reshape(U.gradtarg,nd,3,nt),[3 2 1]);
    F(:,2:1+dim,:) = g(:,1:dim,:);
end
if pgt >= 3
    h = permute(reshape(U.hesstarg,nd,6,nt),[3 2 1]);
    F(:,2+dim:end,:) = h(:,ih,:);
end
F = reshape(F,nt,p.ncomp,nq + ndip,ndens);

pot = zeros(m*nt,ndens);
for e = 1:nent
    i = p.ent(e,1);  l = p.ent(e,2);
    k = find(lq == l);
    if ~isempty(k)
        pot(i:m:end,:) = pot(i:m:end,:) + ...
            reshape(sum(W(:,:,e).*F(:,:,k,:),2),nt,ndens);
    end
    k = find(ld == l);
    if ~isempty(k)
        pot(i:m:end,:) = pot(i:m:end,:) + ...
            reshape(sum(V(:,:,e).*F(:,:,nq+k,:),2),nt,ndens);
    end
end

end
