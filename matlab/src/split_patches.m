function [Sover, ipatchid, uvstot, interp_mat, ivec, jvec, vvec] = ...
    split_patches(S, patch_ids)
% SPLIT_PATCHES  Refine selected patches of a surface discretization.
%
% Syntax:
%   [Sover, ipatchid, uvstot, interp_mat] = split_patches(S, patch_ids)
%   [Sover, ipatchid, uvstot, interp_mat, ivec, jvec, vvec] = split_patches(S, patch_ids)
%
%   Refines a subset of patches in the surfer object S by subdividing each
%   patch listed in patch_ids into four subpatches. The routine constructs
%   the refined surface geometry, returns the uv-coordinates and patch
%   indices of all nodes after refinement, and builds a sparse interpolation
%   matrix mapping data from the original coarse patch layout to the refined
%   one.
%
% Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * patch_ids: indices of patches to refine
%
% Output arguments:
%    * Sover: refined surfer object with updated source geometry
%    * ipatchid: patch index of each evaluation node after refinement
%    * uvstot: uv-coordinates of all nodes on the refined surface
%    * interp_mat: sparse interpolation matrix such that f = interp_mat * c
%                maps values c on coarse patches to values f on the refined grid
%    * ivec,jvec,vvec: sparse matrix triplet representation of interp_mat such
%                    that interp_mat(ivec,jvec) = vvec
%
% Notes:
%   Assumes uniform input order
%   Supported interpolation types are determined by S.iptype:
%        * iptype = 1,  triangular patch discretized using
%                       Rokhlin-Vioreanu nodes
%        * iptype = 11, quadrangular patch discretized using tensor
%                       product Gauss-Legendre nodes
%        * iptype = 12, quadrangular patch discretized using tensor
%                       product Chebyshev

% get uniform patch order
norder = S.norders(1);
% number of patches
npatches = S.npatches;
% number of patches on refined surface
npatches_interp = npatches + 3*numel(patch_ids);
% get patch type
ip0 = S.iptype(1);
% get quadrature nodes and interpolation maps
if (ip0==1)
    % simplex vertices
    v0 = [0 1 0;0 0 1];
    % Vioreanu-Rokhlin nodes on the triangle simplex
    [uvs] = koorn.rv_nodes(norder);
    % number of nodes
    npols = (norder+1)*(norder+2)/2;
    % children vertices
    vcs = gettrichildren(v0);
    % map uvs nodes from simplex to children
    uvssub1 = mapuv_tri(vcs(:,:,1),uvs);
    uvssub2 = mapuv_tri(vcs(:,:,2),uvs);
    uvssub3 = mapuv_tri(vcs(:,:,3),uvs);
    uvssub4 = mapuv_tri(vcs(:,:,4),uvs);
    % map: coefficients in expansion on patch to values on subpatches
    c2vover1 = koorn.coefs2vals(norder,uvssub1);
    c2vover2 = koorn.coefs2vals(norder,uvssub2);
    c2vover3 = koorn.coefs2vals(norder,uvssub3);
    c2vover4 = koorn.coefs2vals(norder,uvssub4);
    % map: values on patch to coefficients in expansion on patch
    v2c = koorn.vals2coefs(norder,uvs);
    % map: coefficients in expansion on patch to values on patch
    c2v = koorn.coefs2vals(norder,uvs);
elseif(ip0 == 11)
    % simplex vertices
    v0 = [-1,1,-1;-1,-1,1];
    % tensor grid of Gauss-Legendre nodes on [-1,1]^2
    [uvs] = polytens.lege.nodes(norder);
    % number of nodes
    npols = (norder+1)^2;
    % children vertices      
    vcs = getquadchildren(v0);    
    % map uvs nodes from simplex to children
    uvssub1 = mapuv_quad(vcs(:,:,1),uvs);
    uvssub2 = mapuv_quad(vcs(:,:,2),uvs);
    uvssub3 = mapuv_quad(vcs(:,:,3),uvs);
    uvssub4 = mapuv_quad(vcs(:,:,4),uvs);
    % map: coefficients in expansion on patch to values on subpatches
    c2vover1 = polytens.lege.coefs2vals(norder,uvssub1);
    c2vover2 = polytens.lege.coefs2vals(norder,uvssub2);
    c2vover3 = polytens.lege.coefs2vals(norder,uvssub3);
    c2vover4 = polytens.lege.coefs2vals(norder,uvssub4);
    % map: values on patch to coefficients in expansion on patch
    v2c = polytens.lege.vals2coefs(norder,uvs);
    % map: coefficients in expansion on patch to values on patch
    c2v = polytens.lege.coefs2vals(norder,uvs);
elseif(ip0 == 12)   
    % simplex vertices
    v0 = [-1,1,-1;-1,-1,1];
    % tensor grid of Chebyshev nodes on [-1,1]^2
    [uvs] = polytens.cheb.nodes(norder);
    % number of nodes
    npols = (norder+1)^2;
    % children vertices      
    vcs = getquadchildren(v0);    
    % map uvs nodes from simplex to children
    uvssub1 = mapuv_quad(vcs(:,:,1),uvs);
    uvssub2 = mapuv_quad(vcs(:,:,2),uvs);
    uvssub3 = mapuv_quad(vcs(:,:,3),uvs);
    uvssub4 = mapuv_quad(vcs(:,:,4),uvs);
    % map: coefficients in expansion on patch to values on subpatches
    c2vover1 = polytens.cheb.coefs2vals(norder,uvssub1);
    c2vover2 = polytens.cheb.coefs2vals(norder,uvssub2);
    c2vover3 = polytens.cheb.coefs2vals(norder,uvssub3);
    c2vover4 = polytens.cheb.coefs2vals(norder,uvssub4);
    % map: values on patch to coefficients in expansion on patch
    v2c = polytens.cheb.vals2coefs(norder,uvs);
    % map: coefficients in expansion on patch to values on patch
    c2v = polytens.cheb.coefs2vals(norder,uvs);
else
    fprintf('Invalid type of patch, exiting\n');
    return
end

% store nodes for children
uvssub = [uvssub1 uvssub2 uvssub3 uvssub4];


% build sparse matrix


% identity block (maps a patch to itself)
id = eye(npols);

% interp_mat_blocks: 1×5 cell array of interpolation matrices, each of size
% npols × npols. Entry 1 is the identity mapping (patch to itself). Entries
% 2–5 map a patch to each of its four subpatches.
interp_mat_blocks = {id,c2vover1*v2c,c2vover2*v2c,c2vover3*v2c,c2vover4*v2c};

% determine the block structure of the global sparse interpolation matrix.
% interp_mat_list(k) specifies which block type (1–5) is used in row-block k.
% rows_idx and cols_idx record the corresponding block-row and block-column.
interp_mat_list = [];
rows_idx = [];
row_ctr = 0;
cols_idx = [];

uvstot = [];
ipatchid = [];

if nargin < 2
    patch_ids = [];
end

srcvals = cell(0);

for ii = 1:npatches
    if ismember(ii,patch_ids)
        % map the geometry of patch ii at the quadrature nodes to
        % of each of its four subpatches.
        panval1 = S.srccoefs{ii}(1:9,:)*c2vover1.';
        panval2 = S.srccoefs{ii}(1:9,:)*c2vover2.';
        panval3 = S.srccoefs{ii}(1:9,:)*c2vover3.';
        panval4 = S.srccoefs{ii}(1:9,:)*c2vover4.';
        srcvals = [srcvals {panval1,panval2,panval3,panval4}];
        % store local uv coordinates
        uvstot = [uvstot uvssub];
        % store id of patch ii new points on subpatches belong to
        ipatchid = [ipatchid; ii*ones(length(uvssub),1)];
        % refined patch: contributes four subpatch block-rows.
        interp_mat_list = [interp_mat_list [2 3 4 5]];
        cols_idx = [cols_idx ii*ones(1,4)];
        rows_idx = [rows_idx row_ctr+(1:4)];
        row_ctr = row_ctr + 4;
    else
        srcvals = [srcvals {S.srccoefs{ii}(1:9,:)*c2v.'}];
        uvstot = [uvstot uvs];
        ipatchid = [ipatchid; ii*ones(length(uvs),1)];
        % unrefined patch: contributes a single identity block-row.
        interp_mat_list = [interp_mat_list 1];
        cols_idx = [cols_idx ii];
        row_ctr = row_ctr + 1;
        rows_idx = [rows_idx row_ctr];
    end
end

srcmat = [];

% compute jacobians
jac = zeros(2,2,4);
for jj = 1:4
jac(:,:,jj) = [vcs(1,2,jj) - vcs(1,1,jj),vcs(1,3,jj) - vcs(1,1,jj);
            vcs(2,2,jj) - vcs(2,1,jj),vcs(2,3,jj) - vcs(2,1,jj)];
end

if (ip0 == 11) || (ip0 == 12)
    jac = jac /2;
end

for ii = 1:length(srcvals)
    tmp = zeros(12,npols);
    if interp_mat_list(ii) == 1
        tmp(1:9,:) = srcvals{ii}(1:9,:);
    else
        tmp(1:3,:) = srcvals{ii}(1:3, :);

        jj = interp_mat_list(ii)-1;

        tmp(4:6,:) = srcvals{ii}(4:6, :) * jac(1,1,jj) + srcvals{ii}(7:9, :)*jac(1,2,jj);
        tmp(7:9,:) = srcvals{ii}(4:6, :) * jac(2,1,jj) + srcvals{ii}(7:9, :)*jac(2,2,jj);

    end
    srcmat = [srcmat tmp];
end

srcmat(10:12,:) = cross(srcmat(4:6,:),srcmat(7:9,:));
srcmat(10:12,:) = srcmat(10:12,:) ./ vecnorm(srcmat(10:12,:));

nptchover = S.npatches + 3*length(patch_ids);

Sover = surfer(nptchover,norder,srcmat,ip0);

% build sparse interpolation matrix interp_mat

% assemble the global sparse interpolation matrix.
nblocks = numel(interp_mat_blocks);
block_i = 1:npatches_interp;
block_j = 1:npatches;
% number of block-rows
nrow_blocks  = max(block_i);
% number of block-cols
ncol_blocks  = max(block_j);
% allocate sparse matrix
interp_mat = sparse(nrow_blocks*npols, ncol_blocks*npols);

for imat = 1:nblocks
    % identify all block locations using the current block type.
    irows = rows_idx(interp_mat_list == imat);
    icols = cols_idx(interp_mat_list == imat);

    % skip if this block type does not occur.
    if isempty(irows), continue; end

    % sparse indicator matrix marking block placement locations.
    Ssp = sparse(irows, icols, 1, nrow_blocks, ncol_blocks);

    % insert blocks using a Kronecker tensor product.
    interp_mat = interp_mat + kron(Ssp, interp_mat_blocks{imat});
end

% get sparse triplet representation
[ivec,jvec,vvec] = find(interp_mat);

return

end


function [vout] = gettrichildren(v0)

%     output
vout = zeros(2,3,4);
% v1 = zeros(2,3);
% v2 = zeros(2,3);
% v3 = zeros(2,3);
% v4 = zeros(2,3);

%     midpoints

vm = zeros(2,3);
vm(:,1) = (v0(:,1)+v0(:,2))/2;
vm(:,2) = (v0(:,3)+v0(:,2))/2;
vm(:,3) = (v0(:,1)+v0(:,3))/2;

%     first triangle
vout(:,1,1) = v0(:,1);
vout(:,2,1) = vm(:,1);
vout(:,3,1) = vm(:,3);

%      second triangle
vout(:,1,2) = vm(:,2);
vout(:,2,2) = vm(:,3);
vout(:,3,2) = vm(:,1);

%      third triangle
vout(:,1,3) = vm(:,1);
vout(:,2,3) = v0(:,2);
vout(:,3,3) = vm(:,2);

%      fourth triangle
vout(:,1,4) = vm(:,3);
vout(:,2,4) = vm(:,2);
vout(:,3,4) = v0(:,3);

return

end

function [uvout] = mapuv_tri(verts,uvs)

dx = verts(1,2)-verts(1,1);
dy = verts(2,3)-verts(2,1) ;
uvout = cat(1,verts(1,1) + dx*uvs(1,:),verts(2,1) + dy*uvs(2,:));
return

end

function [vout] = getquadchildren(v0)

%     output
vout = zeros(2,3,4);
% v1 = zeros(2,3);
% v2 = zeros(2,3);
% v3 = zeros(2,3);
% v4 = zeros(2,3);

%     midpoints
       vm = zeros(2,4);

      vm(1,1) = (v0(1,1)+v0(1,2))/2;
      vm(2,1) = v0(2,1);

      vm(1,2) = (v0(1,3)+v0(1,2))/2;
      vm(2,2) = (v0(2,3)+v0(2,2))/2;

      vm(1,3) = v0(1,1);
      vm(2,3) = (v0(2,1)+v0(2,3))/2;

      vm(1,4) = v0(1,2);
      vm(2,4) = (v0(2,1) + v0(2,3))/2;

      vm(1,5) = (v0(1,1)+v0(1,2))/2;
      vm(2,5) = v0(2,3);

%     first quad

  vout(:,1,1) = v0(:,1);
      v1(1,1) = v0(1,1);
      v1(2,1) = v0(2,1);

  vout(:,2,1) = vm(:,1);
      v1(1,2) = vm(1,1);
      v1(2,2) = vm(2,1);

  vout(:,3,1) = vm(:,3);
      v1(1,3) = vm(1,3);
      v1(2,3) = vm(2,3);

%      second quad

  vout(:,1,2) = vm(:,1);
      v2(1,1) = vm(1,1);
      v2(2,1) = vm(2,1);

  vout(:,2,2) = v0(:,2);
      v2(1,2) = v0(1,2);
      v2(2,2) = v0(2,2);

  vout(:,3,2) = vm(:,2);
      v2(1,3) = vm(1,2);
      v2(2,3) = vm(2,2);

%      third quad
  vout(:,1,3) = vm(:,3);
      v3(1,1) = vm(1,3);
      v3(2,1) = vm(2,3);

  vout(:,2,3) = vm(:,2);
      v3(1,2) = vm(1,2);
      v3(2,2) = vm(2,2);

  vout(:,3,3) = v0(:,3);
      v3(1,3) = v0(1,3);
      v3(2,3) = v0(2,3);

%      fourth quad
  vout(:,1,4) = vm(:,2);
      v4(1,1) = vm(1,2);
      v4(2,1) = vm(2,2);

  vout(:,2,4) = vm(:,4);
      v4(1,2) = vm(1,4);
      v4(2,2) = vm(2,4);

  vout(:,3,4) = vm(:,5);
      v4(1,3) = vm(1,5);
      v4(2,3) = vm(2,5);


return

end

function [uvout] = mapuv_quad(verts,uvs)

dx = verts(1,2)-verts(1,1);
dy = verts(2,3)-verts(2,1) ;
uvout = cat(1,verts(1,1) + dx*(uvs(1,:)+1)/2,verts(2,1) ...
    + dy*(uvs(2,:)+1)/2);
return

end