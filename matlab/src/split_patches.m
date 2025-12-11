function [Sover,ipatchid,uvstot,interp_mat,ivec,jvec,vvec] = split_patches(S,patch_ids)
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
        % Vioreanu-Rokhlin nodes on the triangle simplex
        [uvs] = koorn.rv_nodes(norder);
        % number of nodes
        npols = size(uvs,2);
        % split into 4 subpatches
        uvssub = uvs / 2;
        uvssub1 = uvssub;
        uvssub2 = -uvssub + [1/2;1/2];
        uvssub3 = uvssub + [1/2;0];
        uvssub4 = uvssub + [0;1/2];
        % map: coefficients in expansion on patch to values on subpatches
        c2vover1 = koorn.coefs2vals(norder,uvssub1);
        c2vover2 = koorn.coefs2vals(norder,uvssub2);
        c2vover3 = koorn.coefs2vals(norder,uvssub3);
        c2vover4 = koorn.coefs2vals(norder,uvssub4);
        % map: values on patch to coefficients in expansion on patch
        v2c = koorn.vals2coefs(norder,uvs);
        % map: coefficients in expansion on patch to values on patch
        c2v = koorn.coefs2vals(norder,uvs);
        % map: values of patch basis at uvs, to values of derivative of 
        %      patch basis, wrt u and v, at uvs        
        [~,dersu,dersv] = koorn.ders(norder,uvs);
        val2du = dersu.'*v2c;
        val2dv = dersv.'*v2c;
    elseif(ip0 == 11)
        % tensor grid of Gauss-Legendre nodes on [-1,1]^2
        [uvs] = polytens.lege.nodes(norder);
        % number of nodes
        npols = size(uvs,2);
        % split into 4 subpatches
        uvssub = uvs / 2;
        uvssub1 = uvssub + [-1/2;-1/2];
        uvssub2 = uvssub + [1/2;-1/2];
        uvssub3 = uvssub + [-1/2;1/2];
        uvssub4 = uvssub + [1/2;1/2];
        % map: coefficients in expansion on patch to values on subpatches
        c2vover1 = polytens.lege.coefs2vals(norder,uvssub1);
        c2vover2 = polytens.lege.coefs2vals(norder,uvssub2);
        c2vover3 = polytens.lege.coefs2vals(norder,uvssub3);
        c2vover4 = polytens.lege.coefs2vals(norder,uvssub4);
        % map: values on patch to coefficients in expansion on patch
        v2c = polytens.lege.vals2coefs(norder,uvs);
        % map: coefficients in expansion on patch to values on patch
        c2v = polytens.lege.coefs2vals(norder,uvs);
        % map: values of patch basis at uvs, to values of derivative of 
        %      patch basis, wrt u and v, at uvs        
        [~,dersu,dersv] = polytens.lege.ders(norder,uvs);
        val2du = dersu.'*v2c;
        val2dv = dersv.'*v2c;
    elseif(ip0 == 12)
        % tensor grid of Chebyshev nodes on [-1,1]^2
        [uvs] = polytens.cheb.nodes(norder);
        % number of nodes
        npols = size(uvs,2);
        % split into 4 subpatches
        uvssub = uvs / 2;
        uvssub1 = uvssub + [-1/2;-1/2];
        uvssub2 = uvssub + [1/2;-1/2];
        uvssub3 = uvssub + [-1/2;1/2];
        uvssub4 = uvssub + [1/2;1/2];
        % map: coefficients in expansion on patch to values on subpatches
        c2vover1 = polytens.cheb.coefs2vals(norder,uvssub1);
        c2vover2 = polytens.cheb.coefs2vals(norder,uvssub2);
        c2vover3 = polytens.cheb.coefs2vals(norder,uvssub3);
        c2vover4 = polytens.cheb.coefs2vals(norder,uvssub4);
        % map: values on patch to coefficients in expansion on patch
        v2c = polytens.cheb.vals2coefs(norder,uvs);
        % map: coefficients in expansion on patch to values on patch
        c2v = polytens.cheb.coefs2vals(norder,uvs);
        % map: values of patch basis at uvs, to values of derivative of 
        %      patch basis, wrt u and v, at uvs        
        [~,dersu,dersv] = polytens.cheb.ders(norder,uvs);
        val2du = dersu.'*v2c;
        val2dv = dersv.'*v2c;
    else
        fprintf('Invalid type of patch, exiting\n');
        return
    end
    % store nodes for subpatches
    uvssplit = [uvssub1 uvssub2 uvssub3 uvssub4];

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
            panval1 = S.srccoefs{ii}(1:3,:)*c2vover1.';
            panval2 = S.srccoefs{ii}(1:3,:)*c2vover2.';
            panval3 = S.srccoefs{ii}(1:3,:)*c2vover3.';
            panval4 = S.srccoefs{ii}(1:3,:)*c2vover4.';
            srcvals = [srcvals {panval1,panval2,panval3,panval4}];
            % store local uv coordinates
            uvstot = [uvstot uvssplit];
            % store id of patch ii new points on subpatches belong to
            ipatchid = [ipatchid; ii*ones(length(uvssplit),1)];
            % refined patch: contributes four subpatch block-rows.
            interp_mat_list = [interp_mat_list [2 3 4 5]];
            cols_idx = [cols_idx ii*ones(1,4)];
            rows_idx = [rows_idx row_ctr+(1:4)];
            row_ctr = row_ctr + 4;                    
        else
            srcvals = [srcvals {S.srccoefs{ii}(1:3,:)*c2v.'}];
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

    for ii = 1:length(srcvals)
        tmp = zeros(12,size(uvs,2));
        tmp(1:3,:) = srcvals{ii};
        tmp(4:6,:) = srcvals{ii}*val2du.';
        tmp(7:9,:) = srcvals{ii}*val2dv.';
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


