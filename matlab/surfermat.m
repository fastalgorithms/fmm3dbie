function [sysmat,objover,rfac] = surfermat(surferobj,kern,eps,opts)
%SURFERMAT build matrix for given kernel and surfer description of boundary.
%
% Syntax: [sysmat, objover, rfac] = surfermat(S, kern, eps, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern      - kernel3d object or matrix of kernel3d objects
%   eps       - quadrature tolerance
%
% Optional input:
%   opts  - options structure
%           opts.nonsmoothonly  = (false) return only near-field entries as sparse
%           opts.corrections    = (false) return corrections to smooth rule as sparse
%           opts.l2scale        = (false) scale rows by sqrt(wts), cols by 1/sqrt(wts)
%           opts.unif_nover     = (0) if nonzero, enforce uniform oversampling order
%           opts.ifoversamp     = (1) if 0, skip oversampling (set novers=NaN)
%           opts.ifreturnovers  = (1) if 1, second output is {surfers_over, xinterps}
%                                 (a 2-element cell of nsurfers x nsurfers cells)
%                                 instead of the cell array of oversampling orders
%
% Output:
%   sysmat - system matrix (dense, or sparse if nonsmoothonly/corrections)
%   objover - by default, cell(nsurfers,nsurfers) of per-patch oversampling
%             order vectors; if opts.ifreturnovers=1, a 2-element cell
%             {surfers_over, xinterps} of precomputed oversampled objects.
%   rfac   - cell(nsurfers,nsurfers) of rfac scalars from getquad;
%            NaN for blocks where getquad is not called or does not store rfac
%
% Author: Tristan Goodwill

% Accept either an array of surfers or a cell array of surfers.
if iscell(surferobj)
    surfers = surferobj;
else
    surfers = num2cell(surferobj(:).',1);
end

if ~isa(kern,'kernel3d'), error('SURFERMAT: second input kern not of supported type'); end

if nargin < 4
    opts = [];
end

nonsmoothonly = false;
corrections = false;
l2scale = false;
selfquad = true;

% get opts from struct if available
if isfield(opts,'selfquad'), selfquad = opts.selfquad; end
if isfield(opts,'l2scale'), l2scale = opts.l2scale; end
if isfield(opts,'nonsmoothonly'), nonsmoothonly = opts.nonsmoothonly; end
if isfield(opts,'corrections'), corrections = opts.corrections; end
if corrections, nonsmoothonly = true; end

adaptive_correction = true;
if isfield(opts,'adaptive_correction'), adaptive_correction = opts.adaptive_correction; end

forcesmooth = false;
if isfield(opts,'forcesmooth'), forcesmooth = opts.forcesmooth; end
if forcesmooth, adaptive_correction = false; selfquad = false; end

unif_nover = 0;
if isfield(opts,'unif_nover'), unif_nover = opts.unif_nover; end
if isfield(opts,'unif_novers'), unif_nover = opts.unif_novers; end

ifoversamp = 1;
if isfield(opts,'ifoversamp'), ifoversamp = opts.ifoversamp; end
ifreturnovers = 1;
if isfield(opts,'ifreturnovers'), ifreturnovers = opts.ifreturnovers; end
nsurfers = length(surfers);

lsurfer    = zeros(nsurfers,1);
for i=1:nsurfers
    lsurfer(i) = surfers{i}.npts;
end

if numel(kern) == 1
    opdims_mat = repmat(reshape(kern.opdims, 2, 1, 1), 1, nsurfers, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers, nsurfers);
end

novers       = cell(nsurfers, nsurfers);
surfers_over = cell(nsurfers, nsurfers);
xinterps     = cell(nsurfers, nsurfers);
rfac         = cell(nsurfers, nsurfers);
for i=1:nsurfers
    for j=1:nsurfers
        if numel(kern) == 1
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end
        if ismethod(surfers{j},'oversample')
            novers{i,j} = ktmp.get_overs_orders(surfers{j},surfers{i},eps);
        else
            novers{i,j} = NaN;
        end
        rfac{i,j} = NaN;
    end
end
if unif_nover
    for j=1:nsurfers
        noversj = cell2mat(novers(:,j).');
        noversj = repmat(max(noversj,[],2),1,nsurfers);
        novers(:,j) = mat2cell(noversj,size(noversj,1),ones(nsurfers,1));
    end
end
if ~ifoversamp
    for i=1:nsurfers
        for j=1:nsurfers
            novers{i,j} = NaN*novers{i,j};
        end
    end
end

 for i = 1:nsurfers
    assert(all(opdims_mat(1,i,:) == opdims_mat(1,i,1)), ...
        'SURFERMAT: opdims(1) is not constant across block-row %d', i);
end
for j = 1:nsurfers
    assert(all(opdims_mat(2,:,j) == opdims_mat(2,1,j)), ...
        'SURFERMAT: opdims(2) is not constant across block-column %d', j);
end

irowlocs = zeros(nsurfers+1,1);
icollocs = zeros(nsurfers+1,1);

irowlocs(1) = 1;
icollocs(1) = 1;
for i=1:nsurfers
   icollocs(i+1) = icollocs(i) + lsurfer(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + lsurfer(i)*opdims_mat(1,i,1);
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

if (~nonsmoothonly)
    sysmat = zeros(nrows,ncols);
else 
    sysmat = sparse(nrows,ncols);
    isysmat = [];
    jsysmat = [];
    vsysmat = [];
end

%% Now build quadratures

for i = 1:nsurfers
    surferi = surfers{i};
    for j = 1:nsurfers
        surferj = surfers{j};
        if (surferi.npts < 1 || surferj.npts<1)
            sysmat_tmp = [];
            break
        end

        if (size(kern) == 1)
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end

        % When i~=j, strip patch_id/uvs_targ so evaluators treat targets
        % as off-surface points.
        if i == j
            surferi_targ = surferi;
        else
            surferi_targ = []; surferi_targ.r = surferi.r(:,:);
            for field = ktmp.targ_fields(:).'
                surferi_targ.(field{1}) = surferi.(field{1})(:,:);
            end
        end

        opdims = reshape(opdims_mat(:,i,j),[2,1]);
        wts = repmat(surferj.wts(:).', opdims(2), 1);
        wts = wts(:).';

        noversj = novers{i,j};
        if ismethod(surferj,'oversample') && ~any(isnan(noversj))
            [surferjover, xinterp] = surferj.oversample(noversj);
            xinterp = kron(xinterp, eye(opdims(2)));
        else
            novers{i,j} = NaN;
            surferjover = surferj;
            xinterp = eye(numel(wts));
        end
        surfers_over{i,j} = surferjover;
        xinterps{i,j}     = xinterp;

        wtsover = repmat(surferjover.wts(:).', opdims(2), 1);
        wtsover = wtsover(:).';

        sysmat_tmp = [];
        if (~nonsmoothonly)
            nrows = opdims(1)*size(surferi_targ.r(:,:),2);
            ncols = opdims(2)*size(surferjover.r(:,:),2);
            sysmat_tmp = zeros(nrows, ncols);
            for p = 1:surferjover.npatches
                src_inds = surferjover.ixyzs(p):(surferjover.ixyzs(p+1)-1);
                srcp = [];
                srcp.r = surferjover.r(:,src_inds);
                for field = ktmp.src_fields(:).'
                    srcp.(field{1}) = surferjover.(field{1})(:,src_inds);
                end
                col_inds = (1:opdims(2)).' + opdims(2)*(src_inds-1);
                col_inds = col_inds(:);
                if i == j
                    % Only this patch can have coincident source/target points.
                    % Use eval_mask to zero them out.
                    targ_inds_self = surferi_targ.ixyzs(p):(surferi_targ.ixyzs(p+1)-1);
                    row_inds_self  = (1:opdims(1)).' + opdims(1)*(targ_inds_self-1);
                    row_inds_self  = row_inds_self(:);

                    targp_self = []; targp_self.r = surferi_targ.r(:,targ_inds_self);
                    for field = ktmp.targ_fields(:).'
                        targp_self.(field{1}) = surferi_targ.(field{1})(:,targ_inds_self);
                    end
                    targp_self.uvs_targ = surferi_targ.uvs_targ(:,targ_inds_self);
                    targp_self.patch_id = surferi_targ.patch_id(targ_inds_self);
                    sysmat_tmp(row_inds_self, col_inds) = ktmp.eval_mask(srcp, targp_self);

                    % Off-diagonal target blocks: plain eval
                    off_targ_inds = [1:targ_inds_self(1)-1, targ_inds_self(end)+1:surferi_targ.npts];
                    if ~isempty(off_targ_inds)
                        targp_off = []; targp_off.r = surferi_targ.r(:,off_targ_inds);
                        for field = ktmp.targ_fields(:).'
                            targp_off.(field{1}) = surferi_targ.(field{1})(:,off_targ_inds);
                        end
                        row_inds_off = (1:opdims(1)).' + opdims(1)*(off_targ_inds-1);
                        row_inds_off = row_inds_off(:);
                        sysmat_tmp(row_inds_off, col_inds) = ktmp.eval(srcp, targp_off);
                    end
                else
                    sysmat_tmp(:, col_inds) = ktmp.eval(srcp, surferi_targ);
                end
            end

            sysmat_tmp = (sysmat_tmp.*wtsover)*xinterp;
        end

        if (adaptive_correction && (i~=j)) || ((i==j) && selfquad) && ~isempty(ktmp.getquad)
            sysmat_quad = ktmp.getquad(surferj,eps,surferi_targ);
            if corrections
                [rp, ci] = get_rsc_pattern(surferj, sysmat_quad, opdims);
                sysmat_smth = smooth_sparse_quad(ktmp, surferi_targ, surferj, rp, ci, noversj);
                sysmat_quad = sysmat_quad-sysmat_smth;
            end

            if (~nonsmoothonly)
                inds = abs(sysmat_quad(:))>0;
                sysmat_tmp(inds) = sysmat_quad(inds);
            else
                sysmat_tmp = sysmat_quad;
            end
        end

        if isempty(sysmat_tmp), continue, end

        if l2scale
            wts = sqrt(wts); 
            wtsrow = surferi.wts; wtsrow = sqrt(wtsrow(:))';
            wtsrow = repmat(wtsrow,opdims(1),1); wtsrow = wtsrow(:);
            sysmat_tmp = wtsrow.*sysmat_tmp./wts;
        end
        
        if (~nonsmoothonly)
            irowinds = irowlocs(i):(irowlocs(i+1)-1);
            icolinds = icollocs(j):(icollocs(j+1)-1);
            sysmat(irowinds,icolinds) = sysmat_tmp;
        else
            if adaptive_correction
                [isys,jsys,vsys] = find(sysmat_tmp);
                isysmat = [isysmat;isys+irowlocs(i)-1];
                jsysmat = [jsysmat;jsys+icollocs(j)-1];
                vsysmat = [vsysmat;vsys];
            end
        end
    
    end
end

if (nonsmoothonly)
    sysmat = sparse(isysmat,jsysmat,vsysmat,nrows,ncols);
end

if ifreturnovers
    objover = {surfers_over, xinterps};
else
    objover = novers;
end
end
