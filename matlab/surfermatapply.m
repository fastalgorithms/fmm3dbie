function pot = surfermatapply(surferobj, kern, dens, eps, objover, cors, opts)
%SURFERMATAPPLY apply the integral operator for given kernel and surfer
% description of boundary to a density, using the FMM-accelerated smooth
% quadrature rule and adding a precomputed sparse correction.
%
% Syntax: pot = surfermatapply(surferobj, kern, dens, eps, objover, cors)
%         pot = surfermatapply(surferobj, kern, dens, eps, objover, cors, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern  - kernel3d object or matrix of kernel3d objects
%   dens  - density vector, length ncols (consistent with surferobj and kern)
%   eps   - quadrature tolerance
%   objover - oversampling specification, one of:
%     []       recompute from scratch
%     a vector or cell of oversampling orders, broadcast as needed
%     the precomputed oversampling info returned by surfermat
%
%   cors  - sparse matrix of near-field quadrature corrections,
%           as returned by surfermat(S,kern,opts) with nonsmoothonly=true,
%           or corrections=true. If empty, start by calling surfermat.
%
% Optional input:
%   opts  - options structure
%           opts.usefmm     = (1) use FMM for smooth quadrature
%           opts.unif_nover = (0) if nonzero, enforce uniform oversampling
%
% A single kern uses one FMM for all surfers when each surfer is
% oversampled the same way for all targets (e.g. opts.unif_nover = 1).
%
% Output:
%   pot  - result of applying the operator to dens, length nrows
%
% Author: Tristan Goodwill

if iscell(surferobj)
    surfers = surferobj;
else
    surfers = num2cell(surferobj(:).',1);
end

if ~isa(kern, 'kernel3d')
    error('SURFERMATAPPLY: second input kern not of supported type');
end

if nargin < 5, objover = []; end
if nargin < 6, cors = []; end
if nargin < 7, opts = []; end

usefmm = 1;
if isfield(opts,'usefmm'), usefmm = opts.usefmm; end

unif_nover = 0;
if isfield(opts,'unif_nover'), unif_nover = opts.unif_nover; end
if isfield(opts,'unif_novers'), unif_nover = opts.unif_novers; end

if isempty(cors)
    coropts = opts;
    coropts.corrections = 1;
    [cors,objover] = surfermat(surferobj,kern,eps,coropts);
end


nsurfers = length(surfers);

lsurfer    = zeros(nsurfers, 1);
for i=1:nsurfers
    lsurfer(i) = surfers{i}.npts;
end

if numel(kern) == 1
    opdims_mat = repmat(reshape(kern.opdims, 2, 1, 1), 1, nsurfers, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers, nsurfers);
end

% --- Parse objover into a canonical form ---
%
% objover_mode == false : objover is an nsurfers x nsurfers cell of order vectors
% objover_mode == true  : surfers_over and xinterps are precomputed;
%                         surfers_over{i,j} / xinterps{i,j} used directly.

objover_mode = false;
surfers_over = [];
xinterps     = [];

if iscell(objover) && numel(objover) == 2 && iscell(objover{1})
    % 2-element cell: {surfers_over, xinterps}
    objover_mode = true;
    surfers_over_in = objover{1};
    xinterps_in     = objover{2};

    % Normalise surfers_over to an nsurfers x nsurfers cell
    if ~iscell(surfers_over_in)
        % Array of surfer objects — one per source column
        surfers_over = cell(nsurfers, nsurfers);
        for j = 1:nsurfers
            for i = 1:nsurfers
                surfers_over{i,j} = surfers_over_in(j);
            end
        end
    elseif numel(surfers_over_in) == nsurfers
        % Cell vector of surfer objects — one per source column
        surfers_over = cell(nsurfers, nsurfers);
        for j = 1:nsurfers
            for i = 1:nsurfers
                surfers_over{i,j} = surfers_over_in{j};
            end
        end
    else
        surfers_over = surfers_over_in;
    end

    % Normalise xinterps to an nsurfers x nsurfers cell
    if numel(xinterps_in) == nsurfers && (isvector(xinterps_in) || iscell(xinterps_in))
        xinterps = cell(nsurfers, nsurfers);
        for j = 1:nsurfers
            for i = 1:nsurfers
                xinterps{i,j} = xinterps_in{j};
            end
        end
    else
        xinterps = xinterps_in;
    end

elseif ~isempty(objover) && ~iscell(objover)
    % Plain vector: broadcast to all (i,j) pairs
    novers_vec = objover(:);
    objover = cell(nsurfers, nsurfers);
    for i = 1:nsurfers
        for j = 1:nsurfers
            objover{i,j} = novers_vec;
        end
    end
end

% If still empty (or a standard nsurf x nsurf cell), compute orders as needed
if ~objover_mode
    novers_tmp = cell(nsurfers, nsurfers);
    for i=1:nsurfers
        for j=1:nsurfers
            if numel(kern) == 1
                ktmp = kern;
            else
                ktmp = kern(i,j);
            end
            if isempty(objover)
                if ismethod(surfers{j},'oversample')
                    novers_tmp{i,j} = ktmp.get_overs_orders(surfers{j},surfers{i},eps);
                else
                    novers_tmp{i,j} = NaN;
                end
            end
        end
    end
    if isempty(objover)
        if unif_nover
            for j=1:nsurfers
                noversj = cell2mat(novers_tmp(:,j).');
                noversj = repmat(max(noversj,[],2),1,nsurfers);
                novers_tmp(:,j) = mat2cell(noversj,size(noversj,1),2);
            end
        end
        objover = novers_tmp;
    end
end

for i = 1:nsurfers
    assert(all(opdims_mat(1,i,:) == opdims_mat(1,i,1)), ...
        'SURFERMATAPPLY: opdims(1) is not constant across block-row %d', i);
end
for j = 1:nsurfers
    assert(all(opdims_mat(2,:,j) == opdims_mat(2,1,j)), ...
        'SURFERMATAPPLY: opdims(2) is not constant across block-column %d', j);
end

irowlocs = zeros(nsurfers+1, 1);
icollocs = zeros(nsurfers+1, 1);

irowlocs(1) = 1;
icollocs(1) = 1;
for i = 1:nsurfers
    icollocs(i+1) = icollocs(i) + lsurfer(i)*opdims_mat(2,1,i);
    irowlocs(i+1) = irowlocs(i) + lsurfer(i)*opdims_mat(1,i,1);
end

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

if objover_mode
    srcover = @(i,j,nd,d) over_src(surfers{j}, surfers_over{i,j}, xinterps{i,j}, nd, d);
    ovkey = @(i,j) surfers_over{i,j}.norders;
else
    srcover = @(i,j,nd,d) over_src(surfers{j}, objover{i,j}, [], nd, d);
    ovkey = @(i,j) objover{i,j};
end

% One FMM for all pairs for a single kernel, if each source surfer is
% oversampled the same way for all targets
batch = numel(kern) == 1 && usefmm && ~isempty(kern.fmm) && all(lsurfer > 0);
for j = 1:nsurfers
    for i = 2:nsurfers
        batch = batch && isequaln(ovkey(i,j), ovkey(1,j));
    end
end

pot = zeros(nrows, 1);

if batch
    src = cell(1, nsurfers);
    sig = cell(nsurfers, 1);
    for j = 1:nsurfers
        [src{j}, sig{j}] = srcover(1, j, kern.opdims(2), dens(icollocs(j):(icollocs(j+1)-1)));
    end
    pot = kern.fmm(eps, merge([src{:}]), merge([surfers{:}]), cell2mat(sig));
    pot = pot(:) + cors * dens;
    return
end

for i = 1:nsurfers
    surferi = surfers{i};
    for j = 1:nsurfers
        surferj = surfers{j};
        if (surferi.npts < 1 || surferj.npts < 1)
            continue
        end

        if (numel(kern) == 1)
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end

        surferi_targ = []; surferi_targ.r = surferi.r(:,:);
        for field = ktmp.targ_fields(:).'
            surferi_targ.(field{1}) = slice_field(surferi.(field{1}), 1:surferi.npts);
        end

        if i==j
            surferi_targ.patch_id = surferi.patch_id;
            surferi_targ.uvs_targ = surferi.uvs_targ;
            surferi_targ.npts = surferi.npts;
            surferi_targ.ixyzs = surferi.ixyzs;
        end

        icolinds = icollocs(j):(icollocs(j+1)-1);
        irowinds = irowlocs(i):(irowlocs(i+1)-1);

        [surferjover, dens_j] = srcover(i, j, ktmp.opdims(2), dens(icolinds));
            if usefmm && ~isempty(ktmp.fmm)
                pot_ij = ktmp.fmm(eps,surferjover,surferi_targ,dens_j);
            elseif i == j
                % Loop over source patches: only the self patch needs
                % eval_mask; all other target points use eval.
                pot_ij = zeros(ktmp.opdims(1)*surferi_targ.npts, 1);
                for p = 1:surferjover.npatches
                    src_inds = surferjover.ixyzs(p):(surferjover.ixyzs(p+1)-1);
                    srcp = []; srcp.r = surferjover.r(:,src_inds);
                    for field = ktmp.src_fields(:).'
                        srcp.(field{1}) = slice_field(surferjover.(field{1}), src_inds);
                    end
                    col_inds = (1:ktmp.opdims(2)).' + ktmp.opdims(2)*(src_inds-1);
                    col_inds = col_inds(:);
                    dens_p = dens_j(col_inds);
    
                    targ_inds_self = surferi_targ.ixyzs(p):(surferi_targ.ixyzs(p+1)-1);
                    row_inds_self  = (1:ktmp.opdims(1)).' + ktmp.opdims(1)*(targ_inds_self-1);
                    row_inds_self  = row_inds_self(:);
                    targp_self = []; targp_self.r = surferi_targ.r(:,targ_inds_self);
                    for field = ktmp.targ_fields(:).'
                        targp_self.(field{1}) = slice_field(surferi_targ.(field{1}), targ_inds_self);
                    end
                    pot_ij(row_inds_self) = pot_ij(row_inds_self) + ktmp.eval_mask(srcp,targp_self)*dens_p;
    
                    off_targ_inds = [1:targ_inds_self(1)-1, targ_inds_self(end)+1:surferi_targ.npts];
                    if ~isempty(off_targ_inds)
                        targp_off = []; targp_off.r = surferi_targ.r(:,off_targ_inds);
                        for field = ktmp.targ_fields(:).'
                            targp_off.(field{1}) = slice_field(surferi_targ.(field{1}), off_targ_inds);
                        end
                        row_inds_off = (1:ktmp.opdims(1)).' + ktmp.opdims(1)*(off_targ_inds-1);
                        row_inds_off = row_inds_off(:);
                        pot_ij(row_inds_off) = pot_ij(row_inds_off) + ktmp.eval(srcp,targp_off)*dens_p;
                    end
                end
            else
                pot_ij = ktmp.eval(surferjover,surferi_targ)*dens_j;
            end
        pot(irowinds) = pot(irowinds) + pot_ij(:);
    end
end

%% Add near-field correction
pot = pot + cors * dens;

end

function [sover, sig] = over_src(surf, ov, xinterp, nd, dens)
% Oversampled source and its weighted density; ov is the oversampled
% surfer (with xinterp) or the oversampling orders (xinterp empty)
if isempty(xinterp)
    if ismethod(surf,'oversample') && ~any(isnan(ov))
        [sover, xinterp] = surf.oversample(ov);
        xinterp = kron(xinterp, eye(nd));
    else
        sover = surf;
        xinterp = speye(nd*surf.npts);
    end
else
    sover = ov;
end
wtsover = repmat(sover.wts(:).', nd, 1);
sig = wtsover(:).*(xinterp*dens);
end

function v = slice_field(A, idx)
% Slice a per-point field by point index, tolerating either orientation
if isvector(A) && iscolumn(A)
    v = A(idx).';
else
    v = A(:, idx);
end
end
