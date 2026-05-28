function pot = surfermatapply(surferobj, kern, dens, eps, novers, cors, opts)
%SURFERMATAPPLY apply the integral operator for given kernel and surfer
% description of boundary to a density, using the FMM-accelerated smooth
% quadrature rule and adding a precomputed sparse correction.
%
% Syntax: pot = surfermatapply(S, kern, dens, cors, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern  - kernel3d object or matrix of kernel3d objects
%   dens  - density vector, length ncols (consistent with surferobj and kern)
%   cors  - sparse matrix of near-field quadrature corrections,
%           as returned by surfermat(S,kern,opts) with nonsmoothonly=true,
%           or as assembled by the user.  cors*dens is added to the output.
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

if nargin < 5, novers = []; end
if nargin < 6, cors = []; end
if nargin < 7, opts = []; end

usematlab = 1;
if isfield(opts,'usematlab'), usematlab = opts.usematlab; end
usefmm = 1;
if isfield(opts,'usefmm'), usefmm = opts.usefmm; end

unif_nover = 0;
if isfield(opts,'unif_nover'), unif_nover = opts.unif_nover; end
if isfield(opts,'unif_novers'), unif_nover = opts.unif_novers; end


if isempty(cors)
    % warning('surfermatapply is inefficient with empty corrections')
    coropts = opts;
    if usematlab
        coropts.corrections = 1;
    else
        coropts.nonsmoothonly = 1;
    end
    [cors,novers] = surfermat(surferobj,kern,eps,coropts);
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

novers_tmp = cell(nsurfers, nsurfers);
for i=1:nsurfers
    for j=1:nsurfers
        if numel(kern) == 1
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end
        if isempty(novers)
            if ismethod(surfers{j},'oversample')
                novers_tmp{i,j} = ktmp.get_overs_orders(surfers{j},surfers{i},eps);
            else
                novers_tmp{i,j} = NaN;
            end
        end
    end
end
if isempty(novers)
    if unif_nover
        for j=1:nsurfers
            noversj = cell2mat(novers_tmp(:,j).');
            noversj = repmat(max(noversj,[],2),1,nsurfers);
            novers_tmp(:,j) = mat2cell(noversj,size(noversj,1),2);
        end
    end
    novers = novers_tmp;
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

pot = zeros(nrows, 1);

for i = 1:nsurfers
    surferi = surfers{i};
    for j = 1:nsurfers
        surferj = surfers{j};
        if (surferi.npts < 1 || surferj.npts < 1)
            continue
        end

        surferi_targ = []; surferi_targ.r = surferi.r(:,:);
        for field = ktmp.targ_fields(:).'
            surferi_targ.(field{1}) = surferi.(field{1})(:,:);
        end

        if (numel(kern) == 1)
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end

        icolinds = icollocs(j):(icollocs(j+1)-1);
        irowinds = irowlocs(i):(irowlocs(i+1)-1);

        dens_j = dens(icolinds);

        % call layer eval or handle things in matlab
        if usematlab
            noversj = novers{i,j};

            if ismethod(surferj,'oversample') && ~any(isnan(noversj))
                [surferjover, xinterp] = surferj.oversample(noversj);
                xinterp = kron(xinterp, eye(ktmp.opdims(2)));
            else
                surferjover = surferj;
                xinterp = eye(numel(wts));
            end
            wtsover = repmat(surferjover.wts(:).', ktmp.opdims(2), 1);
            dens_j = wtsover(:).*(xinterp*dens_j);
            if usefmm && ~isempty(ktmp.fmm)
                pot_ij = ktmp.fmm(eps,surferjover,surferi_targ,dens_j);
            elseif i == j
                % Loop over source patches: only the co-located patch needs
                % eval_mask; all other target points use plain eval.
                pot_ij = zeros(opdims(1)*surferi_targ.npts, 1);
                for p = 1:surferjover.npatches
                    src_inds = surferjover.ixyzs(p):(surferjover.ixyzs(p+1)-1);
                    srcp = []; srcp.r = surferjover.r(:,src_inds);
                    for field = ktmp.src_fields(:).'
                        srcp.(field{1}) = surferjover.(field{1})(:,src_inds);
                    end
                    col_inds = (1:opdims(2)).' + opdims(2)*(src_inds-1);
                    col_inds = col_inds(:);
                    dens_p = dens_j(col_inds);

                    targ_inds_self = surferi_targ.ixyzs(p):(surferi_targ.ixyzs(p+1)-1);
                    row_inds_self  = (1:opdims(1)).' + opdims(1)*(targ_inds_self-1);
                    row_inds_self  = row_inds_self(:);
                    targp_self = []; targp_self.r = surferi_targ.r(:,targ_inds_self);
                    for field = ktmp.targ_fields(:).'
                        targp_self.(field{1}) = surferi_targ.(field{1})(:,targ_inds_self);
                    end
                    pot_ij(row_inds_self) = pot_ij(row_inds_self) + ktmp.eval_mask(srcp,targp_self)*dens_p;

                    off_targ_inds = [1:targ_inds_self(1)-1, targ_inds_self(end)+1:surferi_targ.npts];
                    if ~isempty(off_targ_inds)
                        targp_off = []; targp_off.r = surferi_targ.r(:,off_targ_inds);
                        for field = ktmp.targ_fields(:).'
                            targp_off.(field{1}) = surferi_targ.(field{1})(:,off_targ_inds);
                        end
                        row_inds_off = (1:opdims(1)).' + opdims(1)*(off_targ_inds-1);
                        row_inds_off = row_inds_off(:);
                        pot_ij(row_inds_off) = pot_ij(row_inds_off) + ktmp.eval(srcp,targp_off)*dens_p;
                    end
                end
            else
                pot_ij = ktmp.eval(surferjover,surferi_targ)*dens_j;
            end
        else
            % Build a proper Q struct from the (i,j) block of cors
            cors_ij = cors(irowinds, icolinds);
            rsc_ij  = conv_spmat_to_rsc(surferj, cors_ij, ktmp.rsc_to_interleave);

            % TODO we could get rfac faster than this
            rsc_near = getnear(surferj, surferi_targ);

            Q = [];
            Q.targinfo     = surferi_targ;
            Q.format       = 'rsc';  Q.wavenumber = ktmp.zk;
            Q.kernel_order = ktmp.kernel_order; 
            Q.nquad   = rsc_ij.nquad;   Q.row_ptr = rsc_ij.row_ptr;
            Q.col_ind = rsc_ij.col_ind; Q.iquad   = rsc_ij.iquad;
            Q.wnear   = rsc_ij.wnear;   Q.rfac    = rsc_near.rfac;

            opts_eval = [];
            opts_eval.precomp_quadrature = Q;
            pot_ij = ktmp.layer_eval(surferj, dens_j, surferi_targ, eps, opts_eval);
        end
        pot(irowinds) = pot(irowinds) + pot_ij(:);
    end
end

%% Add near-field correction (usematlab path only; layer_eval path folds
%  cors into precomp_quadrature above)

if usematlab
    pot = pot + cors * dens;
end

end
