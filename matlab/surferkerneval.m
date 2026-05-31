function pot = surferkerneval(surferobj, kern, dens, targobj, eps, opts)
%SURFERKERNEVAL evaluate layer potential for given kernel, surfer description
% of boundary, and off-surface target points.
%
% Syntax: pot = surferkerneval(S, kern, dens, targobj, eps, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern      - kernel3d object or vector of kernel3d objects (one per surfer)
%   dens      - density vector, length ncols (consistent with surferobj and kern)
%   targobj   - target point description: surfer object, struct with field
%               .r (and optionally .du, .dv, .n), or numeric (3,:) array
%   eps       - quadrature tolerance
%
% Optional input:
%   opts  - options structure
%           opts.corrections = sparse correction matrix as returned by
%                              surferkernevalmat with nonsmoothonly/corrections;
%                              recomputed if empty
%           opts.novers      = cell(nsurfers,1) of oversampling orders;
%                              recomputed if empty
%           opts.usematlab   = (1) use matlab smooth quadrature + corrections;
%                              0 uses layer_eval (Fortran) path instead
%           opts.usefmm      = (1) use FMM for smooth quadrature (usematlab=1 only)
%           opts.rfac        - rfac scalar or cell(nsurfers,1) as returned by
%                              surferkernevalmat; used when usematlab=0.
%                              Ignored (or treated as NaN) if not provided.
%
% Output:
%   pot  - result of applying the operator to dens, length ntarg*opdims(1)
%
% Author: Tristan Goodwill

if iscell(surferobj)
    surfers = surferobj;
else
    surfers = num2cell(surferobj(:).',1);
end

if ~isa(kern, 'kernel3d')
    error('SURFERKERNEVAL: second input kern not of supported type');
end

if nargin < 6, opts = []; end

cors = [];
if isfield(opts,'corrections'), cors = opts.corrections; end
novers = [];
if isfield(opts,'novers'), novers = opts.novers; end
usematlab = 1;
if isfield(opts,'usematlab'), usematlab = opts.usematlab; end
usefmm = 1;
if isfield(opts,'usefmm'), usefmm = opts.usefmm; end
rfac_in = [];
if isfield(opts,'rfac'), rfac_in = opts.rfac; end


if usematlab && isempty(cors)
    coropts = opts;
    coropts.corrections = 1;
    [cors,novers] = surferkernevalmat(surferobj,kern, targobj, eps,coropts);
end

if isnumeric(targobj)
    targinfo = targobj.r;
else
    targinfo = targobj;
end



nsurfers = length(surfers);

opdims_mat = zeros(2, nsurfers);
lsurfer    = zeros(nsurfers, 1);

for j = 1:nsurfers
    lsurfer(j) = surfers{j}.npts;

    if (numel(kern) == 1)
        opdims_mat(:,j) = kern.opdims;
    else
        opdims_mat(:,j) = kern(j).opdims;
    end
end


if size(dens,2) > 1
    if size(dens,2) > 10
        warning('Use surferkernevalmat for layer potentials with many densities')
    end
    pot = zeros(opdims_mat(1,1)*size(targinfo.r(:,:),2),size(dens,2));
    opts.corrections = cors;
    for i = 1:size(dens,2)
        pot(:,i) = surferkerneval(surferobj, kern, dens(:,i), targinfo, eps, opts);
    end
    return
end
novers_tmp = cell(nsurfers, 1);
for j=1:nsurfers
    lsurfer(j) = surfers{j}.npts;

    if (size(kern) == 1)
        ktmp = kern;
    else
        ktmp = kern(j);
    end
    opdims_mat(:,j) = ktmp.opdims;
    if isempty(novers)
        if ismethod(surfers{j},'oversample') 
            novers_tmp{j} = ktmp.get_overs_orders(surfers{j},targobj,eps);
        else
            novers_tmp{j} = NaN;
        end
    end
end
if isempty(novers)
    novers = novers_tmp;
end

assert(all(opdims_mat(1,:) == opdims_mat(1,1)), ...
    'SURFERKERNEVAL: opdims(1) is not constant across all source blocks');

ntarg  = size(targinfo.r, 2);
opdim1 = opdims_mat(1,1);

icollocs = zeros(nsurfers+1, 1);
icollocs(1) = 1;
for j = 1:nsurfers
    icollocs(j+1) = icollocs(j) + lsurfer(j)*opdims_mat(2,j);
end

nrows = ntarg * opdim1;
ncols = icollocs(end)-1;

pot = zeros(nrows, 1);

%% Apply smooth part via layer_eval with nonsmoothonly = true

for j = 1:nsurfers
    surferj = surfers{j};
    if (surferj.npts < 1 || ntarg < 1)
        continue
    end

    if (numel(kern) == 1)
        ktmp = kern;
    else
        ktmp = kern(j);
    end

    icolinds = icollocs(j):(icollocs(j+1)-1);
    dens_j = dens(icolinds);

    if usematlab
        noversj = novers{j};

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
            pot_j = ktmp.fmm(eps,surferjover,targinfo,dens_j);
        else
            pot_j = ktmp.eval(surferjover,targinfo)*dens_j;
        end
    else
        opts_eval = [];
        if ~isempty(cors)
            targinfosuse = []; targinfosuse.r = targinfo.r(:,:);
            for field = ktmp.targ_fields(:).'
                targinfosuse.(field{1}) = targinfo.(field{1})(:,:);
            end
            
            % Build a proper Q struct from the (i,j) block of cors
            cors_ij = cors(:, icolinds);
            rsc_ij  = conv_spmat_to_rsc(surferj, cors_ij, ktmp.rsc_to_interleave);

            % Use provided rfac if available, otherwise call getnear
            rfac_j = NaN;
            if ~isempty(rfac_in)
                if iscell(rfac_in)
                    rfac_j = rfac_in{j};
                else
                    rfac_j = rfac_in;
                end
            end
            if isscalar(rfac_j) && ~isnan(rfac_j)
                rfac_use = rfac_j;
            else
                rsc_near = getnear(surferj, targinfosuse);
                rfac_use = rsc_near.rfac;
            end

            Q = [];
            Q.targinfo     = targinfosuse;
            Q.format       = 'rsc';  Q.wavenumber = ktmp.zk;
            Q.kernel_order = ktmp.kernel_order;
            Q.nquad   = rsc_ij.nquad;   Q.row_ptr = rsc_ij.row_ptr;
            Q.col_ind = rsc_ij.col_ind; Q.iquad   = rsc_ij.iquad;
            Q.wnear   = rsc_ij.wnear;   Q.rfac    = rfac_use;

            opts_eval.precomp_quadrature = Q;
        end
        pot_j = ktmp.layer_eval(surferj, dens_j, targinfo, eps, opts_eval);
    end
    pot = pot + pot_j(:);
end

%% Add near-field correction
if ~isempty(cors) && usematlab
pot = pot + cors * dens;
end
end
