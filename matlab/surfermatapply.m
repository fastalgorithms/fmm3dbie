function pot = surfermatapply(surferobj, kern, dens, cors, opts)
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
%   opts  - options structure (optional)
%           opts.eps = (1e-10) tolerance for FMM / oversampling
%
% Output:
%   pot  - result of applying the operator to dens, length nrows
%
% Author: Tristan Goodwill

surfers = surferobj;

if ~isa(kern, 'kernel3d')
    error('SURFERMATAPPLY: second input kern not of supported type');
end

if nargin < 5
    opts = [];
end

eps = 1e-10;
if isfield(opts, 'eps')
    eps = opts.eps;
end

nsurfers = length(surfers);

opdims_mat = zeros(2, nsurfers, nsurfers);
lsurfer    = zeros(nsurfers, 1);

for i = 1:nsurfers
    lsurfer(i) = surfers(i).npts;

    for j = 1:nsurfers
        if (numel(kern) == 1)
            opdims_mat(:,i,j) = kern.opdims;
        else
            opdims_mat(:,i,j) = kern(i,j).opdims;
        end
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

pot = zeros(nrows, 1);

%% Apply smooth part via layer_eval with nonsmoothonly = true

opts_eval = [];
opts_eval.nonsmoothonly = true;

for i = 1:nsurfers
    surferi = surfers(i);
    for j = 1:nsurfers
        surferj = surfers(j);
        if (surferi.npts < 1 || surferj.npts < 1)
            continue
        end

        % When i~=j, strip patch_id/uvs_targ so evaluators treat targets
        % as off-surface points.
        if i == j
            surferi_targ = surferi;
        else
            targs = extract_targ_array(surferi);
            surferi_targ = [];
            surferi_targ.r = targs(1:3,:);
            if size(targs,1) == 12
                surferi_targ.du = targs(4:6,:);
                surferi_targ.dv = targs(7:9,:);
                surferi_targ.n  = targs(10:12,:);
            end
        end

        if (numel(kern) == 1)
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end

        icolinds = icollocs(j):(icollocs(j+1)-1);
        irowinds = irowlocs(i):(irowlocs(i+1)-1);

        dens_j = dens(icolinds);

        pot_ij = ktmp.layer_eval(surferj, dens_j, surferi_targ, eps, opts_eval);

        pot(irowinds) = pot(irowinds) + pot_ij(:);
    end
end

%% Add near-field correction

pot = pot + cors * dens;

end
