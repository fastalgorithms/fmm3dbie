function pot = surferkerneval(surferobj, kern, dens, targobj, opts)
%SURFERKERNEVAL apply the evaluation operator for given kernel, surfer
% description of boundary, and off-surface target points to a density,
% using the FMM-accelerated smooth quadrature rule and adding a precomputed
% sparse correction.
%
% Syntax: pot = surferkerneval(S, kern, dens, targobj, cors, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern  - kernel3d object or vector of kernel3d objects (one per surfer)
%   targobj   - target point description: surfer object, struct with field
%               .r (and optionally .du, .dv, .n), or numeric (3,:) array
%   dens  - density vector, length ncols (consistent with surferobj and kern)
%   cors  - sparse matrix of near-field quadrature corrections,
%           as returned by surferkernevalmat(S,kern,targobj,opts) with
%           nonsmoothonly=true, or as assembled by the user.
%           cors*dens is added to the output.
%   opts  - options structure (optional)
%           opts.eps = (1e-10) tolerance for FMM / oversampling
%           opts.corrections = [];
%
% Output:
%   pot  - result of applying the operator to dens, length nrows
%          (= ntarg * opdims(1))
%
% Author: Tristan Goodwill

surfers = surferobj;

if ~isa(kern, 'kernel3d')
    error('SURFERKERNEVAL: second input kern not of supported type');
end

if nargin < 6
    opts = [];
end

eps = 1e-10;
if isfield(opts, 'eps')
    eps = opts.eps;
end
cors = [];
if isfield(opts,'corrections')
    cors = opts.corrections;
end


% Assign appropriate object to targinfo
if isa(targobj, 'surfer') || isstruct(targobj)
    targs_arr = extract_targ_array(targobj);
elseif isnumeric(targobj)
    targs_arr = targobj;
else
    error('SURFERKERNEVAL: input 3 is not a supported type');
end
targinfo = [];
targinfo.r = targs_arr(1:3,:);
if size(targs_arr, 1) == 12
    targinfo.du = targs_arr(4:6,:);
    targinfo.dv = targs_arr(7:9,:);
    targinfo.n  = targs_arr(10:12,:);
end

nsurfers = length(surfers);

opdims_mat = zeros(2, nsurfers);
lsurfer    = zeros(nsurfers, 1);

for j = 1:nsurfers
    lsurfer(j) = surfers(j).npts;

    if (numel(kern) == 1)
        opdims_mat(:,j) = kern.opdims;
    else
        opdims_mat(:,j) = kern(j).opdims;
    end
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

opts_eval = [];
if ~isempty(cors)
    opts_eval.nonsmoothonly = true;
end
for j = 1:nsurfers
    surferj = surfers(j);
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

    pot_j = ktmp.layer_eval(surferj, dens_j, targinfo, eps, opts_eval);

    pot = pot + pot_j(:);
end

%% Add near-field correction
if ~isempty(cors)
pot = pot + cors * dens;
end
end
