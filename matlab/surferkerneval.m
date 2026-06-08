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
%               .r or numeric (3,:) array
%   eps       - quadrature tolerance
%
% Optional input:
%   opts  - options structure
%           opts.corrections = sparse correction matrix as returned by
%                              surferkernevalmat with nonsmoothonly/corrections;
%                              recomputed if empty
%           opts.objover     = oversampling specification; one of:
%                              []               recompute from scratch
%                              vector           oversampling orders, broadcast to all surfers
%                              cell(nsurfers,1) per-surfer order vectors
%                              {surfers_over, xinterps}  precomputed objects,
%                                as returned by surferkernevalmat with
%                                opts.ifreturnovers=1
%           opts.usefmm      = (1) use FMM for smooth quadrature
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
objover = [];
if isfield(opts,'objover'), objover = opts.objover; end
usefmm = 1;
if isfield(opts,'usefmm'), usefmm = opts.usefmm; end

if isempty(cors)
    coropts = opts;
    coropts.corrections = 1;
    [cors,objover] = surferkernevalmat(surferobj,kern, targobj, eps,coropts);
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
    opts.objover     = objover;
    for i = 1:size(dens,2)
        pot(:,i) = surferkerneval(surferobj, kern, dens(:,i), targinfo, eps, opts);
    end
    return
end

% --- Parse objover into a canonical form ---
%
% objover_mode == false : objover is a length-nsurfers cell of order vectors
% objover_mode == true  : surfers_over and xinterps are precomputed length-nsurfers cells

objover_mode = false;
surfers_over = [];
xinterps     = [];

if iscell(objover) && numel(objover) == 2 && iscell(objover{1})
    % 2-element cell: {surfers_over, xinterps}
    objover_mode    = true;
    surfers_over_in = objover{1};
    xinterps_in     = objover{2};

    % Normalise to length-nsurfers cell vectors
    if ~iscell(surfers_over_in)
        surfers_over = num2cell(surfers_over_in(:));
    else
        surfers_over = surfers_over_in(:);
    end
    if ~iscell(xinterps_in)
        xinterps = num2cell(xinterps_in(:));
    else
        xinterps = xinterps_in(:);
    end

elseif ~isempty(objover) && ~iscell(objover)
    % Plain vector: broadcast to all surfers
    novers_vec = objover(:);
    objover = cell(nsurfers, 1);
    for j = 1:nsurfers
        objover{j} = novers_vec;
    end
end

if ~objover_mode
    novers_tmp = cell(nsurfers, 1);
    for j = 1:nsurfers
        lsurfer(j) = surfers{j}.npts;
        if numel(kern) == 1
            ktmp = kern;
        else
            ktmp = kern(j);
        end
        opdims_mat(:,j) = ktmp.opdims;
        if isempty(objover)
            if ismethod(surfers{j},'oversample')
                novers_tmp{j} = ktmp.get_overs_orders(surfers{j},targobj,eps);
            else
                novers_tmp{j} = NaN;
            end
        end
    end
    if isempty(objover)
        objover = novers_tmp;
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

    if objover_mode
            surferjover = surfers_over{j};
            xinterp     = xinterps{j};
        else
            noversj = objover{j};
            if ismethod(surferj,'oversample') && ~any(isnan(noversj))
                [surferjover, xinterp] = surferj.oversample(noversj);
                xinterp = kron(xinterp, eye(ktmp.opdims(2)));
            else
                surferjover = surferj;
                xinterp = eye(numel(wts));
            end
        end
        wtsover = repmat(surferjover.wts(:).', ktmp.opdims(2), 1);
        dens_j = wtsover(:).*(xinterp*dens_j);
        if usefmm && ~isempty(ktmp.fmm)
            pot_j = ktmp.fmm(eps,surferjover,targinfo,dens_j);
        else
            pot_j = ktmp.eval(surferjover,targinfo)*dens_j;
        end
    pot = pot + pot_j(:);
end

%% Add near-field correction
if ~isempty(cors)
    pot = pot + cors * dens;
end
end
