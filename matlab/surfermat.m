function [sysmat,varargout] = surfermat(surferobj,kern,opts)
% [sysmat, opts, novers_out] = surfermat(...)
%   novers_out is a cell(nsurfers,nsurfers) of per-patch oversampling order
%   vectors, one per block, as returned by kern.get_overs_orders.
%SURFERERMAT build matrix for given kernel and surfer description of 
% boundary. This is a wrapper for various quadrature routines. Optionally,
% return only those interactions which do not use the smooth integration
% rule in the sparse matrix format.
%
% Syntax: sysmat = surverermat(S,kern,opts)
%
% Input:
%   Surferobj - array of surfer objects describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.n - unit normals (2,:)
%
% Optional input:
%   opts  - options structure. available options (default settings)%           
%
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactions) and return
%                         in a sparse array.
%           opts.corrections = boolean (false), if true, only compute the
%                         corrections to the smooth quadrature rule and 
%                         return in a sparse array, see opts.nonsmoothonly
%           opts.l2scale = boolean (false), if true scale rows by 
%                           sqrt(whts) and columns by 1/sqrt(whts)
%           opts.quad = quad type 'kern' to call kern.get_quad, or 'native'
%               to use the smooth rule
%           opts.eps = (1e-10) tolerance for adaptive quadrature
%
% Output:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by S
%
%
% Examples:
%   sysmat = surfermat(S,kern); % standard options
%   sysmat = surfermat(S,kern,opts);
%
% Author: Tristan Goodwill

surfers = surferobj;

if ~isa(kern,'kernel3d')
        error('SURFERMAT: second input kern not of supported type');
end

if nargin < 3
    opts = [];
end

nonsmoothonly = false;
corrections = false;
l2scale = false;
selfquad = true;

% get opts from struct if available
if isfield(opts,'selfquad')
    selfquad = opts.selfquad;
end
if isfield(opts,'l2scale')
    l2scale = opts.l2scale;
end
if isfield(opts,'nonsmoothonly')
    nonsmoothonly = opts.nonsmoothonly;
end
if isfield(opts,'corrections')
    corrections = opts.corrections;
end
if corrections
    nonsmoothonly = true;
end
adaptive_correction = true;
if isfield(opts,'adaptive_correction')
    adaptive_correction = opts.adaptive_correction;
end
quad_eps = eps();
if isfield(opts,'eps')
    quad_eps = opts.eps;
end

nsurfers = length(surfers);

opdims_mat = zeros(2,nsurfers,nsurfers);
lsurfer    = zeros(nsurfers,1);

for i=1:nsurfers
    lsurfer(i) = surfers(i).npts;
    
    for j=1:nsurfers
        if (size(kern) == 1)
            opdims_mat(:,i,j) = kern.opdims;
        else
            opdims_mat(:,i,j) = kern(i,j).opdims;
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

novers_out = cell(nsurfers, nsurfers);

for i = 1:nsurfers
    surferi = surfers(i);
    for j = 1:nsurfers
        surferj = surfers(j);
        if (surferi.npts < 1 || surferj.npts<1)
            sysmat_tmp = [];
            break
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

        opdims = reshape(opdims_mat(:,i,j),[2,1]);
        wts = repmat(surferj.wts(:).', opdims(2), 1);
        wts = wts(:).';

        if (size(kern) == 1)
            ktmp = kern;
        else
            ktmp = kern(i,j);
        end

        novers = ktmp.get_overs_orders(surferj,surferi_targ,quad_eps);
        novers_out{i,j} = novers;
        if ismethod(surferj,'oversample')
            [surferjover, xinterp] = surferj.oversample(novers);
            xinterp = kron(xinterp, eye(opdims(2)));
        else
            surferjover = surferj;
            xinterp = eye(numel(wts));
        end

        wtsover = repmat(surferjover.wts(:).', opdims(2), 1);
        wtsover = wtsover(:).';

        if (~nonsmoothonly)
            sysmat_tmp = (ktmp.eval(surferjover,surferi_targ).*wtsover);
            % Zero out entries where oversampled source and target coincide.
            if i == j
                src_norm = max(vecnorm(surferjover.r), 1);
                targ_r   = surferi_targ.r;
                for q = 1:size(surferjover.r, 2)
                    diff_norms = vecnorm(targ_r - surferjover.r(:,q));
                    self_mask  = diff_norms < 1e-14 * src_norm(q);
                    if any(self_mask)
                        row_off = (1:opdims(1)) + opdims(1)*(find(self_mask).' - 1);
                        col_off = (1:opdims(2)) + opdims(2)*(q-1);
                        sysmat_tmp(row_off(:), col_off) = 0;
                    end
                end
            end
            sysmat_tmp = sysmat_tmp*xinterp;
        end

        if (adaptive_correction || (i==j)) && selfquad && ~isempty(ktmp.getquad)
            Qj = ktmp.getquad(surferj,quad_eps,surferi_targ);
            if isfield(Qj, 'row_ptr')
                sysmat_quad = conv_rsc_to_spmat(surferj,Qj.row_ptr,Qj.col_ind,Qj.wnear);
            else
                sysmat_quad = Qj;
            end
            if corrections
                sysmat_smth = smooth_sparse_quad(ktmp, surferi_targ, surferj, Qj.row_ptr, Qj.col_ind, novers);
                sysmat_quad = sysmat_quad-sysmat_smth;
            end

            if (~nonsmoothonly)
                inds = abs(sysmat_quad(:))>0;
                sysmat_tmp(inds) = sysmat_quad(inds);
            else
                sysmat_tmp = sysmat_quad;
            end
        end

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
    % Fix sparse entry format to use rcip matrix entries for repeats
    % instead of using the precomputed self correction
    
    ijind = [isysmat jsysmat];
    [~,idx] = unique(ijind, 'rows','last');
    
    sysmat = sparse(isysmat(idx),jsysmat(idx),vsysmat(idx),nrows,ncols);
end


if (nargout > 1)
    varargout{1} = opts;
end
if (nargout > 2)
    varargout{2} = novers_out;
end

end
