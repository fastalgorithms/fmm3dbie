function mat = kernbyindexeval(i, j, surfers, kern, targobj, eps, novers, Qsparse, opts)
%KERNBYINDEXEVAL  Evaluate smooth-rule off-surface matrix entries by index.
%
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, eps)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, eps, novers)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, eps, novers, Qsparse)
% mat = byindex.kernbyindexeval(i, j, surfers, kern, targobj, eps, novers, Qsparse, opts)
%
% Analogous to byindex.kernbyindex but targets are off-surface points in targobj.
%
% Input:
%   i       - row indices into the flat target DOF array (length ntarg*opdim1)
%   j       - col indices into the source DOF array
%   surfers - surfer object or array of surfer objects (sources)
%   kern    - kernel3d object, scalar or vector of length nsurfers
%   targobj - off-surface targets: struct with .r (3 x ntarg)
%             or numeric (3 x ntarg) array
%   eps     - quadrature tolerance
%   novers  - (optional) cell(nsurfers,1) of per-patch oversampling orders
%   Qsparse - (optional) sparse correction matrix (ntarg*opdim1 x ncols)
%   opts    - options struct
%     opts.replace_quadcorr - logical (default true)
%
% Output:
%   mat - length(i) x length(j) submatrix

if ~isa(kern, 'kernel3d')
    error('BYINDEX.KERNBYINDEXEVAL: kern must be a kernel3d object');
end
if nargin < 7, novers  = []; end
if nargin < 8, Qsparse = []; end
if nargin < 9, opts    = []; end

replace_quadcorr = true;
if isstruct(opts) && isfield(opts, 'replace_quadcorr')
    replace_quadcorr = opts.replace_quadcorr;
end

if isnumeric(targobj)
    targinfo_full.r = targobj;
else
    targinfo_full = targobj;
end

nsurfers = length(surfers);

if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers);
end

assert(all(opdims_mat(1,:) == opdims_mat(1,1)), ...
    'BYINDEX.KERNBYINDEXEVAL: opdims(1) must be constant across all source surfers');

opdim1 = opdims_mat(1,1);

icollocs = zeros(nsurfers+1, 1); icollocs(1) = 1;
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2,k);
end

% Decode row indices -> unique target points
ntarg = size(targinfo_full.r, 2);
ipts  = double(idivide(int64(i(:) - 1), int64(opdim1))) + 1;
[iuni, ~, iiuni] = unique(ipts);
iiuni2 = (iiuni - 1)*opdim1 + mod(i(:) - 1, opdim1) + 1;

targinfo = [];
targinfo.r = targinfo_full.r(:, iuni);
for f = {'n','du','dv'}
    if isfield(targinfo_full, f{1})
        targinfo.(f{1}) = targinfo_full.(f{1})(:, iuni);
    end
end
ntrg_uni = length(iuni);

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j), return; end

for isrc = 1:nsurfers
    srfj = surfers(isrc);

    mat_col_start = icollocs(isrc);
    mat_col_end   = icollocs(isrc+1) - 1;
    flag_j = (j >= mat_col_start) & (j <= mat_col_end);
    mat_col_ind = j(flag_j);
    if isempty(mat_col_ind), continue; end

    if numel(kern) == 1, ktmp = kern; else, ktmp = kern(isrc); end
    opdims_src = opdims_mat(2, isrc);

    jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
    [juni, ~, ijuni] = unique(jpts);
    nsrc_uni = length(juni);
    ijuni2 = (ijuni - 1)*opdims_src + mod(mat_col_ind(:) - mat_col_start, opdims_src) + 1;

    if isempty(novers)
        novers_j = ktmp.get_overs_orders(srfj, targinfo_full, eps);
    else
        novers_j = novers{isrc}(:);
    end
    norders_j = srfj.norders(:);

    if all(novers_j == norders_j) || any(isnan(novers_j))
        srcp   = slice_surfer(srfj, juni, ktmp.src_fields);
        wts    = repmat(srfj.wts(juni(:)).', opdims_src, 1);
        wts    = wts(:).';
        matuni = ktmp.eval(srcp, targinfo) .* wts;
        mat(true(size(i)), flag_j) = matuni(iiuni2, ijuni2);
    else
        matuni = zeros(ntrg_uni * opdim1, nsrc_uni * opdims_src);

        ntmp = [norders_j, novers_j, srfj.iptype(:)];
        [ntmp_uni, ~, intmp] = unique(ntmp, 'rows');

        vmat_cache    = cell(size(ntmp_uni,1), 1);
        xinterp_cache = cell(size(ntmp_uni,1), 1);
        wts_raw_cache = cell(size(ntmp_uni,1), 1);

        for tt = 1:size(ntmp_uni,1)
            norder = ntmp_uni(tt,1); nover = ntmp_uni(tt,2); iptype = ntmp_uni(tt,3);
            if iptype == 1
                rnodes      = koorn.rv_nodes(norder);
                rnodes_over = koorn.rv_nodes(nover);
                vmat = koorn.coefs2vals(norder, rnodes_over);
                umat = koorn.vals2coefs(norder, rnodes);
                wts_raw_cache{tt} = koorn.rv_weights(nover);
            elseif iptype == 11
                rnodes      = polytens.lege.nodes(norder);
                rnodes_over = polytens.lege.nodes(nover);
                vmat = polytens.lege.coefs2vals(norder, rnodes_over);
                umat = polytens.lege.vals2coefs(norder, rnodes);
                wts_raw_cache{tt} = polytens.lege.weights(nover);
            elseif iptype == 12
                rnodes      = polytens.cheb.nodes(norder);
                rnodes_over = polytens.cheb.nodes(nover);
                vmat = polytens.cheb.coefs2vals(norder, rnodes_over);
                umat = polytens.cheb.vals2coefs(norder, rnodes);
                wts_raw_cache{tt} = polytens.cheb.weights(nover);
            end
            vmat_cache{tt} = vmat;
            xi = vmat * umat;
            [lia, iindb] = ismembertol(rnodes', rnodes_over', 1e-7, 'ByRows', true);
            isp = [find(lia), iindb(lia)];
            if ~isempty(isp)
                i1 = isp(:,1); i2 = isp(:,2);
                xi(i2,:) = 0; xi(i2,i1) = 1;
            end
            xinterp_cache{tt} = xi;
        end

        for p = 1:srfj.npatches
            p_inds = srfj.ixyzs(p):(srfj.ixyzs(p+1)-1);
            [tf, loc_in_juni] = ismember(juni, p_inds);
            if ~any(tf), continue; end
            loc_juni_p = find(tf);
            loc_in_p   = loc_in_juni(tf);

            tt   = intmp(p);
            vmat = vmat_cache{tt};

            srcover_vals = srfj.srccoefs{p} * vmat';
            ru   = srcover_vals(4:6,:);
            rv   = srcover_vals(7:9,:);
            rtmp = cross(ru, rv);
            jac  = vecnorm(rtmp, 2);
            nhat = rtmp ./ jac;

            srcp_over = [];
            srcp_over.r  = srcover_vals(1:3,:);
            srcp_over.du = srcover_vals(4:6,:);
            srcp_over.dv = srcover_vals(7:9,:);
            srcp_over.n  = nhat;

            wts_over     = wts_raw_cache{tt}(:) .* jac(:);
            wts_over_rep = repmat(wts_over(:).', opdims_src, 1);
            wts_over_rep = wts_over_rep(:).';

            sub_kern_over = ktmp.eval(srcp_over, targinfo) .* wts_over_rep;
            xinterp_sub   = xinterp_cache{tt}(:, loc_in_p);
            sub_kern_smth = sub_kern_over * kron(xinterp_sub, eye(opdims_src));

            col_inds = expand_cols(loc_juni_p, opdims_src);
            matuni(:, col_inds) = matuni(:, col_inds) + sub_kern_smth;
        end

        mat(true(size(i)), flag_j) = matuni(iiuni2, ijuni2);
    end
end

if ~isempty(Qsparse)
    [isp, jsp, vsp] = find(Qsparse(i(:), j(:)));
    linsp = isp + (jsp-1)*length(i(:));
    if replace_quadcorr
        mat(linsp) = vsp;
    else
        mat(linsp) = mat(linsp) + vsp;
    end
end

end


function srcp = slice_surfer(srfj, pts, src_fields)
    srcp = [];
    srcp.r = srfj.r(:, pts);
    for k = 1:length(src_fields)
        f = src_fields{k};
        srcp.(f) = srfj.(f)(:, pts);
    end
end

function col_inds = expand_cols(loc_juni_p, opdims_src)
    starts   = (loc_juni_p(:) - 1) * opdims_src + 1;
    col_inds = cell2mat(arrayfun(@(s) s:(s+opdims_src-1), starts, 'UniformOutput', false)');
end
