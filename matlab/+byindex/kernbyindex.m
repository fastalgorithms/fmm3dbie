function mat = kernbyindex(i, j, surfers, kern, eps, novers, Qsparse, opts)
%KERNBYINDEX  Evaluate smooth-rule matrix entries by index for surfer arrays.
%
% mat = byindex.kernbyindex(i, j, surfers, kern, eps)
% mat = byindex.kernbyindex(i, j, surfers, kern, eps, novers)
% mat = byindex.kernbyindex(i, j, surfers, kern, eps, novers, Qsparse)
% mat = byindex.kernbyindex(i, j, surfers, kern, eps, novers, Qsparse, opts)
%
% Input:
%   i, j      - row/col indices to evaluate
%   surfers   - surfer object or array of surfer objects
%   kern      - kernel3d object (scalar or nsurfers x nsurfers)
%   eps       - precision requested (e.g. 1e-10)
%   novers    - (optional) cell(nsurfers,nsurfers) of per-patch oversampling
%               order vectors, as returned by surfermat.  novers{itrg,isrc}
%               is a vector of length surfers(isrc).npatches.
%               If omitted or empty, kern.get_overs_orders is called per block.
%   Qsparse   - (optional) sparse matrix of quadrature corrections
%               (same size as the full matrix).  Pass [] to skip.
%               Two modes, selected by opts.replace_quadcorr:
%                 true (default): replacement, i.e. overwrite smooth entries
%                   with Qsparse values where nonzero.  Use with nonsmoothonly=1
%                   output from surfermat, which stores the full near-singular
%                   quadrature values directly.
%                 false: additive, i.e. mat += Qsparse.  Use with corrections=1
%                   output from surfermat, which stores (near-singular - smooth)
%                   so that smooth + correction = full.
%   opts      - options struct
%     opts.replace_quadcorr - logical (default true), see Qsparse above
%
% Output:
%   mat - length(i) x length(j) matrix of smooth-rule entries (with singular
%         corrections applied where Qsparse is provided)

if ~isa(kern, 'kernel3d')
    error('BYINDEX.KERNBYINDEX: kern must be a kernel3d object');
end

if nargin < 5
    error('BYINDEX.KERNBYINDEX: eps is required');
end
if nargin < 6, novers  = []; end
if nargin < 7, Qsparse = []; end
if nargin < 8, opts    = []; end

replace_quadcorr = true;
if isstruct(opts) && isfield(opts, 'replace_quadcorr')
    replace_quadcorr = opts.replace_quadcorr;
end

quad_eps = eps;

nsurfers = length(surfers);

% Build opdims_mat (2 x nsurfers x nsurfers) from kern
if numel(kern) == 1
    opdims_mat = repmat(reshape(kern.opdims, 2, 1, 1), 1, nsurfers, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers, nsurfers);
end

use_get_overs = isempty(novers);

% Build row/col offset tables
irowlocs = zeros(nsurfers+1, 1);
icollocs = zeros(nsurfers+1, 1);
irowlocs(1) = 1;
icollocs(1) = 1;
for k = 1:nsurfers
    irowlocs(k+1) = irowlocs(k) + surfers(k).npts * opdims_mat(1, k, 1);
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end

mat = zeros(length(i), length(j));
if isempty(i) || isempty(j)
    return;
end

for itrg = 1:nsurfers
    srfi = surfers(itrg);

    mat_row_start = irowlocs(itrg);
    mat_row_end   = irowlocs(itrg+1) - 1;
    flag_i = (i >= mat_row_start) & (i <= mat_row_end);
    mat_row_ind = i(flag_i);
    if isempty(mat_row_ind), continue; end

    opdims_row = opdims_mat(1, itrg, 1);

    ipts = idivide(int64(mat_row_ind(:) - mat_row_start), int64(opdims_row)) + 1;
    [iuni, ~, iiuni] = unique(ipts);
    iiuni2 = (iiuni - 1)*opdims_row + mod(mat_row_ind(:) - mat_row_start, opdims_row) + 1;

    targinfo = [];
    targinfo.r  = srfi.r(:, iuni);
    targinfo.n  = srfi.n(:, iuni);
    targinfo.du = srfi.du(:, iuni);
    targinfo.dv = srfi.dv(:, iuni);

    ntrg_uni = length(iuni);

    % Full target struct for get_overs_orders
    targs_full = [];
    targs_full.r  = srfi.r;
    targs_full.n  = srfi.n;
    targs_full.du = srfi.du;
    targs_full.dv = srfi.dv;

    for isrc = 1:nsurfers
        srfj = surfers(isrc);

        mat_col_start = icollocs(isrc);
        mat_col_end   = icollocs(isrc+1) - 1;
        flag_j = (j >= mat_col_start) & (j <= mat_col_end);
        mat_col_ind = j(flag_j);
        if isempty(mat_col_ind), continue; end

        opdims_src = opdims_mat(2, itrg, isrc);

        if numel(kern) == 1
            ktmp = kern;
        else
            ktmp = kern(itrg, isrc);
        end

        jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
        [juni, ~, ijuni] = unique(jpts);
        nsrc_uni = length(juni);
        ijuni2 = (ijuni - 1)*opdims_src + mod(mat_col_ind(:) - mat_col_start, opdims_src) + 1;

        % Get per-patch oversampling orders for this block
        if use_get_overs
            if itrg == isrc
                targ_for_overs = srfi;
            else
                targ_for_overs = targs_full;
            end
            novers_j = ktmp.get_overs_orders(srfj, targ_for_overs, quad_eps);
        else
            novers_j = novers{itrg, isrc}(:);
        end
        norders_j = srfj.norders(:);

        if all(novers_j == norders_j) || any(isnan(novers_j))
            % No oversampling needed: evaluate at original nodes directly
            srcp = slice_surfer(srfj, juni, ktmp.src_fields);
            wts  = repmat(srfj.wts(juni(:)).', opdims_src, 1);
            wts  = wts(:).';
            matuni = ktmp.eval(srcp, targinfo) .* wts;

            if itrg == isrc
                for p = 1:srfj.npatches
                    p_inds = srfj.ixyzs(p):(srfj.ixyzs(p+1)-1);
                    src_in_p = find(ismember(juni, p_inds));
                    trg_in_p = find(ismember(iuni, p_inds));
                    if isempty(src_in_p) || isempty(trg_in_p), continue; end
                    src_r_p = srfj.r(:, juni(src_in_p));
                    trg_r_p = targinfo.r(:, trg_in_p);
                    src_norm_p = max(vecnorm(src_r_p), 1);
                    for q = 1:length(src_in_p)
                        diff_norms = vecnorm(trg_r_p - src_r_p(:,q));
                        self_mask = diff_norms < 1e-14 * src_norm_p(q);
                        if any(self_mask)
                            row_off = (1:opdims_row) + opdims_row*(trg_in_p(self_mask).' - 1);
                            col_off = (1:opdims_src) + opdims_src*(src_in_p(q) - 1);
                            matuni(row_off(:), col_off) = 0;
                        end
                    end
                end
            end

            mat(flag_i, flag_j) = matuni(iiuni2, ijuni2);

        else
            % Oversampling: patchwise accumulation
            matuni = zeros(ntrg_uni * opdims_row, nsrc_uni * opdims_src);

            % Cache interp matrices per unique (norder, nover, iptype)
            ntmp = [norders_j, novers_j, srfj.iptype(:)];
            [ntmp_uni, ~, intmp] = unique(ntmp, 'rows');

            vmat_cache    = cell(size(ntmp_uni, 1), 1);
            xinterp_cache = cell(size(ntmp_uni, 1), 1);
            wts_raw_cache = cell(size(ntmp_uni, 1), 1);

            for tt = 1:size(ntmp_uni, 1)
                norder = ntmp_uni(tt, 1);
                nover  = ntmp_uni(tt, 2);
                iptype = ntmp_uni(tt, 3);

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

                if norder > nover
                    [lia, iindb] = ismembertol(rnodes_over', rnodes', 1e-7, 'ByRows', true);
                    isp = [iindb(lia), find(lia)];
                else
                    [lia, iindb] = ismembertol(rnodes', rnodes_over', 1e-7, 'ByRows', true);
                    isp = [find(lia), iindb(lia)];
                end
                if ~isempty(isp)
                    i1 = isp(:,1); i2 = isp(:,2);
                    xi(i2, :) = 0;
                    xi(i2, i1) = 1;
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

                wts_over = wts_raw_cache{tt}(:) .* jac(:);
                wts_over_rep = repmat(wts_over(:).', opdims_src, 1);
                wts_over_rep = wts_over_rep(:).';

                sub_kern_over = ktmp.eval(srcp_over, targinfo) .* wts_over_rep;

                if itrg == isrc
                    trg_in_p = find(ismember(iuni, p_inds));
                    if ~isempty(trg_in_p)
                        trg_r_p = targinfo.r(:, trg_in_p);
                        src_norm_p = max(vecnorm(srcp_over.r), 1);
                        for q = 1:size(srcp_over.r, 2)
                            diff_norms = vecnorm(trg_r_p - srcp_over.r(:,q));
                            self_mask = diff_norms < 1e-14 * src_norm_p(q);
                            if any(self_mask)
                                row_off = (1:opdims_row) + opdims_row*(trg_in_p(self_mask).' - 1);
                                col_off = (1:opdims_src) + opdims_src*(q-1);
                                sub_kern_over(row_off(:), col_off) = 0;
                            end
                        end
                    end
                end

                xinterp_sub = xinterp_cache{tt}(:, loc_in_p);
                sub_kern_smth = sub_kern_over * kron(xinterp_sub, eye(opdims_src));

                col_inds = expand_cols(loc_juni_p, opdims_src);
                matuni(:, col_inds) = matuni(:, col_inds) + sub_kern_smth;
            end

            mat(flag_i, flag_j) = matuni(iiuni2, ijuni2);
        end
    end
end

% Apply singular quadrature corrections
if ~isempty(Qsparse)
    P = zeros(max(i), 1);
    corr = spget_quadcorr(i(:), j(:), P, Qsparse);
    if replace_quadcorr
        % nonsmoothonly mode: Qsparse holds full near-singular values;
        % replace smooth entries wherever Qsparse is nonzero
        mask = (corr ~= 0);
        mat(mask) = corr(mask);
    else
        % corrections mode: Qsparse holds (near-singular - smooth);
        % add to smooth entries
        mat = mat + corr;
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
    starts = (loc_juni_p(:) - 1) * opdims_src + 1;
    col_inds = cell2mat(arrayfun(@(s) s:(s+opdims_src-1), starts, ...
                                 'UniformOutput', false)');
end
