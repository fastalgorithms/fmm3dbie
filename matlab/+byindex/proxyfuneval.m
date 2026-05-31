function [Kpxy, nbr] = proxyfuneval(slf, nbr, l, ctr, surfers, kern, targr, pr, pn, pw, pin)
%PROXYFUNEVAL  Proxy function for off-surface evaluation (no ifaddtrans).
%
% [Kpxy, nbr] = byindex.proxyfuneval(slf, nbr, l, ctr, surfers, kern,
%                                    targr, pr, pn, pw, pin)
%
% Analogous to byindex.proxyfun but nbr indexes off-surface target points
% (rows of a flat target DOF array) rather than on-surface DOF rows.
% No ifaddtrans: the transpose block is not needed for off-surface eval.
%
% Input:
%   slf, nbr  - source/target DOF indices (as in proxyfun)
%   l, ctr    - box half-length and center
%   surfers   - surfer object or array (sources)
%   kern      - kernel3d, scalar or vector of length nsurfers
%   targr     - (3, ntarg) array of off-surface target positions
%   pr        - (3,npxy) proxy points, unit-scale
%   pn        - (3,npxy) proxy normals
%   pw        - (npxy,1) proxy weights
%   pin       - function handle: logical mask for points inside proxy surface,
%               called on (3,n) unit-scale point array
%
% Output:
%   Kpxy - (npxy*opdim1 x nslf) proxy matrix
%   nbr  - filtered nbr indices (off-surface target DOFs inside proxy removed)

nsurfers = length(surfers);

if numel(kern) == 1
    opdims_mat = repmat(kern.opdims(:), 1, nsurfers);
else
    opdims_mat = reshape([kern.opdims], 2, nsurfers);
end

opdim1 = opdims_mat(1,1);

% Column (source) offsets
icollocs = zeros(nsurfers+1, 1); icollocs(1) = 1;
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2,k);
end

% Scale proxy points/weights
pxy = pr * l + ctr(:);
pw  = pw  * l;
npxy = size(pxy, 2);

targinfo_pxy = []; targinfo_pxy.r = pxy; targinfo_pxy.n = pn;

% Filter nbr: remove off-surface target DOFs inside the proxy surface.
nbr_pts = idivide(int64(nbr(:)) - 1, int64(opdim1)) + 1;
dr  = (targr(:, nbr_pts) - ctr(:)) / l;
nbr = nbr(pin(dr));

Kpxy = zeros(npxy * opdim1, length(slf));

for isrc = 1:nsurfers
    srfj = surfers(isrc);
    mat_col_start = icollocs(isrc);
    mat_col_end   = icollocs(isrc+1) - 1;
    f_col = slf >= mat_col_start & slf <= mat_col_end;
    mat_col_ind = slf(f_col);
    if isempty(mat_col_ind), continue; end

    if numel(kern) == 1, ktmp = kern; else, ktmp = kern(isrc); end
    opdims_src = opdims_mat(2, isrc);

    jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
    [juni, ~, ijuni] = unique(jpts);
    ijuni2 = (ijuni-1)*opdims_src + mod(mat_col_ind(:)-mat_col_start, opdims_src) + 1;

    srcp = slice_surfer(srfj, juni, ktmp.src_fields);
    wsrc = repmat(srfj.wts(juni(:)).', opdims_src, 1);
    wsrc = wsrc(:).';

    matuni = ktmp.eval(srcp, targinfo_pxy) .* wsrc;  % (npxy*opdim1) x nsrc_uni
    Kpxy(:, f_col) = matuni(:, ijuni2);
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
