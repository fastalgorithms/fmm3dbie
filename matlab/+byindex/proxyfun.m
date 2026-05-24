function [Kpxy, nbr] = proxyfun(slf, nbr, l, ctr, surfers, kern, kern2, ...
    opdims_mat, pr, pn, pw, pin, ifaddtrans)
%PROXYFUN  Proxy function for rskelf, for kernels defined on arrays of surfers.
%
% [Kpxy, nbr] = byindex.proxyfun(slf, nbr, l, ctr, surfers, kern, kern2,
%     opdims_mat, pr, pn, pw, pin, ifaddtrans)
%
% pr, pn, pw  - (3,npxy) proxy points/normals and (npxy,1) weights, unit-scale
% pin         - function handle: logical mask for points inside proxy surface,
%               called on (3,n) rescaled/recentered point array
% kern        - kernel3d, sources on surfers -> proxy targets
% kern2       - kernel3d, proxy sources -> targets on surfers (transpose block)
%               pass [] or kern to use kern for both
% ifaddtrans  - if true, append transpose block

if isempty(kern2), kern2 = kern; end

nsurfers = length(surfers);

if ~(ndims(opdims_mat) == 3 && isequal(size(opdims_mat), [2, nsurfers, nsurfers]))
    v = opdims_mat(:);
    tmp = zeros(2, nsurfers, nsurfers);
    for ii_ = 1:nsurfers
        for jj_ = 1:nsurfers
            tmp(1,ii_,jj_) = v(1); tmp(2,ii_,jj_) = v(2);
        end
    end
    opdims_mat = tmp;
end

% Column (source) offsets
icollocs = zeros(nsurfers+1, 1); icollocs(1) = 1;
for k = 1:nsurfers
    icollocs(k+1) = icollocs(k) + surfers(k).npts * opdims_mat(2, 1, k);
end

% Row (target) offsets — used for nbr filtering
irowlocs = zeros(nsurfers+1, 1); irowlocs(1) = 1;
for k = 1:nsurfers
    irowlocs(k+1) = irowlocs(k) + surfers(k).npts * opdims_mat(1, k, 1);
end

% Scale proxy points/weights
pxy = pr * l + ctr(:);
pw  = pw  * l;
npxy = size(pxy, 2);

targinfo_pxy = []; targinfo_pxy.r = pxy; targinfo_pxy.n = pn;
srcinfo_pxy  = []; srcinfo_pxy.r  = pxy; srcinfo_pxy.n  = pn;

% Filter nbr to points inside proxy surface
new_nbr = [];
for k = 1:nsurfers
    opdim = opdims_mat(1, k, 1);
    fnbr = nbr >= irowlocs(k) & nbr < irowlocs(k+1);
    inbr = nbr(fnbr);
    if isempty(inbr), continue; end
    pts = idivide(int64(inbr(:) - irowlocs(k)), int64(opdim)) + 1;
    dr  = (surfers(k).r(:, pts) - ctr(:)) / l;
    new_nbr = [new_nbr, inbr(pin(dr))];
end
nbr = new_nbr;

% Proxy row offsets
npxy_rows = zeros(nsurfers+1, 1); npxy_rows(1) = 1;
for k = 1:nsurfers
    npxy_rows(k+1) = npxy_rows(k) + npxy * opdims_mat(1, k, 1);
end
Kpxy = zeros(npxy_rows(end)-1, length(slf));

if ifaddtrans
    npxy_rows2 = zeros(nsurfers+1,1); npxy_rows2(1) = 1;
    for k = 1:nsurfers
        npxy_rows2(k+1) = npxy_rows2(k) + npxy * opdims_mat(2, k, 1);
    end
    Ktrans = zeros(length(slf), npxy_rows2(end)-1);
end

for isrc = 1:nsurfers
    srfj = surfers(isrc);
    mat_col_start = icollocs(isrc);
    mat_col_end   = icollocs(isrc+1) - 1;
    f_col = slf >= mat_col_start & slf <= mat_col_end;
    mat_col_ind = slf(f_col);
    if isempty(mat_col_ind), continue; end

    opdims_src = opdims_mat(2, 1, isrc);
    jpts = idivide(int64(mat_col_ind(:) - mat_col_start), int64(opdims_src)) + 1;
    [juni, ~, ijuni] = unique(jpts);
    ijuni2 = (ijuni-1)*opdims_src + mod(mat_col_ind(:)-mat_col_start, opdims_src) + 1;

    srcp = slice_surfer(srfj, juni);
    wsrc = repmat(srfj.wts(juni(:)).', opdims_src, 1);
    wsrc = wsrc(:).';

    for itrg = 1:nsurfers
        if numel(kern) == 1, ktmp = kern; else, ktmp = kern(itrg,isrc); end

        matuni = ktmp.eval(srcp, targinfo_pxy) .* wsrc;
        Kpxy(npxy_rows(itrg):npxy_rows(itrg+1)-1, f_col) = matuni(:, ijuni2);

        if ifaddtrans
            if numel(kern2) == 1, ktmp2 = kern2; else, ktmp2 = kern2(isrc,itrg); end
            wpxy = repmat(pw(:).', opdims_mat(2,isrc,itrg), 1);
            wpxy = wpxy(:).';
            matuni2 = ktmp2.eval(srcinfo_pxy, srcp) .* wpxy;
            Ktrans(f_col, npxy_rows2(isrc):npxy_rows2(isrc+1)-1) = matuni2(ijuni2, :);
        end
    end
end

if ifaddtrans
    Kpxy = [Kpxy; Ktrans.'];
end

end


function srcp = slice_surfer(srfj, pts)
    src_field = fieldnames(srfj)';
    srcp = [];
    for i = 1:length(src_field)
        try, val = srfj.(src_field{i}); catch, continue; end
        if isempty(val), continue; end
        if mod(size(val(:,:), 2), srfj.npts) == 0
            srcp.(src_field{i}) = srfj.(src_field{i})(:, pts);
        end
    end
end
