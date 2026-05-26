function U = fmm(eps, src, targ, type, sigma, coefs)
%STOK3D.FMM   Fast multipole method for evaluating Stokes layer
%   potentials in 3D.
%
% Syntax:
%   U  = stok3d.fmm(eps, zk, srcinfo, targinfo, type, sigma, coefs)
%
% Input:
%   eps      - precision requested
%   src      - ptinfo struct: .r (3,:), .n (3,:)
%   targ     - ptinfo struct: .r (3,:), .n (3,:)
%   type     - 's'           single layer S
%              'd'           double layer D
%              'c'           combined layer coefs(1)*S + coefs(2)*D
%   sigma    - density (3 x ns)
%   coefs    - [alpha; beta] for combined layer (type 'c' only)
%
% Output:
%   pot  - potential (3 x nt)

switch lower(type)
    case {'s', 'single'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        srcinfo.stoklet = sigma;
    case {'d', 'double', 'traction'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        srcinfo.strslet = sigma;
        srcinfo.strsvec = src.n;
    case {'c', 'combined'}
        srcinfo = struct();
        srcinfo.sources = src.r;
        if coefs(1) ~= 0
            srcinfo.stoklet = coefs(1) * sigma;
        end
        if coefs(2) ~= 0
            srcinfo.strslet = coefs(2) * sigma;
            srcinfo.strsvec = src.n;
        end
end
U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
U = reshape(U_out.pottarg, 3, []);
end