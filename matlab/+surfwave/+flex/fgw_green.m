function [val, grad, hess, third] = fgw_green(S, which, src, targ)
%SURFWAVE.FLEX.FGW_GREEN  finite depth minus deep water parts of the
%   flexural Green's functions, in the chunkie Green's function calling
%   sequence, so they can be added to the deep-water ones in RADPP_SUM for
%   KERNEL3D.RADPP and KERNEL.RADPP:
%
%     [val, grad, hess, third] = surfwave.flex.fgw_green(S, which, src, targ)
%
%   S = SURFWAVE.FLEX.FGW_SETUP(alpha, gamma, depth) (lengths in the units of
%   the problem), and which is one of
%     'gs'          2*G1_diff, the finite-depth part of GSFLEX (gsflex = 2*G1)
%     'gphi'        G2_diff,   the finite-depth part of GPHIFLEX
%     'bilap_gphi'  Delta^2 G2_diff (value only), the finite-depth part of
%                   the bilaplacian of gphi, so that
%                   (alpha/2) Delta^2 Gphi = S3d_sum - gamma/2 Gphi_deep + (alpha/2) bilap_gphi
%   with G1, G2 the deflection and potential kernels of FGW_SETUP. These are
%   smooth (analytic in a strip of width ~2*depth about the real r axis).
%
%   src, targ - (2+, ns) and (2+, nt) arrays or structs with field r; only
%   the first two rows are used. Derivatives are with respect to the target,
%   as in GSFLEX: val (nt,ns), grad (nt,ns,2), hess (nt,ns,3) = xx,xy,yy,
%   third (nt,ns,4) = xxx,xxy,xyy,yyy.
%
%   Example (finite-depth G_S as one radpp kernel):
%     fgw = surfwave.flex.fgw_setup(alpha, gamma, depth);
%     gs  = @(s,t) surfwave.flex.gsflex(rts, ejs, s, t);
%     gsd = @(s,t) surfwave.flex.fgw_green(fgw, 'gs', s, t);
%     K   = kernel3d.radpp(radpp_sum([1 1], gs, gsd), rmax, 's', ...);
%
% See also SURFWAVE.FLEX.FGW_SETUP, SURFWAVE.FLEX.FGW_EVAL_XY, RADPP_SUM

if isstruct(src),  src  = src.r;  end
if isstruct(targ), targ = targ.r; end
ns = size(src, 2);  nt = size(targ, 2);
rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);

switch lower(which)
    case 'bilap_gphi'
        if nargout > 1
            error('SURFWAVE.FLEX.FGW_GREEN: ''bilap_gphi'' is value only.');
        end
        % Delta^2 g = r^4 Th^4 g + 8 r^2 Th^3 g + 8 Th^2 g,  Th = (1/r) d/dr
        r = hypot(rx(:), ry(:));
        [~, T2] = surfwave.flex.fgw_eval(S, r, 'diff', 'theta');
        val = reshape(r.^4.*T2(:,5) + 8*r.^2.*T2(:,4) + 8*T2(:,3), nt, ns);
        return
    case 'gs'
        [D, ~] = surfwave.flex.fgw_eval_xy(S, rx, ry, 'diff');
        D = 2*D;
    case 'gphi'
        [~, D] = surfwave.flex.fgw_eval_xy(S, rx, ry, 'diff');
    otherwise
        error('SURFWAVE.FLEX.FGW_GREEN: unknown kernel ''%s''.', which);
end
Dab = @(a, b) reshape(D(:, a+1, b+1), nt, ns);       % d^(a+b) G / dx^a dy^b
val = Dab(0, 0);
if nargout > 1, grad  = cat(3, Dab(1,0), Dab(0,1)); end
if nargout > 2, hess  = cat(3, Dab(2,0), Dab(1,1), Dab(0,2)); end
if nargout > 3, third = cat(3, Dab(3,0), Dab(2,1), Dab(1,2), Dab(0,3)); end
end
