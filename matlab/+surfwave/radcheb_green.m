function R = radcheb_green(vgfun, rmax, opts)
%SURFWAVE.RADCHEB_GREEN  Radial Chebyshev representation of a Green's function
%   and its gradient.
%
%   R = surfwave.radcheb_green(vgfun, rmax) or (vgfun, rmax, opts)
%
%   VGFUN(src,targ) is the usual evaluator of a radial Green's function
%   G(|x-y|), e.g. @(s,t) surfwave.flex.gsflex(rts,ejs,s,t), returning
%      val  - (nt,ns) values of G
%      grad - (nt,ns,2) gradient with respect to the TARGET, grad(:,:,k) =
%             dG/dx_k = G'(r) (targ_k - src_k)/r
%   (further outputs are ignored). src/targ are (2,:) or (3,:) arrays or
%   ptinfo structs.
%
%   Two kernel3d.radcheb interpolants are built along the x axis, one for
%   G(r) and one for G'(r) (fit directly, not by differentiating the first,
%   so the gradient is as accurate as the value). Then
%
%      [val, grad] = R.eval(src, targ)
%
%   has the same format as VGFUN, so a kernel that is written in terms of
%   [val, grad] = gsflex(...) (and takes dot products with normals and
%   tangents) can use R.eval instead, see surfwave.flex.kernel_eval.
%
%   RMAX is the largest radius needed, or [rmin, rmax]. OPTS is passed to
%   kernel3d.radcheb (opts.eps, opts.norder, ...); opts.rmin (default 1e-12)
%   is the smallest radius. As in radcheb_eval, radii outside [rmin, rmax]
%   give NaN.
%
%   R is a struct with fields eval, rmin, rmax, breaks_val, coefs_val,
%   breaks_der, coefs_der.
%
% See also KERNEL3D.RADCHEB, RADCHEB_EVAL, SURFWAVE.CAPILLARY.RADCHEB_KERNELS

if nargin < 3 || isempty(opts), opts = struct(); end
rmin = 1e-12;
if isfield(opts, 'rmin'), rmin = opts.rmin; end
if numel(rmax) == 2
    rmin = rmax(1);  rmax = rmax(2);
end

kv = kernel3d.radcheb(@(s,t) rg_probe(vgfun, s, t, 0), [rmin rmax], opts);
kd = kernel3d.radcheb(@(s,t) rg_probe(vgfun, s, t, 1), [rmin rmax], opts);

R = [];
R.rmin       = rmin;
R.rmax       = rmax;
R.breaks_val = kv.params.breaks;   R.coefs_val = kv.params.coefs;
R.breaks_der = kd.params.breaks;   R.coefs_der = kd.params.coefs;
R.err        = [kv.params.err, kd.params.err];
R.eval = @(src, targ) rg_eval(src, targ, R.breaks_val, R.coefs_val, ...
                              R.breaks_der, R.coefs_der, rmin, rmax);

end


function out = rg_probe(vgfun, s, t, which)
%  source at the origin, targets along the x axis: dG/dx = G'(r)
[v, g] = vgfun(s, t);
if which == 0
    out = v(:, 1);
else
    out = g(:, 1, 1);
end
end


function [val, grad] = rg_eval(src, targ, bv, cv, bd, cd, ~, ~)

if isstruct(src),  src  = src.r;  end
if isstruct(targ), targ = targ.r; end

dx = targ(1,:).' - src(1,:);
dy = targ(2,:).' - src(2,:);
r  = sqrt(dx.^2 + dy.^2);

val = radcheb_eval(r, bv, cv);

if nargout > 1
    fac = radcheb_eval(r, bd, cd) ./ r;
    grad = cat(3, fac.*dx, fac.*dy);
end

end
