function [Kgs, Kgp, info] = radcheb_kernels(rts, ejs, rmax, opts)
%SURFWAVE.CAPILLARY.RADCHEB_KERNELS  chunkie kernel objects for the capillary
%   boundary-to-volume operators G_S and G_phi, evaluated through a radial
%   Chebyshev interpolant (kernel3d.radcheb) instead of Struve/Hankel
%   functions.
%
%   [Kgs, Kgp] = surfwave.capillary.radcheb_kernels(rts, ejs, rmax)
%   [Kgs, Kgp, info] = surfwave.capillary.radcheb_kernels(rts, ejs, rmax, opts)
%
%   RMAX is the largest distance between a boundary point and a target, or
%   [rmin, rmax]; radii outside [rmin, rmax] give NaN, as in radcheb_eval, so use a safe bound, e.g. 2*max(vecnorm([targs, chnkr.r(:,:)])).
%   OPTS goes to kernel3d.radcheb (opts.eps, ...; default eps 1e-13).
%
%   Kgs and Kgp are scalar ([1 1]) chunkie kernel objects, drop-in for
%   @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'gs_s') / 'gphi_s' in
%   chunkerkernevalmat.
%
% See also KERNEL3D.RADCHEB, SURFWAVE.RADCHEB_GREEN

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'eps'), opts.eps = 1e-13; end
rmin = 1e-12;
if isfield(opts, 'rmin'), rmin = opts.rmin; end
if numel(rmax) == 2
    rmin = rmax(1);  rmax = rmax(2);
end

[Kgs, i1] = rc_kernel(@(s,t) surfwave.capillary.kern(rts, ejs, s, t, 'gs_s'),   rmin, rmax, opts, 'gs');
[Kgp, i2] = rc_kernel(@(s,t) surfwave.capillary.kern(rts, ejs, s, t, 'gphi_s'), rmin, rmax, opts, 'gphi');
info = struct('gs', i1, 'gphi', i2);

end


function [K, info] = rc_kernel(f, rmin, rmax, opts, name)

kr = kernel3d.radcheb(f, [rmin rmax], opts);
b  = kr.params.breaks;
c  = kr.params.coefs;
info = struct('err', kr.params.err, 'npanels', size(c,1), 'rmin', rmin, 'rmax', rmax);

K = kernel(@(s,t) rc_eval(s, t, b, c));
K.name   = 'capillary radcheb';
K.type   = name;
K.opdims = [1 1];
K.sing   = 'log';

end


function val = rc_eval(s, t, b, c)

dx = t.r(1,:).' - s.r(1,:);
dy = t.r(2,:).' - s.r(2,:);
r  = sqrt(dx.^2 + dy.^2);
val = radcheb_eval(r, b, c);

end
