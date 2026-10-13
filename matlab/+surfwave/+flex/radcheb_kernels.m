function [Kgs, Kgp] = radcheb_kernels(nu, rts, ejs, rmax, opts)
%SURFWAVE.FLEX.RADCHEB_KERNELS  chunkie kernel objects for the polynya
%   (flexural) boundary-to-volume free-plate evaluation, G_S and G_phi,
%   evaluated through radial Chebyshev interpolants of G and its gradient
%   (surfwave.radcheb_green).
%
%   [Kgs, Kgp] = surfwave.flex.radcheb_kernels(nu, rts, ejs, rmax)
%   [Kgs, Kgp] = surfwave.flex.radcheb_kernels(nu, rts, ejs, rmax, opts)
%
%   Kgs, Kgp are [1 3] chunkie kernels (columns G_ny, (1+nu)/2 G_tauy, G),
%   drop-in for surfwave.flex.kernel_eval('free_plate_gs_eval', ...) and
%   'free_plate_gphi_eval' in chunkerkernevalmat. Same calling pattern as
%   surfwave.capillary.radcheb_kernels.
%
%   RMAX is the largest source-target distance, or [rmin, rmax]; radii
%   outside the range give NaN (as radcheb_eval). OPTS goes to
%   kernel3d.radcheb (opts.eps, ...; default eps 1e-13).
%
% See also SURFWAVE.RADCHEB_GREEN, SURFWAVE.FLEX.KERNEL_EVAL,
%   SURFWAVE.CAPILLARY.RADCHEB_KERNELS

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'eps'), opts.eps = 1e-13; end

Kgs = surfwave.flex.kernel_eval('free_plate_gs_eval',   nu, rts, ejs, ...
    struct('rmax', rmax, 'radcheb', opts));
Kgp = surfwave.flex.kernel_eval('free_plate_gphi_eval', nu, rts, ejs, ...
    struct('rmax', rmax, 'radcheb', opts));

end
