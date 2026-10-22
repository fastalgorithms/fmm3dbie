function obj = capillary(type, rts, ejs, green, varargin)
%KERNEL3D.CAPILLARY   Capillary surface wave kernels, for a flat surfer
%   in z = 0.
%
%   K = KERNEL3D.CAPILLARY(type, rts, ejs)
%   K = KERNEL3D.CAPILLARY(type, rts, ejs, green)
%   K = KERNEL3D.CAPILLARY(type, rts, ejs, green, 'rmax', rmax, ...)
%
%   RTS and EJS are the dispersion roots and their partial fraction
%   residues, see SURFWAVE.CAPILLARY.FIND_ROOTS_CAPILLARY. GREEN is a
%   weight vector w on
%
%      [G_S, G_phi, S3d]
%
%   with G_S = SURFWAVE.CAPILLARY.GSHELM, G_phi = SURFWAVE.CAPILLARY.GPHIHELM
%   and S3d = S3d G_phi + the Laplace single layer 1/(4 pi r) (values only,
%   so type 's'), giving the kernel on w(1) G_S + w(2) G_phi + w(3) S3d.
%   Trailing zeros may be left out. 'gs', 'gphi' and 's3d' (default 'gs')
%   are [1 0 0], [0 1 0] and [0 0 1].
%
%   TYPE is any type of KERNEL3D.RADPP with 'dim' = 2 that needs at most
%   two derivatives ('s', 'sp', 'grad', 'hess', ...), or 'lap' (surface
%   Laplacian G_xx + G_yy).
%
%   By default this is a KERNEL3D.RADGREEN kernel (exact eval, near
%   quadrature from a radpp fit); with 'rmax', rmax it is a KERNEL3D.RADPP
%   kernel (eval from a fit on [1e-12 rmax], fast), see SURFWAVE.KERNEL.
%   Further arguments go to KERNEL3D.RADGREEN or KERNEL3D.RADPP, e.g. fit
%   options 'opts', struct('eps', 1e-12).
%
%   Example, the capillary v2b operator cs*S'_S + cp*S'_phi:
%      Kv2b = kernel3d.capillary('sp', rts, ejs, [cs cp], 'rmax', 50);
%
% See also SURFWAVE.CAPILLARY.GREEN, SURFWAVE.KERNEL, KERNEL3D.ICEFLEX, KERNEL3D.GRAVITY

if ( nargin < 3 )
    error('KERNEL3D.CAPILLARY: need type, rts and ejs.');
end
if ( nargin < 4 || isempty(green) ), green = 'gs'; end

rts = rts(:);  ejs = ejs(:);
obj = surfwave.kernel(surfwave.capillary.green(rts, ejs, green), type, ...
                     max(abs(rts)), varargin{:});

end
