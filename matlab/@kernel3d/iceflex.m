function obj = iceflex(type, rts, ejs, green, varargin)
%KERNEL3D.ICEFLEX   Flexural (ice plate) surface wave kernels, for a flat
%   surfer in z = 0.
%
%   K = KERNEL3D.ICEFLEX(type, rts, ejs)
%   K = KERNEL3D.ICEFLEX(type, rts, ejs, green)
%   K = KERNEL3D.ICEFLEX(type, rts, ejs, green, ...)
%
%   RTS and EJS are the dispersion roots and their partial fraction
%   residues, see SURFWAVE.FLEX.FIND_ROOTS_FLEX. GREEN is a weight vector
%   w on
%
%      [G_S, G_phi, S3d, Bilap G_phi]
%
%   with G_S = SURFWAVE.FLEX.GSFLEX, G_phi = SURFWAVE.FLEX.GPHIFLEX,
%   S3d = S3d G_phi + the Laplace single layer 1/(4 pi r) and the
%   bilaplacian of G_phi (the last two values only, so type 's'). Trailing
%   zeros may be left out. 'gs' (default), 'gphi', 's3d' and 'gphi_bilap'
%   are the unit vectors.
%
%   TYPE is any type of KERNEL3D.RADPP with 'dim' = 2 ('s', 'sp', 'grad',
%   'hess', 'clamped_bc', 'supported_bc', 'free_bc', ...), or 'lap', see
%   SURFWAVE.KERNEL. The plate types take the Poisson ratio as 'nu', nu.
%
%   By default this is a KERNEL3D.RADGREEN kernel (exact eval); with
%   'rmax', rmax it is a KERNEL3D.RADPP kernel (eval from a fit on
%   [1e-12 rmax], fast), see SURFWAVE.KERNEL.
%
%   'add', {gs_d, gphi_d, s3d_d, bilap_d} adds smooth corrections to the
%   Green's functions (e.g. finite depth), weighted by w, see
%   SURFWAVE.FLEX.GREEN.
%
%   Further arguments go to KERNEL3D.RADGREEN or KERNEL3D.RADPP.
%
%   Example, the polynya v2b operator:
%      Kv2b = kernel3d.iceflex('free_bc', rts, ejs, [-gamma/4, 1/2], ...
%                              'nu', nu, 'rmax', 50);
%
% See also SURFWAVE.FLEX.GREEN, SURFWAVE.KERNEL, KERNEL3D.CAPILLARY

if ( nargin < 3 )
    error('KERNEL3D.ICEFLEX: need type, rts and ejs.');
end
if ( nargin < 4 || isempty(green) ), green = 'gs'; end

rts = rts(:);  ejs = ejs(:);

add = {};
k = find(strcmpi(varargin(1:2:end), 'add'), 1);
if ( ~isempty(k) )
    add = varargin{2*k};
    varargin(2*k-1:2*k) = [];
end

obj = surfwave.kernel(surfwave.flex.green(rts, ejs, green, add), type, ...
                     max(abs(rts)), varargin{:});

end

