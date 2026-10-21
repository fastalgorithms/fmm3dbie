function obj = kernel2d(type, nu, rts, ejs)
%SURFWAVE.FLEX.KERNEL2D   chunkie kernel objects for the 2D (boundary)
%   flexural / polynya kernels.
%
%   K = SURFWAVE.FLEX.KERNEL2D(type, nu, rts, ejs)
%
%   NU is the Poisson ratio; RTS and EJS are the dispersion roots and their
%   partial fraction residues from SURFWAVE.FLEX.FIND_ROOTS_FLEX.
%
%   K.sing and K.opdims are set, so chunkermat picks the singularity up
%   from the kernel; only the quadrature *rule* is a chunkermat option.
%
%   Boundary to boundary, [2 3], acting on the interleaved density
%   [rho1; H*rho1; rho2] (same layout as the b2v densmap):
%
%      'b2b_log'  - pieces built from gsflex (first part + hilbert),
%                   sing = 'log'.     chunkermat(chnkr, K)
%      'b2b_bh'   - pieces built from bhgreen (first part bh + hilbert bh),
%                   sing = 'smooth'.  chunkermat(chnkr, K, struct('quad','native'))
%                   Self entries are NaN; the caller replaces them with the
%                   diagonal limit. (chunkermat maps 'smooth' to the log GGQ
%                   rule under quad = 'ggq', which these evaluators do not
%                   tolerate, hence the explicit 'native'.)
%
%   The same operator split as chnk.flex2d.kern 'free_plate_log' /
%   'free_plate_smooth' (preferred: the 'b2b_log' pair above puts the full
%   biharmonic part of K21H into the log kernel, which loses digits at the
%   near-singular GGQ nodes):
%
%      'free_plate_log'    - everything but the "bad" K21H terms, at most
%                            log singular, sing = 'log'.   chunkermat(chnkr, K)
%      'free_plate_smooth' - the "bad" K21H terms, sing = 'smooth'.
%                            chunkermat(chnkr, K, struct('quad','native'))
%                            The self entries come from kappa' = d kappa/ds,
%                            which must be in chnkr.data(1,:), e.g.
%                              chnkr = makedatarows(chnkr, 1);
%                              chnkr.data(1,:) = arclengthder(chnkr, signed_curvature(chnkr));
%                            No diagonal fix by the caller.
%
%   Scalar, [1 1]:
%
%      'hilb'     - Hilbert transform kernel, sing = 'pv'
%
%   Boundary to volume, [1 3] (columns G_ny, (1+nu)/2 G_tauy, G):
%
%      'gs_eval'   - SURFWAVE.FLEX.KERNEL_EVAL('free_plate_gs_eval', ...)
%      'gphi_eval' - SURFWAVE.FLEX.KERNEL_EVAL('free_plate_gphi_eval', ...)
%      Combine with kernel arithmetic, e.g. Kb2v = 0.5*Kgs - Kgp.
%
% See also SURFWAVE.FLEX.KERN, SURFWAVE.FLEX.KERNEL_EVAL, KERNEL3D.ICEFLEX

if ( nargin < 4 )
    error('SURFWAVE.FLEX.KERNEL2D: need type, nu, rts and ejs.');
end

switch lower(type)
    case {'gs_eval', 'free_plate_gs_eval'}
        obj = surfwave.flex.kernel_eval('free_plate_gs_eval', nu, rts, ejs);
        return
    case {'gphi_eval', 'free_plate_gphi_eval'}
        obj = surfwave.flex.kernel_eval('free_plate_gphi_eval', nu, rts, ejs);
        return
end

obj = kernel();
obj.name   = 'iceflex';
obj.params = struct('nu', nu, 'rts', rts, 'ejs', ejs);

switch lower(type)
    case {'b2b_log', 'free plate b2b log'}
        obj.type   = 'b2b_log';
        obj.opdims = [2 3];
        obj.sing   = 'log';
        obj.eval   = @(s,t) surfwave.flex.kern(s, t, 'free plate b2b log', nu, rts, ejs);
    case {'b2b_bh', 'free plate b2b bh'}
        obj.type   = 'b2b_bh';
        obj.opdims = [2 3];
        obj.sing   = 'smooth';
        obj.eval   = @(s,t) surfwave.flex.kern(s, t, 'free plate b2b bh', nu, rts, ejs);
    case {'free_plate_log', 'free plate log'}
        obj.type   = 'free_plate_log';
        obj.opdims = [2 3];
        obj.sing   = 'log';
        obj.eval   = @(s,t) surfwave.flex.kern(s, t, 'free plate log', nu, rts, ejs);
    case {'free_plate_smooth', 'free plate smooth'}
        obj.type   = 'free_plate_smooth';
        obj.opdims = [2 3];
        obj.sing   = 'smooth';
        obj.eval   = @(s,t) surfwave.flex.kern(s, t, 'free plate smooth', nu, rts, ejs);
    case {'hilb', 'hilbert'}
        obj.type   = 'hilb';
        obj.opdims = [1 1];
        obj.sing   = 'pv';
        obj.eval   = @(s,t) surfwave.flex.kern(s, t, 'hilb', nu, rts, ejs);
    otherwise
        error('SURFWAVE.FLEX.KERNEL2D: unknown type ''%s''.', type);
end

end
