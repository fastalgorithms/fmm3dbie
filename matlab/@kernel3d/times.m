function f = times(f, g)
% .* Pointwise multiplication for kernel3d class
%
% Currently only supported for scalars: if g is a numeric scalar then
% f.*g or g.*f returns a kernel3d object where eval, fmm, layer_eval,
% and getquad are each scaled by g.

if (~isa(f, 'kernel3d'))
    f = times(g, f);
    return
elseif (isnumeric(g) && isscalar(g))

    if (isa(f.eval, 'function_handle'))
        f.eval = @(varargin) g * f.eval(varargin{:});
    end

    if (isa(f.fmm, 'function_handle'))
        f.fmm = @(varargin) g * f.fmm(varargin{:});
    else
        f.fmm = [];
    end

    if (isa(f.layer_eval, 'function_handle'))
        f.layer_eval = @(varargin) g * f.layer_eval(varargin{:});
    else
        f.layer_eval = [];
    end

    if (isa(f.getquad, 'function_handle'))
        fgetquad = f.getquad;
        f.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
            fgetquad(S, eps, varargin{:}), S, g);
    else
        f.getquad = [];
    end

else
    error('KERNEL3D:times:invalid', ...
        'F or G must be a scalar and the other a kernel3d class object');
end
end
