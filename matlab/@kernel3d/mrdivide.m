function f = mrdivide(f, g)
% / Matrix right division for kernel3d class
%
% Currently only supports scalars: returns F/g for kernel F and scalar g.

if (isnumeric(g) && isscalar(g))

    if (isa(f.eval, 'function_handle'))
        f.eval = @(varargin) f.eval(varargin{:}) / g;
    end

    if (isa(f.fmm, 'function_handle'))
        f.fmm = @(varargin) f.fmm(varargin{:}) / g;
    else
        f.fmm = [];
    end

    if (isa(f.getquad, 'function_handle'))
        fgetquad = f.getquad;
        f.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
            fgetquad(S, eps, varargin{:}), 1/g);
    else
        f.getquad = [];
    end

else
    error('KERNEL3D:mrdivide:invalid', ...
        'F must be a kernel3d class object and G a scalar');
end
end
