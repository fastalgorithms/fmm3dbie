function f = minus(f, g)
% - Pointwise subtraction for kernel3d class
%
% Currently only supported for subtracting two kernel3d class objects.

if (isa(g, 'kernel3d') && isa(f, 'kernel3d'))
    if kernel3d.radgreen('summable', {f, g})
        f = kernel3d.radgreen('sum', f, 1, g, -1);
        return
    end
    if kernel3d.radpp('summable', {f, g})
        f = kernel3d.radpp('sum', {f, g}, [1 -1]);
        return
    end
    assert(all(f.opdims == g.opdims), ...
        'kernel3d dimensions must agree to subtract');

    fdiag = f.diag; gdiag = g.diag;

    f.name = ['custom ', f.name, ' ', g.name];
    f.type = 'custom';

    f.eval = @(varargin) f.eval(varargin{:}) - g.eval(varargin{:});

    if (isa(f.fmm, 'function_handle') && isa(g.fmm, 'function_handle'))
        f.fmm = @(varargin) f.fmm(varargin{:}) - g.fmm(varargin{:});
    else
        f.fmm = [];
    end

    if (isa(f.getquad, 'function_handle') && isa(g.getquad, 'function_handle'))
        fgetquad = f.getquad;
        ggetquad = g.getquad;
        f.getquad = @(S, eps, varargin) kernel3d.addquad( ...
            fgetquad(S, eps, varargin{:}), ggetquad(S, eps, varargin{:}), -1);
    else
        f.getquad = [];
    end
    f.getquad_sub = f.getquad_sub && g.getquad_sub;
    f.kernel_order = max(f.kernel_order, g.kernel_order);

    f.get_overs_orders = @(S,t,eps) max(f.get_overs_orders(S, t, eps),g.get_overs_orders(S, t, eps));

    src  = union(f.src_fields,  g.src_fields);
    targ = union(f.targ_fields, g.targ_fields);
    if isempty(src),  src  = []; end
    if isempty(targ), targ = []; end
    f.src_fields  = src;
    f.targ_fields = targ;

    % identity/self ("delta") terms subtract
    if isempty(fdiag) && ~isempty(gdiag)
        f.diag = @(t) -gdiag(t);
    elseif ~isempty(fdiag) && ~isempty(gdiag)
        f.diag = @(t) fdiag(t) - gdiag(t);
    end

else
    error('KERNEL3D:minus:invalid', ...
        'F and G must be kernel3d class objects');
end
end
