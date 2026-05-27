function f = plus(f, g)
% + Pointwise addition for kernel3d class
%
% Currently only supported for adding two kernel3d class objects.

if (isa(g, 'kernel3d') && isa(f, 'kernel3d'))
    assert(all(f.opdims == g.opdims), ...
        'kernel3d dimensions must agree to add');

    f.name = ['custom ', f.name, ' ', g.name];
    f.type = 'custom';

    f.eval = @(varargin) f.eval(varargin{:}) + g.eval(varargin{:});

    if (isa(f.fmm, 'function_handle') && isa(g.fmm, 'function_handle'))
        f.fmm = @(varargin) f.fmm(varargin{:}) + g.fmm(varargin{:});
    else
        f.fmm = [];
    end

    if (isa(f.layer_eval, 'function_handle') && isa(g.layer_eval, 'function_handle'))
        f.layer_eval = @(varargin) f.layer_eval(varargin{:}) + g.layer_eval(varargin{:});
    else
        f.layer_eval = [];
    end

    if (isa(f.getquad, 'function_handle') && isa(g.getquad, 'function_handle'))
        fgetquad = f.getquad;
        ggetquad = g.getquad;
        fri = f.rsc_to_interleave;
        f.getquad = @(S, eps, varargin) kernel3d.addquad( ...
            fgetquad(S, eps, varargin{:}), ggetquad(S, eps, varargin{:}), S, 1, fri);
    else
        f.getquad = [];
    end

    f.get_overs_orders = @(S,t,eps) max(f.get_overs_orders(S, t, eps),g.get_overs_orders(S, t, eps));

    src  = union(f.src_fields,  g.src_fields);
    targ = union(f.targ_fields, g.targ_fields);
    if isempty(src),  src  = []; end
    if isempty(targ), targ = []; end
    f.src_fields  = src;
    f.targ_fields = targ;

else
    error('KERNEL3D:plus:invalid', ...
        'F and G must be kernel3d class objects');
end
end
