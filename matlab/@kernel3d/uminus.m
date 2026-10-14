function f = uminus(f)

if kernel3d.radgreen('compatible', f)
    f = kernel3d.radgreen('scale', f, -1);
    return
end
% - Unary negation for kernel3d class

f.eval = @(varargin) -f.eval(varargin{:});

if (isa(f.fmm, 'function_handle'))
    f.fmm = @(varargin) -f.fmm(varargin{:});
else
    f.fmm = [];
end

if (isa(f.getquad, 'function_handle'))
    fgetquad = f.getquad;
    f.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
        fgetquad(S, eps, varargin{:}), -1);
else
    f.getquad = [];
end

if isa(f.diag, 'function_handle')
    fdiag = f.diag;
    f.diag = @(t) -fdiag(t);
end

end
