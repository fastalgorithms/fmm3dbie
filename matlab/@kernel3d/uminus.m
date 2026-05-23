function f = uminus(f)
% - Unary negation for kernel3d class

f.eval = @(varargin) -f.eval(varargin{:});

if (isa(f.fmm, 'function_handle'))
    f.fmm = @(varargin) -f.fmm(varargin{:});
else
    f.fmm = [];
end

if (isa(f.layer_eval, 'function_handle'))
    f.layer_eval = @(varargin) -f.layer_eval(varargin{:});
else
    f.layer_eval = [];
end

if (isa(f.getquad, 'function_handle'))
    fgetquad = f.getquad;
    f.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
        fgetquad(S, eps, varargin{:}), S, -1);
else
    f.getquad = [];
end

end
