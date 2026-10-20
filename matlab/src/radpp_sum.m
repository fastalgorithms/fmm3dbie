function g = radpp_sum(coefs, varargin)
%RADPP_SUM linear combination of Green's functions or fmms, for kernel3d.radpp.
%
%  g = radpp_sum(coefs, g1, g2, ...)
%
%  returns the function handle
%
%     [val, grad, hess, third] = g(src, targ)
%
%  with each output sum_k coefs(k)*(that output of gk). Only as many
%  outputs as are requested are computed, so value-only Green's
%  functions may be combined for kernels that need no derivatives.
%
%  Example: Helmholtz minus Laplace, fitted as one profile
%     gh = @(s,t) helm3d.green(1.1, s, t);
%     K = kernel3d.radpp(radpp_sum([1 -1], gh, @lap3d.green), 4, 's');
%
%  Analytically combining kernels is often more stable, see e.g.
%  HELM3D.HELMDIFF_GREEN.
%
%  See also KERNEL3D.RADPP

if numel(coefs) ~= numel(varargin)
    error('RADPP_SUM: need one coefficient per Green''s function.');
end
greens = varargin;
g = @(varargin) green_sum(coefs, greens, varargin{:});

end


function varargout = green_sum(coefs, greens, varargin)
nout = max(nargout, 1);
varargout = cell(1, nout);
for k = 1:numel(greens)
    out = cell(1, nout);
    [out{:}] = greens{k}(varargin{:});
    for i = 1:nout
        varargout{i} = add(varargout{i}, coefs(k), out{i});
    end
end
end


function a = add(a, c, b)
% a + c*b, field by field for the output struct of an fmm
if ~isstruct(b)
    if isempty(a), a = c*b; else, a = a + c*b; end
    return
end
for f = fieldnames(b).'
    v = c*b.(f{1});
    if isfield(a, f{1}), v = a.(f{1}) + v; end
    a.(f{1}) = v;
end
end
