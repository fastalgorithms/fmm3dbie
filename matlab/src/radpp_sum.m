function g = radpp_sum(coefs, varargin)
%RADPP_SUM linear combination of Green's functions in the chunkie calling
%  sequence, for kernel3d.radpp.
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
%  Example:
%     gs   = @(s,t) surfwave.flex.gsflex(rts, ejs, s, t);
%     gphi = @(s,t) surfwave.flex.gphiflex(rts, ejs, s, t);
%     K = kernel3d.radpp(radpp_sum([-gamma/4 1/2], gs, gphi), rmax, ...
%                        'free_bc', 'nu', nu);
%
%  See also KERNEL3D.RADPP

if numel(coefs) ~= numel(varargin)
    error('RADPP_SUM: need one coefficient per Green''s function.');
end
greens = varargin;
g = @(s,t) green_sum(coefs, greens, s, t);

end


function varargout = green_sum(coefs, greens, s, t)
nout = max(nargout, 1);
varargout = cell(1, nout);
for k = 1:numel(greens)
    out = cell(1, nout);
    [out{:}] = greens{k}(s, t);
    for i = 1:nout
        if k == 1
            varargout{i} = coefs(k)*out{i};
        else
            varargout{i} = varargout{i} + coefs(k)*out{i};
        end
    end
end
end
