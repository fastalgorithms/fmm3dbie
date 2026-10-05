function f = mrdivide(f, g)
% / Matrix right division for kernel3d class
%
% Currently only supports scalars: returns F/g for kernel F and scalar g.

if (isa(f, 'kernel3d') && isnumeric(g) && isscalar(g))

    % same as scalar times, so eval, fmm, getquad, diag and the
    % iszero/isnan flags are all handled in one place
    f = times(f, 1/g);

else
    error('KERNEL3D:mrdivide:invalid', ...
        'F must be a kernel3d class object and G a scalar');
end
end
