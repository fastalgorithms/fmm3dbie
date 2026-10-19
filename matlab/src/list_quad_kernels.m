function names = list_quad_kernels()
%LIST_QUAD_KERNELS  Names of the kernels available to GETNEARQUAD_MEX.
%
%  Syntax:
%    names = list_quad_kernels()
%
%  Output arguments:
%    * names: (1,n) cell array of kernel names
%
%  See also GETNEARQUAD_MEX

    n = 0;
    mex_id_ = 'kernel_registry_count(c io int64_t[x])';
[n] = fmm3dbie_routs(mex_id_, n, 1);
    if n < 1, names = {}; return; end

    nchar = 64;
    codes = zeros(nchar, n);
    lens  = zeros(n, 1);
    mex_id_ = 'kernel_registry_names(c i int64_t[x], c i int64_t[x], c io int64_t[xx], c io int64_t[x])';
[codes, lens] = fmm3dbie_routs(mex_id_, n, nchar, codes, lens, 1, 1, nchar, n, n);

    names = cell(1, n);
    for i = 1:n
        names{i} = char(codes(1:lens(i), i)).';
    end

end
