function ri = rsc_interleave_full(m, n, nker, entries)
%RSC_INTERLEAVE_FULL  rsc_to_interleave descriptor for kernel3d.
%
%   ri = kernel3d.rsc_interleave_full(m, n)
%
%     Full m-by-n block: wnear has nker = m*n rows in column-major order,
%     each mapping 1-to-1 to a block entry with coefficient 1.
%
%   ri = kernel3d.rsc_interleave_full(m, n, nker, entries)
%
%     General form: wnear has nker rows; the block is assembled as
%       B(row_id, col_id) += coef * wnear(ker_id, :)
%     where entries is a struct array with fields ker_id, row_id, col_id, coef.
%     Used for kernels with complex coefficients (e.g. Maxwell nrccie-eval)
%     or symmetry (e.g. Stokes symmetric3).
%
%   All ri structs have fields: .m, .n, .nker, .entries

if nargin == 2
    % Full dense block, column-major, coefficient 1.
    nker = m * n;
    [rows, cols] = meshgrid(1:m, 1:n);   % column-major: row varies fastest
    row_ids = rows(:);
    col_ids = cols(:);
    entries = struct('ker_id', num2cell(1:nker), ...
                     'row_id', num2cell(row_ids.'), ...
                     'col_id', num2cell(col_ids.'), ...
                     'coef',   num2cell(ones(1,nker)));
end

ri.m       = m;
ri.n       = n;
ri.nker    = nker;
ri.entries = entries;

end
