function M = tangent_block(src, vds, scalar_rows, opdim, use_pinv)
%TANGENT_BLOCK   Block-diagonal tangent basis matrix for a point set.
%
%   M = TANGENT_BLOCK(SRC, VDS, SCALAR_ROWS, OPDIM, USE_PINV)
%
%   use_pinv = true  : contraction, M is (opdim_tan*n) x (opdim_cart*n)
%                      C = (B'B)\B'  (2x3, pseudoinverse)
%   use_pinv = false : expansion,   M is (opdim_cart*n) x (opdim_tan*n)
%                      C = B        (3x2)
%
%   In both cases OPDIM is the tangent-space opdim per point
%   (numel(scalar_rows) + 2*k).  The Cartesian opdim is numel(scalar_rows) + 3*k.
%
%   See also KERNEL3D.TANGENT_KERN

n         = size(src.r, 2);
k         = size(vds, 2);
n_sc      = numel(scalar_rows);
opdim_tan  = opdim;
opdim_cart = n_sc + 3*k;

n_nz = n * (n_sc + k * 2 * 3);
ii = zeros(n_nz, 1);
jj = zeros(n_nz, 1);
vv = zeros(n_nz, 1);
ptr = 0;

for p = 1:n
    B = [src.du(:,p), src.dv(:,p)];  % 3x2
    if use_pinv
        C        = (B'*B) \ B';               % 2x3
        row_off  = (p-1) * opdim_tan;
        col_off  = (p-1) * opdim_cart;
        row_sc   = @(ri) row_off + ri;
        col_sc   = @(ri) col_off + scalar_rows(ri);
        vec_row0 = @(bi) row_off + n_sc + (bi-1)*2;
        vec_cols = @(bi) col_off + vds(:,bi);  % 3 Cartesian cols
    else
        C        = B;                          % 3x2
        row_off  = (p-1) * opdim_cart;
        col_off  = (p-1) * opdim_tan;
        row_sc   = @(ri) row_off + scalar_rows(ri);
        col_sc   = @(ri) col_off + ri;
        vec_row0 = @(bi) row_off + vds(1,bi) - 1;  % placeholder, overridden below
        vec_cols = @(bi) col_off + n_sc + (bi-1)*2 + 1; % first of 2 tangent cols (1-based)
    end

    for ri = 1:n_sc
        ptr = ptr + 1;
        ii(ptr) = row_sc(ri);
        jj(ptr) = col_sc(ri);
        vv(ptr) = 1;
    end

    for bi = 1:k
        if use_pinv
            in_cols  = vec_cols(bi);   % 3 Cartesian col indices
            out_row0 = vec_row0(bi);
            for ci = 1:3
                for ri = 1:2
                    ptr = ptr + 1;
                    ii(ptr) = out_row0 + ri;
                    jj(ptr) = in_cols(ci);
                    vv(ptr) = C(ri, ci);
                end
            end
        else
            out_rows = vds(:,bi);      % 3 Cartesian row indices
            in_col0  = vec_cols(bi);   % first of 2 tangent cols
            for ci = 1:2
                for ri = 1:3
                    ptr = ptr + 1;
                    ii(ptr) = row_off + out_rows(ri);
                    jj(ptr) = in_col0 + ci - 1;  % in_col0 is 1-based first slot; ci=1..2
                    vv(ptr) = C(ri, ci);
                end
            end
        end
    end
end

if use_pinv
    M = sparse(ii(1:ptr), jj(1:ptr), vv(1:ptr), opdim_tan*n,  opdim_cart*n);
else
    M = sparse(ii(1:ptr), jj(1:ptr), vv(1:ptr), opdim_cart*n, opdim_tan*n);
end

end
