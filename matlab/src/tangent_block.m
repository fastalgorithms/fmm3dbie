function F = tangent_block(src, vds, scalar_rows, opdim, use_pinv)
%TANGENT_BLOCK   Per-point tangent basis pages for use with kernel3d/times.
%
%   F = TANGENT_BLOCK(SRC, VDS, SCALAR_ROWS, OPDIM, USE_PINV)
%
%   Returns F of size (rows x cols x n) where each page F(:,:,p) is the
%   block for point p.
%
%   use_pinv = true  : contraction, F is (opdim_tan x opdim_cart x n)
%                      C = (B'B)\B'  (2x3)
%   use_pinv = false : expansion,   F is (opdim_cart x opdim_tan x n)
%                      C = B        (3x2)
%
%   See also KERNEL3D.TANGENT_KERN

n         = size(src.r, 2);
k         = size(vds, 2);
n_sc      = numel(scalar_rows);
opdim_tan  = opdim;
opdim_cart = n_sc + 3*k;

if use_pinv
    F = zeros(opdim_tan, opdim_cart, n);
else
    F = zeros(opdim_cart, opdim_tan, n);
end

for p = 1:n
    B = [src.du(:,p), src.dv(:,p)];  % 3x2
    if use_pinv
        C = (B'*B) \ B';              % 2x3
    else
        C = B;                        % 3x2
    end

    % Scalar pass-through: identity block at the scalar rows/cols
    for ri = 1:n_sc
        if use_pinv
            F(ri, scalar_rows(ri), p) = 1;
        else
            F(scalar_rows(ri), ri, p) = 1;
        end
    end

    % Vector blocks
    for bi = 1:k
        if use_pinv
            % rows n_sc+(bi-1)*2+(1:2), cols vds(:,bi)
            row0 = n_sc + (bi-1)*2;
            F(row0+(1:2), vds(:,bi), p) = C;
        else
            % rows vds(:,bi), cols n_sc+(bi-1)*2+(1:2)
            col0 = n_sc + (bi-1)*2;
            F(vds(:,bi), col0+(1:2), p) = C;
        end
    end
end

end
