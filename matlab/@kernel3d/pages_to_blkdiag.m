function F = pages_to_blkdiag(F3, p, q, n)
%KERNEL3D.PAGES_TO_BLKDIAG  Build a sparse block-diagonal from a 3D array.
%
%   F = kernel3d.pages_to_blkdiag(F3, p, q, n)
%
%   F3 is (p x q x n).  Returns sparse (p*n x q*n) with
%   F(p*(i-1)+1:p*i, q*(i-1)+1:q*i) = F3(:,:,i) for i = 1..n.

nz  = p * q * n;
ii  = zeros(nz, 1);
jj  = zeros(nz, 1);
vv  = zeros(nz, 1);
ptr = 1;
for i = 1:n
    row0 = (i-1)*p;
    col0 = (i-1)*q;
    for ci = 1:q
        for ri = 1:p
            ii(ptr) = row0 + ri;
            jj(ptr) = col0 + ci;
            vv(ptr) = F3(ri, ci, i);
            ptr = ptr + 1;
        end
    end
end
F = sparse(ii, jj, vv, p*n, q*n);
end
