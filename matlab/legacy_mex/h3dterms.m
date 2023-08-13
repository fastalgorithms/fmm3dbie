function [nterms] = h3dterms(boxsize,zk,eps)
    nterms = 0;
    mex_id_ = 'h3dterms(i double[x], i dcomplex[x], i double[x], io int[x])';
[nterms] = fmm3dbierouts(mex_id_, boxsize, zk, eps, nterms, 1, 1, 1, 1);
end



