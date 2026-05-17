function [xs, ys, ws] = self_quadrature(norder, ipv, verts, x0, y0, dr)
    nmax = 50000;
    xs = zeros(nmax,1);
    ys = zeros(nmax,1);
    ws = zeros(nmax,1);
    [~, nv] = size(verts);
    druse = reshape(dr, [3,2]);
    nquad = 0;
    ier = 0;

    mex_id_ = 'self_quadrature(i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[x], i double[x], i double[xx], io int64_t[x], io double[x], io double[x], io double[x], io int64_t[x])';
[nquad, xs, ys, ws, ier] = kern_routs(mex_id_, norder, ipv, verts, nv, x0, y0, druse, nquad, xs, ys, ws, ier, 1, 1, 2, nv, 1, 1, 1, 3, 2, 1, nmax, nmax, nmax, 1);

    xs = xs(1:nquad);
    ys = ys(1:nquad);
    ws = ws(1:nquad);
end
%
%
%
%
