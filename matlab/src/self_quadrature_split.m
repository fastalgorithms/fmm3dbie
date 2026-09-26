function [xs, ys, ws] = self_quadrature_split(norder, ipv, ifar, verts, x0, y0, dr)
% SELF_QUADRATURE_SPLIT self quadrature as in SELF_QUADRATURE, with
%   the radial rule split so that it also resolves a singularity of the
%   density at the far end of each ray.
%
%   Syntax
%     [xs, ys, ws] = self_quadrature_split(norder, ipv, ifar, verts, x0, y0, dr)
%
%   Input arguments:
%     * norder: order of polynomials for representing the density
%     * ipv: nature of the kernel singularity, 0 or 1
%     * ifar: rule on the far part of each ray
%         - ifar = 0: the near rule reflected to be singular at the
%           far end
%     * verts: (2,nv) coordinates defining the convex polygon patch
%     * x0, y0: coordinates of the target/singularity on the patch
%     * dr: (3,2) surface Jacobian at (x0,y0)
%
%   Output arguments:
%     * xs, ys: quadrature node coordinates on the patch
%     * ws: quadrature weights
%
    nmax = 50000;
    xs = zeros(nmax,1);
    ys = zeros(nmax,1);
    ws = zeros(nmax,1);
    [~, nv] = size(verts);
    druse = reshape(dr, [3,2]);
    nquad = 0;
    ier = 0;

    mex_id_ = 'self_quadrature_split(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c io int64_t[x], c io double[x], c io double[x], c io double[x], c io int64_t[x])';
[nquad, xs, ys, ws, ier] = fmm3dbie_routs(mex_id_, norder, ipv, ifar, verts, nv, x0, y0, druse, nquad, xs, ys, ws, ier, 1, 1, 1, 2, nv, 1, 1, 1, 3, 2, 1, nmax, nmax, nmax, 1);

    if ier == 8
        error('FMM3DBIE:self_quadrature_split:unsupported_ipv', ...
            'self_quadrature_split supports ipv = 0 or 1 (got ipv=%d)', ipv);
    end

    xs = xs(1:nquad);
    ys = ys(1:nquad);
    ws = ws(1:nquad);
end
%
%
%

%
%
%
%
%

%-------------------------------------------------
%
%%
%%   Laplace dirichlet routines
%
%
%-------------------------------------------------

