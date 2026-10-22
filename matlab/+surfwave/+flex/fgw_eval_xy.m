function [D1, D2] = fgw_eval_xy(S, x, y, kind)
%FGW_EVAL_XY  Cartesian partial derivatives of the point-load kernels up to fourth order.
%
%   [D1, D2] = FGW_EVAL_XY(S, x, y, kind)
%
%   Input
%       S     : structure returned by FGW_SETUP
%       x, y  : Cartesian coordinates of the evaluation points (arrays of equal size, any shape),
%               in the same length units as the parameters given to FGW_SETUP
%       kind  : 'full' (default) - the finite-depth kernels G_j(x,y)
%               'deep'           - the deep-water kernels
%               'diff'           - the difference (finite depth minus deep water)
%
%   Output
%       D1, D2 : complex arrays of size numel(x) x 5 x 5 containing
%
%                    D1(p, a+1, b+1) = d^(a+b) G1 / dx^a dy^b   evaluated at (x(p), y(p)),
%
%                for every a, b >= 0 with a + b <= 4; D2 is the same for G2.  G1 is the deflection
%                kernel (transform of kT/((a k^4+gam)kT-1)), G2 the potential kernel (transform of
%                1/((a k^4+gam)kT-1)); see FGW_EVAL for the conventions.  For example
%
%                    D1(:,1,1)  G            D1(:,3,1)  d2G/dx2         D1(:,5,1)  d4G/dx4
%                    D1(:,2,1)  dG/dx        D1(:,2,2)  d2G/dxdy        D1(:,4,2)  d4G/dx3dy
%                    D1(:,1,2)  dG/dy        D1(:,1,3)  d2G/dy2         D1(:,3,3)  d4G/dx2dy2
%                                            D1(:,4,1)  d3G/dx3         D1(:,2,4)  d4G/dxdy3
%                                            D1(:,3,2)  d3G/dx2dy       D1(:,1,5)  d4G/dy4
%                                            D1(:,2,3)  d3G/dxdy2
%                                            D1(:,1,4)  d3G/dy3
%
%                The 10 entries with a + b > 4 (e.g. D1(:,5,5)) are not computed and are set to NaN.
%
%   At the origin the 'deep' and 'full' kernels are not smooth (G1 ~ r^2 log r, G2 contains an r^3
%   term), so the derivatives that do not have a finite limit there are returned as NaN for points
%   with x = y = 0 exactly.  The 'diff' kernels are smooth everywhere.
%
%   Accuracy is that of the radial evaluator (about 1e-14 relative to the natural scale of each
%   derivative), also close to the origin: the derivatives are assembled from the quantities
%   ((1/r) d/dr)^m G, m = 0..4, which FGW_EVAL(..., 'theta') evaluates directly, so no division by
%   r occurs.  For a radial function G(r) the exact relation is
%       d^(a+b) G / dx^a dy^b = sum over i <= a/2, j <= b/2 of
%           a!/(i!(a-2i)! 2^i) * b!/(j!(b-2j)! 2^j) * x^(a-2i) y^(b-2j) * ((1/r) d/dr)^(a+b-i-j) G.
%
%   Example
%       S = fgw_setup(1e5, 10, 50);                       % a [m^5], gam [m], d [m]
%       [X, Y] = meshgrid(-200:5:200, -150:5:150);
%       [D1, D2] = fgw_eval_xy(S, X, Y);                  % 'full'
%       Gxx = reshape(D1(:,3,1), size(X));                % d2G1/dx2 on the grid

if nargin < 4 || isempty(kind), kind = 'full'; end
x = x(:); y = y(:);
if numel(x) ~= numel(y), error('fgw_eval_xy:args', 'x and y must have the same number of elements'); end
r = hypot(x, y);
[T1, T2] = surfwave.flex.fgw_eval(S, r, kind, 'theta');           % numel(x) x 5: ((1/r) d/dr)^m G_j, m = 0..4
M = numel(x);
D1 = nan(M, 5, 5); D2 = nan(M, 5, 5);
c = [1 0 0; 1 0 0; 1 1 0; 1 3 0; 1 6 3];             % c(a+1, i+1) = a!/(i!(a-2i)! 2^i)
xp = x.^(0:4); yp = y.^(0:4);                         % x^0..x^4, y^0..y^4
for a = 0:4
    for b = 0:4-a
        v1 = zeros(M, 1); v2 = zeros(M, 1);
        for i = 0:floor(a/2)
            for j = 0:floor(b/2)
                pf = c(a+1, i+1)*c(b+1, j+1)*xp(:, a-2*i+1).*yp(:, b-2*j+1);
                m = a + b - i - j;
                v1 = v1 + pf.*T1(:, m+1);
                v2 = v2 + pf.*T2(:, m+1);
            end
        end
        D1(:, a+1, b+1) = v1; D2(:, a+1, b+1) = v2;
    end
end
end
