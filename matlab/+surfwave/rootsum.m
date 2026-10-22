function varargout = rootsum(rts, cjs, src, targ, realtol)
%SURFWAVE.ROOTSUM   sum over the dispersion roots of the
%   capillary radial profiles,
%
%      G(x,y) = sum_j cjs(j) K_j(|x-y|),
%
%   with K_j = -S_0(rho_j r) + 2i H_0^(1)(rho_j r) for a real root rho_j
%   and K_j = S_0(-rho_j r) otherwise, where S_0 = H_0 - Y_0 (Struve minus
%   Bessel). The capillary and flexural G_S, G_phi and S3d G_phi are this
%   sum with cjs = ejs.*rts.^2/2, ejs.*rts/4 and ejs/8 (see
%   SURFWAVE.CAPILLARY.GREEN and SURFWAVE.FLEX.GREEN). A root is real if
%   |angle| < realtol (default 1e-14).
%
%   Each K_j is evaluated with its singular parts removed
%   (SURFWAVE.FLEX.STRUVEKDIFFGREEN, SURFWAVE.FLEX.HELMDIFFGREEN), and the
%   removed parts are added back for the whole sum in closed form:
%
%      sum_j cjs(j) K_j = sum_j cjs(j) K_j^diff - 2/pi sum(cjs) log r
%                         - 2/pi sum(cjs.*rts) r     (derivatives only)
%                         + 2/pi sum(cjs.*rts.^3) T  (third and fourth)
%
%   with T_xxx = (x/r)^3 etc. (and its x and y derivatives),
%
%   so there is no cancellation of 1/r^2 terms between the Bessel and
%   Struve parts as r -> 0.
%
%   [val, grad, hess, third, fourth] = SURFWAVE.ROOTSUM(...): val (nt,ns),
%   grad(:,:,1:2) = G_x, G_y, hess(:,:,1:3) = G_xx, G_xy, G_yy,
%   third(:,:,1:4) = G_xxx, G_xxy, G_xyy, G_yyy and fourth(:,:,1:5)
%   likewise (target derivatives). src (2,ns), targ (2,nt), or structs
%   with field r.
%
% See also SURFWAVE.CAPILLARY.GREEN, SURFWAVE.FLEX.GREEN

if ( nargin < 5 || isempty(realtol) ), realtol = 1e-14; end
if ( isstruct(src) ),  src  = src.r;  end
if ( isstruct(targ) ), targ = targ.r; end
src = src(1:2,:);  targ = targ(1:2,:);

nout = max(nargout, 1);
out = cell(1, nout);
for k = 1:nout, out{k} = 0; end

for j = 1:numel(rts)
    rhoj = rts(j);
    P = cell(1, nout);
    if ( abs(angle(rhoj)) < realtol )
        S = cell(1, nout);  H = cell(1, nout);
        [S{:}] = surfwave.flex.struveKdiffgreen(rhoj, src, targ);
        [H{:}] = surfwave.flex.helmdiffgreen(rhoj, src, targ);
        for k = 1:nout, P{k} = -S{k} + 8*H{k}; end     % 2i*(-4i H)
    else
        [P{:}] = surfwave.flex.struveKdiffgreen(-rhoj, src, targ);
    end
    for k = 1:nout, out{k} = out{k} + cjs(j)*P{k}; end
end

% the removed parts, for the whole sum. The moments sum(ejs.*rts.^k)
% vanish for k = 0, 1, so for G_phi and S3d G_phi c0 is zero up to
% rounding, which the 1/r^2 terms would amplify: snap it to zero.
c0 = -2/pi*snap(cjs(:));
c1 = -2/pi*snap(cjs(:).*rts(:));
c3 =  2/pi*snap(cjs(:).*rts(:).^3);
dx = targ(1,:).' - src(1,:);
dy = targ(2,:).' - src(2,:);
r  = sqrt(dx.^2 + dy.^2);
[L, R] = logr_r(dx, dy, r, nout);
varargout = cell(1, nout);
varargout{1} = out{1} + c0*L{1};
for k = 2:nout
    varargout{k} = out{k} + c0*L{k} + c1*R{k};
end
% the third and fourth derivatives of the K_j^diff also leave out a term
% in rho_j^3 (zero for G_S, S3d G_phi and the Laplacian of G_phi, but not
% for G_phi itself)
if ( nout > 3 && c3 ~= 0 )
    ux = dx./r;  uy = dy./r;
    varargout{4} = varargout{4} + c3*cat(3, ux.^3, ux.^2.*uy, ux.*uy.^2, uy.^3);
    if ( nout > 4 )
        T1 = cat(3, ux.^4, ux.^3.*uy, ux.^2.*uy.^2, ux.*uy.^3, uy.^4);
        T2 = cat(3, 6*ux.^2, 3*ux.*uy, ux.^2 + uy.^2, 3*ux.*uy, 6*uy.^2);
        varargout{5} = varargout{5} + c3*(T2 - 7*T1)./r;
    end
end
end


function [L, R] = logr_r(dx, dy, r, nout)
% log r and r with their target derivatives, packed as SURFWAVE.HELMGREEN
x = dx;  y = dy;
L = {log(r), cat(3, x./r.^2, y./r.^2), ...
     cat(3, (r.^2 - 2*x.^2)./r.^4, -2*x.*y./r.^4, (r.^2 - 2*y.^2)./r.^4), ...
     cat(3, (8*x.^3 - 6*x.*r.^2)./r.^6, 2*y.*(4*x.^2 - r.^2)./r.^6, ...
            2*x.*(4*y.^2 - r.^2)./r.^6, (8*y.^3 - 6*y.*r.^2)./r.^6), ...
     cat(3, 6*(-8*x.^4 + 8*x.^2.*r.^2 - r.^4)./r.^8, ...
            24*x.*y.*(r.^2 - 2*x.^2)./r.^8, ...
            2*(-24*x.^2.*y.^2 + 4*(x.^2 + y.^2).*r.^2 - r.^4)./r.^8, ...
            24*x.*y.*(r.^2 - 2*y.^2)./r.^8, ...
            6*(-8*y.^4 + 8*y.^2.*r.^2 - r.^4)./r.^8)};
R = {r, cat(3, x./r, y./r), ...
     cat(3, (r.^2 - x.^2)./r.^3, -x.*y./r.^3, (r.^2 - y.^2)./r.^3), ...
     cat(3, 3*x.*(x.^2 - r.^2)./r.^5, y.*(3*x.^2 - r.^2)./r.^5, ...
            x.*(3*y.^2 - r.^2)./r.^5, 3*y.*(y.^2 - r.^2)./r.^5), ...
     cat(3, 3*(-5*x.^4 + 6*x.^2.*r.^2 - r.^4)./r.^7, ...
            3*x.*y.*(3*r.^2 - 5*x.^2)./r.^7, ...
            (-15*x.^2.*y.^2 + 3*(x.^2 + y.^2).*r.^2 - r.^4)./r.^7, ...
            3*x.*y.*(3*r.^2 - 5*y.^2)./r.^7, ...
            3*(-5*y.^4 + 6*y.^2.*r.^2 - r.^4)./r.^7)};
L = L(1:nout);  R = R(1:nout);
end

function s = snap(c)
s = sum(c);
if ( abs(s) <= 1e-12*sum(abs(c)) ), s = 0; end
end
