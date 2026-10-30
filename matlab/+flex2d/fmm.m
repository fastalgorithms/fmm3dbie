function U = fmm(eps, zk, srcuse, pg, targ, pgt, iflap)
%FLEX2D.FMM   FMM for flex2d.green, or its Laplacian, with the interface
% of HFMM2D. Charges only, pg = 0.
%
% Syntax: U = flex2d.fmm(eps, zk, srcuse, pg, targ, pgt)
%         U = flex2d.fmm(eps, zk, srcuse, pg, targ, pgt, iflap)
%
% With iflap = true the green is H = Lap G instead of G.
%
% Uses G = (G_{zk1} - G_{zk2})/(zk1^2 - zk2^2), with G_k the Helmholtz
% Green's function (the Laplace one if zk2 = 0), so each call is two
% FMM2D calls (hfmm2d, rfmm2d), and Lap G_k = -k^2 G_k for H. The
% biharmonic case (zk = 0) uses
%
%   |x-y|^2 log|x-y| = (|x|^2 - 2 x.y + |y|^2) log|x-y|,
%
% i.e. four Laplace FMMs with charges sigma, sigma*y1, sigma*y2 and
% sigma*|y|^2, done in one rfmm2d call, with the origin shifted to the
% centroid of the sources to limit cancellation, and H = (log r + 1)/(2 pi).
%
% see also FLEX2D.GREEN, FMM2DTO3D

if ( nargin < 7 ), iflap = false; end
if ( pg ~= 0 || isfield(srcuse, 'dipoles') )
    error('FLEX2D:fmm:pg', 'Only charges and target outputs.');
end
nd = 1;
if ( isfield(srcuse, 'nd') ), nd = srcuse.nd; end
ns = size(srcuse.sources, 2);
nt = size(targ, 2);
q = reshape(srcuse.charges, nd, ns);
src = srcuse.sources(1:2,:);
t = targ(1:2,:);

U = [];
U.pottarg = zeros(nd, nt);
if ( pgt > 1 ), U.gradtarg = zeros(nd, 2, nt); end
if ( pgt > 2 ), U.hesstarg = zeros(nd, 3, nt); end
for j = 1:nd
    [val, grad, hess] = green_fmm(eps, zk, src, t, q(j,:), pgt, iflap);
    U.pottarg(j,:) = val;
    if ( pgt > 1 ), U.gradtarg(j,:,:) = reshape(grad, 1, 2, nt); end
    if ( pgt > 2 ), U.hesstarg(j,:,:) = reshape(hess, 1, 3, nt); end
end

end


function [val, grad, hess] = green_fmm(eps, zk, src, targ, sigma, pgt, iflap)
% G (or H = Lap G) with charges sigma, its gradient (2,nt) and Hessian
% (3,nt: xx, xy, yy)

zk = zk(:).';
if ( isscalar(zk) ), zk = [zk, 1i*zk]; end
if ( any(abs(zk) < 1e-6) )
    zk1 = zk(abs(zk) >= 1e-6);
    zk2 = 0;
    isbh = isempty(zk1);
else
    zk1 = zk(1);
    zk2 = zk(2);
    isbh = false;
end

grad = []; hess = [];
if ( isbh && iflap )
    srcuse = [];
    srcuse.sources = src;
    srcuse.nd = 2;
    srcuse.charges = [real(sigma); imag(sigma)]/(2*pi);
    U = rfmm2d(eps, srcuse, 0, targ, pgt);
    nt = size(targ, 2);
    P = reshape(U.pottarg, 2, nt);
    val = P(1,:) + 1i*P(2,:) + sum(sigma)/(2*pi);
    if ( pgt > 1 )
        G = reshape(U.gradtarg, 2, 2, nt);
        grad = reshape(G(1,:,:) + 1i*G(2,:,:), 2, nt);
    end
    if ( pgt > 2 )
        H = reshape(U.hesstarg, 2, 3, nt);
        hess = reshape(H(1,:,:) + 1i*H(2,:,:), 3, nt);
    end
    return
end
if ( isbh )
    [val, grad, hess] = bh_fmm(eps, src, targ, sigma, pgt);
    return
end

a1 = 1;  a2 = 1;
if ( iflap ), a1 = -zk1^2;  a2 = -zk2^2; end

srcuse = [];
srcuse.sources = src;
srcuse.charges = sigma/(zk1^2 - zk2^2);

% G_{zk1} part
U = hfmm2d(eps, zk1, srcuse, 0, targ, pgt);

% minus G_{zk2} part
if ( zk2 ~= 0 )
    U2 = hfmm2d(eps, zk2, srcuse, 0, targ, pgt);
else
    % Laplace, G_0 = -log|x-y|/(2 pi), real charges only
    c = srcuse.charges;
    srcuse.charges = -real(c)/(2*pi);
    U2 = rfmm2d(eps, srcuse, 0, targ, pgt);
    srcuse.charges = -imag(c)/(2*pi);
    U2i = rfmm2d(eps, srcuse, 0, targ, pgt);
    U2.pottarg = U2.pottarg + 1i*U2i.pottarg;
    if ( pgt > 1 ), U2.gradtarg = U2.gradtarg + 1i*U2i.gradtarg; end
    if ( pgt > 2 ), U2.hesstarg = U2.hesstarg + 1i*U2i.hesstarg; end
end

val = a1*U.pottarg - a2*U2.pottarg;
if ( pgt > 1 ), grad = a1*U.gradtarg - a2*U2.gradtarg; end
if ( pgt > 2 ), hess = a1*U.hesstarg - a2*U2.hesstarg; end

end


function [val, grad, hess] = bh_fmm(eps, src, targ, sigma, pgt)
%BH_FMM  sum_j sigma_j |x-y_j|^2 log|x-y_j| / (8 pi) and its gradient and
% Hessian (xx, xy, yy) at the targets, from Laplace FMMs.
%
% With A, B = (B1,B2), C the log potentials of the charges sigma,
% sigma*y, sigma*|y|^2,
%
%   8 pi u = |x|^2 A - 2 x.B + C

nt = size(targ, 2);

% shift the origin to limit cancellation
cen = mean(src, 2);
y = src - cen;
x = targ - cen;
x1 = x(1,:); x2 = x(2,:); x2n = x1.^2 + x2.^2;

sig = sigma(:).';
q = [sig; sig.*y(1,:); sig.*y(2,:); sig.*(y(1,:).^2 + y(2,:).^2)];

% rfmm2d takes real charges, kernel log|x-y|
srcuse = [];
srcuse.sources = src;
srcuse.nd = 8;
srcuse.charges = [real(q); imag(q)];
U = rfmm2d(eps, srcuse, 0, targ, pgt);

P = reshape(U.pottarg, 8, nt);
P = P(1:4,:) + 1i*P(5:8,:);
A = P(1,:); B1 = P(2,:); B2 = P(3,:); C = P(4,:);
val = (x2n.*A - 2*(x1.*B1 + x2.*B2) + C)/(8*pi);

grad = []; hess = [];
if ( pgt > 1 )
    G = reshape(U.gradtarg, 8, 2, nt);
    G = G(1:4,:,:) + 1i*G(5:8,:,:);
    Ax  = reshape(G(1,1,:),1,nt); Ay  = reshape(G(1,2,:),1,nt);
    B1x = reshape(G(2,1,:),1,nt); B1y = reshape(G(2,2,:),1,nt);
    B2x = reshape(G(3,1,:),1,nt); B2y = reshape(G(3,2,:),1,nt);
    Cx  = reshape(G(4,1,:),1,nt); Cy  = reshape(G(4,2,:),1,nt);
    grad = zeros(2, nt, 'like', val);
    grad(1,:) = 2*x1.*A + x2n.*Ax - 2*B1 - 2*(x1.*B1x + x2.*B2x) + Cx;
    grad(2,:) = 2*x2.*A + x2n.*Ay - 2*B2 - 2*(x1.*B1y + x2.*B2y) + Cy;
    grad = grad/(8*pi);
end
if ( pgt > 2 )
    H = reshape(U.hesstarg, 8, 3, nt);
    H = H(1:4,:,:) + 1i*H(5:8,:,:);
    h = @(k,l) reshape(H(k,l,:),1,nt);
    hess = zeros(3, nt, 'like', val);
    hess(1,:) = 2*A + 4*x1.*Ax + x2n.*h(1,1) - 4*B1x ...
        - 2*(x1.*h(2,1) + x2.*h(3,1)) + h(4,1);
    hess(2,:) = 2*x1.*Ay + 2*x2.*Ax + x2n.*h(1,2) - 2*B1y - 2*B2x ...
        - 2*(x1.*h(2,2) + x2.*h(3,2)) + h(4,2);
    hess(3,:) = 2*A + 4*x2.*Ay + x2n.*h(1,3) - 4*B2y ...
        - 2*(x1.*h(2,3) + x2.*h(3,3)) + h(4,3);
    hess = hess/(8*pi);
end

end
