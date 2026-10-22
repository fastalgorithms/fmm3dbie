function varargout = gsflex(rts,ejs,src,targ,ifr2logr)
%SURFWAVE.FLEX.GSFLEX   flexural Green's function G_S and its target
%   derivatives.
%
%   [val,grad,hess,third,fourth] = SURFWAVE.FLEX.GSFLEX(rts,ejs,src,targ)
%
%   for the dispersion roots rts and residues ejs
%   (SURFWAVE.FLEX.FIND_ROOTS_FLEX): val (nt,ns), grad(:,:,1:2) = G_x,
%   G_y, hess(:,:,1:3) = G_xx, G_xy, G_yy, third(:,:,1:4) and
%   fourth(:,:,1:5) likewise.
%
%   SURFWAVE.FLEX.GSFLEX(rts,ejs,src,targ,true) returns G_S minus its
%   r^2 log r part, sum(ejs.*rts.^4)/(4 pi) r^2 log r.
%
% See also SURFWAVE.FLEX.GREEN, SURFWAVE.ROOTSUM

if ( nargin < 5 ), ifr2logr = false; end
rts = rts(:);  ejs = ejs(:);
nout = max(nargout, 1);
[varargout{1:nout}] = surfwave.rootsum(rts, ejs.*rts.^2/2, src, targ, 1e-8);
if ( ifr2logr )
    D = r2logr(src, targ, nout);
    c = sum(ejs.*rts.^4)/(4*pi);
    for k = 1:nout, varargout{k} = varargout{k} - c*D{k}; end
end
end


function D = r2logr(src, targ, nout)
% r^2 log r and its target derivatives
if ( isstruct(src) ),  src  = src.r;  end
if ( isstruct(targ) ), targ = targ.r; end
x = targ(1,:).' - src(1,:);
y = targ(2,:).' - src(2,:);
r2 = x.^2 + y.^2;  r = sqrt(r2);  lr = log(r);
D = {r2.*lr, cat(3, x.*(2*lr + 1), y.*(2*lr + 1)), ...
     cat(3, 2*x.^2./r2 + 2*lr + 1, 2*x.*y./r2, 2*y.^2./r2 + 2*lr + 1), ...
     cat(3, 2*x.*(x.^2 + 3*y.^2), 2*y.*(y.^2 - x.^2), ...
            2*x.*(x.^2 - y.^2), 2*y.*(3*x.^2 + y.^2))./r2.^2, ...
     cat(3, 2*(-x.^4 - 6*x.^2.*y.^2 + 3*y.^4), 4*x.*y.*(x.^2 - 3*y.^2), ...
            2*(-x.^4 + 6*x.^2.*y.^2 - y.^4), 4*x.*y.*(y.^2 - 3*x.^2), ...
            2*(3*x.^4 - 6*x.^2.*y.^2 - y.^4))./r2.^3};
D = D(1:nout);
end
