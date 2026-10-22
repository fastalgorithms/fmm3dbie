function [val,grad,hess,third,bilap] = gphiflex(rts,ejs,src,targ)
%SURFWAVE.FLEX.GPHIFLEX   flexural Green's function G_phi, its target
%   derivatives and its bilaplacian.
%
%   [val,grad,hess,third,bilap] = SURFWAVE.FLEX.GPHIFLEX(rts,ejs,src,targ)
%
%   for the dispersion roots rts and residues ejs
%   (SURFWAVE.FLEX.FIND_ROOTS_FLEX): val (nt,ns), grad(:,:,1:2) = G_x,
%   G_y, hess(:,:,1:3) = G_xx, G_xy, G_yy, third(:,:,1:4) and bilap
%   = Lap^2 G_phi (away from r = 0).
%
% See also SURFWAVE.FLEX.GREEN, SURFWAVE.ROOTSUM

rts = rts(:);  ejs = ejs(:);
nout = max(nargout, 1);
out = cell(1, min(nout, 4) + (nout > 4));
[out{:}] = surfwave.rootsum(rts, ejs.*rts/4, src, targ, 1e-8);
val = out{1};
if ( nout > 1 ), grad  = out{2}; end
if ( nout > 2 ), hess  = out{3}; end
if ( nout > 3 ), third = out{4}; end
if ( nout > 4 )
    f = out{5};
    bilap = f(:,:,1) + 2*f(:,:,3) + f(:,:,5);
end
end
