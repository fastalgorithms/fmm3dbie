function varargout = green(zk,src,targ)
%FLEX2D.GREEN  2D flexural (thin-plate) Green's function and derivatives.
%
%  [val,grad,hess,der3,der4,der5] = flex2d.green(zk,src,targ)
%
%  Evaluates the Green's function used by FLEX2D.KERN,
%
%     G(x,y) = 1/(zk1^2 - zk2^2) * (G_{zk1}(x,y) - G_{zk2}(x,y)),
%
%  where G_k is the (Laplace-subtracted) 2D Helmholtz Green's function
%  of SURFWAVE.FLEX.HELMDIFFGREEN. The wavenumbers are chosen from zk
%  exactly as in FLEX2D.KERN:
%
%     scalar zk, |zk| < 1e-6 : biharmonic, G = |x-y|^2 log|x-y| / (8 pi)
%     scalar zk              : flexural,   [zk1, zk2] = [zk, 1i*zk]
%     [zk1, zk2], one ~ 0    : Stokes-like, zk2 = 0
%     [zk1, zk2]             : general pair
%
%  src (2,ns) and targ (2,nt) are points in the plane (further rows are
%  ignored). All outputs are (nt,ns,:) with derivatives taken in the
%  target variables, ordered as in SURFWAVE.FLEX.HELMDIFFGREEN:
%     grad: G_x, G_y
%     hess: G_xx, G_xy, G_yy
%     der3: G_xxx, G_xxy, G_xyy, G_yyy
%     der4, der5 likewise.

nout = max(nargout,1);
src  = src(1:2,:);
targ = targ(1:2,:);

zk = zk(:).';

if isscalar(zk)
    if abs(zk) < 1e-6
        % biharmonic: no difference formula needed
        [varargout{1:nout}] = flex2d.bhgreen(src,targ);
        return
    end
    zk1 = zk;
    zk2 = 1i*zk;
elseif numel(zk) == 2
    if any(abs(zk) < 1e-6)
        zk1 = zk(abs(zk) >= 1e-6);
        zk2 = 0;
    else
        zk1 = zk(1);
        zk2 = zk(2);
    end
else
    error('FLEX2D.GREEN: zk must be a scalar or a pair [zk1, zk2].');
end

out1 = cell(1,nout);
[out1{:}] = flex2d.helmdiffgreen(zk1,src,targ);

scale = 1/(zk1^2 - zk2^2);

if zk2 == 0
    varargout = cellfun(@(a) scale*a, out1, 'UniformOutput', false);
else
    out2 = cell(1,nout);
    [out2{:}] = flex2d.helmdiffgreen(zk2,src,targ);
    varargout = cellfun(@(a,b) scale*(a-b), out1, out2, 'UniformOutput', false);
end

end
