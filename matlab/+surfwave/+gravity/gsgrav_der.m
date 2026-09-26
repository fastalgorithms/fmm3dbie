function [val,grad,hess,third,fourth] = gsgrav_der(rts,ejs,src,targ)
% GSGRAV_DER  Gravity Green's function G_S (see GSGRAV) and its target
% derivatives. Evaluated without singularity subtraction, so accurate
% only away from r = 0.
%
%   Output ordering (derivatives on TARGET coords):
%     val   : nt x ns ; grad nt x ns x 2 ; hess x3 ; third x4 ; fourth x5

    rslf = 1e-14;
    if isstruct(src),  src  = src.r;  end
    if isstruct(targ), targ = targ.r; end
    src = src(1:2,:);  targ = targ(1:2,:);
    rts = rts(:).';  ejs = ejs(:).';

    [~,ns] = size(src);
    [~,nt] = size(targ);
    xs = repmat(src(1,:),nt,1);    ys = repmat(src(2,:),nt,1);
    xt = repmat(targ(1,:).',1,ns); yt = repmat(targ(2,:).',1,ns);
    dx = xt - xs;  dy = yt - ys;
    r2 = dx.*dx + dy.*dy;  r = sqrt(r2);
    z  = r < rslf;
    ir = 1./r;

    % algebraic coefficient ej*rho/(2 pi) summed over roots
    calg = sum(ejs.*rts)/2;

    %% ---- special-function radial profile (sum over all roots) ----
    f = 0; fr = 0; frr = 0; frrr = 0; frrrr = 0;
    for iroot = 1:numel(rts)
        rhoj = rts(iroot);  ej = ejs(iroot);

        if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)
        % ---- real-positive (propagating) root: sk0 and h0 = H0^(1) ----
        zt = rhoj*r;

        % Struve R0,R1 raw:  R0'(z) = -R1 + 2i/pi ; R1'(z) = R0 - R1/z.
        [R0,R1] = surfwave.struveR(zt);
        twoi_pi = 2i/pi;
        R0p  = -R1 + twoi_pi;
        R1p  = R0 - R1./zt;
        R0pp = -R1p;
        R1pp = R0p - R1p./zt + R1./zt.^2;
        R0ppp = -R1pp;
        R1ppp = R0pp - R1pp./zt + 2*R1p./zt.^2 - 2*R1./zt.^3;
        R0pppp = -R1ppp;

        % Hankel H0^(1): H0'=-H1 ; H1'=H0 - H1/z.
        H0 = besselh(0,1,zt);  H1 = besselh(1,1,zt);
        H0p  = -H1;
        H1p  = H0 - H1./zt;
        H0pp = -H1p;
        H1pp = H0p - H1p./zt + H1./zt.^2;
        H0ppp= -H1pp;
        H1ppp = H0pp - H1pp./zt + 2*H1p./zt.^2 - 2*H1./zt.^3;
        H0pppp = -H1ppp;

        % per-root: (1/4) ej rho^2 (-sk0 + 2i h0) = (1/4) ej rho^2 * i (R0 + H0)
        cc = 0.25i*ej*rhoj^2;
        q    = cc*(R0    + H0);
        qz   = cc*(R0p   + H0p);
        qzz  = cc*(R0pp  + H0pp);
        qzzz = cc*(R0ppp + H0ppp);
        qzzzz= cc*(R0pppp+ H0pppp);

        % chain to r:  d/dr = rhoj d/dz  (accumulate over roots)
        f     = f     + q;
        fr    = fr    + rhoj  *qz;
        frr   = frr   + rhoj^2*qzz;
        frrr  = frrr  + rhoj^3*qzzz;
        frrrr = frrrr + rhoj^4*qzzzz;

        elseif abs(rhoj) > rslf
        % ---- complex/evanescent root: sk0 only, no outgoing H0^(1) piece ----
        % zt = -rho*r, conjugated for stability off the real axis
        ilow = (imag(rhoj) > 0);
        rr = rhoj;
        if ilow, rr = conj(rhoj); end
        zt = -rr*r;

        [R0,R1] = surfwave.struveR(zt);
        twoi_pi = 2i/pi;
        R0p  = -R1 + twoi_pi;
        R1p  = R0 - R1./zt;
        R0pp = -R1p;
        R1pp = R0p - R1p./zt + R1./zt.^2;
        R0ppp = -R1pp;
        R1ppp = R0pp - R1pp./zt + 2*R1p./zt.^2 - 2*R1./zt.^3;
        R0pppp = -R1ppp;

        H0 = besselh(0,1,zt);  H1 = besselh(1,1,zt);
        H0p  = -H1;
        H1p  = H0 - H1./zt;
        H0pp = -H1p;
        H1pp = H0p - H1p./zt + H1./zt.^2;
        H0ppp= -H1pp;
        H1ppp = H0pp - H1pp./zt + 2*H1p./zt.^2 - 2*H1./zt.^3;
        H0pppp = -H1ppp;

        % per-root: (1/4) ej rho^2 sk0 = (1/4) ej rho^2 * i (H0 - R0)
        cc = 0.25i*ej*rr^2;
        q    = cc*(H0     - R0);
        qz   = cc*(H0p    - R0p);
        qzz  = cc*(H0pp   - R0pp);
        qzzz = cc*(H0ppp  - R0ppp);
        qzzzz= cc*(H0pppp - R0pppp);

        if ilow
            q=conj(q); qz=conj(qz); qzz=conj(qzz); qzzz=conj(qzzz); qzzzz=conj(qzzzz);
        end

        % chain to r:  d/dr = -rhoj d/dz
        f     = f     + q;
        fr    = fr    - rhoj  *qz;
        frr   = frr   + rhoj^2*qzz;
        frrr  = frrr  - rhoj^3*qzzz;
        frrrr = frrrr + rhoj^4*qzzzz;
        end
    end

    %% ---- algebraic term  calg/(pi r) : radial profile and r-derivs ----
    a = calg/pi;
    g0   =  a.*ir;
    gr   = -a.*ir.^2;
    grr  =  2*a.*ir.^3;
    grrr = -6*a.*ir.^4;
    grrrr=  24*a.*ir.^5;

    % total radial profile F(r) = f + g0 and r-derivatives
    F    = f    + g0;
    Fr   = fr   + gr;
    Frr  = frr  + grr;
    Frrr = frrr + grrr;
    Frrrr= frrrr+ grrrr;

    %% ---- chain radial -> Cartesian (target) derivatives ----
    % For a radial function F(r), with r = sqrt(dx^2+dy^2):
    %   F_x   = Fr * (dx/r)
    %   F_xx  = Frr*(dx^2/r^2) + Fr*(dy^2/r^3)
    %   F_xy  = (Frr - Fr/r) * dx*dy/r^2
    %   F_xxx = Frrr*(dx^3/r^3) + 3*(Frr/r - Fr/r^2)*(dx*dy^2/r^2)
    ux = dx.*ir;  uy = dy.*ir;            % unit radial components
    ir2 = ir.^2;  ir3 = ir.^3;

    Fx = Fr.*ux;
    Fy = Fr.*uy;

    Fxx = Frr.*ux.^2 + Fr.*(dy.^2).*ir3;
    Fyy = Frr.*uy.^2 + Fr.*(dx.^2).*ir3;
    Fxy = (Frr - Fr.*ir).*(dx.*dy).*ir2;

    % third derivatives via explicit closed forms (ddx helper)
    Fxxx = ddx(Frr,Fr,Frrr, dx,dy,r, 'xx_x');
    Fxxy = ddx(Frr,Fr,Frrr, dx,dy,r, 'xx_y');
    Fxyy = ddx(Frr,Fr,Frrr, dx,dy,r, 'yy_x');
    Fyyy = ddx(Frr,Fr,Frrr, dx,dy,r, 'yy_y');

    %% zero self entries
    zz = @(M) M.*(~z);
    val   = zz(F);
    grad  = cat(3, zz(Fx),  zz(Fy));
    hess  = cat(3, zz(Fxx), zz(Fxy), zz(Fyy));
    third = cat(3, zz(Fxxx),zz(Fxxy),zz(Fxyy),zz(Fyyy));

    if nargout > 4
        % fourth Cartesian derivatives (ddx4 helper), order:
        %   [G_xxxx, G_xxxy, G_xxyy, G_xyyy, G_yyyy]
        Fxxxx = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, 'xxxx');
        Fxxxy = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, 'xxxy');
        Fxxyy = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, 'xxyy');
        Fxyyy = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, 'xyyy');
        Fyyyy = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, 'yyyy');
        fourth = cat(3, zz(Fxxxx),zz(Fxxxy),zz(Fxxyy),zz(Fxyyy),zz(Fyyyy));
    end
end

function out = ddx4(Fr,Frr,Frrr,Frrrr, dx,dy,r, which)
% Fourth Cartesian derivatives of a radial function F(r) (r>0).
    ir=1./r; ir2=ir.^2; ir3=ir.^3; ir4=ir.^4;
    ux=dx.*ir; uy=dy.*ir;
    % coefficients of the radial-derivative chain for 4th order:
    %  F_ijkl = a4 u_i u_j u_k u_l
    %         + a3 ( sum over pairs  delta_(ij) u_k u_l )
    %         + a2 ( sum  delta_(ij) delta_(kl) )
    % with
    a4 = Frrrr - 6*Frrr.*ir + 15*Frr.*ir2 - 15*Fr.*ir3;
    a3 = (Frrr.*ir - 3*Frr.*ir2 + 3*Fr.*ir3);
    a2 = (Frr.*ir2 - Fr.*ir3);
    % symmetric tensor assembly for 2D indices (x=1,y=2)
    switch which
        case 'xxxx'
            out = a4.*ux.^4 + a3.*(6*ux.^2) + a2.*3;
        case 'yyyy'
            out = a4.*uy.^4 + a3.*(6*uy.^2) + a2.*3;
        case 'xxxy'
            out = a4.*ux.^3.*uy + a3.*(3*ux.*uy) ;
        case 'xyyy'
            out = a4.*ux.*uy.^3 + a3.*(3*ux.*uy) ;
        case 'xxyy'
            out = a4.*ux.^2.*uy.^2 + a3.*(ux.^2+uy.^2) + a2 ;
    end
end

function out = ddx(Frr,Fr,Frrr, dx,dy,r, which)
% Third Cartesian derivatives of a radial function F(r), given Fr,Frr,Frrr.
% Exact closed forms (r>0):
    ir=1./r; ir2=ir.^2; ir3=ir.^3; ir4=ir.^4; ir5=ir.^5;
    dx2=dx.^2; dy2=dy.^2;
    switch which
        case 'xx_x'   % G_xxx
            out = Frrr.*dx.^3.*ir3 ...
                + Frr.*(3*dx.*dy2).*ir4 ...
                - Fr .*(3*dx.*dy2).*ir5;
        case 'xx_y'   % G_xxy = d/dy G_xx
            out = Frrr.*dx2.*dy.*ir3 ...
                + Frr.*(dy.^3 - 2*dx2.*dy).*ir4 ...
                - Fr .*(dy.^3 - 2*dx2.*dy).*ir5;
        case 'yy_x'   % G_xyy = d/dx G_yy
            out = Frrr.*dy2.*dx.*ir3 ...
                + Frr.*(dx.^3 - 2*dy2.*dx).*ir4 ...
                - Fr .*(dx.^3 - 2*dy2.*dx).*ir5;
        case 'yy_y'   % G_yyy
            out = Frrr.*dy.^3.*ir3 ...
                + Frr.*(3*dy.*dx2).*ir4 ...
                - Fr .*(3*dy.*dx2).*ir5;
    end
end
