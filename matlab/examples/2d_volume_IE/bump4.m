function [icoefs,pcoefs,H,alpha] = bump4(targinfo,ifbdry)
    % Coefficients for the variable-thickness plate problem.
    % H = H0 + amp*exp(-(x^2+y^2)/(2*width^2))

    w     = 1;
    E     = 1;
    nu    = 0.33;
    H0    = 1;
    rhoi  = 1;
    amp   = 1;
    width = 1;

    if nargin < 2
        ifbdry = 0;
    end

    X = targinfo.r(1,:);
    Y = targinfo.r(2,:);

    a0   = E*H0^3/(12*(1-nu^2));
    b0   = rhoi*H0*w^2;
    Hbar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));
    H    = H0 + Hbar;
    beta = rhoi*H*w^2;

    alpha = E*H.^3/(12*(1-nu^2));

    Hx = -X.*Hbar/width^2;
    Hy = -Y.*Hbar/width^2;

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    if ifbdry
        dx = targinfo.d(1,:);
        dy = targinfo.d(2,:);
        ds = sqrt(dx.*dx + dy.*dy);

        taux = dx./ds;
        tauy = dy./ds;
        nx   = targinfo.n(1,:);
        ny   = targinfo.n(2,:);

        dadn = ax.*nx  + ay.*ny;
        dadt = ax.*taux + ay.*tauy;

        icoefs = {alpha, dadn, dadt};
    else
        abar = alpha - a0;
        bbar = beta  - b0;

        Hxx = -Hbar/width^2 - X.*Hx/width^2;
        Hxy = -X.*Hy/width^2;
        Hyy = -Hbar/width^2 - Y.*Hy/width^2;

        axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
        axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
        ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

        pcoefs = {a0/E, abar/E, b0/E, bbar/E, ax/E, ay/E, axx/E, axy/E, ayy/E, nu};

        icoefs = zeros(7,length(X));
        icoefs(1,:) =  2*ax;
        icoefs(2,:) =  2*ay;
        icoefs(3,:) =  axx + ayy;
        icoefs(4,:) = -(1-nu)*axx;
        icoefs(5,:) = -(1-nu)*ayy;
        icoefs(6,:) =  2*(1-nu)*axy;
        icoefs(7,:) =  (alpha*b0 - a0*beta) / a0;

        icoefs = icoefs*a0 ./ alpha;
        icoefs = icoefs/E;

        alpha = alpha/E;
    end
end
