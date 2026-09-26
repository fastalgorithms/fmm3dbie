function [val] = gsgrav(rho, ejs, src, targ)
% GSGRAV  Gravity Green's function
%
%   G_S(x,y) = ej*rho/(2 pi |x-y|)
%              + (ej*rho^2/4) [ -K0(rho|x-y|) + 2i H0^(1)(rho|x-y|) ]
%
% G_phi = G_S/g with g = 2*rho, i.e. gsgrav(rho, 1/g, ...).
%
% Inputs:
%   rts  - pole location rho (= g/2); may be a vector (sums over roots)
%   ejs  - residue ej at each pole (matching rts elementwise)
%   src  - (2,ns) or struct with .r
%   targ - (2,nt) or struct with .r

    rslf = 1e-14;
    if isstruct(src),  src  = src.r;  end
    if isstruct(targ), targ = targ.r; end
    src = src(1:2,:);  targ = targ(1:2,:);
    eulergamma = 0.57721566490153286060651209008240243;

    [~,ns] = size(src);  [~,nt] = size(targ);
    xs = repmat(src(1,:),nt,1);    ys = repmat(src(2,:),nt,1);
    xt = repmat(targ(1,:).',1,ns); yt = repmat(targ(2,:).',1,ns);
    r  = sqrt((xt-xs).^2 + (yt-ys).^2);

    rho = rho(:).';  ejs = ejs(:).';
    val = 0;  algsum = 0;
    for i = 1:numel(rho)
        rhoj = rho(i);  ej = ejs(i);
        algsum = algsum + ej*rhoj;
        if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)
            zt = r*rhoj;
            cr0 = surfwave.struveR(zt);
            h0  = green2d.helm(rhoj, src, targ);
            h0(r < rslf) = 1/(2*pi)*(1i*pi/2 - eulergamma + log(2/rhoj));
            h0  = -4i*h0;
            sk0 = -1i*cr0 + 1i*h0;
            val = val + ej*rhoj^2*(-sk0 + 2i*h0);
        elseif abs(rhoj) > rslf
            ilow = (imag(rhoj) > 0);
            if ilow, rhoj = conj(rhoj); end
            zt = -r*rhoj;
            cr0 = surfwave.struveR(zt);
            h0  = green2d.helm(-rhoj, src, targ);
            h0(r < rslf) = 1/(2*pi)*(1i*pi/2 - eulergamma + log(2/rhoj));
            h0  = -4i*h0;
            sk0 = -1i*cr0 + 1i*h0;
            if ilow, sk0 = conj(sk0); rhoj = conj(rhoj); end
            val = val + ej*rhoj^2*sk0;
        end
    end
    val = val/4;
    % algebraic 1/r term, ej*rho/(2 pi r) summed over all roots
    alg = algsum ./ (2*pi*r);
    alg(r < rslf) = 0;
    val = val + alg;

    val(r < rslf) = 0;
end
