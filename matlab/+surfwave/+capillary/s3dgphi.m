function val = s3dgphi(rts,ejs,src,targ)
%
% computes S3d Gphi
%
% outputs are:
% - val is the value of S3d Gphi centered at zero and
%   evaluated at (x,y)
%
% input:
%
% x - x-coordinates array
% y - y-coordinates array
% rts - roots of cubic polynomial
% ejs - residues (see notes)
%

% src = src.r;
% targ = targ.r;

rslf = 1e-14;

if isstruct(src)
    src = src.r;
end

if isstruct(targ)
    targ = targ.r;
end

eulergamma = 0.57721566490153286060651209008240243;

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

r = sqrt((xt-xs).^2 + (yt-ys).^2);

% clear xs ys xt yt

val = 0;

for i = 1:3
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)

       zt = r*rhoj;
       [cr0] = surfwave.capillary.struveR(zt);

       h0 = chnk.helm2d.green(rhoj,src,targ);
       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       h0 = -4i*h0;
        
       sk0 = -1i*cr0+1i*h0;
       
       val = val + ej*(-sk0 + 2i*h0);

    elseif abs(rhoj) > rslf

       ilow = (imag(rhoj) > 0);

        if ilow
            rhoj = conj(rhoj); % flip -rhoj up into upper half plane
        end

       zt = -r*rhoj;
       [cr0] = surfwave.capillary.struveR(zt);

       h0 = chnk.helm2d.green(-rhoj,src,targ);

       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
       h0 = -4i*h0;

       sk0 = -1i*cr0+1i*h0;

       if ilow
           sk0 = conj(sk0);
           rhoj = conj(rhoj); % flip -rhoj up into upper half plane
       end

       val = val + ej*sk0;

    end

end

val = 1/8*val;

end