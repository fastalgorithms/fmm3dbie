function val = s3dgphiflex(rts,ejs,src,targ)
%
% Evaluates the 3D-component contribution to the flexural-gravity wave
% Green's function, G_{3d,phi} = S_{3d}[G_phi], built from the 5 roots
% of the dispersion polynomial  alpha*z^5 + gamma*z - 1 = 0.
% This is the (1/8)*sum_j ej*(-sk0_j + 2i*h0_j) combination used in the
% post-processing evaluation of phi.
%
% Output:
%   val  - (nt,ns) values of G_{3d,phi}(targ,src)
%
% Input:
%   rts  - (5,1) roots of the dispersion polynomial
%   ejs  - (5,1) partial-fraction residues corresponding to rts
%   src  - (2,ns) source positions, or struct with field .r
%   targ - (2,nt) target positions, or struct with field .r
%

rslf = 1e-14;

if isstruct(src)
    src = src.r;
    targ = targ.r;
end

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

r = sqrt(dx.*dx + dy.*dy);

val = 0;
gradx = 0;
grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;
thirdxxx = 0;
thirdxxy = 0;
thirdxyy = 0;
thirdyyy = 0;
bilap = 0;

oneoveralpha = sum(ejs.*rts.^4);

eulergamma = 0.5772156649015328606065120900824;


for i = 1:5
    
    rhoj = rts(i);
    ej = ejs(i);

    if (abs(angle(rhoj)) < 1e-8) && (abs(rhoj) > rslf)

       sk0 = surfwave.flex.struveKdiffgreen(rhoj,src,targ);
       h0 = surfwave.flex.helmdiffgreen(rhoj,src,targ);

       h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       
       val = val + ej*(-sk0 + 2i*h0);


    elseif abs(rhoj) > rslf

       sk0 = surfwave.flex.struveKdiffgreen(-rhoj,src,targ);

       val = val + ej*sk0;

    end

end

val = 1/8*val;

end