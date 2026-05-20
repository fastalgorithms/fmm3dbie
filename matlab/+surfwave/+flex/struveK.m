function [ck0,ck1] = struveK(rhoj, z)
    % where rhoj is a complex number and z is an array of real numbers
    % the hankel part of this function has log subtracted

    rslf = 1e-14;
    sz = size(z);
    z = z(:).';
    
    ilow = (imag(rhoj) < 0);

    if ilow
        rhoj = conj(rhoj); % flip rhoj up into upper half plane
    end

    eulergamma = 0.5772156649015328606065120900824;

    zt = z*rhoj;

    [cr0,cr1] = struveR(zt);

    src = [0; 0]; targ = [ z ; 0*z];
    [h0,h1] = surfwave.flex.helmdiffgreen(rhoj,src,targ);
    h0(zt < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
    h0 = -4i*h0.' ;
    h1 = 4i*h1(:,:,1).' ./ rhoj; 
    h1(zt < rslf) = 0; 

    ck0 = -1i*cr0+1i*h0;
    ck1 = -1i*cr1+1i*h1;

    if ilow
        ck0 = conj(ck0);
        ck1 = conj(ck1);
    end

    ck0 = reshape(ck0,sz);
    ck1 = reshape(ck1,sz);

end
