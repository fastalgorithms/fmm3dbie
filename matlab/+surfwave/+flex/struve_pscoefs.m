function cf1 = struve_pscoefs(nterms)
% powerseries coefficients for the struve function H
%
% these are series with the terms z,z^3,z^5,..z^(2*nterms-1)

pow1 = 0.5;

cf1 = zeros(nterms,1);
sgn = 1;

for k = 0:nterms-1
    cf1(k+1) = sgn*pow1/gamma(k+3/2)^2;
    pow1 = pow1*0.25;
    sgn = -sgn;
end

end