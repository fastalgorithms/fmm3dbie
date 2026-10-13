function comps = radpp_tensors(dx,dy,dz,breaks,coefs,nderiv)
%RADPP_TENSORS derivative tensors of a radial function from its on-axis
%  profiles, the matlab twin of radpp_tensors in src/kernels/radpp.f90.
%
%  comps = radpp_tensors(dx,dy,dz,breaks,coefs,nderiv)
%
%  dx, dy, dz are displacements target - source (any shape, dz may be 0).
%  breaks, coefs are a radpp_fit of the on-axis profiles
%    nderiv = 0:  G
%    nderiv = 1:  + G_x
%    nderiv = 2:  + G_xx, G_yy
%    nderiv = 3:  + G_xxx, G_xyy
%  comps is (numel(dx), ncomp), ncomp = 1, 3, 6, 10, in chunkie order
%    [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy]
%  with derivatives taken with respect to the target.
%
%  See also RADPP_FIT, KERNEL3D.RADPP

dx = dx(:);  dy = dy(:);  dz = dz(:);
r = sqrt(dx.^2 + dy.^2 + dz.^2);
c = dx./r;   s = dy./r;
c(r == 0) = 1;  s(r == 0) = 0;

prof = radpp_eval(r,breaks,coefs);

comps = prof(:,1);
if nderiv >= 1
    gr = prof(:,2);
    comps = [comps, gr.*c, gr.*s];
end
if nderiv >= 2
    p = prof(:,3);  q = prof(:,4);
    comps = [comps, p.*c.^2 + q.*s.^2, (p-q).*c.*s, p.*s.^2 + q.*c.^2];
end
if nderiv >= 3
    a = prof(:,5);  b = prof(:,6);
    comps = [comps, a.*c.^3 + 3*b.*c.*s.^2, ...
                    a.*c.^2.*s + b.*(s.^3 - 2*c.^2.*s), ...
                    a.*c.*s.^2 + b.*(c.^3 - 2*c.*s.^2), ...
                    a.*s.^3 + 3*b.*s.*c.^2];
end

end
