function val = radcheb_eval(r,breaks,coefs)
%RADCHEB_EVAL evaluate a radcheb interpolant, vectorized in r
%
%  val = radcheb_eval(r,breaks,coefs), with breaks and coefs as returned
%  by radcheb_fit. Radii outside [breaks(1), breaks(end)] return NaN.

n = size(coefs,2);

sz = size(r);
r  = r(:);

ib   = discretize(r,breaks(:));
ibad = isnan(ib);
if any(ibad)
    ib(ibad) = 1;
end

x   = r - breaks(ib);
val = coefs(ib,1);
for i = 2:n
    val = coefs(ib,i) + x.*val;
end

if any(ibad)
    val(ibad) = NaN;
end

val = reshape(val,sz);

end
