function val = radpp_eval(r,breaks,coefs)
%RADPP_EVAL evaluate a radpp_fit expansion
%
%  val = radpp_eval(r,breaks,coefs), with breaks (nbin+1,1) and coefs
%  (nbin,n,m) from the info output of RADPP_FIT. r is any array; val is
%  (numel(r),m). Radii outside [breaks(1),breaks(end)] give NaN.
%
%  See also RADPP_FIT, RADPP_ROTATE

[~,n,m] = size(coefs);
r = r(:);

ib   = discretize(r,breaks);
ibad = isnan(ib);
ib(ibad) = 1;

x = r - breaks(ib);
val = reshape(coefs(ib,1,:),numel(r),m);
for i = 2:n
    val = reshape(coefs(ib,i,:),numel(r),m) + x.*val;
end
val(ibad,:) = NaN;

end
