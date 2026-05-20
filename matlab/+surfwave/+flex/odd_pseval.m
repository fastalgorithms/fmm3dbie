function [v] = odd_pseval(cf,x)
% evaluates an odd power series starting at n = 1 
% with coefficients given in cf. 
%
% x can be an array

x2 = x.*x;

xp = x;

v = zeros(size(x));

nterms = length(cf);
for j = 1:nterms
    v = v + cf(j)*xp;
    xp = xp.*x2;
end
