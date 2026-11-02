function U = fmm(eps, srcuse, pg, targ, pgt)
%LAP2D.FMM   FMM for lap2d.green with the interface of RFMM2D, for
% complex charges.
%
% Syntax: U = lap2d.fmm(eps, srcuse, pg, targ, pgt)
%
% see also FMM2DTO3D

nd = 1;
if ( isfield(srcuse, 'nd') ), nd = srcuse.nd; end
ns = size(srcuse.sources, 2);
nt = size(targ, 2);
q = reshape(srcuse.charges, nd, ns);

src = [];
src.sources = srcuse.sources;
if ( isreal(q) )
    src.nd = nd;
    src.charges = -q/(2*pi);
else
    src.nd = 2*nd;
    src.charges = -[real(q); imag(q)]/(2*pi);
end
V = rfmm2d(eps, src, pg, targ, pgt);

U = [];
U.pottarg = combine(V.pottarg, nd, isreal(q), [1 nt]);
if ( pgt > 1 ), U.gradtarg = combine(V.gradtarg, nd, isreal(q), [2 nt]); end
if ( pgt > 2 ), U.hesstarg = combine(V.hesstarg, nd, isreal(q), [3 nt]); end

end


function W = combine(V, nd, isr, sz)
if ( isr )
    W = reshape(V, [nd sz]);
else
    V = reshape(V, [2*nd sz]);
    W = V(1:nd,:,:) + 1i*V(nd+1:end,:,:);
end
if ( sz(1) == 1 ), W = reshape(W, nd, sz(2)); end
end
