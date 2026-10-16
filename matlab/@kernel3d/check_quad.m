function [err, Q, Qs] = check_quad(K, eps, opts)
%KERNEL3D.CHECK_QUAD check the near quadrature of a kernel, e.g. one
%   built with kernel3d(@green, type, ...), against a smooth integral.
%
%   [err, Q, Qs] = kernel3d.check_quad(K)
%   [err, Q, Qs] = kernel3d.check_quad(K, eps, opts)
%
%   Builds one curved quad patch z = (a x^2 + b y^2)/2, |x|,|y| <= h and a
%   target at opts.dist*h (default 1) in the normal direction, with fields
%   n and mean_curv there. Q = K.getquad on the patch is compared with the
%   integral of K on the same patch at order opts.nover (default 30)
%   interpolated to the nodes.
% 
%   Reports the err = |Q - Qs|/|Qs|.
%
%   opts: norder (default 8), h (0.5), a, b (1, 3), dist, nover.

if nargin < 2 || isempty(eps), eps = 1e-7; end
if nargin < 3, opts = struct(); end
norder = getopt(opts, 'norder', 8);
nover  = getopt(opts, 'nover', 30);
h      = getopt(opts, 'h', 0.5);
a      = getopt(opts, 'a', 1);
b      = getopt(opts, 'b', 3);
dist   = getopt(opts, 'dist', 1);

c1 = cos(0.7);  s1 = sin(0.7);  c2 = cos(1.1);  s2 = sin(1.1);
R = [1 0 0; 0 c1 -s1; 0 s1 c1]*[c2 0 s2; 0 1 0; -s2 0 c2];

S  = patch(norder, h, a, b, R);
So = patch(nover, h, a, b, R);

targ.r = dist*h*R(:,3);
targ.n = R(:,3);
targ.mean_curv = -(a + b)/2;            % surfer convention, at the centre
targ.patch_id = -1;
targ.uvs_targ = [0; 0];

m = K.opdims(1);  n = K.opdims(2);
np = S.npts;

Q = full(K.getquad(S, eps, targ));                 % (m, n*np)

% smooth integral at order nover, interpolated to the patch nodes
x  = polytens.lege.nodes(norder);
xo = polytens.lege.nodes(nover);
I  = polytens.lege.pols(norder, xo).'*polytens.lege.vals2coefs(norder, x);
[~, ~, ~, ~, ~, wo] = extract_arrays(So);
src.r = So.r;  src.n = So.n;
A = K.eval(src, targ);                             % (m, n*no)
Qs = zeros(m, n*np);
for k = 1:n
    Qs(:, k:n:end) = (A(:, k:n:end).*wo(:).')*I;
end

err = norm(Q - Qs, 'fro')/norm(Qs, 'fro');
fprintf('%s %s: dist %g h, getquad vs smooth %.2e\n', K.name, K.type, ...
        dist, err);
end


function S = patch(norder, h, a, b, R)
% one tensor Legendre patch (iptype 11) of z = (a x^2 + b y^2)/2
uv = polytens.lege.nodes(norder);
x = h*uv(1,:);  y = h*uv(2,:);
r  = [x; y; (a*x.^2 + b*y.^2)/2];
du = h*[ones(size(x)); zeros(size(x)); a*x];
dv = h*[zeros(size(x)); ones(size(x)); b*y];
nn = cross(du, dv);  nn = nn./vecnorm(nn);
S = surfer(1, norder, [R*r; R*du; R*dv; R*nn], 11);
end


function v = getopt(opts, name, v)
if isfield(opts, name), v = opts.(name); end
end
