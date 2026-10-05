function [ipars,dpars,info] = radcheb_fit(f,rmin,rmax,opts)
%RADCHEB_FIT build an adaptive piecewise Chebyshev interpolant of f(r)
%
%  Syntax:
%    [ipars,dpars,info] = radcheb_fit(f,rmin,rmax)
%    [ipars,dpars,info] = radcheb_fit(f,rmin,rmax,opts)
%
%  f must be vectorized in r and may be real or complex. The interpolant
%  is returned in the vpp format of src/kernels/vpp.f90.
%
%  Input arguments:
%    * f: function handle of one variable
%    * rmin, rmax: interval on which to build the interpolant, rmin > 0
%    * opts: options struct
%        opts.norder: terms per panel (default 6)
%        opts.eps: tolerance (default 1e-12)
%        opts.nlevmax: max subdivision levels (default 50)
%        opts.maxsub: max number of panels (default 10000)
%
%  Output arguments:
%    * ipars: (5,1) [nbin, norder, ibins, icfs, nv]
%    * dpars: panel endpoints followed by scaled monomial coefficients
%    * info: struct with fields nbin, nlev, err, breaks, coefs
%        breaks (nbin+1,1) and coefs (nbin,norder) are the same data in
%        the form used by radcheb_eval
%
%  See also RADCHEB_EVAL, KERNEL3D.RADCHEB
%

if nargin < 4 || isempty(opts), opts = struct(); end

norder  = 6;
tol     = 1e-12;
nlevmax = 50;
maxsub  = 20000;
if isfield(opts,'norder'),  norder  = opts.norder;  end
if isfield(opts,'eps'),     tol     = opts.eps;     end
if isfield(opts,'nlevmax'), nlevmax = opts.nlevmax; end
if isfield(opts,'maxsub'),  maxsub  = opts.maxsub;  end

if ~(rmin > 0) || ~(rmax > rmin)
    error('RADCHEB_FIT: need 0 < rmin < rmax.');
end

n = norder;

% first kind Chebyshev nodes
x = cos((2*(1:n).'-1)*pi/(2*n));

% the same nodes on each half of [-1,1]
x2 = [-1 + (x+1)/2; (x+1)/2];

amat = radcheb_interpmat(x,x2);

a = rmin;
b = rmax;
fv = f(a + (x+1)*(b-a)/2);
fv = fv(:);

ad = []; bd = []; fd = zeros(n,0);   % resolved panels
err = 0;

for ilev = 1:nlevmax

    na = numel(a);
    t  = a.' + (x2+1).*(b-a).'/2;
    fo = f(t(:));
    fo = reshape(fo,2*n,na);
    fi = amat*fv;

    d    = fo - fi;
    derr = max(max(abs(real(d)),[],1), max(abs(imag(d)),[],1));
    drel = max(max(abs(real(fo)),[],1), max(abs(imag(fo)),[],1));
    e    = derr./max(drel,1);
    ikeep = e(:) <= tol;
    if any(ikeep), err = max(err,max(e(ikeep))); end

    ad = [ad; a(ikeep)];       %#ok<AGROW>
    bd = [bd; b(ikeep)];       %#ok<AGROW>
    fd = [fd, fv(:,ikeep)];    %#ok<AGROW>

    ibad = find(~ikeep);
    if isempty(ibad), break; end

    if numel(ad) + 2*numel(ibad) > maxsub
        error('RADCHEB_FIT: maxsub exceeded, %d panels.', ...
              numel(ad) + 2*numel(ibad));
    end

    amid = (a(ibad) + b(ibad))/2;
    fv   = [fo(1:n,ibad), fo(n+1:2*n,ibad)];
    a    = [a(ibad); amid];
    b    = [amid; b(ibad)];

end

if ~isempty(ibad)
    error('RADCHEB_FIT: nlevmax exceeded, %d panels unresolved.', ...
          numel(ibad));
end

[ad,isort] = sort(ad);
bd = bd(isort);
fd = fd(:,isort);

nbin   = numel(ad);
breaks = [ad; bd(end)];

% monomial coefficients in r-a on each panel, highest degree first
tt = (x+1)/2;
vm = tt.^(n-1:-1:0);
cf = (vm\fd).';
cf = cf.*((1./(bd-ad)).^(n-1:-1:0));

ibins = 1;
icfs  = nbin + 2;
if isreal(cf)
    nv = 1;
    cfv = cf.';
else
    nv = 2;
    cfv = zeros(2,n,nbin);
    cfv(1,:,:) = reshape(real(cf).',[1,n,nbin]);
    cfv(2,:,:) = reshape(imag(cf).',[1,n,nbin]);
end

ipars = [nbin; n; ibins; icfs; nv];
dpars = [breaks(:); cfv(:)];

info = struct('nbin',nbin,'nlev',ilev,'err',err, ...
              'breaks',breaks,'coefs',cf);

end


function amat = radcheb_interpmat(x,x2)
%  barycentric interpolation matrix from the nodes x to the nodes x2

n  = numel(x);
bw = ((-1).^((1:n).'-1)).*sin((2*(1:n).'-1)*pi/(2*n));

dd = x2 - x.';
w  = bw.'./dd;
amat = w./sum(w,2);

[i0,j0] = find(dd == 0);
for k = 1:numel(i0)
    amat(i0(k),:)     = 0;
    amat(i0(k),j0(k)) = 1;
end

end
