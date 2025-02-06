%
%  This file illustrates how to decouple quadrature generation
%  from its use with an FMM using the normal derivative of 
%  a single layer potential for the Helmholtz equation
%
%  NOTE: this file is dependent on an independent compilation
%  of the MATLAB wrappers for FMM3D and is compatible with
%  the version of FMM3D that comes included as a submodule with
%  this package. 
%
r = 1;
na = 4;

zk = 1.1;

tol = 1e-5;
S = geometries.sphere(r, na);

Q = helm3d.neumann.get_quadrature_correction(S, tol, zk, 0.0);
%%
Q2 = Q;
Q2.wnear = Q.wnear(1,:).';
Q2.spmat = conv_rsc_to_spmat(S, Q2.row_ptr, Q2.col_ind, Q2.wnear);
%%
jn = @(n,z) sqrt(pi/2/z)*besselj(n+0.5,z);
hn = @(n,z) sqrt(pi/2/z)*besselh(n+0.5,1,z);

jnp = @(n,z) 0.5*(jn(n-1,z) - (jn(n,z) + z*jn(n+1,z))/z);
hnp = @(n,z) 0.5*(hn(n-1,z) - (hn(n,z) + z*hn(n+1,z))/z);

zfac = 1j*zk*zk*(jn(1,zk)*hnp(1,zk) + jnp(1,zk)*hn(1,zk))/2;
%%
rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
sig = S.r(3,:).'./rr.';

u = eval_sprime(S, zk, Q2, sig, tol);

%%
uex = sig.*zfac;
err = norm((u-uex).*sqrt(S.wts(:)))/norm(uex.*sqrt(S.wts(:)));
fprintf('error in sprime=%d\n',err);

function [u] = eval_sprime(S, zk, Q, sig, tol)
    srcuse = [];
    srcuse.sources = S.r(:,:);
    srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
    pg = 2;
    U = hfmm3d(tol, zk, srcuse, pg);

    u = U.grad(1,:).*S.n(1,:) + U.grad(2,:).*S.n(2,:) + ...
        U.grad(3,:).*S.n(3,:);
    u = u(:);

    ixyzs = S.ixyzs(:);
    npatches = S.npatches;

    istarts = Q.row_ptr(1:end-1);
    iends = Q.row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    for i=1:S.npts
        iinds = horzcat(isrcinds{Q.col_ind(istarts(i):iends(i))});
        srcinfo = [];
        srcinfo.r = S.r(:,iinds);
        srcinfo.n = S.n(:,iinds);
        targinfo = [];
        targinfo.r = S.r(:,i);
        targinfo.n = S.n(:,i);
        submat = helm3d.kern(zk, srcinfo, targinfo, 'sprime');
        submat(isnan(submat)) = 0;
        u(i) = u(i) - submat*(srcuse.charges(iinds).')*4*pi;
    end

    u = u + Q.spmat*sig;

end
