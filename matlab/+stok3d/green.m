function [val,grad,hess] = green(src,targ)
%STOKES.GREEN evaluate the Stokes green's function
% for the given sources and targets
%
% Syntax: [val, grad, hess] = stok3d.green(src, targ)
%
% Input:
%   src  - (3, ns) source locations
%   targ - (3, nt) target locations
%
% Output:
%   val  - (3, nt, 3, ns)  G_{ij}(x,y),  Stokeslet kernel
%   grad - (3, 3, nt, 3, ns)  grad(k,i,it,j,is) = dG_{ij}/dx_k
%   hess - (3, 3, nt, 3, 3, ns)  hess(k,l,it,i,j,is) = d^2 G_{ij}/dx_k dx_l

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);
zs = repmat(src(3,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);
zt = repmat(targ(3,:).',1,ns);

[m,n] = size(xs);


rx = reshape(xt-xs,[1,m,1,n]);
ry = reshape(yt-ys,[1,m,1,n]);
rz = reshape(zt-zs,[1,m,1,n]);

rx2 = rx.*rx;
ry2 = ry.*ry;
rz2 = rz.*rz;

r2 = rx2+ry2+rz2;

r = sqrt(r2);

gmat = zeros(3,m,3,n);
gmat(1,:,1,:) = rx2;
gmat(2,:,2,:) = ry2;
gmat(3,:,3,:) = rz2;

gmat(2,:,1,:) = rx.*ry;
gmat(3,:,1,:) = rx.*rz;
gmat(3,:,2,:) = ry.*rz;


gmat(1,:,2,:) = gmat(2,:,1,:);
gmat(1,:,3,:) = gmat(3,:,1,:);
gmat(2,:,3,:) = gmat(3,:,2,:);

fact = 1/8/pi;

if nargout > 0
    val = gmat;
    val(1,:,1,:) = val(1,:,1,:) + r2;
    val(2,:,2,:) = val(2,:,2,:) + r2;
    val(3,:,3,:) = val(3,:,3,:) + r2;
    val = fact*val./(r.^3);
end


gfact = 3/4/pi;
if nargout > 1
    grad = zeros(3,3,m,3,n);
    gmatg = gfact*reshape(gmat./(r.^5),[1,3,m,3,n]);
    grad(1,:,:,:,:) = gmatg.*reshape(rx,[1,1,m,1,n]);
    grad(2,:,:,:,:) = gmatg.*reshape(ry,[1,1,m,1,n]);
    grad(3,:,:,:,:) = gmatg.*reshape(rz,[1,1,m,1,n]);
end

if nargout > 2
    % hess(k,l,it,i,j,is) = d^2 G_{ij} / dx_k dx_l
    %
    % Full formula (r = targ - src):
    %   H_{kl,ij} = (1/8pi) [
    %     -delta_{ij}(delta_{kl}/r^3 - 3 r_k r_l/r^5)
    %     + delta_{ik}(delta_{jl}/r^3 - 3 r_j r_l/r^5)
    %     + delta_{jk}(delta_{il}/r^3 - 3 r_i r_l/r^5)
    %     - 3(delta_{il} r_j r_k + delta_{jl} r_i r_k + delta_{kl} r_i r_j)/r^5
    %     + 15 r_i r_j r_k r_l / r^7
    %   ]
    %
    % Ordering: off-diagonal kl pairs first, then diagonal, matching
    % helm3d.green convention for hess entries.
    %
    % Note: rx,ry,rz here have shape (1,m,1,n) = (1,nt,1,ns).
    % We work with flat (m,n) = (nt,ns) arrays for simplicity, then reshape.

    rflat  = reshape(r,  [m,n]);
    rxf    = reshape(rx, [m,n]);
    ryf    = reshape(ry, [m,n]);
    rzf    = reshape(rz, [m,n]);
    rx2f   = rxf.*rxf;
    ry2f   = ryf.*ryf;
    rz2f   = rzf.*rzf;

    ir3h = 1./rflat.^3;   % (nt,ns)
    ir5h = 1./rflat.^5;
    ir7h = 1./rflat.^7;

    % r-components indexed by spatial direction (cell array for convenience)
    rc = {rxf, ryf, rzf};

    hess = zeros(3,3,m,3,3,n);

    for kk = 1:3
      for ll = 1:3
        for ii = 1:3
          for jj = 1:3
            val_h = 0;
            % Term 1: -delta_{ij}(delta_{kl}/r^3 - 3 r_k r_l/r^5)
            if ii==jj
              val_h = val_h - ((kk==ll)*ir3h - 3*rc{kk}.*rc{ll}.*ir5h);
            end
            % Term 2: delta_{ik}(delta_{jl}/r^3 - 3 r_j r_l/r^5)
            if ii==kk
              val_h = val_h + ((jj==ll)*ir3h - 3*rc{jj}.*rc{ll}.*ir5h);
            end
            % Term 3: delta_{jk}(delta_{il}/r^3 - 3 r_i r_l/r^5)
            if jj==kk
              val_h = val_h + ((ii==ll)*ir3h - 3*rc{ii}.*rc{ll}.*ir5h);
            end
            % Term 4: -3(delta_{il} r_j r_k + delta_{jl} r_i r_k + delta_{kl} r_i r_j)/r^5
            val_h = val_h - 3*((ii==ll)*rc{jj}.*rc{kk} + ...
                                (jj==ll)*rc{ii}.*rc{kk} + ...
                                (kk==ll)*rc{ii}.*rc{jj}).*ir5h;
            % Term 5: 15 r_i r_j r_k r_l / r^7
            val_h = val_h + 15*rc{ii}.*rc{jj}.*rc{kk}.*rc{ll}.*ir7h;

            hess(kk,ll,:,ii,jj,:) = fact*reshape(val_h,[1,m,1,1,1,n]);
          end
        end
      end
    end
end
