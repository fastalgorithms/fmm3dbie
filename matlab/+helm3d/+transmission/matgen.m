function A = matgen(iin, jin, S, zpars,P, Q,opts)

% Convert system indices into original indices
ipts = idivide(int64(iin(:)-1),int64(2))+1;
jpts = idivide(int64(jin(:)-1),int64(2))+1;

[i,~,iiuni] = unique(ipts);
[j,~,ijuni] = unique(jpts);

if isempty(i) || isempty(j)
    A = zeros(2*length(i), 2*length(j));
    % Extract relevant rows and columns in A_uni
    iiuni2 = (iiuni-1)*2 + mod(iin(:)-1, 2)+1;
    ijuni2 = (ijuni-1)*2 + mod(jin(:)-1, 2)+1;

    A = A(iiuni2, ijuni2);
    return
end

[I, J] = ndgrid(i, j);
Idiag = find(I==J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if(nargin <=6)
       opts = [];
     end
     
     ifdiag = 1;
     if(isfield(opts,'ifdiag'))
        ifdiag = opts.ifdiag;
     end


     l2scale = 0;
     if(isfield(opts,'l2scale'))
        l2scale = opts.l2scale;
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assemble (block) system matrix for transmission problem
A = zeros(2*length(i), 2*length(j));

% Compute block operators for transmission problem
area = S.wts.';


z_k0 = zpars(1)*sqrt(zpars(2)*zpars(3));
z_k=   zpars(1)*sqrt(zpars(4)*zpars(5));
deps0  = zpars(2);
deps   = zpars(4);

srcinfo = [];
srcinfo.r = S.r(:,j);
srcinfo.n = S.n(:,j);
targinfo = [];
targinfo.r = S.r(:,i);
targinfo.n = S.n(:,i);



%K_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k, nx), area(j));
%K = bsxfun(@times, helm_dirichlet_kernel(targets, sources, zdtmp, ny), area(j));
%Sl = bsxfun(@times, helm_dirichlet_kernel(targets, sources, zstmp, ny), area(j));
%T = bsxfun(@times, helm_hypersingular_kernel(targets, sources, z_k, nx, ny), area(j));
K_prime = bsxfun(@times,helm3d.kern(z_k,srcinfo,targinfo,'sprime'),area(j));
K = bsxfun(@times,helm3d.kern(z_k,srcinfo,targinfo,'d'),area(j));
Sl= bsxfun(@times,helm3d.kern(z_k,srcinfo,targinfo,'s'),area(j));
T = bsxfun(@times,helm3d.kern(z_k,srcinfo,targinfo,'dprime'),area(j));
K_prime(Idiag) = 0;
K(Idiag)       = 0;
Sl(Idiag)      = 0;
T(Idiag)       = 0;

%K0_prime = bsxfun(@times, helm_sound_hard_kernel(targets, sources, z_k0, nx), area(j));
%K0 = bsxfun(@times, helm_dirichlet_kernel(targets, sources, z0dtmp, ny), area(j));
%Sl0 = bsxfun(@times, helm_dirichlet_kernel(targets, sources, z0stmp, ny), area(j));
%T0 = bsxfun(@times, helm_hypersingular_kernel(targets, sources, z_k0, nx, ny), area(j));
K0_prime = bsxfun(@times,helm3d.kern(z_k0,srcinfo,targinfo,'sprime'),area(j));
K0 = bsxfun(@times,helm3d.kern(z_k0,srcinfo,targinfo,'d'),area(j));
Sl0= bsxfun(@times,helm3d.kern(z_k0,srcinfo,targinfo,'s'),area(j));
T0 = bsxfun(@times,helm3d.kern(z_k0,srcinfo,targinfo,'dprime'),area(j));
K0_prime(Idiag) = 0;
K0(Idiag)       = 0;
Sl0(Idiag)      = 0;
T0(Idiag)       = 0;

% A11 
A(1:2:end, 1:2:end) = -(deps*K-deps0*K0);
% A12
A(1:2:end, 2:2:end) = -(deps^2*Sl-deps0^2*Sl0);
% A21
A(2:2:end, 1:2:end) = T-T0;
% A22
A(2:2:end, 2:2:end) = deps*K_prime-deps0*K0_prime;


if(nargin > 4)
    % Compute quadrature corrections
    M1 = spget_quadcorr(i, j, P, Q{2});
    idx1 = abs(M1) ~= 0;

    M2 = spget_quadcorr(i, j, P, Q{1});
    idx2 = abs(M2) ~= 0;

    M3 = -spget_quadcorr(i, j, P, Q{4});
    idx3 = abs(M3) ~= 0;

    M4 = -spget_quadcorr(i, j, P, Q{3});
    idx4 = abs(M4) ~= 0;


    % Add quadrature corrections to each block
    tmp = A(1:2:end, 1:2:end);
    tmp(idx1) = M1(idx1);
    A(1:2:end, 1:2:end) = tmp;

    tmp = A(1:2:end, 2:2:end);
    tmp(idx2) = M2(idx2);
    A(1:2:end, 2:2:end) = tmp;

    tmp = A(2:2:end, 1:2:end);
    tmp(idx3) = M3(idx3);
    A(2:2:end, 1:2:end) = tmp;

    tmp = A(2:2:end, 2:2:end);
    tmp(idx4) = M4(idx4);
    A(2:2:end, 2:2:end) = tmp;
end


if(ifdiag)
    % Add identity parts to A11 and A22 blocks
    fac = (deps+deps0)/2;

    tmp = A(1:2:end, 1:2:end);
    tmp(I == J) = tmp(I == J) + fac;
    A(1:2:end, 1:2:end) = tmp;

    tmp = A(2:2:end, 2:2:end);
    tmp(I == J) = tmp(I == J) + fac;
    A(2:2:end, 2:2:end) = tmp;
end


if (l2scale)
    % Compute weights for each block
    tmp = A(1:2:end, 1:2:end);
    tmp = bsxfun(@times, sqrt(area(i)).',tmp);
    A(1:2:end, 1:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

    tmp = A(1:2:end, 2:2:end);
    tmp = bsxfun(@times, sqrt(area(i)).',tmp);
    A(1:2:end, 2:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

    tmp = A(2:2:end, 1:2:end);
    tmp = bsxfun(@times, sqrt(area(i)).',tmp);
    A(2:2:end, 1:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));

    tmp = A(2:2:end, 2:2:end);
    tmp = bsxfun(@times, sqrt(area(i)).',tmp);
    A(2:2:end, 2:2:end) = bsxfun(@times, tmp,1.0./sqrt(area(j)));
end


% Extract relevant rows and columns in A_uni
iiuni2 = (iiuni-1)*2 + mod(iin(:)-1, 2)+1;
ijuni2 = (ijuni-1)*2 + mod(jin(:)-1, 2)+1;

A = A(iiuni2, ijuni2);

end
