S = geometries.disk([],[],[4 4 4],6);

% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

V = eval_gauss(S.r);
%rhs = eval_gauss(xs);
zk = pi;

%% v2v


% A = helm2d.slp_matgen(S,zk,1e-9);
% lhs_11 = -eye(S.npts) + V.*A;

tic;
[v2v_cor,nover] = helm2d.get_quad_cor_sub(S,zk, 1e-8);
toc;

v2v_apply = @(mu) helm2d.apply_v2v(S,zk,mu,v2v_cor,nover,1e-14);
lhs_11 = @(mu) -mu + V.*v2v_apply(mu);


%% b2v


% lhs_12 = V.*chunkerkernevalmat(chnkr,h2d_s,S.r(1:2,:));

h2d_s = kernel('h','s',zk);

opts = []; opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr,h2d_s,S.r(1:2,:),opts);

apply_b2v = @(mu) helm2d.apply_b2v_neu(S,zk,chnkr,mu,Ab2v_cor,1e-12);
lhs_12 = @(mu) V.*apply_b2v(mu);


%% v2b


targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
% lhs_21 = helm2d.v2b_neu(S,zk,targinfo,1e-8);

tic;
[Av2b_cor,nover] = helm2d.get_quad_cor_v2b_neu(S, zk,targinfo, 1e-8);
toc;

v2b_apply = @(mu) helm2d.apply_v2b_neu(S,zk,targinfo,mu,Av2b_cor,nover,1e-14);
lhs_21 = @(mu) v2b_apply(mu);


%% b2b


h2d_sp = kernel.helm2d('sp',zk);
lhs_22 = 0.5*eye(chnkr.npt)+chunkermat(chnkr,h2d_sp); %0.5*eye(n)... -

%%

% lhs = [lhs_11, lhs_12; lhs_21, lhs_22];
lhs = @(x) [lhs_11(x(1:S.npts)) + lhs_12(x(S.npts+1:end)); lhs_21(x(1:S.npts)) + lhs_22*x(S.npts+1:end)];


% analytic solution : u = sin(x)sin(y)
rhs_1 = (-2+zk^2+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_2 = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; % sin(xs(1,:))


rhs = [rhs_1 rhs_2].';

% dens = lhs\rhs;
dens = gmres(lhs,rhs,[],1e-12,2000);
mu = dens(1:S.npts);
rho = dens(S.npts+1:end);
% fkern_d = kernel('l','d');
% u = A*mu + chunkerkerneval(chnkr,fkern_d,rho,S.r(1:2,:));
u = v2v_apply(mu) + chunkerkerneval(chnkr,h2d_s,rho,S.r(1:2,:));

ref_u = sin(S.r(1,:)).*sin(S.r(2,:));
err = abs(u - ref_u(:)) / max(abs(u));

figure(1); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); colorbar

fprintf('relative L^2/L^2 error: %5.5e\n', vecnorm(err .* S.wts(:)) / vecnorm(u .* S.wts(:)) );

return


%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end