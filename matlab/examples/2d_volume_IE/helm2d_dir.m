% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')

% run('/Users/yuguan/software/chunkie/startup.m')
% run('/Users/yuguan/Dropbox/fmm3dbie/matlab/startup.m')
% addpath '/Users/yuguan/Dropbox/fmm3dbie/src'

S = geometries.disk([],[],[4 4 4],8);

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

V = 0*eval_gauss(S.r);
%rhs = eval_gauss(xs);
zk = pi;
zk = sqrt(2);

%% v2v


A = helm2d.slp_matgen(S,zk,1e-9);
lhs_11 = -eye(S.npts) + V.*A;



%% b2v


fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'d');
lhs_12 = V.*chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));



%% v2b


targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
lhs_21 = helm2d.v2b_dir(zk,S,targinfo,1e-8);



%% b2b


h2d_d = kernel.helm2d('d',zk);lhs_22 = -0.5*eye(chnkr.npt)+chunkermat(chnkr,h2d_d); %0.5*eye(n)... -

%%

lhs = [lhs_11, lhs_12; lhs_21, lhs_22];
% lhs = @(x) [lhs_11(x(1:S.npts)) + lhs_12(x(S.npts+1:end)); lhs_21(x(1:S.npts)) + lhs_22*x(S.npts+1:end)];


% analytic solution : u = sin(x)sin(y)
rhs_1 = (-2+zk^2+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_2 = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:)); % sin(xs(1,:))

rhs = [rhs_1 rhs_2].';

dens = lhs\rhs;
mu = dens(1:S.npts);
rho = dens(S.npts+1:end);

% fkern_d = kernel('l','d');
% u = A*mu + chunkerkerneval(chnkr,fkern_d,rho,S.r(1:2,:));
u = A*mu + chunkerkerneval(chnkr,h2d_d,rho,S.r(1:2,:));

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