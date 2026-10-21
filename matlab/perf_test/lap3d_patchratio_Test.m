%
% This file tests the Laplace single layer
% potenial on the sphere
%
%
run ../startup.m

errs = zeros(4,1);

for ii = 1:4

% For triangles
S = geometries.stellarator(ii*[5,2],8,1);

figure(1); clf
plot(S,rand(S.npatches,1),'FaceAlpha',0.8)
drawnow

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
rr2 = sum(wts);
rr3 = norm(S.r-S.n);

% tic, plot(S); toc;

dpars = [1.0, 0.0];
ndeg = 1;


rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

rhs = rhs(:);
eps = 1e-9;


% p = lap3d.dirichlet.eval(S,rhs,S,eps,dpars);
% 
% p_ex = rhs/(2*ndeg+1);
% 
% err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));
% 
% fprintf('Error in single layer potential=%d\n',err1);
% 
% %% Now test eval routine with precomputed quadrature corrections
% opts_quad = [];
% opts_quad.format='rsc';
% Q = lap3d.dirichlet.get_quadrature_correction(S, ...
%    eps,dpars,S,opts_quad);
% opts_use = [];
% opts_use.precomp_quadrature = Q;
% p = lap3d.dirichlet.eval(S,rhs,S,eps,dpars,opts_use);
% err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));
% 
% fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [4.6;0;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = lap3d.kern(src_info,S,'s');


dpars = [1,1];
sig = lap3d.dirichlet.solver(S,rhs,eps,dpars);

targ_info = [];
targ_info.r = xyz_out;

dat = lap3d.kern(S,targ_info,'c',dpars(1),dpars(2));

pot = dat*(sig.*wts);
pot_ex = lap3d.kern(src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));

errs(ii) = abs(pot-pot_ex)/abs(pot_ex);

end


%%

plot(log10(1:4),log10(errs))
hold on
plot(log10((1:4)),-3.5+log10((1:4))*-8)

