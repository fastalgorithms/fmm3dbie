%
% This file tests the convergence of Laplace single layer
% potenial on the sphere
%
%

% For triangles

norders = 1:16;
nas = [1 2 3 4 5];

convmat = zeros(length(norders),length(nas));
patch_diams = zeros(length(nas),1);

for ii = 1:length(norders)

    norders(ii)

    for jj = 1:length(nas) 

        S = geometries.sphere(1, nas(jj), [0;0;0], norders(ii), 1);

        patch_diams(jj) = mean(S.rads);
        
        eps = 1e-13;
        
        tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
        
        xyz_in = [0.3;0.5;0.1];
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

        convmat(ii,jj) = abs(pot-pot_ex)/abs(pot_ex);

    end
    
end

%%

figure(1);clf
plot(log10(patch_diams(1:end-1)),log10(convmat(:,1:end-1)),'x-')
legend("n = " + norders)
