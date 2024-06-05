function [p,dersu,dersv] = cheb_pols(norder,uv)
   [~,m] = size(uv);
   chebpols_u = zeros(m,norder+1);
   chebpols_v = zeros(m,norder+1);

   chebders_u = zeros(m,norder+1);
   chebders_v = zeros(m,norder+1);

   u = uv(1,:);
   v = uv(2,:);
   u = u(:);
   v = v(:);
   chebpols_u(:,1+0) = 1;
   chebpols_v(:,1+0) = 1;
   
   chebders_u(:,1+0) = 0;
   chebders_v(:,1+0) = 0;
   
   chebpols_u(:,1+1) = u;
   chebpols_v(:,1+1) = v;
   
   chebders_u(:,1+1) = 1;
   chebders_v(:,1+1) = 1;
   
   pkp1_u = chebpols_u(:,2);
   pk_u = chebpols_u(:,1);
   
   pkdp1_u = chebders_u(:,2);
   pkd_u = chebders_u(:,1);
   
   pkp1_v = chebpols_v(:,2);
   pk_v = chebpols_v(:,1);
   
   pkdp1_v = chebders_v(:,2);
   pkd_v = chebders_v(:,1);
   
   for k=1:norder-1
       pkm1_u = pk_u;
       pkm1_v = pk_v;

       pkdm1_u = pkd_u;
       pkdm1_v = pkd_v;
       
       pk_u = pkp1_u;
       pk_v = pkp1_v;

       pkd_u = pkdp1_u;
       pkd_v = pkdp1_v;
       
       pkp1_u = (2*u.*pk_u -pkm1_u);
       pkp1_v = (2*v.*pk_v -pkm1_v);
       
       pkdp1_u = 2*(pk_u + u.*pkd_u) -pkdm1_u;
       pkdp1_v = 2*(pk_v + v.*pkd_v) -pkdm1_v;
       
       chebpols_u(:,k+2) = pkp1_u;
       chebpols_v(:,k+2) = pkp1_v;

       chebders_u(:,k+2) = pkdp1_u;
       chebders_v(:,k+2) = pkdp1_v;

   end
   
   npols = (norder+1)*(norder+1);
   p = zeros(m,npols);
   dersu = zeros(m,npols);
   dersv = zeros(m,npols);
   iii = 0;
   for i=0:norder
       for j=0:norder
           iii = iii + 1;
           p(:,iii) = chebpols_u(:,j+1).*chebpols_v(:,i+1);
           dersu(:,iii) = chebders_u(:,j+1).*chebpols_v(:,i+1);
           dersv(:,iii) = chebpols_u(:,j+1).*chebders_v(:,i+1);
       end
   end
   
   p = reshape(p.',[npols,m]);
   dersu = reshape(dersu.',[npols,m]);
   dersv = reshape(dersv.',[npols,m]);
end
