function [p, dersu, dersv] = lege_ders(norder,uv)
   [~,m] = size(uv);

   legpols_u = zeros(m,norder+1);
   legpols_v = zeros(m,norder+1);

   legders_u = zeros(m,norder+1);
   legders_v = zeros(m,norder+1);

   u = uv(1,:);
   v = uv(2,:);
   u = u(:);
   v = v(:);
   legpols_u(:,1+0) = 1;
   legpols_v(:,1+0) = 1;

   legders_u(:,1+0) = 0;
   legders_v(:,1+0) = 0;
   
   legpols_u(:,1+1) = u;
   legpols_v(:,1+1) = v;
   
   legders_u(:,1+1) = 1;
   legders_v(:,1+1) = 1;

   pkp1_u = legpols_u(:,2);
   pk_u = legpols_u(:,1);
   
   pkp1_v = legpols_v(:,2);
   pk_v = legpols_v(:,1);


   pkdp1_u = legders_u(:,2);
   pkd_u = legders_u(:,1);

   pkdp1_v = legders_v(:,2);
   pkd_v = legders_v(:,1);
   
   for k=1:norder-1
       pkm1_u = pk_u;
       pkm1_v = pk_v;

       pkdm1_u = pkd_u;
       pkdm1_v = pkd_v;
       
       pk_u = pkp1_u;
       pk_v = pkp1_v;

       pkd_u = pkdp1_u;
       pkd_v = pkdp1_v;
       
       pkp1_u = ((2*k+1)*u.*pk_u -k*pkm1_u)/(k+1);
       pkp1_v = ((2*k+1)*v.*pk_v -k*pkm1_v)/(k+1);

       pkdp1_u = ((2*k+1)*(pk_u + u.*pkd_u) -k*pkdm1_u)/(k+1);
       pkdp1_v = ((2*k+1)*(pk_v + v.*pkd_v) -k*pkdm1_v)/(k+1);

       
       legpols_u(:,k+2) = pkp1_u;
       legpols_v(:,k+2) = pkp1_v; 

       legders_u(:,k+2) = pkdp1_u;
       legders_v(:,k+2) = pkdp1_v;
   end
   
   npols = (norder+1)*(norder+1);
   p = zeros(m,npols);
   dersu = zeros(m,npols);
   dersv = zeros(m,npols);
   iii = 0;
   for i=0:norder
       for j=0:norder
           iii = iii + 1;
           p(:,iii) = legpols_u(:,j+1).*legpols_v(:,i+1);
           dersu(:,iii) = legders_u(:,j+1).*legpols_v(:,i+1);
           dersv(:,iii) = legpols_u(:,j+1).*legders_v(:,i+1);
       end
   end
   
   p = reshape(p.',[npols,m]);
   dersu = reshape(dersu.',[npols,m]);
   dersv = reshape(dersv.',[npols,m]);
end
