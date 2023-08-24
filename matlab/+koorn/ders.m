function [pols,ders]=ders(nmax,uv)
%   !
%   ! This subroutine evalutes a bunch of orthogonal polynomials (and
%   ! their first partial derivatives) on the simplex with
%   ! vertices (0,0), (1,0), (0,1). The polynomials computed by this
%   ! routine are normalized and rotated Koornwinder polynomials, which
%   ! are classically defined on the triangle with vertices (0,0),
%   ! (1,0), (1,1), and given analytically as:
%   !
%   !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
%   !
%   ! After mapping to the uv-simplex via the change of variables
%   !
%   !      u = y,  v = 1-x
%   !
%   ! these polynomials are given as
%   !
%   !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
%   !        P_k((2u+v-1)/(1-v))
%   !
%   ! See Koornwinder 1975 for details, or the NIST handbook.
%   !
%   ! Input:
%   !   uv - a point in the simplex (0,0), (1,0), (0,1)
%   !   nmax - maximum degree polynomials to compute, a total of
%   !       (nmax+1)(nmax+2)/2 values are returned
%   !
%   ! Output:
%   !   pols - values of all the polynomials, ordered in the following
%   !       (n,k) manner:
%   !
%   !            (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), etc.
%   !
%   !   ders - partial derivatives with respect to u and v
%   !
%   !


  [~,m] = size(uv);
  legpols=zeros(m,nmax+1);
  legu   =zeros(m,nmax+1);
  legv   =zeros(m,nmax+1);
  jacpols=zeros(m,nmax+1,nmax+1);
  jacv   =zeros(m,nmax+1,nmax+1);

  u = uv(1,:);
  v = uv(2,:);
  u = u(:);
  v = v(:);
  z = 2*u+v-1;
  y = 1-v;

  legpols(:,1+0) = 1;
  legpols(:,1+1) = z;
  legu(:,1+0) = 0;
  legu(:,1+1) = 2;
  legv(:,1+0) = 0;
  legv(:,1+1) = 1;

  pkp1 = legpols(:,2);
  pk = legpols(:,1);

% Compute legendre polynomials
  for k=1:nmax
    pkm1 = pk;
    pk = pkp1;
    pkp1 = ((2*k+1)*z.*pk - k*pkm1.*y.*y)/(k+1);
    legpols(:,k+2) = pkp1;
    legu(:,k+2) = ((2*k+1)*(2*legpols(:,k+1) + z.*legu(:,k+1)) - ...
        k*legu(:,k).*y.*y)/(k+1);
    legv(:,k+2) = ((2*k+1)*(legpols(:,k+1) + z.*legv(:,1+k)) - ...
        k*(legv(:,k).*y.*y - 2*legpols(:,k).*y ))/(k+1);
  end

  x = 1-2*v;

% Compute Jacobi polynomials
  for k = 0:nmax
    beta = 2*k+1;
    jacpols(:,1+0,1+k) = 1;
    jacpols(:,1+1,1+k) = (-beta + (2+beta)*x)/2;

    jacv(:,1+0,1+k) = 0;
    jacv(:,1+1,1+k) = -(2+beta);
    
    for n = 1:nmax-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1);
      bn = (-beta^2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta);
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta);
      jacpols(:,1+n+1,1+k) = (an*x + bn).*jacpols(:,1+n,1+k) - ...
            cn*jacpols(:,1+n-1,1+k);
      jacv(:,1+n+1,1+k) = -2*an*jacpols(:,1+n,1+k) + ...
          (an.*x + bn).*jacv(:,1+n,1+k)- cn.*jacv(:,1+n-1,1+k);
    end

  end

  npols = (nmax+1)*(nmax+2)/2;
  pols = zeros(m,npols);
  deru = zeros(m,npols);
  derv = zeros(m,npols);
  iii = 0;
  for n = 0:nmax
    for k = 0:n
      sc = sqrt(1.0/(2*k+1)/(2*n+2));
      iii = iii + 1;
      pols(:,iii) = legpols(:,1+k).*jacpols(:,1+n-k,1+k)/sc;
      deru(:,iii) = legu(:,1+k).*jacpols(:,1+n-k,1+k)/sc;
      derv(:,iii) = (legv(:,1+k).*jacpols(:,1+n-k,1+k) + ...
          legpols(:,k+1).*jacv(:,1+n-k,1+k))/sc;
    end
  end

  pols = reshape(pols.',[npols,m]);
  deru = reshape(deru.',[npols,m]);
  derv = reshape(derv.',[npols,m]);
  ders = zeros([npols,m,2]);
  ders(:,:,1) = deru;
  ders(:,:,2) = derv;



  
end 










