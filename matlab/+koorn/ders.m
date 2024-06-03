function [pols,dersu,dersv] = ders(norder, uv)
  %
  % This subroutine evalutes a bunch of orthogonal polynomials (and
  % their first partial derivatives) on the simplex with
  % vertices (0,0), (1,0), (0,1). The polynomials computed by this
  % routine are normalized and rotated Koornwinder polynomials, which
  % are classically defined on the triangle with vertices (0,0),
  % (1,0), (1,1), and given analytically as:
  %
  %   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
  %
  % After mapping to the uv-simplex via the change of variables
  %
  %      u = y,  v = 1-x
  %
  % these polynomials are given as
  %
  %   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
  %        P_k((2u+v-1)/(1-v))
  %
  % See Koornwinder 1975 for details, or the NIST handbook.
  %
  % Input:
  %   uv - a point in the simplex (0,0), (1,0), (0,1)
  %   nmax - maximum degree polynomials to compute, a total of
  %       (nmax+1)(nmax+2)/2 values are returned
  %
  % Output:
  %   pols - values of all the polynomials, ordered in the following
  %       (n,k) manner:
  %
  %            (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), etc.
  %
  %   ders - partial derivatives with respect to u and v
  %
  %

  npols = (norder+1)*(norder+2)/2;

  pols = zeros(npols,size(uv,2));
  dersu = zeros(npols,size(uv,2));
  dersv = zeros(npols,size(uv,2));

  legpols = zeros(1,norder+1);
  jacpols = zeros(norder+1,norder+1);
  legu = zeros(1,norder+1);
  legv = zeros(1,norder+1);
  jacv = zeros(norder+1,norder+1);



  done = 1;
  if (norder >= 100) 
    fprintf("norder too large, nmae = %i\n", norder);
    return
  end
  for i = 1:size(uv,2)

  %
  % first run the recurrence for P_k(z/y)*y^k
  %
  u = uv(1,i);
  v = uv(2,i);
  z = 2*u+v-1;
  y = 1-v;

  legpols(1) = 1;
  legpols(2) = z;

  legu(1) = 0;
  legu(2) = 2;

  legv(1) = 0;
  legv(2) = 1;

  for k=1:norder
    legpols(k+2) = ((2*k+1)*z*legpols(k+1) - k*legpols(k)*y*y)/(k+1);
    legu(k+2) = ((2*k+1)*(2*legpols(k+1) + z*legu(k+1)) - ...
        k*legu(k)*y*y)/(k+1);
    legv(k+2) = ((2*k+1)*(legpols(k+1) + z*legv(k+1)) - ...
        k*(legv(k)*y*y - 2*legpols(k)*y ))/(k+1);
  end


  %
  % now loop over degrees, in reverse order,
  % and for each run the jacobi recurrence
  %
  x = 1-2*v;
  
  for k = 0:norder
    beta = 2*k+1;
    jacpols(1,k+1) = 1;
    jacpols(2,k+1) = (-beta + (2+beta)*x)/2;

    jacv(1,k+1) = 0;
    jacv(2,k+1) = -(2+beta);

    for n = 1:norder-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1);
      bn = (-beta^2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta);
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta);
      jacpols(n+2,k+1) = (an*x + bn)*jacpols(n+1,k+1) - cn*jacpols(n,k+1);
      jacv(n+2,k+1) = -2*an*jacpols(n+1,k+1) + (an*x + bn)*jacv(n+1,k+1) ...
          - cn*jacv(n,k+1);
    end

  end


  %
  % now assemble the ORTHONORMAL koornwinder polynomials
  %
  iii = 0;
  for n = 0:norder
    for k = 0:n
      sc = sqrt(1.0/(2*k+1)/(2*n+2)) ;
      iii = iii + 1;
      pols(iii,i) = legpols(k+1)*jacpols(n-k+1,k+1)/sc ;
      dersu(iii,i) = legu(k+1)*jacpols(n-k+1,k+1)/sc ;
      dersv(iii,i) = (legv(k+1)*jacpols(n-k+1,k+1) + ...
           legpols(k+1)*jacv(n-k+1,k+1))/sc;
    end
  end
  end
end 
