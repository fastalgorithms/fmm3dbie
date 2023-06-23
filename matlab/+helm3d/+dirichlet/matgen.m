function A = matgen(i,j,S,zpars,P,Q,opts)
%
% 
%  helm3d.dirichlet.matgen
%    This subroutine returns the matrix block A(I,J) corresponding 
%    to the discretization points indexed by i\in I and j\in J. 
%  
%  Syntax: 
%    A = helm3d.dirichlet.matgen(i,j,S,zpars)
%    A = helm3d.dirichlet.matgen(i,j,S,zpars,P,Q)
%    A = helm3d.dirichlet.matgen(i,j,S,zpars,P,Q,opts)
%
%  Integral representation
%     u = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%  
%  S_{k}, D_{k}: helmholtz single and double layer potentials respectively
%
%  k, \alpha, beta = zpars(1:3)
%
%  Note that in this subroutine, the default functionality it to return
%  the matrix entries with the identitiy term included
%
%  Input arguments:
%    * i: row indices
%    * j: column indices
%    * S: surfer object
%    * zpars: kernel parameters
%        zpars(1) - wave number
%        zpars(2) - single layer strength
%        zpars(3) - double layer strength
%    * P: temporary array needed for efficient extraction of 
%         quadrature correction
%    * Q: quadrature correction stored as a sparse matrix
%    * opts: options struct
%        opts.diag (true), whether to include the diagonal correction or not
%     
     if isempty(i) || isempty(j)
        A = zeros(length(i),length(j));
        return;
     end
     [I,J] = ndgrid(i,j);
     srcinfo = [];
     srcinfo.r = S.r(:,i);
     srcinfo.n = S.n(:,i);
     targinfo = [];
     targinfo.r = S.r(:,j);
     targinfo.n = S.n(:,j);
     if(nargin <=6)
       opts = [];
     end

     ifdiag = 1;
     if(isfield(opts,'ifdiag'))
        ifdiag = opts.ifdiag;
     end

     ifout = 1;
     if(isfield(opts,'ifout'))
        ifout = opts.ifout;
     end

     ifl2scale = 0;
     if(isfield(opts,'l2scale'))
        ifl2scale = opts.l2scale;
     end
     
     A = bsxfun(@times,helm3d.dirichlet.kern(zpars(1),srcinfo,targinfo,'c',zpars(2),zpars(3)), ...
            S.wts(j));

     if(nargin > 4)
        M = spget_quadcorr(i,j,P,Q.spmat);
        idx = abs(M) ~=0;
        A(idx) = M(idx); 
     end

     if(ifdiag)
        A(I==J) = A(I==J) - (-1)**ifout*0.5*zpars(3) 
     end

     if(l2scale)
        A = bsxfun(@times,sqrt(S.wts(i)).',A);
        A = bsxfun(@times,A,1.0/sqrt(S.wts(j)));
     end
end
