function A = slp_matgen(S, eps)
%
%  lap2d.slp_matgen
%
%  Syntax
%   A = lap2d.slp_matgen(S,eps)
%
%  Integral representation
%     pot = (Update representation here)
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * eps: precision requested
%


    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;
    npp1 = npatches+1;
    n3 = 3;
    row_ptr = 1:npatches:(npatches*npts+1);
    col_ind = repmat(1:npatches,[1,npts]);
    row_ptr = row_ptr(:);
    col_ind = col_ind(:);
    zpuse = 1j;
    nnz = npatches*npts;
    nnzp1 = nnz + 1;
    iquadtype = 1;
    ntp1 = npts + 1;
    npols = ixyzs(2) - ixyzs(1);
    iquad = 1:npols:(npts*npts+1);
    A = zeros(npts,npts);
    rfac0 = 1.25;
    nquad = npts*npts;

    mex_id_ = 'getnearquad_lap2d_gv2v(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 1, 1, ntp1, nnz, nnzp1, 1, 1, nquad);
    
    A = reshape(A,[S.npts,S.npts]).';

end
%
%
rior square in the discretization
%    * npars(3): integer (optional, [3,3,3])
%        number of quad patches in u, and v direction. if iptype is 1,
%        then each quad patch is further divided into two triangular
%        patches
%    * norder: (optional, 4)
%        order of discretization
%    * iptype: (optional, 1)
%        type of patch to be used in the discretization
%        * iptype = 1, triangular patch discretized using
%                      Rokhlin-Vioreanu nodes
%        * iptype = 11, quadrangular patch discretized using tensor
%                       product Gauss-Legendre nodes
%        * iptype = 12, quadrangular patch discretized using tensor
%                       product Chebyshev
%    * iort: (optional, 1)
%         orientation vector, normals point in the unbounded region
%         if iort = 1, and in the interior of the startorus otherwise
%
% Example
%   % create 8th-order GL-quad patch for an ellipse with semi-major
%   axes (1,0.5):
%   S = geometries.disk([1,0.5], [], [], [], 8, 11)
%
  if nargin < 1 || isempty(scales)
    scales = [1;1];
  end

  if nargin < 2 || isempty(rmid)
    rmid = 0.5;
  end

  if nargin < 3 || isempty(npars)
    npars = [3;3;3];
  end

  if nargin < 4 || isempty(norder)
    norder = 4;
  end

  if nargin < 5 || isempty(iptype)
    iptype = 1;
  end

  if nargin < 6 || isempty(iort)
    iort = -1;
  end


  if iptype == 1
    npatches = 2*(npars(1)*npars(1) + 4*npars(2)*npars(3))
    npts = npatches*(norder+1)*(norder+2)/2
  elseif iptype == 11 || iptype == 12
    npatches = npars(1)*npars(1) + 4*npars(2)*npars(3)
    npts = npatches*(norder+1)*(norder+1)
  end

  npp1 = npatches+1;
  ixys = zeros(npp1,1);
  ptinfo = zeros(6,npts);

  mex_id_ = 'mesh_circle_pts(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io double[xx])';
[ixys, ptinfo] = kern_routs(mex_id_, rmid, npars, iort, norder, iptype, npatches, npts, ixys, ptinfo, 1, 3, 1, 1, 1, 1, 1, npp1, 6, npts);

  norders = norder*ones(npatches,1);
  srcvals = zeros(12,npts);
  srcvals(1,:) = scales(1)*ptinfo(1,:);
  srcvals(2,:) = scales(2)*ptinfo(2,:);
  srcvals(4,:) = scales(1)*ptinfo(3,:);
  srcvals(5,:) = scales(2)*ptinfo(4,:);
  srcvals(7,:) = scales(1)*ptinfo(5,:);
  srcvals(8,:) = scales(2)*ptinfo(6,:);
  srcvals(12,:) = iort/abs(iort);

  iptype_all = iptype*ones(npatches,1);

  S = surfer(npatches, norders, srcvals, iptype_all);



end
%
%
%
%
%
