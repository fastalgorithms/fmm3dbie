function xmat = get_quad_cor_sub(S, type, zpars, eps, gs_kern, ivpp)
%SURFWAVE.FLEX.GET_QUAD_COR_SUB near-field quadrature correction matrix for
%  a single flexural-gravity wave kernel by name, on-surface targets.
%
% Syntax:
%   xmat = surfwave.flex.get_quad_cor_sub(S, type, zpars, eps, gs_kern)
%   xmat = surfwave.flex.get_quad_cor_sub(S, type, zpars, eps, gs_kern, ivpp)
%
% Computes the sparse near-field quadrature correction matrix for the
% named kernel type on-surface (targets = source nodes). 
%
% Input:
%   S       - surfer object describing the surface discretization
%   type    - string kernel name:
%               'gs'        -> iker=1  (G_s)
%               'gphi'      -> iker=2  (G_phi)
%               'bilapgs'   -> iker=3  (Delta^2 G_s)
%               'bilapgphi' -> iker=4  (Delta^2 G_phi)
%               's3dgphi'   -> iker=5  (S_{3d} G_phi)
%               's3d'       -> iker=6  (Laplace S_{3d})
%   zpars   - complex array of kernel parameters:
%               zpars(1:5)  = dispersion roots
%               zpars(6:10) = partial-fraction residues
%   eps     - requested quadrature precision
%   gs_kern - kernel function handle (used by smooth_sparse_quad; may be
%             empty for on-surface-only usage)
%   ivpp    - (optional) flag: 1 = Chebyshev series MEX path
%             (default 1), 0 = standard MEX path
%
% Output:
%   xmat - (S.npts, S.npts) sparse complex near-field quadrature
%          correction matrix

    if strcmp(type,'gs')
        iker = 1;
    elseif strcmp(type,'gphi')
        iker = 2;
    elseif strcmp(type,'bilapgs')
        iker = 3;
    elseif strcmp(type,'bilapgphi')
        iker = 4;
    elseif strcmp(type,'s3dgphi')
        iker = 5;
    elseif strcmp(type,'s3d')
        iker = 6;
    else
        error('kernel name not recognized') 
    end 

    if nargin < 6
        ivpp = 1;
    end

    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npp1 = npatches+1;
    npts = S.npts;
    ntp1 = npts+1;
    n9 = 9;
    n12 = 12;
    ndtarg = 3;

    % [patch_id, uvs_targ] = get_patch_id_uvs(S);

    %this might need fixing

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % 
    % find near quadrature correction interactions
    % 

    targs = srcvals(1:3,:);
    ntarg = npts;

    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, npts);
    nnzp1 = nnz+1;

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,npts,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = cell(length(iker),1);
    for i = 1:length(iker)
    if iker < 6
        A = zeros(nquad,1,'like',1i);
        if ivpp
            maxdist = 1.5*sqrt((max(S.r(1,:)) - min(S.r(1,:))).^2 + (max(S.r(2,:)) - min(S.r(2,:))).^2);
            mex_id_ = 'getnearquad_flex_all_vpp(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, maxdist, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, 1, nquad);
        else
            mex_id_ = 'getnearquad_flex_all(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io dcomplex[x])';
[A] = kern_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, iker, A, 1, npatches, npp1, npatches, 1, n9, npts, n12, npts, 1, 13, 1, 1, ntp1, nnz, nnzp1, 1, 1, 1, nquad);
        end
        wnear{i} = A;
    elseif iker == 6
        opts.rep = 'eval';
        opts.format = 'rsc';
        wstruc = lap3d.neumann.get_quadrature_correction(S,eps,S,opts);
        wnear{i} = wstruc.wnear;
        row_ptr = wstruc.row_ptr;
        iquad = wstruc.iquad;
        col_ind = wstruc.col_ind;
        nnz = length(iquad)-1;
        nquad = iquad(nnz+1)-1;
    end

    end

    if length(iker) == 1
        wnear = wnear{1};
    end

   xmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);

end
%
%
%%
