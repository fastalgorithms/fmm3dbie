function xmat = get_quad_cor(S, eps, zpars, iker)
%SURFWAVE.FLEX.GET_QUAD_COR near-field quadrature correction matrix for
%  flexural-gravity wave kernels, on-surface targets.
%
% Syntax:
%   xmat = surfwave.flex.get_quad_cor(S, eps, zpars, iker)
%
% Computes the sparse near-field quadrature correction matrix for one or
% more flexural-gravity wave kernels evaluated on-surface (targets = source
% nodes).
%
% Input:
%   S     - surfer object describing the surface discretization
%   eps   - requested quadrature precision
%   zpars - complex array of kernel parameters:
%             zpars(1:5)  = dispersion roots
%             zpars(6:10) = partial-fraction residues
%   iker  - kernel selector, scalar or vector of integers:
%             1 = G_s,          2 = G_phi,
%             3 = Delta^2 G_s,  4 = Delta^2 G_phi,
%             5 = S_{3d} G_phi, 6 = Laplace S_{3d},
%             7 = S_{3d} G_phi + Laplace S_{3d}
%
% Output:
%   xmat - if iker is scalar: (S.npts, S.npts) sparse complex correction
%          matrix.  If iker is a vector: cell array of sparse matrices,
%          one per entry of iker.


    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npts = S.npts;
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
    wnear{i} = surfwave.flex.getnearquad_flexural(npatches,norders,ixyzs, ...
          iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
          row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker(i)) ;

    end

    if length(iker) == 1
        wnear = wnear{1};
    end

   xmat = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
   
end
%
%
%
