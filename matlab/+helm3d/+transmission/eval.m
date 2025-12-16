function p = eval(S,densities,targinfo,eps,zks,rep_params,varargin)
%
%  helm3d.transmission.eval
%    Evaluates the helmholtz combined layer potential at a collection 
%    of targets either in the interior/exterior
%
%  Syntax
%   pot = helm3d.transmission.eval(S,densities,targinfo,eps,zks,rep_params)
%
%  PDE:
%      (\Delta + k0^2)u0 = 0   in the exterior
%      (\Delta + k1^2)u1 = 0   in the interior
%
%      alpha0 u0 - alpha1 u1 = f 
%      beta0 dudn0 - beta1 dudn1 = g
%
%  Integral representation
%     pot0 = 1/beta0 (ik_{0} S_{k0} [\sigma] + D_{k0} [\rho]) (exterior representation)
%     pot1 = 1/beta1 (ik_{1} S_{k1} [\sigma] + D_{k1} [\rho]) (interior representation)
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  densities(2,1:npts) for helm_comb, with densities(1,:) = rho, 
%  and densities(2,:) = sigma
%
%  rep_params(1) = alpha0 
%  rep_params(2) = beta0
%  rep_params(3) = alpha1 
%  rep_params(4) = beta1
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * densities: layer potential densities
%    * targinfo: target info 
%       targinfo.r = (3,nt) target locations
%       targinfo.n = normal info
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * eps: precision requested
%    * zks: wave numbers (k)
%    * rep_params: [alpha0; beta0; alpha1; beta1] above 
%
% todo: add options



% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    if(isfield(targinfo,'patch_id') || isprop(targinfo,'patch_id'))
      patch_id = targinfo.patch_id;
    else
      patch_id = zeros(ntarg,1);
    end

    if(isfield(targinfo,'uvs_targ') || isprop(targinfo,'uvs_targ'))
      uvs_targ = targinfo.uvs_targ;
    else
      uvs_targ = zeros(2,ntarg);
    end

    if(length(patch_id)~=ntarg)
      fprintf('Incorrect size of patch id in target info struct. Aborting! \n');
    end

    [n1,n2] = size(uvs_targ);
    if(n1 ~=2 && n2 ~=ntarg)
      fprintf('Incorrect size of uvs_targ array in targinfo struct. Aborting! \n');
    end


    p = complex(zeros(ntarg,1));

    ndz = 6;
    zpars = complex(zeros(6,1));
    zpars(1) = zks(1);
    zpars(2) = rep_params(1);
    zpars(3) = rep_params(2);
    zpars(4) = zks(2);
    zpars(5) = rep_params(3);
    zpars(6) = rep_params(4);
    ndim_s = 2;
% Call the layer potential evaluator
    mex_id_ = 'helm_comb_trans_eval(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i double[xx], i double[x], i dcomplex[x], i dcomplex[xx], io dcomplex[x])';
[p] = fmm3dbie_routs(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, patch_id, uvs_targ, eps, zpars, densities, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 6, 2, npts, ntarg);
end    
%
%
%
%----------------------------------
%
