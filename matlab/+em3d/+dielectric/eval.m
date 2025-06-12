function [E, H] = eval(S, densities, targinfo, eps, om, rep_params, varargin)
%
%  em3d.dielectric.eval
%
%    This subroutine evaluates the electric and magnetic
%    field at a colelction of targets given the solution
%    to the corresponding integral equation
%
%
%  Notes for this routine:
%  The PDE takes the form
%  1v) \nabla \times E =  i\om \mu H
%  2v) \nabla \cdot  E =     0
%  3v) \nabla \times H = -i\om \ep E
%  4v) \nabla \cdot  H =     0
%  
%  where E is the electric field, H is the magnetic field, 
%  and \om is the wavenumber, \ep is the permittivity, 
%  and \mu is the permeability
%
%  The dielectric boundary conditions are given by
%  1b) n \times (E0 + E_in) = n \times E1
%  2b) n \times (H0 + H_in) = n \times H1
%
%  where (E_in, H_in) are the incoming electric and magnetic 
%  fields, E0, H0 are the fields in the exterior, and 
%  E1, H1 are the fields in the interior
% 
%  This routine will support the following representations:
%  * muller   (Muller integral equation)
%
%  For notes on the specific representations, boundary integral equations,
%  and order of kernels returned by this routine, checkout
%  em3d.dielectric.Contents.m
%
%  Syntax
%   [E, H] = em3d.dielectric.eval(S, densities, targinfo, eps, zk, rep_params)
%   [E, H] = em3d.dielectric.eval(S, densities, targinfo, eps, zk, rep_params, opts)
%
%  Note: No quadrature corrections are currently used 
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * densities: layer potential densities, of size (ndim, npts)
%        where ndim depends on the integral representation used
%    * targinfo: target info 
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * eps: precision requested
%    * om : wave number
%    * rep_params: parameters for integral representation 
%                  for muller, it should be a 4 vector
%                  consisting of [\ep0, \mu0, \ep1, \mu1]
%    * opts: options struct
%        opts.rep - integral representation being used
%                         Supported representations
%                         'muller'
%        opts.nonsmoothonly - use smooth quadrature rule for
%                             evaluating layer potential (true)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format (currently not supported) 
%        opts.in - boolean index of which points are inside the domain
%           targ(:,in) are the interior targets, and targ(:,~in) are
%           the exterior targets
%
%
%

    if(nargin < 7) 
      opts = [];
    else
      opts = varargin{1};
    end

    nonsmoothonly = true;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    isprecompq = false;
    if isfield(opts, 'precomp_quadrature')
      Q = opts.precomp_quadrature;
      isprecompq = true;
    end

    rep = 'muller';

    if isfield(opts, 'rep')
      rep = opts.rep;
    end

    if isfield(opts, 'in')
      in = opts.in;
    else
      sigma = ones(S.npts,1);
      flagext = lap3d.eval(S, 'double', sigma, targinfo, eps);
      in = find(flagext<=-0.5);
    end


% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    ff = 'rsc';

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

    targin = targs(:,in);
    targout = targs(:,~in);
    
    [~, ntin] = size(targin);
    [~, ntout] = size(targout);

    E = complex(zeros(3,ntarg));
    H = complex(zeros(3,ntarg));
    
    Ein = complex(zeros(3,ntin));
    Hin = complex(zeros(3,ntin));

    Eout = complex(zeros(3,ntout));
    Hout = complex(zeros(3,ntout));
    zpars = complex([om rep_params(:).']);

    
    if (ntin > 0)
      iside = 1;
      mex_id_ = 'em_muller_trans_eval_oneside(i double[x], i dcomplex[xx], i dcomplex[x], i int[x], i double[x], i double[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx])';
[Ein, Hin] = fmm3dbie_routs(mex_id_, eps, densities, zpars, npts, wts, srcvals, ntin, targin, iside, Ein, Hin, 1, npts, 4, 5, 1, npts, 12, npts, 1, 3, ntin, 1, 3, ntin, 3, ntin);
      E(:,in) = Ein;
      H(:,in) = Hin;
    end

    if(ntout > 0)
      iside = 0;
      mex_id_ = 'em_muller_trans_eval_oneside(i double[x], i dcomplex[xx], i dcomplex[x], i int[x], i double[x], i double[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx])';
[Eout, Hout] = fmm3dbie_routs(mex_id_, eps, densities, zpars, npts, wts, srcvals, ntout, targout, iside, Eout, Hout, 1, npts, 4, 5, 1, npts, 12, npts, 1, 3, ntout, 1, 3, ntout, 3, ntout);
      E(:,out) = Eout;
      H(:,out) = Hout;
    end

end    
%
%
%
%----------------------------------
%
