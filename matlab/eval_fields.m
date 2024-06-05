function [varargout] = eval_fields(S, pde, bc, densities, targinfo, eps, varargin)
%
%EVALUATE_FIELDS: Evaluate the fields for solution from integral equation 
%  
%  Syntax
%    p = eval_fields(S, pde, bc, densities, eps, targinfo) 
%    p = eval_fields(S, pde, bc, densities, eps, targinfo, opts) 
%    [E, H] = eval_fields(S, pde, bc, densities, eps, targinfo) 
%    [E, H] = eval_fields(S, pde, bc, densities, eps, targinfo, opts) 
%
%  Currently supported
%  combination solvers for pde and bc are
%
%      pde                               bc
%      ---                               --
%      
%      'laplace' ('lap', 'l')            'dir', 'neu'
%      'helmholtz' ('helm', 'h')         'dir', 'neu', 'imp'
%      'stokes' ('stok', 's')            'vel'
%      'maxwell' ('em', 'm')             'pec'
%
%  The boundary conditions may also be written in longer form, e.g.
%    'dirichlet', 'neumann', 'impedance', 'velocity',
%    'perfect electric conductor'. 
%  
%  Different solvers require different parameters to be specified. 
%  Such parameters should be passed as trailing arguments.
%  For example, the Laplace Dirichlet, and Stokes velocity boundary
%  value problems require a double array of size 2 which determines
%  the strength of the single and double layer potentials in the
%  representation. 
%
%  Similarly, the helmholtz, and maxwell representations, require
%  a complex array, which includes the wavenumber as the first
%  parameter, and other parameters specific to the representation.
%  
%  The size of densities for the solver is also dependent on
%  the pde, bc, and representation being used. It should be 
%  an array of size (nd, S.npts) 
%  which correspond to the densities on surface if nd > 1, and
%  of size (S.npts,1) if nd == 1.
%
%  And finally, for Maxwell boundary value problems, the data
%  must be a complex array of size (6, S.npts), where
%  (1:3, :) are the cartesian components of the incident electric
%  field, and (4:6, :) are the cartesian components of the incident
%  magnetic fields.
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * pde, bc: strings for determining pde being solved and representation
%               being used
%    * densities: layer potential densities, of size (ndim, npts)
%    * eps: tolerance requested
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%
%  SEE ALSO: lap3d.eval, helm3d.eval, stok3d.eval, and em3d.eval 
%
     switch lower(pde)
       case {'laplace', 'lap', 'l'}
         [varargout{1:nargout}] = lap3d.eval(S, bc, densities, ...
                                               targinfo, eps, varargin{:});
       case {'helmholtz', 'helm', 'h'}
         [varargout{1:nargout}] = helm3d.eval(S, bc, densities, ...
                                               targinfo, eps, varargin{:});
       case {'stokes', 'stok', 's'}
         [varargout{1:nargout}] = stok3d.eval(S, bc, densities, ...
                                               targinfo, eps, varargin{:});  
       case {'maxwell', 'em', 'm'}
         [varargout{1:nargout}] = em3d.eval(S, bc, densities, ...
                                               targinfo, eps, varargin{:});  
       otherwise
         error('PDE ''%s'' not found.', pde);
    end
end
