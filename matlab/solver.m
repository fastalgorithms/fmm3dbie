function [densities, varargout] = solver(S, pde, bc, rhs, eps, varargin)
%
%SOLVER: Solve elliptic PDE with integral equations
%
%  densities = solver(S, pde, bc, rhs, eps, opts). Currently supported
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
%    'perfect electric conductor'
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
%  The right hand side (rhs) for the solver is also dependent on
%  the pde being solved. For Laplace problems, the right hand side
%  must be a double array of size (S.npts,1) 
%
%  For Stokes boundary value problems, the data must be (3, S.npts)
%  i.e. the cartesian components of the velocity on the boundary
%
%  For Helmholtz Dirichlet, Neumann, and Impedance problems, the
%  data must be a complex array of size (S.npts,1)
%
%  And finally, for Maxwell boundary value problems, the data
%  must be a complex array of size (6, S.npts), where
%  (1:3, :) are the cartesian components of the incident electric
%  field, and (4:6, :) are the cartesian components of the incident
%  magnetic fields
%
%  The output densities is an array of size (nd, S.npts) 
%  which correspond to the densities on surface if nd > 1, and
%  of size (S.npts,1) if nd == 1.
%
%  SEE ALSO: lap3d.solver, helm3d.solver, stok3d.solver, and em3d.solver 
%
     switch lower(pde)
       case {'laplace', 'lap', 'l'}
         [densities, varargout{2:nargout}] = lap3d.solver(S, bc, rhs, ...
                                               eps, varargin{:});
       case {'helmholtz', 'helm', 'h'}
         [densities, varargout{2:nargout}] = helm3d.solver(S, bc, rhs, ...
                                               eps, varargin{:});
       case {'stokes', 'stok', 's'}
         [densities, varargout{2:nargout}] = stok3d.solver(S, bc, rhs, ...
                                               eps, varargin{:});  
       case {'maxwell', 'em', 'm'}
         [densities, varargout{2:nargout}] = em3d.solver(S, bc, rhs, ...
                                               eps, varargin{:});  

       otherwise
         error('PDE ''%s'' not found.', pde);
    end
end
