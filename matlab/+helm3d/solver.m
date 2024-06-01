function [densities, varargout] = solver(S, bc, rhs, eps, varargin)
%LAP3D.SOLVER: solve Helmholtz boundary value problems using integral 
% equations
%
% Syntax for Dirichlet problems
%   [densities] = helm3d.solver(S, 'dir', rhs, eps, zk);
%   [densities] = helm3d.solver(S, 'dir', rhs, eps, zk, rep_params);
%   [densities] = helm3d.solver(S, 'dir', rhs, eps, zk, rep_params, opts);
%
% Syntax for Neumann problems
%   [densities] = helm3d.solver(S, 'neu', rhs, eps, zk);
%   [densities] = helm3d.solver(S, 'neu', rhs, eps, zk, alpha);
%   [densities] = helm3d.solver(S, 'neu', rhs, eps, zk, alpha, opts);
%
% Syntax for Impedance problems
%   [densities] = helm3d.solver(S, 'imp', rhs, eps, zlams, zk);
%   [densities] = helm3d.solver(S, 'imp', rhs, eps, zlams, zk, alpha);
%   [densities] = helm3d.solver(S, 'imp', rhs, eps, zlams, zk, alpha, opts);
%
    switch lower(bc)
      case {'dir', 'd', 'dirichlet'}
        if nargin < 5
          error('HELM3D.SOLVER: not enough input arguments for dirichlet bc');
        elseif nargin < 6
          zk = varargin{1};
          rep_params = [-1j*zk; 1];
          opts = [];
        elseif nargin < 7
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
          opts = [];
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
          opts = varargin{3};
        end
        [densities, varargout{2:nargout}] = helm3d.dirichlet.solver(S, ...
            rhs, eps, zk, rep_params, opts);

      case {'neu', 'n', 'neumann'}
        if nargin < 5
          error('HELM3D.SOLVER: not enough input arguments for neumann bc');
        elseif nargin < 6
          zk = varargin{1};
          alpha = 1;
          opts = [];
        elseif nargin < 7
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.SOLVER: alpha, should be a complex number');
          end
          opts = [];
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.SOLVER: alpha should be a complex number');
          end
          opts = varargin{3};
        end
        [densities, varargout{2:nargout}] = helm3d.neumann.solver(S, ...
            rhs, eps, zk, alpha, opts);

      case {'imp', 'i', 'impedance'}
        if nargin < 6
          error('HELM3D.SOLVER: not enough input arguments for impedance bc');
        elseif nargin < 7
          if isa(varargin{1}, 'double')
              zlams = varargin{1};
          else
              error('HELM3D.SOLVER: zlams must be a complex array');
          end
          zk = varargin{2};
          alpha = 1;
          opts = [];
        elseif nargin < 8
          if isa(varargin{1}, 'double')
              zlams = varargin{1};
          else
              error('HELM3D.SOLVER: zlams must be a complex array');
          end
          zk = varargin{2};
          if isa(varargin{3}, 'double')
            alpha = varargin{3};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
          opts = [];
        else 
          if isa(varargin{1}, 'double')
              zlams = varargin{1};
          else
              error('HELM3D.SOLVER: zlams must be a complex array');
          end
          zk = varargin{2};
          if isa(varargin{3}, 'double')
            alpha = varargin{3};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
          opts = varargin{4};
        end
        [densities, varargout{2:nargout}] = helm3d.impedance.solver(S, ...
            zlams, rhs, eps, zk, alpha, opts);


    end
    
end
