function [densities, varargout] = solver(S, bc, rhs, eps, varargin)
%HELM3D.SOLVER  Solve Helmholtz boundary value problems using BIE
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
% Syntax for Transmission problems
%   [densities] = helm3d.solver(S, 'trans', rhs, eps, zk, rep_params);
%   [densities] = helm3d.solver(S, 'trans', rhs, eps, zk, rep_params, opts);
%
% 
% SEE ALSO: helm3d.eval, helm3d.dirichlet.solver, helm3d.neumann.solver, helm3d.impendance.solver, 
%           helm3d.transmission.solver

%--------------------
%
    switch lower(bc)
      case {'dir', 'dirichlet'}
        if nargin < 5
          error('HELM3D.SOLVER: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        rep_params = [-1j*zk; 1];
        opts = [];

        if nargin > 5
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
        end

        if nargin > 6
          opts = varargin{3};
        end

        [densities, varargout{1:nargout-1}] = helm3d.dirichlet.solver(S, ...
            rhs, eps, zk, rep_params, opts);

      case {'neu', 'neumann'}
        if nargin < 5
          error('HELM3D.SOLVER: not enough input arguments for neumann bc');
        end
        zk = varargin{1};
        alpha = 1;
        opts = [];

        if nargin > 5
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.SOLVER: alpha, should be a complex number');
          end
        end

        if nargin > 6
          opts = varargin{3};
        end
        [densities, varargout{1:nargout-1}] = helm3d.neumann.solver(S, ...
            rhs, eps, zk, alpha, opts);

      case {'imp', 'impedance'}
        if nargin < 6
          error('HELM3D.SOLVER: not enough input arguments for impedance bc');
        end
        
        if isa(varargin{1}, 'double')
          zlams = varargin{1};
        else
          error('HELM3D.SOLVER: zlams must be a complex array');
        end
        zk = varargin{2};
        alpha = 1;
        opts = [];

        if nargin > 6
          if isa(varargin{3}, 'double')
            alpha = varargin{3};
          else
            error('HELM3D.SOLVER: rep_params, should be a double array');
          end
        end

        if nargin > 7
          opts = varargin{4};
        end
        [densities, varargout{1:nargout-1}] = helm3d.impedance.solver(S, ...
            zlams, rhs, eps, zk, alpha, opts);

      case {'trans', 'transmission'}
        if nargin < 6
          error('HELM3D.SOLVER: not enough input args for helm transmission solver\n');
        end
        opts = [];
        if nargin > 6 
          opts = varargin{3};
        end

        rep_params = varargin{2};
        if length(rep_params) ~=4
          error('HELM3D.SOLVER: transmission problem rep_params must be of length 4');
        end
        rep_params = rep_params(:);

        zks = varargin{1};
        if length(zks) ~=2
          error('HELM3D.SOLVER: transmission problem requires 2 wavenumbers');
        end
        zks = zks(:);
        [densities, varargout{1:nargout-1}] = helm3d.transmission.solver(S, ...
            rhs, eps, zks, rep_params, opts);

        otherwise
            error('HELM3D.SOLVER: boundary condition %s not found\n',bc);
    end
    
end
