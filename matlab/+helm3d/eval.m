function [p, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%HELM3D.eval: solve Helmholtz boundary value problems using integral 
% equations
%
% Syntax for Dirichlet problems
%   [p] = helm3d.eval(S, 'dir', densities, targinfo, eps, zk);
%   [p] = helm3d.eval(S, 'dir', densities, targinfo, eps, zk, rep_params);
%   [p] = helm3d.eval(S, 'dir', densities, targinfo, eps, zk, rep_params, opts);
%
% Syntax for Neumann problems
%   [p] = helm3d.eval(S, 'neu', densities, targinfo, eps, zk);
%   [p] = helm3d.eval(S, 'neu', densities, targinfo, eps, zk, alpha);
%   [p] = helm3d.eval(S, 'neu', densities, targinfo, eps, zk, alpha, opts);
%
% Syntax for Impedance problems
%   [p] = helm3d.eval(S, 'imp', densities, targinfo, eps, zk);
%   [p] = helm3d.eval(S, 'imp', densities, targinfo, eps, zk, alpha);
%   [p] = helm3d.eval(S, 'imp', densities, targinfo, eps, zk, alpha, opts);
%
    switch lower(bc)
      case {'dir', 'd', 'dirichlet'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        elseif nargin < 7
          zk = varargin{1};
          rep_params = [-1j*zk; 1];
          opts = [];
        elseif nargin < 8
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
          opts = [];
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
          opts = varargin{3};
        end
        [p, varargout{2:nargout}] = helm3d.dirichlet.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);

      case {'neu', 'n', 'neumann'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for neumann bc');
        elseif nargin < 7
          zk = varargin{1};
          alpha = 1;
          opts = [];
        elseif nargin < 8
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: alpha, should be a complex number');
          end
          opts = [];
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: alpha should be a complex number');
          end
          opts = varargin{3};
        end
        [p, varargout{2:nargout}] = helm3d.neumann.eval(S, ...
            densities, targinfo, eps, zk, alpha, opts);

      case {'imp', 'i', 'impedance'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for impedance bc');
        elseif nargin < 7
          zk = varargin{1};
          alpha = 1;
          opts = [];
        elseif nargin < 8
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
          opts = [];
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
          opts = varargin{3};
        end
        [p, varargout{2:nargout}] = helm3d.impedance.eval(S, ...
           densities, targinfo, eps, zk, alpha, opts);


    end
    
end
