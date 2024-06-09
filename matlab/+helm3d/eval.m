function [p, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%HELM3D.eval: eval Helmholtz BVP boundary integral equation solution
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
      case {'dir', 'dirichlet'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        rep_params = [-1j*zk; 1];
        opts = [];

        if nargin > 6
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
        end 
        if nargin > 7
          opts = varargin{3};
        end
        [p, varargout{1:nargout-1}] = helm3d.dirichlet.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);

      case {'s', 'single'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        rep_params = [1.0; 0.0];
        opts = [];

        if nargin > 6
          opts = varargin{2};
        end
        [p, varargout{1:nargout-1}] = helm3d.dirichlet.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);

      case {'d', 'double'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        rep_params = [0.0; 1.0];
        opts = [];

        if nargin > 6
          opts = varargin{2};
        end
        [p, varargout{1:nargout-1}] = helm3d.dirichlet.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);

      case {'c', 'comb'}
        if nargin < 7
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        if isa(varargin{2}, 'double')
          rep_params = varargin{2};
        else
          error('HELM3D.eval: rep_params, should be a double array');
        end
        opts = [];

        if nargin > 7
          opts = varargin{3};
        end
        [p, varargout{1:nargout-1}] = helm3d.dirichlet.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);

      case {'neu', 'neumann'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for neumann bc');
        end
        zk = varargin{1};
        alpha = 1;
        opts = [];

        if nargin > 6
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: alpha, should be a complex number');
          end
        end

        if nargin > 7
          opts = varargin{3};
        end
        [p, varargout{1:nargout-1}] = helm3d.neumann.eval(S, ...
            densities, targinfo, eps, zk, alpha, opts);

      case {'imp', 'impedance'}
        if nargin < 6
          error('HELM3D.eval: not enough input arguments for impedance bc');
        end

        zk = varargin{1};
        alpha = 1;
        opts = [];
        if nargin > 6
          if isa(varargin{2}, 'double')
            alpha = varargin{2};
          else
            error('HELM3D.eval: rep_params, should be a double array');
          end
        end

        if nargin > 7
          opts = varargin{3};
        end
        [p, varargout{1:nargout-1}] = helm3d.impedance.eval(S, ...
           densities, targinfo, eps, zk, alpha, opts);

      case {'trans', 'transmission'}
        if nargin < 7
          error('HELM3D.eval: not enough input arguments for dirichlet bc');
        end
        zks = varargin{1};
        if length(zks) ~=2
          error('HELM3D.SOLVER: transmission problem requires 2 wavenumbers');
        end

        zks = zks(:);

        if isa(varargin{2}, 'double')
          rep_params = varargin{2};
        else
          error('HELM3D.eval: rep_params, should be a double array');
        end
        if length(rep_params) ~=4
          error('HELM3D.SOLVER: transmission problem rep_params must be of length 4');
        end

        rep_params = rep_params(:);
        opts = [];

        if nargin > 7
          opts = varargin{3};
        end
        [p, varargout{1:nargout-1}] = helm3d.transmission.eval(S, ...
            densities, targinfo, eps, zks, rep_params, opts);

        otherwise
            error('HELM3D.EVAL: representation %s not found\n',bc);

    end
    
end
