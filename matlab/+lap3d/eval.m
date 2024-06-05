function [p, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%LAP3D.eval: solve Laplace boundary value problems using integral 
% equations
%
% Syntax for Dirichlet problems
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps);
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps, rep_pars);
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps, rep_pars, opts);
%
%   [p] = lap3d.eval(S, 'single', densities, targinfo, eps); 
%   [p] = lap3d.eval(S, 's', densities, targinfo, eps, opts); 
%   [p] = lap3d.eval(S, 'double', densities, targinfo, eps);
%   [p] = lap3d.eval(S, 'd', densities, targinfo, eps, opts);
%   [p] = lap3d.eval(S, 'comb', densities, targinfo, eps, rep_params);
%   [p] = lap3d.eval(S, 'c', densities, targinfo, eps, rep_params,opts);
%
% Syntax for Neumann problems
%   [densities] = lap3d.eval(S, 'neu', densities, targinfo, eps);
%   [densities] = lap3d.eval(S, 'neu', densities, targinfo, eps, opts);
%
%
%  SEE ALSO: lap3d.solver, helm3d.eval, stok3d.eval, and em3d.eval  
%
    switch lower(bc)
      case {'dir', 'dirichlet'}
        rep_params = [1; 1];
        opts = [];

        if nargin > 5
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('LAP3D.eval: rep_params, should be a double array');
          end
        end

        if nargin > 6
          opts = varargin{2};
        end
        [p, varargout{2:nargout}] = lap3d.dirichlet.eval(S, densities, targinfo, eps, ... 
                                    rep_params, opts);

      case {'s', 'single'}
        opts = [];
        if nargin > 5
          opts = varargin{1};
        end
        rep_params = [1.0; 0.0];
        [p, varargout{2:nargout}] = lap3d.dirichlet.eval(S, densities, targinfo, eps, ...
                                      rep_params, opts);

      case {'d', 'double'}
        opts = [];
        if nargin > 5
          opts = varargin{1};
        end
        rep_params = [0.0; 1.0];
        [p, varargout{2:nargout}] = lap3d.dirichlet.eval(S, densities, targinfo, eps, ...
                                      rep_params, opts);

      case {'c', 'comb'}
        if nargin < 6
          error('LAP3D.eval: not enough params for comb');
        end
        rep_params = varargin{1};
        opts = [];
        if nargin > 6
          opts = varargin{2};
        end
        [p, varargout{2:nargout}] = lap3d.dirichlet.eval(S, densities, targinfo, eps, ...
                                      rep_params, opts);
      case {'neu', 'neumann'}
        opts = [];
        if nargin > 5
          opts = varargin{1};
        end
        [p, varargout{2:nargout}] = lap3d.neumann.eval(S, densities, targinfo, eps, opts);
    end
    
end
