function [p, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%LAP3D.eval: solve Laplace boundary value problems using integral 
% equations
%
% Syntax for Dirichlet problems
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps);
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps, rep_pars);
%   [p] = lap3d.eval(S, 'dir', densities, targinfo, eps, rep_pars, opts);
%
% Syntax for Neumann problems
%   [densities] = lap3d.eval(S, 'neu', densities, targinfo, eps);
%   [densities] = lap3d.eval(S, 'neu', densities, targinfo, eps, opts);
%
%
%  SEE ALSO: lap3d.solver, helm3d.eval, stok3d.eval, and em3d.eval  
%
    switch lower(bc)
      case {'dir', 'd', 'dirichlet'}
        if nargin < 6
          rep_params = [1; 1];
          opts = [];
        elseif nargin < 7
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('LAP3D.eval: rep_params, should be a double array');
          end
          opts = [];
        else 
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('LAP3D.eval: rep_params, should be a double array');
          end
          opts = varargin{2};
        end
        [p, varargout{2:nargout}] = lap3d.dirichlet.eval(S, densities, targinfo, eps, ... 
                                    rep_params, opts);

      case {'neu', 'n', 'neumann'}
        if nargin < 6
          opts = [];
        else
          opts = varargin{1};
        end
        [p, varargout{2:nargout}] = lap3d.neumann.eval(S, densities, targinfo, eps, opts);
    end
    
end
