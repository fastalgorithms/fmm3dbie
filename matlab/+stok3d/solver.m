function [densities, varargout] = solver(S, bc, rhs, eps, varargin)
%stok3d.SOLVER: solve Stokes boundary value problems using integral 
% equations
%
% Syntax for Velocity problems
%   [densities] = stok3d.solver(S, 'vel', rhs, eps);
%   [densities] = stok3d.solver(S, 'vel', rhs, eps, rep_pars);
%   [densities] = stok3d.solver(S, 'vel', rhs, eps, rep_pars, opts);
% 
% SEE ALSO: stok3d.eval, stok3d.velocity.solver
%
    switch lower(bc)
      case {'vel', 'velocity'}
        if nargin < 5
          rep_params = [1; 1];
          opts = [];
        elseif nargin < 6
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('stok3d.SOLVER: rep_params, should be a double array');
          end
          opts = [];
        else
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('stok3d.SOLVER: rep_params, should be a double array');
          end
          opts = varargin{2};
        end
        [densities, varargout{1:nargout-1}] = stok3d.velocity.solver(S, rhs, eps, ... 
                                    rep_params, opts);

        otherwise
            error('STOK3D.SOLVER: boundary condition %s not found\n',bc);
    end

    
end
