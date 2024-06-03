function [vel, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%stok3d.eval: solve Stokes boundary value problems using integral 
% equations
%
% Syntax for Velocity problems
%   [vel] = stok3d.eval(S, 'vel', densities, targinfo, eps);
%   [vel] = stok3d.eval(S, 'vel', densities, targinfo, eps, rep_pars);
%   [vel] = stok3d.eval(S, 'vel', densities, targinfo, eps, rep_pars, opts);
%
%
    switch lower(bc)
      case {'v', 'vel', 'velocity'}
        if nargin < 6
          rep_params = [1; 1];
          opts = [];
        elseif nargin < 7
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('stok3d.eval: rep_params, should be a double array');
          end
          opts = [];
        else
          if isa(varargin{1}, 'double')
            rep_params = varargin{1};
          else
            error('stok3d.eval: rep_params, should be a double array');
          end
          opts = varargin{2};
        end
        [vel, varargout{2:nargout}] = stok3d.velocity.eval(S, densities, targinfo, eps, ... 
                                    rep_params, opts);
    end

    
end
