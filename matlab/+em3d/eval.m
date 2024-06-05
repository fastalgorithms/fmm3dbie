function [E, H, varargout] = eval(S, bc, densities, targinfo, eps, varargin)
%EM3D.eval: solve Maxwell boundary value problems using integral 
% equations
%
% Syntax for perfect conductor problems
%   [E, H] = em3d.eval(S, 'pec', densities, targinfo, eps, zk);
%   [E, H] = em3d.eval(S, 'pec', densities, targinfo, eps, zk, rep_params);
%   [E, H] = em3d.eval(S, 'pec', densities, targinfo, eps, zk, rep_params, opts);
%
%
    switch lower(bc)
      case {'pec', 'perfect electric conductor'}
        if nargin < 6
          error('em3D.eval: not enough input arguments for dirichlet bc');
        end
        zk = varargin{1};
        opts = [];
        opts.rep = 'nrccie';
        rep_params = 1.0;

        if nargin < 6
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('em3D.eval: rep_params, should be a double array');
          end
        end

        if nargin > 7
          opts = varargin{3};
        end
        
        [E, H, varargout{3:nargout}] = em3d.pec.eval(S, ...
            densities, targinfo, eps, zk, rep_params, opts);
    end
    
end
