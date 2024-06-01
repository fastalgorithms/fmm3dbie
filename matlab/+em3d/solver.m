function [densities, varargout] = solver(S, bc, rhs, eps, varargin)
%EM3D.SOLVER: solve Maxwell boundary value problems using integral 
% equations
%
% Syntax for perfect conductor problems
%   [densities] = em3d.solver(S, 'pec', rhs, eps, zk);
%   [densities] = em3d.solver(S, 'pec', rhs, eps, zk, rep_params);
%   [densities] = em3d.solver(S, 'pec', rhs, eps, zk, rep_params, opts);
%
%
    switch lower(bc)
      case {'pec', 'p', 'perfect electric conductor'}
        if nargin < 5
          error('em3D.SOLVER: not enough input arguments for dirichlet bc');
        elseif nargin < 6
          zk = varargin{1};
          opts = [];
          opts.rep = 'nrccie';
          rep_params = 1;
        elseif nargin < 7
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('em3D.SOLVER: rep_params, should be a double array');
          end
          opts = [];
          opts.rep = 'nrccie';
        else 
          zk = varargin{1};
          if isa(varargin{2}, 'double')
            rep_params = varargin{2};
          else
            error('em3D.SOLVER: rep_params, should be a double array');
          end
          opts = varargin{3};
        end
        
        if isa(rhs, 'double')
          [m, n] = size(rhs);
          if m~=6 && n~= S.npts
             error('em3d.SOLVER: rhs should be of size (6, npts)');
          end
          einc = rhs(1:3,:);
          hinc = rhs(4:6,:);
        elseif isa(rhs, 'struct')
          einc = rhs.einc;
          hinc = rhs.hinc;
        else
          error('em3d.SOLVER: rhs must be a struct or an array');
        end
        [densities, varargout{2:nargout}] = em3d.pec.solver(S, ...
            einc, hinc, eps, zk, rep_params, opts);
    end
    
end
