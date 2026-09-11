function [sol, flag, relres, iter, resvec, times] = surfermatgmres(surferobj, kern, rhs, eps, tol, maxit, objover, cors, opts)
%SURFERMATGMRES Solve the integral equation system for given kernel and
% surfer array using unrestarted GMRES. 
%
% Syntax: sol = surfermatgmres(surferobj, kern, rhs)
%         [sol,flag,relres,iter,resvec,times] = surfermatgmres(surferobj, kern, rhs, ...
%             eps, tol, maxit, objover, cors, opts)
%
% Input:
%   surferobj - array of surfer objects describing boundary
%   kern      - kernel3d object or matrix of kernel3d objects
%   rhs       - right-hand side, size nrows x nrhs
%
% Optional input:
%   eps     - (1e-6) quadrature tolerance
%   tol     - (eps) gmres relative residual tolerance
%   maxit   - (min(nrows,200)) maximum number of gmres iterations
%   objover - oversampling specification, see surfermat
%   cors    - sparse matrix of near-field quadrature corrections, from
%             surfermat with corrections=true.
%   opts    - see surfermatapply for supported options
%
% Output:
%   sol    - solution, nrows x nrhs
%   flag   - nrhs x 1, gmres convergence flag for each right-hand side
%   relres - nrhs x 1, gmres relative residual for each right-hand side
%   iter   - nrhs x 1, gmres iteration count for each right-hand side
%   resvec - 1 x nrhs cell array of gmres residual history vectors
%   times  - (nrhs+1) x 1, times(1) is the precomputation time, 
%            times(1+i) is the gmres solve time for right-hand side i
%
% See also SURFERMAT, SURFERMATAPPLY.
%
% Author: Tristan Goodwill

if nargin < 4, eps = 1e-6; end
if nargin < 5, tol = eps; end
if nargin < 6, maxit = []; end
if nargin < 7, objover = []; end
if nargin < 8, cors = []; end
if nargin < 9, opts = []; end

tprecomp = tic;
if isempty(cors)
    coropts = opts;
    coropts.corrections = 1;
    [cors,objover] = surfermat(surferobj,kern,eps,coropts);
end

nrows = size(cors,1);
sysapply = @(x) surfermatapply(surferobj, kern, x, eps, objover, cors, opts);

if isempty(maxit), maxit = min(nrows,200); end

nrhs = size(rhs,2);

sol    = zeros(nrows, nrhs);
flag   = zeros(nrhs, 1);
relres = zeros(nrhs, 1);
iter   = zeros(nrhs, 2);
resvec = cell(1, nrhs);
times  = zeros(nrhs+1, 1);
times(1) = toc(tprecomp);

for i = 1:nrhs
    tsolve = tic;
    if nargout < 2
        sol(:,i) = gmres(sysapply, rhs(:,i), [], tol, maxit);
    else
        [sol(:,i), flag(i), relres(i), iter(i,:), resvec{i}] = ...
            gmres(sysapply, rhs(:,i), [], tol, maxit);
    end
    times(1+i) = toc(tsolve);
end

iter = iter(:,2);
end
