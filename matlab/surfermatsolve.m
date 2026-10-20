function [sol, flag, relres, iter, resvec, times] = surfermatsolve(surferobj, kern, rhs, eps, tol, maxit, objover, cors, opts)
%SURFERMATSOLVE Solve the integral equation system for given kernel and
% surfer array, using the solver selected by opts.solver.
%
% Syntax: sol = surfermatsolve(surferobj, kern, rhs)
%         [sol,flag,relres,iter,resvec,times] = surfermatsolve(surferobj, kern, rhs, ...
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
%   opts    - options structure, also passed to surfermatapply
%       opts.solver = ('gmres') solver to use. Supported:
%                     'gmres' - unrestarted gmres with surfermatapply
%                     'dense' - direct solve with the dense matrix
%       opts.checkcors = (true) warn if a user-supplied cors does not
%                     appear to be a correction
%
% Output:
%   sol    - solution, nrows x nrhs
%   flag   - nrhs x 1, solver convergence flag for each right-hand side
%   relres - nrhs x 1, relative residual for each right-hand side
%   iter   - nrhs x 1, gmres iteration count (0 for dense)
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

solver = 'gmres';
if isfield(opts,'solver'), solver = lower(opts.solver); end
if ~any(strcmp(solver,{'gmres','dense'}))
    error('surfermatsolve: unsupported solver ''%s''', solver);
end

checkcors = ~isempty(cors);
if isfield(opts,'checkcors'), checkcors = checkcors && opts.checkcors; end

surfers = surferobj;
if ~iscell(surfers), surfers = num2cell(surferobj(:).',1); end

nrhs   = size(rhs,2);
flag   = zeros(nrhs, 1);
relres = zeros(nrhs, 1);
iter   = zeros(nrhs, 2);
resvec = cell(1, nrhs);
times  = zeros(nrhs+1, 1);

tprecomp = tic;
switch solver
case 'gmres'
    if isempty(cors)
        coropts = opts;
        coropts.corrections = 1;
        [cors,objover] = surfermat(surferobj,kern,eps,coropts);
    end
    sysapply = @(x) surfermatapply(surferobj, kern, x, eps, objover, cors, opts);
    if checkcors, check_cors(surfers, kern, sysapply); end

    nrows = size(cors,1);
    if isempty(maxit), maxit = min(nrows,200); end
    sol = zeros(nrows, nrhs);
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

case 'dense'
    matopts = opts;
    matopts.corrections = 0;
    matopts.nonsmoothonly = 0;
    if isempty(cors)
        A = surfermat(surferobj, kern, eps, matopts);
    else
        % smooth rule + supplied corrections
        matopts.forcesmooth = 1;
        matopts.ifreturnovers = 0;
        [A, novers] = surfermat(surferobj, kern, eps, matopts);
        if ~isempty(objover) && ~overs_match(objover, novers, surfers)
            error('surfermatsolve: objover orders do not match surfermat');
        end
        A = A + cors;
        if checkcors, check_cors(surfers, kern, @(x) A*x); end
    end
    dA = decomposition(A);
    times(1) = toc(tprecomp);

    sol = zeros(size(A,1), nrhs);
    for i = 1:nrhs
        tsolve = tic;
        sol(:,i) = dA \ rhs(:,i);
        times(1+i) = toc(tsolve);
        relres(i) = norm(A*sol(:,i) - rhs(:,i)) / norm(rhs(:,i));
        resvec{i} = relres(i);
    end
end

iter = iter(:,2);
end


function ok = overs_match(objover, novers, surfers)
% Compare oversampling orders in objover to novers from surfermat.
ns = numel(surfers);
ok = true;
for i = 1:ns
    for j = 1:ns
        if iscell(objover) && numel(objover) == 2 && iscell(objover{1})
            so = objover{1};
            if ~iscell(so), so = num2cell(so); end
            if numel(so) == ns, so = so{j}; else, so = so{i,j}; end
            if all(isnan(novers{i,j}))
                nij = surfers{j}.norders;
            else
                nij = novers{i,j};
            end
            ok = isequal(so.norders(:), nij(:));
        elseif iscell(objover)
            ok = isequaln(objover{i,j}(:), novers{i,j}(:));
        else
            ok = isequaln(objover(:), novers{i,j}(:));
        end
        if ~ok, return; end
    end
end
end


function check_cors(surfers, kern, sysapply)
% Apply the operator to a smooth density and warn if the result is not
% smooth, which indicates cors is not a correction to the smooth rule.

ns = numel(surfers);

if numel(kern) == 1
    opdims = repmat(kern.opdims(:), 1, ns, ns);
else
    opdims = reshape([kern.opdims], 2, ns, ns);
end

dens = [];
for j = 1:ns
    r = surfers{j}.r;
    f = cos(r(1,:) + 2*r(2,:) - r(3,:));
    dens = [dens; reshape(repmat(f, opdims(2,1,j), 1), [], 1)];
end

pot = sysapply(dens);

tol = 5e-4;
iloc = 0;
for i = 1:ns
    nd = opdims(1,i,1);
    npts = surfers{i}.npts;
    poti = reshape(pot(iloc+(1:nd*npts)), nd, npts);
    iloc = iloc + nd*npts;
    scale = max(abs(poti(:)));
    if scale == 0, continue; end
    err = max(surf_fun_error(surfers{i}, poti), [], 'all') / scale;
    if err > tol
        warning('surfermatsolve:cors', ['Applying the operator gives a non-smooth result %d ' ...
            '(rel. tail %.1e). Was cors built with corrections=true?'], i, err);
        return
    end
end
end
