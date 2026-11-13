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
%                     'gmres'  - unrestarted gmres with surfermatapply
%                     'dense'  - direct solve with the dense matrix
%                     'rskelf' - direct solve with FLAM's rskelf to eps,
%                                no oversampling
%                     'ifmm'   - gmres with FLAM's ifmm to eps, with
%                                oversampling
%                   The FLAM solvers need FLAM on the path, and do not
%                   take objover, cors or opts.l2scale. Quadrature options
%                   go in opts.flamfuns.quadopts.
%       opts.precond = ('none') preconditioner for 'gmres' and 'ifmm',
%                     built from the matrix without oversampling:
%                     'rskelf'      - rskelf factorization
%                     'blocklu'     - LU of each surfer's diagonal block
%                     'blockrskelf' - rskelf of each diagonal block
%                     or a function handle M, with M(x) approximating
%                     the inverse of the matrix applied to x
%                   The rskelf preconditioners need FLAM on the path.
%       opts.precond_tol = (sqrt(eps)) rank tolerance of the rskelf
%                     preconditioners
%       opts.occ    = (64) FLAM leaf occupancy
%       opts.flamopts = options for rskelf/ifmm, e.g. symm or
%                     proxybylevel. For 'ifmm', store defaults to 'a'.
%       opts.flamfuns = options for byindex.flamfuns
%       opts.unif_nover = for 'gmres', see surfermat. Defaults to 1 for a
%                     single kern, so each apply uses one FMM.
%       opts.checkcors = (true) warn if a user-supplied cors does not
%                     appear to be a correction
%
% Output:
%   sol    - solution, nrows x nrhs
%   flag   - nrhs x 1, solver convergence flag for each right-hand side
%   relres - nrhs x 1, relative residual for each right-hand side
%   iter   - nrhs x 1, gmres iteration count (0 for dense and rskelf)
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
if ~any(strcmp(solver,{'gmres','dense','rskelf','ifmm'}))
    error('surfermatsolve: unsupported solver ''%s''', solver);
end
precond = get_opt(opts, 'precond', 'none');
if isa(precond, 'function_handle')
    precname = 'function handle';
else
    precond = lower(precond);
    precname = precond;
    if ~any(strcmp(precond,{'none','rskelf','blocklu','blockrskelf'}))
        error('surfermatsolve: unsupported preconditioner ''%s''', precond);
    end
end
if ~strcmp(precname,'none') && ~any(strcmp(solver,{'gmres','ifmm'}))
    error('surfermatsolve: solver ''%s'' does not take a preconditioner', solver);
end

% FLAM routines needed by the solver and preconditioner
flamneeds = {};
if strcmp(solver,'rskelf'), flamneeds = {'rskelf','rskelf_sv','rskelf_mv'}; end
if strcmp(solver,'ifmm'),   flamneeds = {'ifmm','ifmm_mv'}; end
if any(strcmp(precname,{'rskelf','blockrskelf'}))
    flamneeds = [flamneeds, {'rskelf','rskelf_sv'}];
end
for f = flamneeds
    if exist(f{1}, 'file') ~= 2
        error('surfermatsolve: solver ''%s'' with preconditioner ''%s'' needs FLAM on the path', ...
            solver, precname);
    end
end
if any(strcmp(solver,{'rskelf','ifmm'}))
    if ~isempty(objover) || ~isempty(cors)
        error('surfermatsolve: solver ''%s'' does not use objover or cors', solver);
    end
    if get_opt(opts, 'l2scale', false)
        error('surfermatsolve: solver ''%s'' does not support opts.l2scale', solver);
    end
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
    % With a single kernel, uniform oversampling lets each apply use one FMM
    if numel(kern) == 1 && ~isfield(opts, 'unif_nover'), opts.unif_nover = 1; end
    if isempty(cors)
        coropts = opts;
        coropts.corrections = 1;
        [cors,objover] = surfermat(surferobj,kern,eps,coropts);
    end
    sysapply = @(x) surfermatapply(surferobj, kern, x, eps, objover, cors, opts);
    if checkcors, check_cors(surfers, kern, sysapply); end
    M = make_precond(precond, surfers, kern, eps, opts);

    nrows = size(cors,1);
    if isempty(maxit), maxit = min(nrows,200); end
    sol = zeros(nrows, nrhs);
    times(1) = toc(tprecomp);

    for i = 1:nrhs
        tsolve = tic;
        if nargout < 2
            sol(:,i) = gmres(sysapply, rhs(:,i), [], tol, maxit, M);
        else
            [sol(:,i), flag(i), relres(i), iter(i,:), resvec{i}] = ...
                gmres(sysapply, rhs(:,i), [], tol, maxit, M);
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

case 'rskelf'
    [F, funs] = flam_rskelf(surfers, kern, eps, eps, opts);
    times(1) = toc(tprecomp);

    sol = zeros(size(funs.x, 2), nrhs);
    for i = 1:nrhs
        tsolve = tic;
        sol(:,i) = rskelf_sv(F, rhs(:,i));
        times(1+i) = toc(tsolve);
        % residual of the compressed matrix
        relres(i) = norm(rskelf_mv(F, sol(:,i)) - rhs(:,i)) / norm(rhs(:,i));
        resvec{i} = relres(i);
    end

case 'ifmm'
    fopts = get_opt(opts, 'flamfuns', []);
    fopts.type = 'ifmm';
    flamopts = get_opt(opts, 'flamopts', []);
    if ~isfield(flamopts, 'store'), flamopts.store = 'a'; end
    if isfield(flamopts, 'symm'), fopts.symm = flamopts.symm; end
    funs = byindex.flamfuns([surfers{:}], kern, eps, [], fopts);
    F = ifmm(funs.Afun, funs.rx, funs.cx, get_opt(opts, 'occ', 64), eps, ...
        funs.pxyfun, flamopts);
    sysapply = @(x) ifmm_mv(F, funs.xinterp*x, funs.Afun) + funs.Q*x;

    M = make_precond(precond, surfers, kern, eps, opts);

    nrows = size(funs.rx, 2);
    if isempty(maxit), maxit = min(nrows, 200); end
    sol = zeros(size(funs.Q, 2), nrhs);
    times(1) = toc(tprecomp);

    for i = 1:nrhs
        tsolve = tic;
        [sol(:,i), flag(i), relres(i), iter(i,:), resvec{i}] = ...
            gmres(sysapply, rhs(:,i), [], tol, maxit, M);
        times(1+i) = toc(tsolve);
    end
end

iter = iter(:,2);
end


function [F, funs] = flam_rskelf(surfers, kern, eps, tol, opts)
% rskelf factorization of the matrix without oversampling, to rank
% tolerance tol
fopts = get_opt(opts, 'flamfuns', []);
fopts.type = 'rskelf';
flamopts = get_opt(opts, 'flamopts', []);
if isfield(flamopts, 'symm'), fopts.symm = flamopts.symm; end
funs = byindex.flamfuns([surfers{:}], kern, eps, [], fopts);
F = rskelf(funs.Afun, funs.x, get_opt(opts, 'occ', 64), tol, funs.pxyfun, ...
    flamopts);
end


function M = make_precond(precond, surfers, kern, eps, opts)
% Preconditioner from the matrix without oversampling, [] for none
if isa(precond, 'function_handle')
    M = precond;
    return
end
ptol = get_opt(opts, 'precond_tol', sqrt(eps));
switch precond
case 'none'
    M = [];
case 'rskelf'
    F = flam_rskelf(surfers, kern, eps, ptol, opts);
    M = @(x) rskelf_sv(F, x);
case {'blocklu', 'blockrskelf'}
    % One solve per surfer's diagonal block
    ns = numel(surfers);
    solves = cell(ns, 1);
    inds   = cell(ns, 1);
    iloc = 0;
    for i = 1:ns
        if numel(kern) == 1, kii = kern; else, kii = kern(i,i); end
        n = kii.opdims(1)*surfers{i}.npts;
        inds{i} = iloc + (1:n);
        iloc = iloc + n;
        if strcmp(precond, 'blocklu')
            dA = decomposition(surfermat(surfers{i}, kii, eps, struct('ifoversamp', 0)));
            solves{i} = @(x) dA \ x;
        else
            F = flam_rskelf(surfers(i), kii, eps, ptol, opts);
            solves{i} = @(x) rskelf_sv(F, x);
        end
    end
    M = @(x) block_solve(x, solves, inds);
end
end


function y = block_solve(x, solves, inds)
y = zeros(size(x), 'like', x);
for i = 1:numel(solves)
    y(inds{i},:) = solves{i}(x(inds{i},:));
end
end


function v = get_opt(opts, name, default)
if isfield(opts, name), v = opts.(name); else, v = default; end
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
