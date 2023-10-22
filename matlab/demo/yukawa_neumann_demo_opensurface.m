% Add fmm3dbie path
addpath(genpath('~/git/fmm3dbie/matlab'))
run ~/git/surface-hps/setup.m

rng(0)
n = 10;                     % Number of nodes per dimension for each mesh element
nref = 4;                  % Number of levels of uniform mesh refinement

m = 18.0;                     % mass term
E = 15.0;                   % energy term
om = sqrt(m^2 - E^2); 

xyzsrc = [0.2; 0.3; -0.2]; % Source location for evaluating Green's function

xm = max(2,norm(xyzsrc));  % Bounding box for perturbation and b/c
tol = 1e-7;
dm = xm + log(100/tol)/om; % Bounding box size including buffer region

rect = [-xm xm -xm xm];    % Bounding box for HPS
zk = E;                    % Helmholtz wavenumber
eta = zk;                  % Impedance parameter

% Incoming plane wave
uincf    = @(x,y,z) exp(1i*zk*x);
uincf_dx = @(x,y,z) 1i*zk*exp(1i*zk*x);
uincf_dy = @(x,y,z) 0*x;

%% Make a compactly supported perturbation of a flat surface

ifflat = 0;

% Build a uniformly refined flat square
dom = surfacemesh.square(n, nref, rect);
x = dom.x;
y = dom.y;
z = dom.z;
if ~ifflat
    f = @(x,y) 0.5+0*x;
    b = @(x,y,z) f(x,y).*exp(-20*(x.^2+y.^2));  
    for k = 1:length(dom)
        z{k} = b(x{k},y{k},z{k});
    end
    dom = surfacemesh(x, y, z);
end




% Collapse the arrays
npols = n*n;
npatches = size(x,1);
npt = npatches*npols;

[S, dom] = surfer.surfacemesh_to_surfer(dom);

x = dom.x;
y = dom.y;
z = dom.z;



xx = zeros(npt,1);
yy = zeros(npt,1);
zz = zeros(npt,1);
isrcinds = cell(npatches,1);

for i=1:npatches
    istart = (i-1)*npols + 1;
    iend = i*npols;
    iind = istart:iend;
    
    xx(iind) = x{i}(:);
    yy(iind) = y{i}(:);
    zz(iind) = z{i}(:);
    isrcinds{i} = iind;
end

% Find interior and exterior indices, truncate at patch boundaries to
% avoid any interpolation issues

iind_out = find(xx > xm | xx < -xm | yy > xm | yy < -xm);
ipat = unique(ceil(iind_out/npols));
iind_out = horzcat(isrcinds{ipat});

figure
clf
plot(xx(iind_out),yy(iind_out),'k.'); hold on;

iind_all = 1:npt;

iind_in = setdiff(iind_all,iind_out);

plot(xx(iind_in), yy(iind_in),'r.');

%% get layer potential correction
zpars = 1j*om;
opts = [];
opts.int_rep = 's';

tic, Q = helm3d.neumann.get_quadrature_correction(S, zpars, tol, S, opts); toc

%% Now test accuracy of computed quadrature

% generate boundary data
srcinfo = [];
srcinfo.r = xyzsrc;
ubdry = real(helm3d.kern(1j*om, srcinfo, S, 'sprime'));

% Iterate

matvec = @(x) eval_yukawa_sp(S, om, Q, x, tol);
[sol, ~, relres, iter] = gmres(matvec, ubdry, [], 1e-8, 100);

%% Test solution in the upper half plane
xyztarg = [0.21; 0.22; 0.14];
targinfo = [];
targinfo.r = xyztarg;
potex = real(helm3d.kern(1j*om, srcinfo, targinfo, 's'));

udata = real(helm3d.kern(1j*om, S, targinfo, 's')).';
pot = sum(udata.*sol.*S.wts);

err = norm(pot-potex)/abs(potex);
fprintf('error in potential = %d\n', err)



function [u] = eval_yukawa_sp(S, om, Q, sig, tol)
    zk = 1j*om;
    siguse = complex(sig);
    opts = [];
    opts.int_rep = 's';
    u = helm3d.neumann.eval(S, zk, siguse, tol, S, Q, opts);
    u = real(u(:));
    u = u - sig(:)/2;
end



function [u] = eval_yukawa_part(S, om, m, Q, sig, tol)
    zk = 1j*om;
    zpars = [zk; 1; 0];
    siguse = complex(sig);
    u = helm3d.dirichlet.eval(S, zpars, siguse, tol, S, Q, opts);
    u = -2*m*real(u(:));
    u = u + sig(:);
end



