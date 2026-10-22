% SURFWAVE_KERNELSTEST  Surface wave kernels (kernel3d.capillary, iceflex,
% gravity) on a flat disk.
%
% Each kernel is a kernel3d.radgreen kernel on a SURFWAVE Green's function.
% Checks:
%   - K.eval agrees with the Green's function it was built on;
%   - sums of kernels of one type are one kernel, with the summed eval,
%     and match the weight vector form;
%   - the 'rmax' (radpp fit) kernels agree with the exact ones;
%   - the volume potential of a smooth bump at on-surface and boundary
%     targets agrees with a reference at a higher order and tolerance.

% run ../startup.m

rng(42);

nu = 0.33;

[rts_c, ejs_c] = surfwave.capillary.find_roots_capillary(1.3, -0.7);
[rts_f, ejs_f] = surfwave.flex.find_roots_flex(23, 0.5);
g = 1.7;

kerns = { ...
  'cap s gs',          kernel3d('capillary', 's',   rts_c, ejs_c, 'gs'); ...
  'cap s gphi',        kernel3d('capillary', 's',   rts_c, ejs_c, 'gphi'); ...
  'cap s s3d',         kernel3d('capillary', 's',   rts_c, ejs_c, 's3d'); ...
  'cap sp gs',         kernel3d('capillary', 'sp',  rts_c, ejs_c, 'gs'); ...
  'cap sp gphi',       kernel3d('capillary', 'sp',  rts_c, ejs_c, 'gphi'); ...
  'cap lap gphi',      kernel3d('capillary', 'lap', rts_c, ejs_c, 'gphi'); ...
  'cap lap gs+gphi',   kernel3d('capillary', 'lap', rts_c, ejs_c, [0.3 0.5]); ...
  'flex s gs',         kernel3d('iceflex', 's', rts_f, ejs_f, 'gs'); ...
  'flex s gphi',       kernel3d('iceflex', 's', rts_f, ejs_f, 'gphi'); ...
  'flex s s3d',        kernel3d('iceflex', 's', rts_f, ejs_f, 's3d'); ...
  'flex s gphi_bilap', kernel3d('iceflex', 's', rts_f, ejs_f, 'gphi_bilap'); ...
  'flex free gs',      kernel3d('iceflex', 'free_bc', rts_f, ejs_f, 'gs', 'nu', nu); ...
  'flex free gphi',    kernel3d('iceflex', 'free_bc', rts_f, ejs_f, 'gphi', 'nu', nu); ...
  'grav s gs',         kernel3d('gravity', 's',    g, 'gs'); ...
  'grav s gphi',       kernel3d('gravity', 's',    g, 'gphi'); ...
  'grav grad gs',      kernel3d('gravity', 'grad', g, 'gs'); ...
};

test_eval(kerns, rts_c, ejs_c, rts_f, ejs_f, g);
test_sums(rts_c, ejs_c, rts_f, ejs_f, nu);
test_rmax(rts_c, ejs_c, rts_f, ejs_f, nu, g);
test_volume(kerns);
test_gravity_vol(g);

fprintf('surfwave_kernelsTest: all checks passed.\n');


function check(name, err, tol)
fprintf('%-34s rel err %5.2e\n', name, err);
assert(err < tol, '%s: error %5.2e exceeds %5.2e', name, err, tol);
end


function [s, t] = points()
% sources and targets in the plane, with target normals and tangents
s = []; s.r = [randn(2, 6); zeros(1, 6)];
th = 2*pi*rand(1, 5);
t = []; t.r = [3*cos(th); 3*sin(th); 0*th];
t.n = [cos(th); sin(th); 0*th];
t.d = [-3*sin(th); 3*cos(th); 0*th];
t.d2 = -t.r;
end


function test_eval(kerns, rts_c, ejs_c, rts_f, ejs_f, g)
% K.eval against the Green's functions
[s, t] = points();
ref = containers.Map();
xs = s.r(1:2,:);  xt = t.r(1:2,:);
nx = repmat(t.n(1,:).', 1, size(xs, 2));
ny = repmat(t.n(2,:).', 1, size(xs, 2));
dn = @(gr) gr(:,:,1).*nx + gr(:,:,2).*ny;

[v, gr, h] = surfwave.capillary.gshelm(rts_c, ejs_c, xs, xt);
ref('cap s gs') = v;
ref('cap sp gs') = dn(gr);
[v, gr, h] = surfwave.capillary.gphihelm(rts_c, ejs_c, xs, xt);
ref('cap s gphi') = v;
ref('cap sp gphi') = dn(gr);
ref('cap lap gphi') = h(:,:,1) + h(:,:,3);
ref('cap s s3d') = surfwave.flex.lap3dkern(xs, xt) + ...
                    surfwave.capillary.s3dgphi(rts_c, ejs_c, xs, xt);
ref('flex s gs') = surfwave.flex.gsflex(rts_f, ejs_f, xs, xt);
ref('flex s gphi') = surfwave.flex.gphiflex(rts_f, ejs_f, xs, xt);
ref('flex s s3d') = surfwave.flex.lap3dkern(xs, xt) + ...
                      surfwave.flex.s3dgphiflex(rts_f, ejs_f, xs, xt);
[~, ~, ~, ~, bl] = surfwave.flex.gphiflex(rts_f, ejs_f, xs, xt);
ref('flex s gphi_bilap') = bl;
ref('grav s gs') = surfwave.gravity.gsgrav(g/2, 1, xs, xt);
ref('grav s gphi') = surfwave.gravity.gsgrav(g/2, 1/g, xs, xt);
[~, gr] = surfwave.gravity.gsgrav_der(g/2, 1, xs, xt);
G = zeros(2*size(xt, 2), size(xs, 2));
G(1:2:end,:) = gr(:,:,1);  G(2:2:end,:) = gr(:,:,2);
ref('grav grad gs') = G;

for k = 1:size(kerns, 1)
    [name, K] = kerns{k,:};
    if ~isKey(ref, name), continue; end
    r = ref(name);
    A = K.eval(s, t);
    check(['eval ' name], norm(A - r, 'fro')/norm(r, 'fro'), 1e-12);
end
end


function test_sums(rts_c, ejs_c, rts_f, ejs_f, nu)
% a*K1 + b*K2 of one type is one radgreen kernel with the summed eval
[s, t] = points();
a = 0.3 - 0.2i;  b = 0.5;
K1 = kernel3d.capillary('sp', rts_c, ejs_c, 'gs');
K2 = kernel3d.capillary('sp', rts_c, ejs_c, 'gphi');
K  = a*K1 + b*K2;
assert(strcmp(K.name, 'radgreen'), 'capillary sum is not a radgreen kernel');
A = K.eval(s, t);  B = a*K1.eval(s, t) + b*K2.eval(s, t);
check('sum cap sp', norm(A - B, 'fro')/norm(B, 'fro'), 1e-13);
C = kernel3d.capillary('sp', rts_c, ejs_c, [a b]).eval(s, t);
check('weights cap sp', norm(C - B, 'fro')/norm(B, 'fro'), 1e-13);

K1 = kernel3d.iceflex('free_bc', rts_f, ejs_f, 'gs', 'nu', nu);
K2 = kernel3d.iceflex('free_bc', rts_f, ejs_f, 'gphi', 'nu', nu);
K  = a*K1 - b*K2;
assert(strcmp(K.name, 'radgreen'), 'iceflex sum is not a radgreen kernel');
A = K.eval(s, t);  B = a*K1.eval(s, t) - b*K2.eval(s, t);
check('sum flex free_bc', norm(A - B, 'fro')/norm(B, 'fro'), 1e-13);
C = kernel3d.iceflex('free_bc', rts_f, ejs_f, [a -b], 'nu', nu).eval(s, t);
check('weights flex free_bc', norm(C - B, 'fro')/norm(B, 'fro'), 1e-13);
end


function test_rmax(rts_c, ejs_c, rts_f, ejs_f, nu, g)
% 'rmax' kernels (eval from a radpp fit) against the exact (radgreen) ones,
% singly and summed
[s, t] = points();
ro = {'opts', struct('eps', 1e-12)};
cases = { ...
  'cap s gs',      @(varargin) kernel3d.capillary('s', rts_c, ejs_c, 'gs', varargin{:}); ...
  'cap sp gphi',   @(varargin) kernel3d.capillary('sp', rts_c, ejs_c, 'gphi', varargin{:}); ...
  'cap s s3d',     @(varargin) kernel3d.capillary('s', rts_c, ejs_c, 's3d', varargin{:}); ...
  'flex free gs',  @(varargin) kernel3d.iceflex('free_bc', rts_f, ejs_f, 'gs', 'nu', nu, varargin{:}); ...
  'grav s gphi',   @(varargin) kernel3d.gravity('s', g, 'gphi', varargin{:}); ...
};
for k = 1:size(cases, 1)
    [name, mk] = cases{k,:};
    A = mk('rmax', 20, ro{:}).eval(s, t);
    B = mk().eval(s, t);
    check(['rmax ' name], norm(A - B, 'fro')/norm(B, 'fro'), 1e-10);
end
a = 0.3 - 0.2i;  b = 0.5;
K1 = kernel3d.capillary('s', rts_c, ejs_c, 'gs',  'rmax', 20, ro{:});
K2 = kernel3d.capillary('s', rts_c, ejs_c, 'gphi', 'rmax', 20, ro{:});
K  = a*K1 + b*K2;
assert(strcmp(K.name, 'radpp'), 'capillary rmax sum is not a radpp kernel');
A = K.eval(s, t);  B = a*K1.eval(s, t) + b*K2.eval(s, t);
check('rmax sum cap s', norm(A - B, 'fro')/norm(B, 'fro'), 1e-13);
end


function test_volume(kerns)
% volume potential of a smooth bump on the unit disk, on-surface and
% boundary targets, against a reference at a higher order and tolerance
S  = geometries.disk([], [], [1,1,1], 8);
Sr = geometries.disk([], [], [1,1,1], 12);
bump = @(r) exp(-2*sum(r(1:2,:).^2, 1)).';

th = linspace(0, 2*pi, 9); th = th(1:end-1);
tb = []; tb.r = [cos(th); sin(th); 0*th];
tb.n = tb.r;  tb.d = [-sin(th); cos(th); 0*th];  tb.d2 = -tb.r;
[rr, tt] = ndgrid([0.1 0.45 0.8], [0.3 2.4 4.4]);
tv = []; tv.r = [rr(:).'.*cos(tt(:).'); rr(:).'.*sin(tt(:).'); 0*rr(:).'];
tv.n = [cos(tt(:).'); sin(tt(:).'); 0*rr(:).'];
tv.d = [-sin(tt(:).'); cos(tt(:).'); 0*rr(:).'];  tv.d2 = -tv.r;
targs = {'v2v', tv; 'v2b', tb};

for k = 1:size(kerns, 1)
    [name, K] = kerns{k,:};
    for j = 1:size(targs, 1)
        % grad G_S ~ 1/r^2: the gradient of the potential is log
        % singular at the edge of the disk, so no boundary targets
        if strcmp(name, 'grav grad gs') && j == 2, continue; end
        T  = add_patch_ids(S,  targs{j,2});
        Tr = add_patch_ids(Sr, targs{j,2});
        u  = surferkerneval(S,  K, bump(S.r),  T,  1e-8);
        ur = surferkerneval(Sr, K, bump(Sr.r), Tr, 1e-11);
        check(sprintf('%s (%s)', name, targs{j,1}), ...
              norm(u - ur, inf)/norm(ur, inf), 1e-5);
    end
end
end


function t = add_patch_ids(S, t)
[~, pid, uvs, ~, flags] = get_closest_pts(S, t);
assert(all(flags == 0), 'closest point search failed');
t.patch_id = pid(:);
t.uvs_targ = uvs;
end


function test_gravity_vol(g)
% the eval only 'vol' kernel
[s, t] = points();
K = kernel3d.gravity('vol', g, [], 'coefs', [0.7 1.1]);
A = K.eval(s, t);
B = surfwave.gravity.kern(g, s, t, 'vol', [0.7 1.1]);
check('grav vol eval', norm(A - B, 'fro')/norm(B, 'fro'), 1e-14);
assert(isempty(K.getquad), 'grav vol: no near quadrature expected');
end
