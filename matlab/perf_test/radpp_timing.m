% Time helm3d against kernel3d.radpp for 's' and 'cp': the fit (build),
% the near quadrature and eval.

% run ../startup.m

zk = 1.1;  eps = 1e-7;  coefs = [1 1];
S = geometries.sphere(1, 4, [0;0;0], 6, 1);
g = @(s,t) helm3d.green(zk, s, t);

src.r  = S.r;  src.n  = S.n;
targ.r = 1.5*S.r(:,1:2:end);  targ.n = S.n(:,1:2:end);   % off surface
opts.rsc = getnear(S, S);

fprintf('\nradpp timing, %d points\n', S.npts);
for type = {'s', 'cp'}
    type = type{1};
    Ke = kernel3d.helm3d(type, zk, coefs(:));

    tic;  Kr = kernel3d.radpp(g, [1e-9 3], type, 'coefs', coefs);
    tb = toc;

    tic;  Qe = Ke.getquad(S, eps, S, opts);  tqe = toc;
    tic;  Qr = Kr.getquad(S, eps, S, opts);  tqr = toc;

    tic;  Ae = Ke.eval(src, targ);  tee = toc;
    tic;  Ar = Kr.eval(src, targ);  ter = toc;

    fprintf(['%-3s build %6.3fs | getquad helm3d %7.3fs radpp %7.3fs ' ...
             '(err %.1e) | eval helm3d %6.3fs radpp %6.3fs (err %.1e)\n'], ...
            type, tb, tqe, tqr, norm(Qr-Qe,'fro')/norm(Qe,'fro'), ...
            tee, ter, norm(Ar-Ae,'fro')/norm(Ae,'fro'));
end
