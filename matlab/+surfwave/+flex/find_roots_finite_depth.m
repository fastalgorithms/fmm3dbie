function [rts,ejs] = find_roots_finite_depth(alpha,gamma,depth,nkappa)
%%%%%
%
% Finite-depth version of find_roots. Finds roots of the dispersion function
%
%       D(k) = (alpha*k^4 + gamma)*k*tanh(k*depth) - 1
%
% (depth -> inf gives alpha*k^5 + gamma*k - 1, the polynomial of find_roots).
% D is even in k and has infinitely many roots; returned are those with
% Re(k) > 0 plus the first nkappa roots on the positive imaginary axis:
%
%   rts(1)       = k0            real, propagating wavenumber (> deep-water root)
%   rts(2:3)     = k1, conj(k1)  complex pair (-> q1, conj(q1) of find_roots as depth -> inf)
%   rts(4:end)   = 1i*kappa_n    evanescent modes, kappa_n in ((n-1/2)pi/depth, n pi/depth)
%
% (-rts are roots as well.)  Sorted by abs(angle(rts)) like find_roots.
%
% Returns:
% - rts: roots (column vector)
% - ejs: 1./D'(rts), the analogue of the partial fraction coefficients
%        ejs = 1./f'(rts) of find_roots (they agree as depth -> inf)
%
% nkappa (optional, default 10): number of imaginary roots returned
%
% gamma < 0 (plate inertia exceeds buoyancy) is allowed as long as the root structure is unchanged
% (single real deep-water root, one imaginary root per interval); an error is raised otherwise.
%
%%%%%

    if nargin < 4
        nkappa = 10;
    end

    P   = @(k) (alpha*k.^4 + gamma).*k;
    dP  = @(k) 5*alpha*k.^4 + gamma;
    D   = @(k) P(k).*tanh(k*depth) - 1;
    dD  = @(k) dP(k).*tanh(k*depth) + P(k)*depth.*sech(k*depth).^2;

    % deep-water roots (seeds)
    q = roots([alpha 0 0 0 gamma -1]);
    kinf = real(q(abs(imag(q)) < 1e-12*abs(q) & real(q) > 0));
    q1 = q(real(q) > 0 & imag(q) > 0);
    if numel(kinf) ~= 1 || numel(q1) ~= 1
        error('find_roots_finite_depth: unexpected deep-water root structure (gamma = %g too negative?)', gamma);
    end
    if gamma < 0 && (-gamma/alpha)^(1/4) >= pi/(2*depth)
        error(['find_roots_finite_depth: gamma = %g too negative for depth %g ((alpha k^4 + gamma) k < 0 ' ...
               'beyond the first evanescent-root interval)'], gamma, depth);
    end

    % real root k0 > kinf (tanh < 1): bracket and solve
    lo = kinf*(1 - 1e-10);                 % D(lo) < 0 even when tanh(kinf*depth) rounds to 1
    hi = 2*kinf;
    while D(hi) < 0, hi = 2*hi; end
    k0 = fzero(D, [lo hi]);
    for it = 1:2, k0 = k0 - D(k0)/dD(k0); end

    % first-quadrant complex root: Newton from the deep-water root and a few other seeds
    seeds = [q1; 1.2*kinf*exp(1i*pi/3); 2*kinf*exp(1i*pi/4); 0.7*kinf*exp(1.2i)];
    k1 = NaN;
    for s = 1:numel(seeds)
        z = seeds(s);
        for it = 1:100
            dz = D(z)/dD(z); z = z - dz;
            if abs(dz) < 1e-15*abs(z), break; end
        end
        if abs(D(z)) < 1e-10 && real(z) > 0 && imag(z) > 0
            k1 = z; break
        end
    end
    if isnan(k1)
        error('find_roots_finite_depth: no first-quadrant complex root found');
    end

    % imaginary roots i*kappa: tan(kappa*depth) = -1/Pt(kappa), Pt = (alpha kappa^4 + gamma) kappa.
    % write kappa*depth = n*pi + t, t in (-pi/2, 0):  t + atan(1/Pt(kappa)) = 0  (one root per n)
    kappa = zeros(nkappa,1);
    for n = 1:nkappa
        g = @(t) t + atan(1./P((n*pi + t)/depth));
        t = fzero(g, [-pi/2*(1 - 1e-15), 0]);
        kappa(n) = (n*pi + t)/depth;
    end

    rts = [k0; k1; conj(k1); 1i*kappa];
    ejs = 1./dD(rts);

    [~, sort_id] = sort(abs(angle(rts)));
    rts = rts(sort_id);
    ejs = ejs(sort_id);

end
