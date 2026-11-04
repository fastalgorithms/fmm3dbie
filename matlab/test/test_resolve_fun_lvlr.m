% Test spatial level restriction and global refinement limits.
% Assumes pwd is the directory this script is in.

run('../startup.m');

% Graded and sharply refined meshes, including stretched triangles

nrefmax = 4;
iptypes = [1,11,12];
for itype = 1:numel(iptypes)
    iptype = iptypes(itype);
    Sgraded = geometries.square([2,2], 4, iptype);
    for k = 1:3
        bound = -1 + 2^(1-k);
        isplit = find(Sgraded.cms(1,:) < bound & Sgraded.cms(2,:) < bound);
        Sgraded = split_patches(Sgraded, isplit);
    end
    Ssharp = geometries.square([2,2], 4, iptype);
    for k = 1:3
        isplit = find(Ssharp.cms(1,:) < 0 & Ssharp.cms(2,:) < 0);
        Ssharp = split_patches(Ssharp, isplit);
    end

    for stretch = [1,8]
        theta = 0.37; phi = 0.61;
        Rz = [cos(theta),-sin(theta),0; sin(theta),cos(theta),0; 0,0,1];
        Rx = [1,0,0; 0,cos(phi),-sin(phi); 0,sin(phi),cos(phi)];
        mat = Rx*Rz*diag([1,stretch,1]);
        G = affine_transf(Sgraded, mat);
        B = resolve_fun(G, 1e-6, [], 6, nrefmax);
        assert(B.npatches == G.npatches, ...
            'resolve_fun: a graded, stretched mesh should not grow');
        A = affine_transf(Ssharp, mat);
        C = resolve_fun(A, 1e-6, [], 6, nrefmax);
        nexpected = 91;
        if iptype == 1; nexpected = 176; end
        assert(C.npatches == nexpected, ...
            'resolve_fun: sharp-jump refinement disagrees with the flat polygon check');
        assert(abs(area(C) - area(A))/area(A) < 1e-11, ...
            'resolve_fun: balancing changed surface area');
        fprintf('iptype %d, stretch %d: graded %d -> %d, sharp %d -> %d patches\n', ...
            iptype, stretch, G.npatches, B.npatches, A.npatches, C.npatches);

    end
end

% Nearby parallel sheets should not count as touching patches

Scoarse = geometries.square([1,1], 4, 11);
Ssmall = affine_transf(Scoarse, 0.125*eye(3), [0;0;0.01]);
Ssep = merge([Scoarse,Ssmall]);
Ssep_bal = resolve_fun(Ssep, 1e-6, [], 6, nrefmax);
assert(Ssep_bal.npatches == Ssep.npatches, ...
    'resolve_fun: separated parallel sheets should not trigger balancing');

% Function refinement followed by balancing, and disabling balancing

Sfun = geometries.square([2,2],4,11);
f = @(T) (T.r(1,:) < 0 & T.r(2,:) < 0).*T.r(1,:).^3;
ftol = 2e-4;
Su = resolve_fun(Sfun,ftol,f,2,4,Inf,struct('lvlr','n'));
assert(Su.npatches == 67, 'resolve_fun: disabling balancing changed the refinement');
Sb = resolve_fun(Sfun,ftol,f,2,4,Inf);
assert(Sb.npatches == 91, 'resolve_fun: function refinement should be followed by balancing');
err = surf_fun_error(Sb,f(Sb),Inf,2);
assert(max(err(:)) < ftol, 'resolve_fun: balancing lost function resolution');
Sloose = resolve_fun(Su,ftol,[],2,4,Inf,struct('lvlrfac',4.1));
assert(Sloose.npatches == 76, 'resolve_fun: a looser size ratio should refine fewer patches');

% Two mixed patches, already resolved, with complex multi-row data

Stri = slicesurfer(geometries.square([1,1],4,1),1);
Squad = affine_transf(geometries.square([1,1],5,12),eye(3),[3;0;0]);
Smix = merge([Squad,Stri]);
fmix = @(T) [T.r(1,:) + 2*T.r(2,:); (1+2i)*(T.r(1,:) - T.r(2,:))];
Smix_out = resolve_fun(Smix,1e-11,fmix,1,1,Inf);
assert(Smix_out.npatches == 2, 'resolve_fun: resolved mixed patches should not split');
assert(isequal(sortrows([Smix_out.norders(:),Smix_out.iptype(:)]), ...
    sortrows([Smix.norders(:),Smix.iptype(:)])), ...
    'resolve_fun: mixed patch orders and types were changed');
err = surf_fun_error(Smix_out,fmix(Smix_out),Inf,1);
assert(max(err(:)) < 1e-11, 'resolve_fun: complex multi-row data lost resolution');

% A curved 2:1 transition should tolerate physical radius variation

Scurve = split_patches(geometries.square([2,1],4,11),2);
r = Scurve.r; du = Scurve.du; dv = Scurve.dv;
a = 1.25;
r(3,:) = a*r(1,:).^2;
du(3,:) = 2*a*r(1,:).*du(1,:);
dv(3,:) = 2*a*r(1,:).*dv(1,:);
normals = cross(du,dv,1);
normals = normals./vecnorm(normals,2,1);
Scurve = surfer(Scurve.npatches,4,[r;du;dv;normals],11);
ratio = max(Scurve.rads)/min(Scurve.rads);
assert(ratio > 2.1 && ratio < 2.5, 'Curved regression needs a radius ratio between 2.1 and 2.5');
Bcurve = resolve_fun(Scurve,1e-6,[],2,2,Inf);
assert(Bcurve.npatches == 5, 'resolve_fun: the default should preserve this curved transition');
Bstrict = resolve_fun(Scurve,1e-6,[],2,2,Inf,struct('lvlrfac',2.1));
assert(Bstrict.npatches == 8, 'resolve_fun: the stricter ratio should split the coarse curved patch');
assert(max(abs(Bstrict.r(3,:) - a*Bstrict.r(1,:).^2)) < 1e-12, ...
    'resolve_fun: balancing changed the quadratic surface');

% Pass and global patch-count limits

Szero = resolve_fun(Ssharp, 1e-6, @(T) error('test:UnexpectedEvaluation', ...
    'Zero passes should not evaluate f.'), 6, 0);
assert(isequal(Szero.r, Ssharp.r), 'resolve_fun: zero passes should return unchanged');
Slimit = Smix;
fretired = @(T) (T.patch_id.' == 1).*T.uvs_targ(1,:).^3;
Sretired = geometries.square([2,1],4,11);
checks = {@() resolve_fun(Slimit,1e-12,[],2,0,Inf,struct('npatchmax',1)), ...
    @() resolve_fun(Slimit,1e-12,@(T) exp(T.r(1,:)),2,2,Inf,struct('npatchmax',5)), ...
    @() resolve_fun(Sretired,1e-12,fretired,2,2,Inf,struct('lvlr','n','npatchmax',7)), ...
    @() resolve_fun(Ssharp,1e-6,[],6,4,Inf,struct('npatchmax',Ssharp.npatches)), ...
    @() resolve_fun(Ssharp,1e-6,[],6,1)};
expected = {'surfer:resolve_fun:PatchLimit', 'surfer:resolve_fun:PatchLimit', ...
    'surfer:resolve_fun:PatchLimit', 'surfer:resolve_fun:PatchLimit', ...
    'surfer:resolve_fun:LevelRestrictionLimit'};
for k = 1:numel(checks)
    caught = false;
    try
        checks{k}();
    catch exception
        if ~strcmp(exception.identifier, expected{k}); rethrow(exception); end
        caught = true;
    end
    assert(caught, 'resolve_fun: expected limit error was not raised');
end

Sone = geometries.square([1,1],4,11);
Sfour = resolve_fun(Sone,1e-12,@(T) exp(T.r(1,:)),2,1,Inf, ...
    struct('lvlr','n','npatchmax',4));
assert(Sfour.npatches == 4, 'resolve_fun: exactly npatchmax patches should be allowed');

fprintf('test_resolve_fun_lvlr passed\n');
