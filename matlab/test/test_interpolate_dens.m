test_interpolate_dens0();

function test_interpolate_dens0()

iseed = 312;
rng(iseed);

S = geometries.disk([],[],[3 3 3],8,1);
dens = eval_gauss(S.r(1,:),S.r(2,:));

ipatchids = randsample(S.npatches,4);
uvs_targ = [0.1 0.8; 0.3 0.5; 0.7 0.1; 0.24 0.31].';

xs = interpolate_dens(S,S.r(1,:),ipatchids,uvs_targ);
ys = interpolate_dens(S,S.r(2,:),ipatchids,uvs_targ);
vals = interpolate_dens(S,dens,ipatchids,uvs_targ);

vals2 = eval_gauss(xs,ys);

assert(vecnorm( vals - vals2 ) < 1e-10)

S = geometries.disk([],[],[3 3 3],8,11);

iseed = 312;
rng(iseed);
dens = eval_gauss(S.r(1,:),S.r(2,:));

ipatchids = randsample(S.npatches,4);
uvs_targ = [0.1 0.8; 0.3 0.5; 0.7 0.1; 0.24 0.31].';

xs = interpolate_dens(S,S.r(1,:),ipatchids,uvs_targ);
ys = interpolate_dens(S,S.r(2,:),ipatchids,uvs_targ);
vals = interpolate_dens(S,dens,ipatchids,uvs_targ);

vals2 = eval_gauss(xs,ys);

assert(vecnorm( vals - vals2 ) < 1e-10)


S = geometries.disk([],[],[3 3 3],8,12);

iseed = 312;
rng(iseed);
dens = eval_gauss(S.r(1,:),S.r(2,:));

ipatchids = randsample(S.npatches,4);
uvs_targ = [0.1 0.8; 0.3 0.5; 0.7 0.1; 0.24 0.31].';

xs = interpolate_dens(S,S.r(1,:),ipatchids,uvs_targ);
ys = interpolate_dens(S,S.r(2,:),ipatchids,uvs_targ);
vals = interpolate_dens(S,dens,ipatchids,uvs_targ);

vals2 = eval_gauss(xs,ys);

assert(vecnorm( vals - vals2 ) < 1e-10)


end

function val = eval_gauss(x,y)
    val = exp(-x.^2 - y.^2);
end