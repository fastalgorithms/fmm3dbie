test_interpolate_data0();

function test_interpolate_data0()

iseed = 319;
rng(iseed);

S = geometries.disk([],[],[3 3 3],8,1);
dens = [eval_gauss(S.r(1,:),S.r(2,:)) eval_gauss(S.r(2,:),S.r(1,:)).^2];

ipatchids = randi(S.npatches,[4,1]);
uvs_targ = [0.1 0.8; 0.3 0.5; 0.7 0.1; 0.24 0.31].';

xs = S.interpolate_data(S.r(1,:),ipatchids,uvs_targ);
ys = S.interpolate_data(S.r(2,:),ipatchids,uvs_targ);
vals = S.interpolate_data(dens,ipatchids,uvs_targ);

vals2 = [eval_gauss(xs,ys) eval_gauss(xs,ys).^2];

assert(vecnorm( vals(:) - vals2(:) ) < 1e-9)

S = geometries.disk([],[],[3 3 3],8,11);

iseed = 312;
rng(iseed);
dens = [eval_gauss(S.r(1,:),S.r(2,:)) eval_gauss(S.r(2,:),S.r(1,:)).^2].';

ipatchids = randi(S.npatches,[4,1]);

xs = S.interpolate_data(S.r(1,:),ipatchids,uvs_targ);
ys = S.interpolate_data(S.r(2,:),ipatchids,uvs_targ);
vals = S.interpolate_data(dens,ipatchids,uvs_targ);

vals2 = [eval_gauss(xs,ys) eval_gauss(xs,ys).^2].';

assert(vecnorm( vals(:) - vals2(:) ) < 1e-9)

S = geometries.disk([],[],[3 3 3],8,12);

iseed = 312;
rng(iseed);
dens = [eval_gauss(S.r(1,:),S.r(2,:)) eval_gauss(S.r(2,:),S.r(1,:)).^2];

xs = S.interpolate_data(S.r(1,:),ipatchids,uvs_targ);
ys = S.interpolate_data(S.r(2,:),ipatchids,uvs_targ);
vals = S.interpolate_data(dens,ipatchids,uvs_targ);

vals2 = [eval_gauss(xs,ys) eval_gauss(xs,ys).^2];

assert(vecnorm( vals(:) - vals2(:) ) < 1e-9)


end

function val = eval_gauss(x,y)
    val = exp(-x(:).^2 - y(:).^2);
end