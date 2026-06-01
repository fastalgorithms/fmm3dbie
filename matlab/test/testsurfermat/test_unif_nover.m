% Test opts.unif_nover=1 in surfermat.

run ../../startup.m

%% Now run the tests

test_uniform_orders();
test_monotone_on_add();
test_nonuniform_without_flag();


function [S1, S2, S3, kern, eps] = make_setup()
S1 = geometries.ellipsoid([1.5,1,1], 2*[1,1,1], [-0.5;0;0], 4);
S2 = geometries.ellipsoid([1,1,1.5], 2*[1,1,1], [10;0;0],   4);
S3 = geometries.ellipsoid([1,1,1.5], [1,1,1],   [8;0;0],    4);
kern = kernel3d('z', [1,1]);
kern.kernel_order = -1;
kern.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, kern.zk, kern.kernel_order);
eps = 1e-8;
end


function test_uniform_orders()
% With unif_nover=1, novers{:,j} must all equal novers{j,j}.

[S1, S2, ~, kern, eps] = make_setup();
opts = struct('corrections', 1, 'unif_nover', 1);
[~, novers] = surfermat([S1, S2], kern, eps, opts);
for j = 1:size(novers,2)
    assert(all(all([novers{:,j}] == [novers{j,j}], 2)), ...
        'unif_nover: novers{:,%d} not uniform', j);
end

end


function test_monotone_on_add()
% Adding a third surfer must not decrease any existing oversampling order.

[S1, S2, S3, kern, eps] = make_setup();
opts = struct('corrections', 1, 'unif_nover', 1);
[~, novers]  = surfermat([S1, S2],     kern, eps, opts);
[~, novers2] = surfermat([S1, S2, S3], kern, eps, opts);
for j = 1:size(novers,2)
    assert(all(all([novers{:,j}] <= [novers2{j,j}], 2)), ...
        'unif_nover: novers decreased after adding surfer for column %d', j);
end

end


function test_nonuniform_without_flag()
% Without unif_nover, at least one column must have non-uniform orders.

[S1, S2, ~, kern, eps] = make_setup();
[~, novers] = surfermat([S1, S2], kern, eps, struct('corrections', 1));
any_nonuniform = false;
for j = 1:size(novers,2)
    if ~all(all([novers{:,j}] == [novers{j,j}], 2))
        any_nonuniform = true;
    end
end
assert(any_nonuniform, 'expected non-uniform novers without unif_nover flag');

end
