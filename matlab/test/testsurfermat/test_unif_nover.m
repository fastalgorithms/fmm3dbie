%TEST_UNIF_NOVER   Verify that opts.unif_nover=1 in surfermat produces
%   oversampling orders that are constant across target surfaces for each
%   source surface j: novers{i,j} must be the same for all i.

run ../../startup.m

S1 = geometries.ellipsoid([1.5,1,1], 2*[1,1,1], [-0.5;0;0], 4);
S2 = geometries.ellipsoid([1,1,1.5], 2*[1,1,1], [10;0;0], 4);
S3 = geometries.ellipsoid([1,1,1.5], [1,1,1], [8;0;0], 4);

kern = kernel3d('z', [1,1]);
kern.kernel_order = -1;
kern.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, kern.zk, kern.kernel_order);

% test uniform
eps  = 1e-8;
opts = struct('corrections', 1, 'unif_nover', 1);
[~, novers] = surfermat([S1, S2], kern, eps, opts);

pass = true;
for j = 1:size(novers,2)
    if ~all(all([novers{:,j}]==[novers{j,j}],2))
        fprintf('FAIL: novers{:,%d} ~= novers{%d,%d}\n', j, j,j);
    end
end

if pass
    fprintf('PASS\n');
end

% test monotonic when adding a new surfer
[~, novers2] = surfermat([S1, S2, S3], kern, eps, opts);

pass = true;
for j = 1:size(novers,2)
    if ~all(all([novers{:,j}]<=[novers2{j,j}],2))
        fprintf('FAIL: novers{:,%d} > novers2{%d,%d}\n', j,j,j);
    end
end

if pass
    fprintf('PASS\n');
end

% test not equal with unif_nover off
opts = struct('corrections', 1);
[~, novers_nonuni] = surfermat([S1, S2], kern, eps, opts);

pass = false;
for j = 1:size(novers_nonuni,2)
    if ~all(all([novers_nonuni{:,j}]==[novers_nonuni{j,j}],2))
        fprintf('Pass: novers{:,%d} ~= novers{%d,%d}\n', j, j,j);
        pass = true;
    end
end

if pass
    fprintf('PASS\n');
end