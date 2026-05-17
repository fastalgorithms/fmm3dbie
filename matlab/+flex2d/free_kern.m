function u = free_kern(src,targ,nu,zk)
ndt = size(targ,1);
dpars = nu;
zpars = zk;
ndd = 1;
ndz = 2;

ndi = 0;
ipars = [];

u1 = zeros(1,1,'like',1i);
mex_id_ = 'flex2d_gsupp2(i double[x], i int64_t[x], i double[x], i int64_t[x], i double[x], i int64_t[x], i dcomplex[x], i int64_t[x], i int64_t[x], io dcomplex[x])';
[u1] = kern_routs(mex_id_, src, ndt, targ, ndd, dpars, ndz, zpars, ndi, ipars, u1, 3, 1, ndt, 1, ndd, 1, ndz, 1, ndi, 1);
u2 = zeros(1,1,'like',1i);
mex_id_ = 'flex2d_gfree2(i double[x], i int64_t[x], i double[x], i int64_t[x], i double[x], i int64_t[x], i dcomplex[x], i int64_t[x], i int64_t[x], io dcomplex[x])';
[u2] = kern_routs(mex_id_, src, ndt, targ, ndd, dpars, ndz, zpars, ndi, ipars, u2, 3, 1, ndt, 1, ndd, 1, ndz, 1, ndi, 1);
u = [u1;u2];
end
%
%
%
%
%
%
