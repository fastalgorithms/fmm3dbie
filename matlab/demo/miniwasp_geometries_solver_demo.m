geom = 'rhab';
iref = 2;
S = surfer.load_from_file(['~/git/miniwasp/geometries/' geom '_aspfixed_iref' ...
       int2str(iref) '_ippw05_norder05_irlam3.go3']);
plot(S);      
zk = 0.009;


switch lower(geom)
    case 'lens'
        xyz_in = [101.1;-301.0;1.69e4];
    case 'cone'
        xyz_in = [101.1;-301.0;1.69e4];
    case 'rhab'
        xyz_in = [101.1;-301.0;1.69e4];
end


plot(S.r(1,1:42),S.r(2,1:42),'b.','MarkerSize',20)
axis('equal')
return
    
xyz_out = [1.3;-5.2;2.7e4];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');


zpars = [zk,-1j*zk,1];
sig = helm3d.dirichlet.solver(S,zpars,rhs,eps);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',zpars(2),zpars(3));

wts = S.wts;
pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));
