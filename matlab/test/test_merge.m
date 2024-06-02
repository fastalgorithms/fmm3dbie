


% For quads
S = geometries.sphere(1, 2, [0;0;0]);

rad = 3;
nsp = 6;


shifts = zeros(3,nsp);
for ii=1:nsp
    th = (ii-1)*(2*pi)/nsp;
    shift = rad*[cos(th),sin(th),0].';
    shifts(:,ii) = shift;
end

 [S] = translate(S,shifts);
tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;

 zk = 1.1;
rep_pars = [1.0, 0.0];
ndeg = 1;


xyz_in = [rad;0;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');

rep_pars = [-1j*zk, 1];
eps = 1E-4;
sig = helm3d.dirichlet.solver(S,rhs,eps,zk,rep_pars);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',rep_pars(1),rep_pars(2));

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));


euls = [0,pi/2,0];
p1 = euls(1);
p2 = euls(2);
p3 = euls(3);
t1 = [cos(p3),sin(p3),0;-sin(p3),cos(p3),0;0,0,1];
t2 = [1,0,0;0,cos(p2),sin(p2);0,-sin(p2),cos(p2)];
t3 = [cos(p1),sin(p1),0;-sin(p1),cos(p1),0;0,0,1];


xyz_in = t1*t2*t3*xyz_in;
xyz_out = t1*t2*t3*xyz_out;
[S] = rotate(S,euls);

src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');

rep_pars = [-1j*zk, 1];
eps = 1E-4;
sig = helm3d.dirichlet.solver(S,rhs,eps,zk,rep_pars);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',rep_pars(1),rep_pars(2));

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));