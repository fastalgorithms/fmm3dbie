function [varargout] =  test1(ik,npu,norder,iker)
radii = [1.0;2.0;0.25];
scales = [1.2;1.0;1.7];
nu = npu;
nv = npu;
nosc = 5;
if(ik == 1)
    zk = 0.97;
else
    zk = 0.97*nu/10;
end
fname = ['diary_ik' int2str(ik) '_np' int2str(npu) '_norder' int2str(norder) '_qtime.dat'];
diary(fname);
S = wtorus(radii,scales,nosc,nu,nv,norder);

if(iker == 1)
    zpars = complex([zk; 1.0; 0.0]);
elseif(iker==2)
    zpars = complex([zk; 0.0; 1.0]);
elseif(iker==3)
    disp('Here3')
    zpars = complex([zk; 1j*zk; 1.0]);
elseif(iker==4)
    disp('Here4')
    zpars = complex([zk; -1j*zk; 1.0]);
end
fprintf('zk=: %d \n',zk);
disp(zpars)
eps = 0.51e-6;
tic, spmat = helm_near_corr(S,zpars,eps); toc;
diary('off');

zpars2 = complex([zk;1]);
tic, spmat_neu = helm_neu_near_corr(S,zpars2,eps); toc;

varargout{1} = S;
varargout{2} = spmat;
varargout{3} = spmat_neu;
%exit;
end
