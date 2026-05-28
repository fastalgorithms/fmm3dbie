function ri = rsc_interleave_nrccie_eval(zk)
%RSC_INTERLEAVE_NRCCIE_EVAL  rsc_to_interleave for the nrccie-eval kernel.
%
%   ri = kernel3d.rsc_interleave_nrccie_eval(zk)
%
%   The nrccie-eval wnear stores 4 scalar basis kernels:
%     wnear(1,:) = S_k
%     wnear(2,:) = d/dx S_k
%     wnear(3,:) = d/dy S_k
%     wnear(4,:) = d/dz S_k
%
%   These combine to fill the 6x4 block B (rows=E/H components, cols=Jx,Jy,Jz,rho):
%
%     E = ik*S_k[J] - grad S_k[rho]
%     H = curl S_k[J]
%
%   Block layout (row: Ex=1,Ey=2,Ez=3,Hx=4,Hy=5,Hz=6; col: Jx=1,Jy=2,Jz=3,rho=4):
%     B(1,1) += ik * w1      B(1,4) += -w2
%     B(2,2) += ik * w1      B(2,4) += -w3
%     B(3,3) += ik * w1      B(3,4) += -w4
%     B(4,2) += -w4          B(4,3) += +w3
%     B(5,1) += +w4          B(5,3) += -w2
%     B(6,1) += -w3          B(6,2) += +w2

ik = 1i * zk;

e = struct('ker_id', {}, 'row_id', {}, 'col_id', {}, 'coef', {});

% wnear(1) = S_k: diagonal J blocks
e(end+1) = struct('ker_id',1,'row_id',1,'col_id',1,'coef', ik);
e(end+1) = struct('ker_id',1,'row_id',2,'col_id',2,'coef', ik);
e(end+1) = struct('ker_id',1,'row_id',3,'col_id',3,'coef', ik);

% wnear(2) = d/dx S_k
e(end+1) = struct('ker_id',2,'row_id',1,'col_id',4,'coef',-1);   % E_x from rho
e(end+1) = struct('ker_id',2,'row_id',5,'col_id',3,'coef',-1);   % H_y from J_z: -gx
e(end+1) = struct('ker_id',2,'row_id',6,'col_id',2,'coef', 1);   % H_z from J_y: +gx

% wnear(3) = d/dy S_k
e(end+1) = struct('ker_id',3,'row_id',2,'col_id',4,'coef',-1);   % E_y from rho
e(end+1) = struct('ker_id',3,'row_id',4,'col_id',3,'coef', 1);   % H_x from J_z: +gy
e(end+1) = struct('ker_id',3,'row_id',6,'col_id',1,'coef',-1);   % H_z from J_x: -gy

% wnear(4) = d/dz S_k
e(end+1) = struct('ker_id',4,'row_id',3,'col_id',4,'coef',-1);   % E_z from rho
e(end+1) = struct('ker_id',4,'row_id',4,'col_id',2,'coef',-1);   % H_x from J_y: -gz
e(end+1) = struct('ker_id',4,'row_id',5,'col_id',1,'coef', 1);   % H_y from J_x: +gz

ri = kernel3d.rsc_interleave_full(6, 4, 4, e);

end
