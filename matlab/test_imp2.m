clear;
clc;

addpath('./src')

%% Generate geometry and right hand side for impendence problem

n = 200;
src = [0.01;-0.07];
targ = [12.1;5.2];

% Generate geometry
a = 1.1; b=1.3;
[srcinfo,h] = ellipse(a,b,n);


zk = complex(1.1,0.1);
uex = helm_c_p(zk,src,targ);
uin = helm_c_p(zk,src,srcinfo);
dudnin = helm_c_gn(zk,src,srcinfo);

alpha = 2;
norder = 16;

%% Neumann test
eta = complex(1.2);
zpars = complex(zeros(3,1));
zpars(1) = zk;
zpars(2) = alpha;
zpars(3) = eta;

rhs = 1j*zk*alpha*uin + dudnin;

Der = specdiffmat(n,srcinfo);
D = dlp_ext_mat(zk,norder,h,srcinfo);
S = slp_mat(zk,norder,h,srcinfo);
Sp = sprime_ext_mat(zk,norder,h,srcinfo);

rnx = srcinfo(3,:);
rny = srcinfo(4,:);

T = (Der*S*Der + zk.^2.*(diag(rnx)*S*diag(rnx) + ... 
  diag(rny)*S*diag(rny)));
xmat = (T + 1j*eta*Sp) + 1j*zk*alpha*(D + 1j*eta*S);
soln = xmat\rhs;


% Test resulting solution
zs = helm_c_p(zk,srcinfo,targ);
zd = helm_d_p(zk,srcinfo,targ);
ucomp = 1j*eta*(zs.*srcinfo(5,:))*soln*h + (zd.*srcinfo(5,:))*soln*h;

err = abs(ucomp - uex);
fprintf("Error in Impedence test = %d\n\n\n",err);


%% 
% Now resample the curve and perform the same test
n_bd = 3;
var_up        = zeros(1,2*n_bd+1); 
L0 = 2*pi;
[srcout,hout,Lout,~,tt] = resample_curve(srcinfo,L0,n_bd,var_up',n);

plot(srcout(1,:),srcout(2,:),'x'), hold on, axis equal

uin = helm_c_p(zk,src,srcout);
dudnin = helm_c_gn(zk,src,srcout);

rhs = 1j*zk*alpha*uin + dudnin;


Der = specdiffmat(n,srcout);
rsc = 2*pi/Lout;
Der = Der*rsc;
D = dlp_ext_mat(zk,norder,hout,srcout)*rsc;
S = slp_mat(zk,norder,hout,srcout)*rsc;
Sp = sprime_ext_mat(zk,norder,hout,srcout)*rsc;

rnx = srcout(3,:);
rny = srcout(4,:);

T = (Der*S*Der + zk.^2.*(diag(rnx)*S*diag(rnx) + ... 
  diag(rny)*S*diag(rny)));
xmat = (T + 1j*eta*Sp) + 1j*zk*alpha*(D + 1j*eta*S);
soln = xmat\rhs;



% Test resulting solution
zs = helm_c_p(zk,srcout,targ);
zd = helm_d_p(zk,srcout,targ);
ucomp = 1j*eta*(zs.*srcout(5,:))*soln*h + (zd.*srcout(5,:))*soln*h;

err = abs(ucomp - uex);
fprintf("Error in Impedence test = %d\n\n\n",err);

