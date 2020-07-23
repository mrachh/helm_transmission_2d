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


zk = complex(1.1);
uex = helm_c_p(zk,src,targ);
uin = helm_c_p(zk,src,srcinfo);
dudnin = helm_c_gn(zk,src,srcinfo);

alpha = 2;
zk2 = 1j*zk;
norder = 16;
sikmat = slp_mat(zk2,norder,h,srcinfo);


%% Dirichlet test

eta = complex(0,0);
zpars = complex(zeros(3,1));
zpars(1) = zk;
zpars(2) = alpha;
zpars(3) = eta;

rhs = uin + 1j*zk*eta*dudnin;

xmat = rpcomb_dir_ext_mat(zpars,norder,h,srcinfo);
soln = xmat\rhs;
soln2 = sikmat*soln;


% Test resulting solution
zs = helm_c_p(zk,srcinfo,targ);
zd = helm_d_p(zk,srcinfo,targ);
ucomp = (zs.*srcinfo(5,:))*soln*h + 1j*alpha*(zd.*srcinfo(5,:))*soln2*h;

err = abs(ucomp - uex);
fprintf("Error in Dirichlet test = %d\n\n\n",err);

%% Neumann test

zpars = complex(zeros(3,1));
zpars(1) = zk;
zpars(2) = alpha;
zpars(3) = eta;

rhs = dudnin;

xmat = rpcomb_neu_ext_mat(zpars,norder,h,srcinfo);
soln = xmat\rhs;
soln2 = sikmat*soln;


% Test resulting solution
zs = helm_c_p(zk,srcinfo,targ);
zd = helm_d_p(zk,srcinfo,targ);
ucomp = (zs.*srcinfo(5,:))*soln*h + 1j*alpha*(zd.*srcinfo(5,:))*soln2*h;

err = abs(ucomp - uex);
fprintf("Error in Neumann test = %d\n\n\n",err);

%% Neumann test
eta = complex(1.2);
zpars = complex(zeros(3,1));
zpars(1) = zk;
zpars(2) = alpha;
zpars(3) = eta;

rhs = uin + 1j*zk*eta*dudnin;

xmat = rpcomb_imp_ext_mat(zpars,norder,h,srcinfo);
soln = xmat\rhs;
soln2 = sikmat*soln;


% Test resulting solution
zs = helm_c_p(zk,srcinfo,targ);
zd = helm_d_p(zk,srcinfo,targ);
ucomp = (zs.*srcinfo(5,:))*soln*h + 1j*alpha*(zd.*srcinfo(5,:))*soln2*h;

err = abs(ucomp - uex);
fprintf("Error in Impedence test = %d\n\n\n",err);
