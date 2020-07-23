clear;
clc;

addpath('./src')

%% Generate geometry and right hand side for dirichlet problem
n = 200;
src = [0.01;-0.07];
targ = [12.1;5.2];

% Generate geometry
a = 1.1; b=1.3;
[srcinfo,h] = ellipse(a,b,n);


zk = complex(1.1);
uex = helm_c_p(zk,src,targ);
rhs = helm_c_p(zk,src,srcinfo);


%% Single layer test

% Generate single layer potential matrix
norder = 16;
xmat = slp_mat(zk,norder,h,srcinfo);
soln = xmat\rhs;

% Test resulting solution
z = helm_c_p(zk,srcinfo,targ);
ucomp = (z.*srcinfo(5,:))*soln*h;
err = abs(ucomp - uex);
fprintf("Error in single layer test = %d\n\n\n",err);

%% Double layer test

% Generate double layer potential matrix
norder = 16;
xmat = dlp_ext_mat(zk,norder,h,srcinfo);
soln = xmat\rhs;

% Test resulting solution
z = helm_d_p(zk,srcinfo,targ);
ucomp = (z.*srcinfo(5,:))*soln*h;
err = abs(ucomp - uex);
fprintf("Error in double layer test = %d\n\n\n",err);

%% combined layer test

% Generate comb layer potential matrix
norder = 16;
zpars = complex(zeros(3,1));
zpars(1) = zk;
zpars(2) = 1j*zk;
zpars(3) = 1;
xmat = comb_ext_mat(zpars,norder,h,srcinfo);
soln = xmat\rhs;

% Test resulting solution
zs = helm_c_p(zk,srcinfo,targ);
zd = helm_d_p(zk,srcinfo,targ);
z = zpars(2)*zs + zpars(3)*zd;
ucomp = (z.*srcinfo(5,:))*soln*h;
err = abs(ucomp - uex);
fprintf("Error in combined layer test = %d\n\n\n",err);
