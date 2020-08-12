clear;
clc;

addpath('./src')

%% Generate geometry and right hand side for dirichlet problem
n = 200;
targ = [0.01;-0.07];
src = [12.1;5.2];

% Generate geometry
a = 1.1; b=1.3;
[srcinfo,h] = ellipse(a,b,n);

uex = lap_c_p(src,targ);
rhs = lap_c_p(src,srcinfo);


%% Double layer test

% Generate double layer potential matrix
norder = 16;
xmat = lap_dlp_mat(norder,h,srcinfo);

% make it interior matrix by subtracting identity
xmat = xmat - eye(n);
soln = xmat\rhs;

% Test resulting solution
z = lap_d_p(srcinfo,targ);
ucomp = (z.*srcinfo(5,:))*soln*h;
err = abs(ucomp - uex);
fprintf("Error in double layer test = %d\n\n\n",err);