clear;
clc;

addpath('./src')
%
%  This is a test script for solving Helmholtz transmission 
% problem with boundary conditions
% [au]/q = f, [b dudn] = g
% where q = 0.5*(a(1)/b(1) + a(2)/b(2))
%
%  We create an artificial solution by placing a charge
% outside to generate the helmholz solution in the interior
% and a charge inside to generate a helmholtz solution
% in the exterior and read off the jumps between these two 
% solutions. 
%   
%

% Generate geometry and right hand side for dirichlet problem
n = 200;
src = [0.01 12.1;-0.07 5.2];
targ = [-0.02 13.1;0.03 4.7];

a0 = 1.1; b0=1.3;
[srcinfo,h] = ellipse(a0,b0,n);

zks = [complex(1.1);complex(1.7)];
a = [complex(1.0); complex(1.5)];
b = [complex(2.0); complex(4.0)];


% Compute exact solutions
uex_in = helm_c_p(zks(1),src(:,2),targ(:,1));
uex_out = helm_c_p(zks(2),src(:,1),targ(:,2));

nsys = 2*n;


% Generate right hand side
q = (a(1)/b(1) + a(2)/b(2))*0.5;
rhs = complex(zeros(nsys,1));
uin = helm_c_p(zks(1),src(:,2),srcinfo);
uout = helm_c_p(zks(2),src(:,1),srcinfo);

dudnin = helm_c_gn(zks(1),src(:,2),srcinfo);
dudnout = helm_c_gn(zks(2),src(:,1),srcinfo);
rhs(1:2:nsys) = (a(2)*uout - a(1)*uin)/q;
rhs(2:2:nsys) = b(2)*dudnout - b(1)*dudnin;




% Generate transmission matrix
norder = 16;
xmat = transmission_mat(zks,a,b,norder,h,srcinfo);
soln = xmat\rhs;


% test solution at interior target
zs = helm_c_p(zks(1),srcinfo,targ(:,1)).*srcinfo(5,:)*h;
zd = helm_d_p(zks(1),srcinfo,targ(:,1)).*srcinfo(5,:)*h;

ucomp_in = (-zs*soln(1:2:nsys)+zd*soln(2:2:nsys))/b(1);

% test solution at exterior target
zs = helm_c_p(zks(2),srcinfo,targ(:,2)).*srcinfo(5,:)*h;
zd = helm_d_p(zks(2),srcinfo,targ(:,2)).*srcinfo(5,:)*h;

ucomp_out = (-zs*soln(1:2:nsys)+zd*soln(2:2:nsys))/b(2);


err = abs(ucomp_in-uex_in) + abs(ucomp_out-uex_out);
fprintf("Error in transmission test = %d\n\n\n",err);
