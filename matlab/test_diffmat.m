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

D = specdiffmat(n,srcinfo);
rhsdiff = D*rhs;

[ux,uy] = helm_c_g(zk,src,srcinfo);
rhs_diff_exact =  -ux.*srcinfo(4,:)' + uy.*srcinfo(3,:)';
err = norm(rhsdiff-rhs_diff_exact);
fprintf("Error in diff mat = %d\n\n\n",err);


%%%%%
% Now test dk1' - dk2' by constructing dk1' and dk2' 
norder = 16;
zk1 = complex(1.1,0.1);
zk2 = complex(2.1,0.0);
d0mat = d0_neu_mat(norder,h,srcinfo);
dk1mat = dk_neu_mat(zk1,norder,h,srcinfo);
dk2mat = dk_neu_mat(zk2,norder,h,srcinfo,d0mat); % use precomputed d0mat

zpars = complex(zeros(2,1));
zpars(1) = zk1;
zpars(2) = zk2;

dkdiffmat = ddiff_neu_mat(zpars,norder,h,srcinfo);
dkdiffmat2 = dk1mat-dk2mat;

err = norm(dkdiffmat(:)-dkdiffmat2(:));
fprintf("Error in dk diff mat = %d\n\n\n",err);
