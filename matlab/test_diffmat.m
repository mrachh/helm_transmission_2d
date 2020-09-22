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
zk1 = complex(1.1,0.0);
zk2 = complex(2.1,0.0);

%% Compute Dk' using Maue's identity

rnx = srcinfo(3,:);
rny = srcinfo(4,:);

S1 = slp_mat(zk1,norder,h,srcinfo);
dm1 = (D*S1*D + zk1.^2.*(diag(rnx)*S1*diag(rnx) + ... 
  diag(rny)*S1*diag(rny)));
S2 = slp_mat(zk2,norder,h,srcinfo);
dm2 = (D*S2*D + zk2.^2.*(diag(rnx)*S2*diag(rnx) + ... 
  diag(rny)*S2*diag(rny)));


zpars = complex(zeros(2,1));
zpars(1) = zk1+0.0j;
zpars(2) = zk2+0.0j;
zpars = complex(zpars);
dkdiffmat = ddiff_neu_mat(zpars,norder,h,srcinfo);

dkdiffmat3 = (dm1-dm2);
err = norm(dkdiffmat(:)-dkdiffmat3(:))/norm(dkdiffmat(:));
fprintf("Error in dk diff mat = %d\n\n\n",err);

z1 = dkdiffmat*rhs;
z2 = dkdiffmat3*rhs;

err = norm(z1-z2)/norm(z1);
fprintf("Error in dk diff mat apply = %d\n\n\n",err);