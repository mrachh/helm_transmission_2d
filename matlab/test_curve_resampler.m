clear;
clc;
addpath('./src');

x = importdata('bunny10out');
x = x';
nb = 200;
% simple interface
tic, [srcinfo, h] = resample_curve(x,nb); toc;

% Length of curve
[~,nn] = size(srcinfo);
L = nn*h;
fprintf('Length of curve=%d\n',L);

% Full interface
eps = 1e-12;
nmax = 1;
tic, [srcinfo,h,wsave,ts] = resample_curve(x,nb,eps,nmax); toc;




figure(1);
clf();
plot(x(1,:),x(2,:),'rx'), hold on;
plot(srcinfo(1,:),srcinfo(2,:),'k-');

% Routine to resample curve at a given number of nodes
% currently being used to test interpolation error
tic, binfo = eval_curve(ts,wsave); toc;
err = norm(binfo(1,:) - x(1,:)).^2 + norm(binfo(2,:)-x(2,:)).^2;
rr = norm(x(:))^2;
err = sqrt(err/rr);
fprintf('error in interpolation=%d\n',err);