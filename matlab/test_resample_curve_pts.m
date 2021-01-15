clear;
clc;
addpath('./src');

x = importdata('bunny10out');
x = x';
nb = 150;
% simple interface
tic, [rltot, wsave] = resample_curve_pts(x,nb); toc;

fprintf('Length of curve=%d\n',rltot);

% Full interface
eps = 1e-9;
tic, [rltot, wsave, tts] = resample_curve_pts(x,nb,eps,4); toc;


% Evalute curve at user specified points
tic, binfo = eval_curve(tts,wsave); toc;

figure(1);
clf();
plot(x(1,:),x(2,:),'rx'), hold on;
plot(binfo(1,:),binfo(2,:),'k-');

% Routine to resample curve at a given number of nodes
% currently being used to test interpolation error
err = norm(binfo(1,:) - x(1,:)).^2 + norm(binfo(2,:)-x(2,:)).^2;
rr = norm(x(:))^2;
err = sqrt(err/rr);
fprintf('error in interpolation=%d\n',err);

