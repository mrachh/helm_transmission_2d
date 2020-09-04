clear;
clc;

addpath('./src')

%% Generate geometry and right hand side for dirichlet problem
n = 400;
src = [0.01;-0.07];
targ = [12.1;5.2];

% Generate geometry
a = 3.1; b=1.3;
[srcinfo,h] = ellipse(a,b,n);
nh = 15;
hcoefs = rand(2*nh+1,1)*0.1;

rlin = 2*pi;
nout = n;
[srcinfoout,hout,rltot,ier,tts] = resample_curve(srcinfo, ...
  h,rlin,nh,hcoefs,nout,eps);
figure(1)
plot(srcinfoout(1,:),srcinfoout(2,:),'kx'); axis equal;

sfft = fftshift(fft(srcinfo(5,:)));

xfft = fftshift(fft(srcinfoout(1,:)));
figure(2)
semilogy(abs(sfft),'kx')
figure(3)
semilogy(abs(xfft),'kx')