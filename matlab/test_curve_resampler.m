clear;
clc;
addpath('./src');

x = importdata('bunny10out');
x = x';
nb = 400;
% simple interface
tic, [srcinfo, h] = resample_curve(x,nb); toc;

% Length of curve
[~,nn] = size(srcinfo);
L = nn*h;
fprintf('Length of curve=%d\n',L);

% Full interface
eps = 0;
nmax = 4;
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

[~,n] = size(srcinfo);
nfac = 1;
nover = nfac*n;

x = srcinfo(1,:);
y = srcinfo(2,:);

L = n*h;
ts2 = (0:(nover-1))*h/nfac;
tic, srcinfo = eval_curve(ts2,wsave); toc,
[~,n] = size(srcinfo);
L = n*h/nfac;
dxdt1 = perispecdiff(srcinfo(1,:))*2*pi/L;
dydt1 = perispecdiff(srcinfo(2,:))*2*pi/L;

dxdt = -srcinfo(4,:);
dydt = srcinfo(3,:);

err1 = norm(dxdt1-dxdt)/norm(dxdt1)
err2 = norm(dydt1-dydt)/norm(dydt1)


function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = PERISPECDIFF(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points).
%
% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));
end
