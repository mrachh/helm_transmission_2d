
%  this subroutine returns the truncated 2d fourier series expanison
%  of the indicator set of the region whose boundary is defined
%  by the src_info struct.
%
%  The fourier coeffs are organized in the standard ordering

addpath('~/git/finufft/matlab/'); % set path here for finufft directory
n = 100;
t = 0:2*pi/n:2*pi-2*pi/n;
t = t(:);
rt = 1+0.3*cos(3*t);
drdt = -0.9*sin(3*t);
xs = rt.*cos(t);
ys = rt.*sin(t);
dxdt = -rt.*sin(t) + drdt.*cos(t);
dydt = rt.*cos(t) + drdt.*sin(t);
dsdt = sqrt(dxdt.^2 + dydt.^2);
rnx = dydt./dsdt;
rny = -dxdt./dsdt;
h = 2*pi/n;

nmax = 40;
xtest = -pi/2:0.01:pi/2;
ntest=length(xtest);
[xx,yy]=meshgrid(xtest);
x=xx(:);
y=yy(:);

c = compute_trunc_fourier(xs,ys,rnx,rny,dsdt,h,nmax,x,y);
c = reshape(c,[ntest,ntest]);
fh = pcolor(xx,yy,c);
colorbar();
set(fh,'EdgeColor','None');
