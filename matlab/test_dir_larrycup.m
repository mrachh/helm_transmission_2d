clear;
clc;

addpath('./src')

%% Generate geometry and right hand side for dirichlet problem
n = 300;
src = [0.01;-1.1];
targ = [12.1;5.2];

% Generate geometry

a = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/6;  % controls approx opening angle in radians (keep small for resonant)


nhalf = ceil(n/2);
s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
%th = b + sqrt(log(exp((2*(1-b/pi)*(s-pi/2)).^2) + exp(c^2))); % theta

% coords in (rho,z) plane:
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords

z = z*1.2;  % vert stretch! makes ellipse cavity

Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve

% (appropriate for half-integer offset
figure; semilogy(abs(fft(Z))); title('Fourier coeff decay, to close to emach?')
%Z = Z(end:-1:1);
zhat = fft(Z(:))/n;
t1 = (0:(n-1))/n;
h = 1.0/n;
xy = fourierZ(zhat,t1);
dxydt = fourierZp(zhat,t1);
d2xydt2 = fourierZpp(zhat,t1);
srcinfo = zeros(6,n);
srcinfo(1,:) = real(xy);
srcinfo(2,:) = imag(xy);
srcinfo(5,:) = abs(dxydt);
srcinfo(3,:) = imag(dxydt)./srcinfo(5,:);
srcinfo(4,:) = -real(dxydt)./srcinfo(5,:);

figure(2)
plot(srcinfo(1,:),srcinfo(2,:),'k.')

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




% analytic formulae for a Fourier segment --------------

function z = fourierZ(zhat,t)     % must work on vector of t's
t = 2*pi*t;
N =numel(zhat);  % even
z = 0*t;
for k=0:N/2
  z = z + zhat(k+1)*exp(1i*k*t);
end
for k=-N/2:-1
  z = z + zhat(k+1+N)*exp(1i*k*t);
end
end

function zp = fourierZp(zhat,t);  % deriv func Z'
N = numel(zhat);
zp = 2*pi*fourierZ(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);
end

function zpp = fourierZpp(zhat,t);  % deriv func Z''
N = numel(zhat);
zpp = 2*pi*fourierZp(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);
end

% ---------------------

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
