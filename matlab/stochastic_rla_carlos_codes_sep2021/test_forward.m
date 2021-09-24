
% Larry cup parameters

a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/3;  % controls approx opening angle in radians (keep small for resonant)

kh = 25;

% Estimate length of curve
n = 300;
nhalf = ceil(n/2);
s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
z = z*1.2;  % vert stretch! makes ellipse cavity
Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve

% % (appropriate for half-integer offset
% figure; semilogy(abs(fft(Z))); title('Fourier coeff decay, to close to emach?')
% %Z = Z(end:-1:1);
zhat = fft(Z(:))/n;
t1 = (0:(n-1))/n;
h = 1.0/n;
xy = fourierZ(zhat,t1);
dxydt = fourierZp(zhat,t1);
Lorig = sum(abs(dxydt))*h;    


% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];

% Set interior point
src0 = [0.01;-1.1];%horseshoe


% incidence directions
n_dir = 100;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
x_dir = cos(t_dir);
y_dir = sin(t_dir);
dir =[ x_dir; y_dir ];

%generating the boundary
wl = 2*pi/kh;
Nw = Lorig/wl;
n_bd  = ceil(100*Nw);
if (n_bd < 600)
    n_bd = 600;
end
if mod(n_bd,2)
    n_bd = n_bd+1;
end
n = n_bd;
a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/3;  % controls approx opening angle in radians (keep small for resonant)
nhalf = ceil(n/2);
s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
z = z*1.2;  % vert stretch! makes ellipse cavity
Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve
zhat = fft(Z(:))/n;
t1 = (0:(n-1))/n;
h = 1.0/n;
xy = fourierZ(zhat,t1);
dxydt = fourierZp(zhat,t1);
d2xydt2 = fourierZpp(zhat,t1);
src_info = zeros(6,n);
src_info(1,:) = real(xy);
src_info(2,:) = imag(xy);
src_info(5,:) = abs(dxydt);
src_info(3,:) = imag(dxydt)./src_info(5,:);
src_info(4,:) = -real(dxydt)./src_info(5,:);
t_bd  = t1;
h_bd  = h;
L     = sum(src_info(5,:))*h;
src_old = src_info;    

%calculating field measurements  
tic,
[umeas.data,~,~,~,~,~,~,~] = forward_dirichlet(kh,src_info,L,h_bd,dir,tgt,src0);
toc

[work,lw] = hank106_workarray();
tic,
[umeas.data,~,~,~,~,~,~,~] = forward_dirichlet_fast(kh,src_info,L,h_bd,dir,tgt,src0,work,lw);
toc