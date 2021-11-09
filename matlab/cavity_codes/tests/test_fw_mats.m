
n  = 300;

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
src_info = [];
src_info.xs = real(xy);
src_info.ys = imag(xy);
src_info.ds = abs(dxydt);
src_info.dxs = real(dxydt);
src_info.dys = imag(dxydt);
src_info.h = h;
% src_info.H = get_curvature(src_info);
src_info.L = sum(src_info.ds)*h;
plot(src_info.xs,src_info.ys,'b.');

src0 = [0.01;-1.2];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;

% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];


sensor_info = [];
sensor_info.xtgt = x_t;
sensor_info.ytgt = y_t;


kh = 2.2;

bc = [];
bc.type = 'Dirichlet';


[mats,erra] = get_fw_mats(kh,src_info,bc,sensor_info,opts);
fprintf('Error in dirichlet problem: %d\n',erra);


bc = [];
bc.type = 'Neumann';
[mats,erra] = get_fw_mats(kh,src_info,bc,sensor_info,opts);
fprintf('Error in Neumann problem: %d\n',erra);


src_info.lambda = ones(n,1);
bc = [];
bc.type = 'Impedance';
[mats,erra] = get_fw_mats(kh,src_info,bc,sensor_info,opts);
fprintf('Error in Impedance problem: %d\n',erra);
