%driver 
%global variables for old code
clear
global type_obj
global radius
type_obj = 20;
radius = 0.3;

%variables
n  = 64;%half the number of points in the domain
t  = 0:pi/n:2*pi-pi/n; %discretization of boundary
kh = 1; % wavenumber
d  = [1 0]; %incident direction of a plan e wave

%source and target
source = [ 5 5 ]; %source point of the incident field
% target = [ 0 0 ]; %target point inside the domain
% res = (1i/4)*besselh(0,1,kh*norm(target-source)); %value to check

%boundary
% c_t    = cos(t);
% s_t    = sin(t);
% a      = radius;
% p_t    = ( 1 + a * cos(4*t) );
% dp_t   = -4 * a * sin(4*t);
% ddp_t  = -16 * a * cos(4*t);
% dddp_t = 64 * a * sin(4*t);
% xs     = p_t .* c_t;
% ys     = p_t .* s_t;
% dxs    = dp_t .* c_t - p_t .* s_t;
% dys    = dp_t .* s_t + p_t .* c_t;
% ddxs   = ddp_t .* c_t - 2 * dp_t .* s_t - p_t .* c_t;
% d3xs   = dddp_t .* c_t + ddp_t .*(-s_t) - 2 * (ddp_t .* s_t + dp_t .* c_t) - (dp_t .*c_t + p_t .*(-s_t));
% ddys   = ddp_t .* s_t + 2 * dp_t .* c_t - p_t .* s_t;
% d3ys   = dddp_t .* s_t + ddp_t .* c_t + 2 * (ddp_t .* c_t - dp_t .* s_t) - (dp_t .* s_t + p_t .* c_t);
% ds     = sqrt(dxs.^2 + dys.^2);

%boundary
xs   = cos(t) + 0.65 * cos(2*t) - 0.65;
ys   = 1.5 * sin(t);
dxs  = -sin(t) - 1.3 * sin(2*t);
dys  = 1.5 * cos(t);
ddxs = -cos(t) - 2.6 * cos(2*t);
d3xs = sin(t) + 5.2 * sin(2*t);
ddys = -1.5 * sin(t);
d3ys = -1.5 * cos(t);
ds   = sqrt(dxs.^2 + dys.^2);

%target point
x_t = [0 0.1 -0.1  0.1 -0.1];
y_t = [0 0.1 -0.1 -0.1  0.1];

res = (1i/4)*besselh(0,1,kh*sqrt(bsxfun(@minus,x_t',source(1)).^2+bsxfun(@minus,y_t',source(2)).^2));

%incident field
rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
u_inc = (1i/4)*besselh(0,1,kh*rr);
du_inc = -1i*kh./4*besselh(1,1,kh*rr)./rr.*(dys.*(xs-source(1))- dxs.*(ys-source(2)));
% u_plane = exp(1i*kh*xs);

%setting domain
src = zeros(8,2*n);
src(1,:) = xs;
src(2,:) = ys;
src(3,:) = dxs;
src(4,:) = dys;
src(5,:) = ddxs;
src(6,:) = ddys;
src(7,:) = d3xs;
src(8,:) = d3ys;


%setting target
tgt = zeros(2,length(x_t));
tgt(1,:) = x_t;
tgt(2,:) = y_t;

%settiing xhat
t1 = 0:2*pi/2:2*pi-2*pi/2;
xhat = cos(t1);
yhat = sin(t1);

%creating potential
S  = slmat(kh,src,t);
D  = dlmat(kh,src,t);
% S1 = single2D(kh,t,t);
% D1 = double2D(kh,t,t);
Sp = sprimelmat(kh,src,t);
T  = dprimelmat(kh,src,t);

% max(max(abs(D-D1/2))) 
% max(max(abs(S-S1/2)))

%Single layer
%solving for the field 
pot_s = S \ transpose(u_inc);

%Checking the solution in the target point
% rr = sqrt((xs-target(1)).^2 + (ys-target(2)).^2);
% sol_s = (pi/n)*(1i/4)*sum(besselh(0,1,kh*rr).*ds.*transpose(pot_s));

sol_s = slmat_out(kh,src,tgt)*pot_s;

fprintf('For single:%d\n',norm(res-sol_s)/norm(res))

%For Sprime layer
% sol_g = (eye(2*n)/2 + Sp -1i*S) \ transpose(du_inc -1i * u_inc);

%Double layer
%solving for the field 
pot_d = (-eye(2*n)/2 + D) \ transpose(u_inc);

sol_d = dlmat_out(kh,src,tgt)*pot_d;

fprintf('For double:%d\n',norm(res-sol_d)/norm(res))

%Combined layer
%solving for the field 
eta   = kh;
pot_c = (-eye(2*n)/2 + D + 1i*eta*S) \ transpose(u_inc);

%Checking the solution in the target point
sol_c = (dlmat_out(kh,src,tgt)+1i*eta*slmat_out(kh,src,tgt))*pot_c;

fprintf('For combined:%d\n',norm(res-sol_c)/norm(res))

%With plane waves
uinc  = exp(1i*kh*xs);
duinc = 1i * kh * dys .* exp(1i*kh*xs)./ds;

%Combined layer vs double layer vs single layer
tgt = zeros(2,4);
tgt(1,:) = [5 -5 -5  5];
tgt(2,:) = [5  5 -5 -5];
eta = 1;
pot_s   = S \ transpose(- uinc) ;
uscat_s = slmat_out(kh,src,tgt)*pot_s;

% uscat_s

uinf_s = exp(1i*pi/4)/sqrt(8*pi*kh)*2*pi/length(xs)*...
    exp(-1i*kh*(bsxfun(@times,xhat',xs)+bsxfun(@times,yhat',ys)))*(pot_s.*(ds'));

% uinf_s

pot_d   = (eye(2*n)/2 + D ) \ transpose(- uinc) ;
uscat_d = dlmat_out(kh,src,tgt)*pot_d;

% uscat_d

uinf_d = exp(-1i*pi/4)/sqrt(8*pi*kh)*2*pi/length(xs)*...
    exp(-1i*kh*(bsxfun(@times,xhat',xs)+bsxfun(@times,yhat',ys))).*...
    kh .* (bsxfun(@times,xhat',dys)-bsxfun(@times,yhat',dxs))*...
    pot_d;

% uinf_d

pot_c   = (eye(2*n)/2 + D - 1i * eta * S) \ transpose(- uinc) ;
uscat_c = (dlmat_out(kh,src,tgt)-1i*eta*slmat_out(kh,src,tgt))*pot_c;

% uscat_c

uinf_c = exp(-1i*pi/4)/sqrt(8*pi*kh)*2*pi/length(xs)*...
    exp(-1i*kh*(bsxfun(@times,xhat',xs)+bsxfun(@times,yhat',ys))).* ...
    (kh * (bsxfun(@times,xhat',dys)-bsxfun(@times,yhat',dxs)) + eta * repmat(ds,length(xhat),1)) *...
    pot_c;

% uinf_c


fprintf('Single - Double=%d\n',norm(uscat_s-uscat_d))
fprintf('Single - Combined=%d\n',norm(uscat_s-uscat_c))
fprintf('Double - Combined=%d\n',norm(uscat_c-uscat_d))

%Sprime layer
eta = kh;
pot_s = S \ transpose(uinc); %this is dudn

pot_sp = (eye(2*n)/2 + Sp ) \ transpose(duinc);

pot_c  =  (eye(2*n)/2 + Sp - 1i * eta * S) \ transpose(duinc - 1i * eta * uinc);

fprintf('For Sprime:%d\n',norm(pot_sp - pot_s)/norm(pot_s))
fprintf('For Sprime Combined:%d\n',norm(pot_c - pot_s)/norm(pot_s))

%Dprimelayer
fprintf('For Derivative of D!\n')
eta = 1;
% pot_t = (eye(2*n)/2 + D + 1i * eta * S) \ transpose(- uinc) ;
% 
% sol_t0 = ( T + 1i * eta * (Sp - eye(2*n)/2)) * pot_t;
% 
% sol_dp = sol_t0 + transpose(duinc);

pot_t = ( S + 1i * eta * ( D + eye(2*n)/2 )) \ transpose(- uinc) ;

sol_t0 = ( Sp - eye(2*n)/2 + 1i * eta * T) * pot_t;

sol_dp = sol_t0 + transpose(duinc);

fprintf('For Dprime:%d\n',norm(pot_s-sol_dp)/norm(pot_s))

pot_sp = (Sp - eye(2*n)/2) \ transpose(-duinc);

uscat_s1 = slmat_out(kh,src,tgt)*pot_sp

uinf_sp = exp(1i*pi/4)/sqrt(8*pi*kh)*2*pi/length(xs)*...
    exp(-1i*kh*(bsxfun(@times,xhat',xs)+bsxfun(@times,yhat',ys)))*(pot_sp.*(ds'));

uinf_sp

pot_cp1  = (T -1i*eta* (Sp - eye(2*n)/2)) \ transpose(-duinc);

uscat_s2 = (dlmat_out(kh,src,tgt)-1i*eta*slmat_out(kh,src,tgt))*pot_cp1

uinf_dp = exp(-1i*pi/4)/sqrt(8*pi*kh)*2*pi/length(xs)*...
    exp(-1i*kh*(bsxfun(@times,xhat',xs)+bsxfun(@times,yhat',ys))).* ...
    (kh * (bsxfun(@times,xhat',dys)-bsxfun(@times,yhat',dxs)) + eta * repmat(ds,length(xhat),1)) *...
    pot_cp1;

uinf_dp

pot_cp2  = (-1i*eta * T + (Sp - eye(2*n)/2)) \ transpose(-duinc);

uscat_s3 = (-1i*eta * dlmat_out(kh,src,tgt) + slmat_out(kh,src,tgt))*pot_cp2

pot_t = ( T ) \ transpose(-duinc);

uscat_s4 = dlmat_out(kh,src,tgt) * pot_t







