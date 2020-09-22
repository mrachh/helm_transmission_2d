%driver 
clear

addpath('./src')
%boundaydx
N_bd             = 3;
coefs_bd         = zeros(1,2*N_bd+1);
coefs_bd(1)      = 2.;
coefs_bd(N_bd+1) = 0.;


n_bd = 200
t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
h_bd  = 2*pi/n_bd    
c_t   = cos(t_bd);
s_t   = sin(t_bd);
p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
    bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
    bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
d3p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd).*(1:N_bd))',sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
    bsxfun(@times,((1:N_bd).*(1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
xs    = p_t .* c_t;
ys    = p_t .* s_t;
dxs   = dp_t .* c_t + p_t .* (-s_t);
dys   = dp_t .* s_t + p_t .* c_t;
d2xs  = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
d2ys  = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);    
d3xs  = d3p_t .* c_t + d2p_t .*(-s_t) - 2 * (d2p_t .* s_t + dp_t .* c_t) - (dp_t .*c_t + p_t .*(-s_t));
d3ys  = d3p_t .* s_t + d2p_t .* c_t + 2 * (d2p_t .* c_t - dp_t .* s_t) - (dp_t .* s_t + p_t .* c_t);
ds    = sqrt(dxs.^2 + dys.^2);
H     = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
L     = length(ds)*h_bd  

var_up        = zeros(1,2*N_bd+1);    
var_up(1) = 0;
var_up(4) = 0.1;
src_info      = zeros(5,n_bd);
src_info(1,:) = xs;
src_info(2,:) = ys;
src_info(3,:) = dys./ds;
src_info(4,:) = -dxs./ds;
src_info(5,:) = ds;    
% src_old = src_info;
[srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_bd,var_up',n_bd);   
size(srcout,2)
h_bd = hout
L    = Lout
xs1   = srcout(1,:);
ys1   = srcout(2,:);
ds1   = srcout(5,:);
dxs1  = -srcout(4,:).*srcout(5,:);
dys1  = srcout(3,:).*srcout(5,:);
src_info1     = srcout;
err = max(abs(srcout(1,:) - (2+var_up(1)+var_up(4)*cos(3*tt')).*cos(tt')));
disp(err)
err = max(abs(srcout(2,:) - (2+var_up(1)+var_up(4)*cos(3*tt')).*sin(tt')));
disp(err)
dht = -3*var_up(4)*sin(3*tt');
dxs = dht .*cos(tt')-(2+var_up(1)+var_up(4)*cos(3*tt')).*sin(tt');
dys = dht .*sin(tt')+(2+var_up(1)+var_up(4)*cos(3*tt')).*cos(tt');
ds  = sqrt(dxs.^2+dys.^2);
err = max(abs(srcout(3,:) - dys./ds));
disp(err)
err = max(abs(srcout(4,:) - (-dxs./ds)));
disp(err)
return

var_up(1) = 0;
var_up(4) = 0.5;
[srcout1,hout1,Lout1,~,tt] = resample_curve(src_info1,Lout,N_bd,var_up',n_bd);    
xs2   = srcout1(1,:);
ys2   = srcout1(2,:);
ds2   = srcout1(5,:);
dxs2  = -srcout1(4,:).*srcout1(5,:);
dys2  = srcout1(3,:).*srcout1(5,:);
size(srcout1,2)
hout1
Lout1
max(abs(ys1-ys2))
max(abs(xs1-xs2))
max(abs(dys1-dys2))
max(abs(dxs1-dxs2))

% Need to rescale tt between (0,2*pi) because that was the new
% arclength parameterization of src_info1
tt = 2*pi*tt/Lout; 
rex = 2+var_up(1) + var_up(4)*cos(3*tt');
drdtex = -3*var_up(4)*sin(3*tt');
dxsex = -rex.*sin(tt') + drdtex.*cos(tt');
dysex = rex.*cos(tt') + drdtex.*sin(tt');
dsdtex = sqrt(dxsex.^2 + dysex.^2);
rnxex = dysex./dsdtex;
rnxerr = max(abs(srcout1(3,:)-rnxex));
xerr = max(abs(rex.*cos(tt')-srcout1(1,:)));
disp(rnxerr)
disp(xerr)

src = [0.01;-0.07];
targ = [12.1;5.2];
zk = complex(1.1);
uex = helm_c_p(zk,src,targ);
rhs = helm_c_p(zk,src,srcout1);

D = specdiffmat(n_bd,srcout1);
rhsdiff = D*rhs*2*pi/Lout1;

[ux,uy] = helm_c_g(zk,src,srcout1);
rhs_diff_exact =  -ux.*srcout1(4,:)' + uy.*srcout1(3,:)';
err = norm(rhsdiff-rhs_diff_exact);
fprintf("Error in diff mat = %d\n\n\n",err);

