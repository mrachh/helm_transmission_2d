clear

addpath('../src')
addpath('..')

rng('shuffle');

%calculating L
%set up bounday
% Generate geometry
n = 300;
a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/6;  % controls approx opening angle in radians (keep small for resonant)
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
d2xydt2 = fourierZpp(zhat,t1);
srcinfo = zeros(6,n);
srcinfo(1,:) = real(xy);
srcinfo(2,:) = imag(xy);
srcinfo(5,:) = abs(dxydt);
srcinfo(3,:) = imag(dxydt)./srcinfo(5,:);
srcinfo(4,:) = -real(dxydt)./srcinfo(5,:);
xs_orig = real(xy);
ys_orig = imag(xy);
Lorig = sum(srcinfo(5,:))*h;     

%source point for texting the operator
%src0 = [0.01;-0.07];%other domains
src0 = [0.01;-1.1];%horseshoe
% src0 = [-1.0;0.0];%horseshoe

%incident field frequency
k0    = 1.;
dk    = 0.25;
n_kh  = 97;
khv   = k0:dk:k0+(n_kh-1)*dk;

% incidence directions
n_dir = 100;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
x_dir = cos(t_dir);
y_dir = sin(t_dir);
dir =[ x_dir; y_dir ];
    
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];

% Newton stopping criterion
eps_step    = 1e-5;
eps_res     = 1e-5;
max_it      = 200;

%choose to add noise
ifnoise   = 0;
noise_lvl = 0.02;

%generating data
generate = 0;

if generate 

    for ik = 1 : length(khv)    
        %incident data
        kh = khv(ik);

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
        L     = sum(srcinfo(5,:))*h;
        src_old = src_info;    

        %calculating field measurements  
        [umeas(ik).data,~,~,~,~,~,~,~] = forward_dirichlet(kh,src_info,L,h_bd,dir,tgt,src0);
        
        %add noise
        if (ifnoise == 1)
            noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
            umeas(ik).data = umeas(ik).data + noise_lvl ...
                * abs(umeas(ik).data).* noise./abs(noise);
        end

    end
    save('./larry_01_pi6_nk097_nsource100_ntarget100.mat')
    %stop
    %save('./larry_01_pi3_nk197.mat','umeas','khv');
else
    load('./larry_01_pi6_nk097_nsource100_ntarget100.mat')
    %load('./larry_01_pi3_nk197.mat')
    %load('./larry_025_pi2_nk197.mat')
end

%regularization constants
cb = 2;

%parameter stochastic
par_kh = 0.4;

% number of runnings
N_stoc = 5;

%pause
%apply Newton method several times

tic,
for istoc =1:N_stoc

%Set up Newton domain
ik=1;

% Initialize domain to unit circle
n_bd = 300;
if mod(n_bd,2)
    n_bd = n_bd+1;
end        
t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
h_bd     = 2*pi/n_bd;

c_t   = cos(t_bd);
s_t   = sin(t_bd);
p_t   = 1;
dp_t  = 0;
d2p_t = 0;
d3p_t = 0;
xs    = p_t .* c_t;
ys    = p_t .* s_t;
dxs   = dp_t .* c_t + p_t .* (-s_t);
dys   = dp_t .* s_t + p_t .* c_t;
d2xs  = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
d2ys  = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);    
ds    = sqrt(dxs.^2 + dys.^2);
H     = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
L    = length(ds)*h_bd;  

src_info      = zeros(6,n_bd);
src_info(1,:) = xs;
src_info(2,:) = ys;
src_info(3,:) = dys./ds;
src_info(4,:) = -dxs./ds;
src_info(5,:) = ds;    
src_info(6,:) = H;

ik = 1;%turn into iteration value
flag_ik = 1;

while flag_ik %length(kh_problems)
    
    %incident data
    if ik == 1
        kh = 1;
    else    
        par_choice = rand();
        if (kh == 1)
            kh = kh + dk;       
        else
            if par_choice > par_kh
                kh = kh + dk;
            else
                kh = kh - dk;
            end
        end
    end
    ind_khv = find(kh == khv);            

    eta = kh;
    
    kh_used(ik) = kh;
    
    fprintf('Wavenumber=%d\n',kh)
        
    %modes to search for
    N_var      = floor(cb*kh);    
    
    %inverse solver
   [src_info,L,t_bd,h_bd,iesc,it_newton,rhs,step_flag(ik).iteration, ...
    filtering(ik).iteration,ratio_gn(ik).iteration,ratio_sd(ik).iteration, ...
    error_sd(ik).iteration,error_gn(ik).iteration,itdomain(ik).domain] = ...
    inverse_dirichlet_v1(umeas(ind_khv).data,kh,src_info,L,t_bd,dir,tgt,...
               N_var,max_it,eps_step,eps_res,src0);	       
	       
    fprintf('After newton RHS =%d\n',norm(rhs)/norm(umeas(ind_khv).data(:))) 
    bd_sols(ik).xs  = src_info(1,:);
    bd_sols(ik).ys  = src_info(2,:);
    bd_sols(ik).rnx = src_info(3,:);
    bd_sols(ik).rny = src_info(4,:);
    bd_sols(ik).ds  = src_info(5,:);
    bd_sols(ik).H   = src_info(6,:);
    bd_sols(ik).h_bd = h_bd;
    bd_sols(ik).L = L;
    iesc_flag(ik) = iesc;
    it_newtons(ik) = it_newton;
    rhs_mags(ik) = norm(rhs)/norm(umeas(ind_khv).data);
    ik = ik + 1;
    

    %stopping frequency
    if kh >= 25
        flag_ik = 0;
    end

end

stoc(istoc).bd     = bd_sols;
stoc(istoc).newton = it_newtons;
stoc(istoc).rhs    = rhs_mags;
stoc(istoc).iesc   = iesc_flag;
stoc(istoc).step   = step_flag;
stoc(istoc).filtering = filtering;
stoc(istoc).ratio_gn  = ratio_gn;
stoc(istoc).ratio_sd  = ratio_sd;
stoc(istoc).error_sd  = error_sd;
stoc(istoc).error_gn  = error_gn;
stoc(istoc).kh     = kh_used;

clear bd_sols it_newtons rhs_mags
clear iesc_flag step_flag filtering
clear ratio_gn ration_sd error_sd error_gn
clear kh_used

save('multi_larry_01_pi6_k25_cb2_ns1_nostoc46-50.mat')
end
toc

save('multi_larry_01_pi6_k25_cb2_ns1_nostoc46-50.mat')
exit;

