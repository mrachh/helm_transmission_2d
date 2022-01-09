%driver for inverse shape
clear

addpath('./src')

%incident field frequency
k0    = 1.;
dk    = 0.25;
%n_kh  = 117;%->30
n_kh  = 197;%->50
khv   = k0:dk:k0+(n_kh-1)*dk;

%domain jump
rho_d = 1.;

%ratio type for contrast
% q_plus_1 = 0.1, 0.33,  0.5, 0.9, 0.99, 1.01, 1.1, 2, 3, 10
q_plus_1 = 0.99;%need dom2 with inc5 and q10

%setting domain
% dom_type = 1 circle radius 1 with center (0.3,0.3)
%            2 simple plane
%            3 cavity 1
%            4 cavity 2
dom_type = 3;

%setting incident waves and receptors
% inc_type = 1 10 inc dir 200 recep
%            2 50 inc dir 200 recep
%            3 100 inc dir 200 recep
%            4 N inc dir 4N recep
%            5 N inc dir 10N recep
inc_type = 1;

%choose to add noise
% noise_type = 0; no noise
%              1; additive 
%              2; multiplicative
noise_type =0;
noise_lvl = 0.02;

%file to save data
gen_data_file = strcat('../../test-data/data_k',int2str(ceil(khv(n_kh))),...
    '_dom',int2str(dom_type),'_inctype',int2str(inc_type),'_q_func',...
     num2str(q_plus_1),'.mat');
 
res_file = strcat('../../results/transmission/sol_k',int2str(ceil(khv(n_kh))),...
    '_dom',int2str(dom_type),'_inctype',int2str(inc_type),'_q_func',...
     num2str(q_plus_1),'_noise',int2str(noise_type),'_noiselvl',num2str(noise_lvl),'.mat');

%generating data
generate = 1;

if generate

    for ik = 1 : length(khv)
        %incident data
        kh = khv(ik);

        % incidence directions
        if inc_type == 1
            n_dir = 10;
            n_tgt = 200;        
        end
        if inc_type == 2
            n_dir = 50;
            n_tgt = 200;        
        end
        if inc_type == 3
            n_dir = 100;
            n_tgt = 200;        
        end
        if inc_type == 4
            n_dir = floor(2*kh);
            n_tgt = 4*floor(2*kh);        
        end
        if inc_type == 5
            n_dir = floor(2*kh);
            n_tgt = 10*floor(2*kh);        
        end
        
        %incident waves
        t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
        x_dir = cos(t_dir);
        y_dir = sin(t_dir);
        dir =[ x_dir; y_dir ];
        
        %receptors
        r_tgt = 10;
        t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
        x_t   = r_tgt * cos(t_tgt);
        y_t   = r_tgt * sin(t_tgt);
        tgt   = [ x_t; y_t];
        
        %frequency domain
        khi = kh*sqrt(q_plus_1); %interior
        
        %calculating length of domain
        %wl = 2*pi/kh;
        wl = min(2*pi/kh,2*pi/khi);
        if dom_type == 1
            N_bd          = 8;
            coefs_bd      = zeros(1,2*N_bd+1);
            coefs_bd(1)   = 1.;
            center_x      = 0.3;
            center_y      = 0.3;    
            n_bd  = 1000;
            t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
            h_bd  = 2*pi/n_bd;
            c_t   = cos(t_bd);
            s_t   = sin(t_bd);
            p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
                    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
            xs    = center_x + p_t .* c_t;
            ys    = center_y + p_t .* s_t;            
            xo    = xs;
            yo    = ys;
            Ldom = 2*pi;
        end        
        if dom_type == 2
            N_bd          = 8;
            coefs_bd      = zeros(1,2*N_bd+1);
            coefs_bd(1)   = 1.;
            center_x      = 0.;
            center_y      = 0.;
            coefs_bd(1)   = 1.;
            coefs_bd(4)   = 0.2;
            coefs_bd(5)   = 0.02;
            coefs_bd(7)   = 0.1;
            coefs_bd(9)   = 0.1;
            coefs_bd      = 0.9*coefs_bd;
            n_bd  = 1000;
            t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
            h_bd  = 2*pi/n_bd;
            c_t   = cos(t_bd);
            s_t   = sin(t_bd);
            p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
                    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
            dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                    bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
            xs    = center_x + p_t .* c_t;
            ys    = center_y + p_t .* s_t;
            xo    = xs;
            yo    = ys;
            dxs   = dp_t .* c_t + p_t .* (-s_t);
            dys   = dp_t .* s_t + p_t .* c_t;
            ds    = sqrt(dxs.^2+dys.^2);
            Ldom  = sum(ds)*h_bd;
        end
        if dom_type == 3
            n = 1000;
            a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
            b = pi/12;  % controls approx opening angle in radians (keep small for resonant)
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
            xs = real(xy);
            ys = imag(xy);
            xo    = xs;
            yo    = ys;
            ds = srcinfo(5,:);            
            Ldom = sum(srcinfo(5,:))*h;
        end
        if dom_type == 4
            n = 1000;
            a = 0.25;  % thickness: cannot be >1/3 otherwise not smooth to emach
            b = pi/2;  % controls approx opening angle in radians (keep small for resonant)
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
            xs = real(xy);
            ys = imag(xy);
            ds = srcinfo(5,:);            
            Ldom = sum(srcinfo(5,:))*h;
        end
        Nw = Ldom/wl;
        n_bd  = ceil(100*Nw);
        if (n_bd < 600)
            n_bd = 600;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        t_bd  = 0:Ldom/n_bd:Ldom-Ldom/n_bd;
        h_bd  = Ldom/n_bd;
        if dom_type == 1       
            t_bd = 2*pi*t_bd/Ldom;
            h_bd = 2*pi*h_bd/Ldom;
            N_bd          = 8;
            coefs_bd      = zeros(1,2*N_bd+1);
            coefs_bd(1)   = 1.;
            center_x      = 0.3;
            center_y      = 0.3;            
            c_t   = cos(t_bd);
            s_t   = sin(t_bd);
            p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
                    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
            dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                    bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
            d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                    bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';	
            xs   = center_x + p_t .* c_t;
            ys   = center_y + p_t .* s_t;
            dxs  = dp_t .* c_t + p_t .* (-s_t);
            dys  = dp_t .* s_t + p_t .* c_t;
            d2xs = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
            d2ys = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);
            ds   = sqrt(dxs.^2 + dys.^2);
            H    = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );            
            L    = length(ds)*h_bd;
        end        
        if dom_type == 2
            t_bd = 2*pi*t_bd/Ldom;
            h_bd = 2*pi*h_bd/Ldom;
            N_bd          = 8;
            coefs_bd      = zeros(1,2*N_bd+1);
            coefs_bd(1)   = 1.;
            center_x      = 0.;
            center_y      = 0.;
            coefs_bd(1)   = 1.;
            coefs_bd(4)   = 0.2;
            coefs_bd(5)   = 0.02;
            coefs_bd(7)   = 0.1;
            coefs_bd(9)   = 0.1;
            coefs_bd      = 0.9*coefs_bd;
            c_t   = cos(t_bd);
            s_t   = sin(t_bd);
            p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
                    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
            dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                    bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
            d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                    bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';	
            xs   = center_x + p_t .* c_t;
            ys   = center_y + p_t .* s_t;
            dxs  = dp_t .* c_t + p_t .* (-s_t);
            dys  = dp_t .* s_t + p_t .* c_t;
            d2xs = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
            d2ys = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);
            ds   = sqrt(dxs.^2 + dys.^2);
            H    = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
            L    = length(ds)*h_bd;
        end
        if dom_type == 3
            n = n_bd;
            a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
            b = pi/12;  % controls approx opening angle in radians (keep small for resonant)
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
            xs  = real(xy);
            ys  = imag(xy);
            ds  = abs(dxydt);
            dxs = real(dxydt);
            dys = imag(dxydt);
            L   = length(ds)*h;
            h_bd = h;
        end
        if dom_type == 4
            n = n_bd;
            a = 0.25;  % thickness: cannot be >1/3 otherwise not smooth to emach
            b = pi/2;  % controls approx opening angle in radians (keep small for resonant)
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
            xs  = real(xy);
            ys  = imag(xy);
            ds  = abs(dxydt);
            dxs = real(dxydt);
            dys = imag(dxydt);
            L   = length(ds)*h;
            h_bd = h;
        end
        
        src_info      = zeros(6,n_bd);
        src_info(1,:) = xs;
        src_info(2,:) = ys;
        src_info(3,:) = dys./ds;
        src_info(4,:) = -dxs./ds;
        src_info(5,:) = ds;
        
        zks = [complex(khi);complex(kh)];
        a = [complex(1.0); complex(1.0)];
        b = [complex(1.0); complex(rho_d)];
        q = (a(1)/b(1) + a(2)/b(2))*0.5;        
        
        %caulating field
        [uscat,~,~,~,~,~,inv_Fw,D_tgt,S_tgt] = forward_transmission(kh,khi,rho_d,src_info,L,h_bd,dir,tgt);
        umeas(ik).data = uscat;
                 
        % Analytic solution test to make sure correct data is generated
        %transmission_out
        if ((dom_type == 1) || (dom_type == 2))
            src0 = [0.01;-0.07];
        else
            src0 = [0.01;-1.2];
        end
        uex_out = helm_c_p(zks(2),src0,tgt);
        uout = helm_c_p(zks(2),src0,src_info);
        dudnout = helm_c_gn(zks(2),src0,src_info);
        rhsa =[ (uout)/q; rho_d*dudnout ];
        sola = inv_Fw * rhsa;
        utestout = (D_tgt *sola(1:n_bd) - S_tgt *  sola(n_bd+1:end))/rho_d;
        errs(ik) = norm(utestout-uex_out)/norm(uex_out);
        fprintf('For kh=%d with n_bd=%d, error=%d\n',kh,n_bd,errs(ik));
        inv_Fw = 0; D_tgt = 0; S_tgt = 0;
        if mod(ik,15) == 1
            save(gen_data_file);                        
        end
    end
    save(gen_data_file);                        
else
    load(gen_data_file);    
end

if (noise_type == 1)
   for ik=1:length(khv)
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                * abs(umeas(ik).data).* noise./abs(noise);
   end
end

if (noise_type == 2)
   for ik=1:length(khv)
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                .* noise./abs(noise);
   end
end

%Inverse problem
ik=1;

while ik <= length(khv)
    
    %incident data
    kh = khv(ik);
    khi = kh*sqrt(q_plus_1);
    rho = rho_d;
   
    if inc_type == 1
            n_dir = 10;
            n_tgt = 200;
    end
    if inc_type == 2
            n_dir = 50;
            n_tgt = 200;
    end
    if inc_type == 3
            n_dir = 100;
            n_tgt = 200;
    end
    if inc_type == 4
            n_dir = floor(2*kh);
            n_tgt = 4*floor(2*kh);
    end
    if inc_type == 5
            n_dir = floor(2*kh);
            n_tgt = 10*floor(2*kh);
    end

    %incident waves
    t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
    x_dir = cos(t_dir);
    y_dir = sin(t_dir);
    dir =[ x_dir; y_dir ];

    %receptors
    r_tgt = 10;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
    x_t   = r_tgt * cos(t_tgt);
    y_t   = r_tgt * sin(t_tgt);
    tgt   = [ x_t; y_t];

    fprintf('Wavenumber=%d\n',kh)
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    eps_step    = 1e-5;
    eps_res     = 1e-5;
    max_it      = 200;
        
    %set up initial domain        
    if ik == 1 && it_newton == 1
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
    end

    %modes to search for
    N_var = floor(3*kh);    

    %%%%%%%%%%%Shape%%%%%%%%%%%%%%
    %generating the boundary
    %wl = 2*pi/kh;
    wl = min(2*pi/kh,2*pi/khi);
    Nw = L/wl;

    n_bd  = ceil(60*Nw);
    %n_bd  = ceil(100*Nw);
    if (n_bd < 300)
        n_bd = 300;
    end
    if mod(n_bd,2)
        n_bd = n_bd+1;
    end
    %resample
    [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,zeros(2*N_var+1,1),n_bd);
    n_bd = size(srcout,2);
    h_bd = hout;
    L    = Lout;
    t_bd = 0:h_bd:L-h_bd;

    a = [complex(1.0); complex(1.0)];
    b = [complex(1.0); complex(rho)];
    q = (a(1)/b(1) + a(2)/b(2))*0.5;
    
    src_info = srcout;
            
    [uscat,u,uin,dudn,dudnin,Der,inv_Fw,D_tgt,S_tgt] = forward_transmission(kh,khi,rho,src_info,L,h_bd,dir,tgt);        
    
     %right hand side
     rhs      = umeas(ik).data-uscat;     
     rhs      = rhs(:);
     rhs_old  = rhs;

     while flag_newton
        
        if it_newton == 1
            src_info_init = src_info;            
        end

        fprintf('Iteration =%d\n',it_newton)
                        
        %delta shape
        DFw_var    = zeros(length(rhs),2*N_var+1);
        for ivar = 1 : (2*N_var+1)
            
            %basis delta shape            
            delta_var  = zeros(1,2*N_var+1);    
            delta_var(ivar) = 1;     
            t_h      = t_bd*2*pi/L;
            h_t      = (delta_var(1)+cos(bsxfun(@times,t_h',1:N_var))*delta_var(2:N_var+1)'+...
                       sin(bsxfun(@times,t_h',1:N_var))*delta_var(N_var+2:end)')';                                                 
            h_nu     = repmat(h_t',1,n_dir); 

            
            %derivative u        
            der_rhs_u  = h_nu.*(dudn - dudnin);

            %derivative dudn        
            part1      = h_nu.*(khi^2*uin -kh^2*rho*u);    
            part2      = Der*(h_nu.*(Der*(uin-rho*u)));
            der_rhs_du = part1 + part2;

            der_rhs = -[der_rhs_u/q; der_rhs_du];

            soln = inv_Fw * der_rhs;

            % measure delta and constructing the matrix
            DFw_col = (D_tgt *soln(1:n_bd,:) - S_tgt *  soln(n_bd+1:end,:))/rho;                        
            DFw_var(:,ivar) = DFw_col(:);
                    
        end
                
        Minv = [real(DFw_var); imag(DFw_var)];
        fprintf('Condition number for matrix %d\n',cond(Minv))        
        delta_gn =  Minv \ [ real(rhs); imag(rhs) ];
        delta_sd =  Minv' * [ real(rhs); imag(rhs) ];
        
        scale_sd = norm(delta_sd);
        scale_gn = norm(delta_gn);
        
        t = delta_sd' * delta_sd/norm(Minv*delta_sd)^2;               
        
        fprintf('|gn|=%d, |t*sd|=%d\n',scale_gn,t*scale_sd)

        %old values for domain
        src_info_old   = src_info;  
        L_old = L;	
        n_bd_old = n_bd;
        h_bd_old = h_bd;
        t_bd_old = t_bd;

        %parameteres
        %wl = 2*pi/kh;
        wl = min(2*pi/kh,2*pi/khi);
        Nw = L/wl;
        n_bd  = ceil(60*Nw);
        %n_bd  = ceil(100*Nw);
        if (n_bd < 300)
            n_bd = 300;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        
        sigma_bd  = 1;
        it_filtering  = 0;
        max_filtering = 10;
        
        while it_filtering < max_filtering
            
            it_filtering = it_filtering + 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%Checking domains%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % setting domain for gauss-newton
            delta_var = delta_gn';

            %resample
            [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,delta_var',n_bd);
            n_bd_gn = size(srcout,2);
            h_bd_gn = hout;
            L_gn    = Lout;
            t_bd_gn = 0:h_bd_gn:L_gn-h_bd_gn;
            src_info_gn = srcout;

            %calculating wiggling curvature
            xs  = srcout(1,:);
            ys  = srcout(2,:);
            ds  = srcout(5,:);
            dxs = -srcout(4,:).*srcout(5,:);
            dys = srcout(3,:).*srcout(5,:);
            rsc = 2*pi/L;
            Der = specdiffmat(n_bd,srcout)*rsc;
            H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );
            Hhat = fft(H);
            NH = floor(N_var);
            Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];                
            ratio_gn(ik).it(it_newton) = 1 - norm(Hhat_bounded)/norm(Hhat);
            fprintf('Gauss-Newton ratio of high frequency modes is %d\n',ratio_gn(ik).it(it_newton))

            % setting domain for s-descent
            t = delta_sd' * delta_sd/norm(Minv*delta_sd)^2;
            delta_var = t * delta_sd';% is this the right one? That is the step size that they use to compare in doleg
            
            %resample
            [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,delta_var',n_bd);
            n_bd_sd = size(srcout,2);
            h_bd_sd = hout;
            L_sd   = Lout;
            t_bd_sd = 0:h_bd_sd:L_sd-h_bd_sd;
            src_info_sd = srcout;

            %calculating wiggling curvature
            xs  = srcout(1,:);
            ys  = srcout(2,:);
            ds  = srcout(5,:);
            dxs = -srcout(4,:).*srcout(5,:);
            dys = srcout(3,:).*srcout(5,:);
            rsc = 2*pi/L;
            Der = specdiffmat(n_bd,srcout)*rsc;
            H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );
            Hhat = fft(H);
            NH = floor(N_var);
            Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];                
            ratio_sd(ik).it(it_newton) = 1 - norm(Hhat_bounded)/norm(Hhat);
            fprintf('SD ratio of high frequency modes is %d\n',ratio_sd(ik).it(it_newton))

            if ((ratio_gn(ik).it(it_newton) < 10^-3) && (ratio_sd(ik).it(it_newton) < 10^-3))

                % Newton method
                [uscat_gn,u_gn,uin_gn,dudn_gn,dudnin_gn,Der_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                    forward_transmission(kh,khi,rho,src_info_gn,L_gn,h_bd_gn,dir,tgt);        

                rhs_gn = umeas(ik).data-uscat_gn;            
                rhs_gn = rhs_gn(:);
                error_gn(ik).it(it_newton) = norm(rhs_gn)/norm(umeas(ik).data(:));
                fprintf('Error GN=%d\n',norm(rhs_gn)/norm(umeas(ik).data(:)))

                %S-descent
                [uscat_sd,u_sd,uin_sd,dudn_sd,dudnin_sd,Der_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                        forward_transmission(kh,khi,rho,src_info_sd,L_sd,h_bd_sd,dir,tgt);        

                rhs_sd = umeas(ik).data-uscat_sd;
                rhs_sd = rhs_sd(:);
                error_sd(ik).it(it_newton) = norm(rhs_sd)/norm(umeas(ik).data(:));
                fprintf('Error SD=%d\n',norm(rhs_sd)/norm(umeas(ik).data(:)))
                
                if norm(rhs_sd)<norm(rhs_gn)
                    flag_domain_step = 0;
                    fprintf('SD step\n')

                    % domain
                    src_info   = src_info_sd;  
                    L = L_sd;	
                    n_bd = n_bd_sd;
                    h_bd = h_bd_sd;
                    t_bd = 0:h_bd:L-h_bd;

                    %field and operators
                    uscat = uscat_sd;
                    u = u_sd;
                    uin = uin_sd;
                    dudn = dudn_sd;
                    dudnin = dudnin_sd;
                    Der = Der_sd;
                    inv_Fw = inv_Fw_sd;
                    D_tgt = D_tgt_sd;
                    S_tgt = S_tgt_sd;

                    %residual
                    rhs = rhs_sd(:);

                    step_flag(ik).iteration(it_newton) = 0;

                else
                    flag_domain_step = 0;
                    fprintf('Newton step\n')

                    % domain
                    src_info   = src_info_gn;  
                    L = L_gn;	
                    n_bd = n_bd_gn;
                    h_bd = h_bd_gn;
                    t_bd = 0:h_bd:L-h_bd;

                    %field and operators            
                    uscat = uscat_gn;
                    u = u_gn;
                    uin = uin_gn;
                    dudn = dudn_gn;
                    dudnin = dudnin_gn;
                    Der = Der_gn;
                    inv_Fw = inv_Fw_gn;
                    D_tgt = D_tgt_gn;
                    S_tgt = S_tgt_gn;

                    %residual
                    rhs = rhs_gn(:);

                    step_flag(ik).iteration(it_newton) = 1;

                end
                filtering(ik).sigma(it_newton) = sigma_bd;
                break;
                   
            else
		error_gn(ik).it(it_newton) = 0;
		error_sd(ik).it(it_newton) = 0;
                if ratio_gn(ik).it(it_newton) < 10^-3
                    [uscat_gn,u_gn,uin_gn,dudn_gn,dudnin_gn,Der_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                    forward_transmission(kh,khi,rho,src_info_gn,L_gn,h_bd_gn,dir,tgt);        

                    rhs_gn = umeas(ik).data-uscat_gn;            
                    rhs_gn = rhs_gn(:);
                    error_gn(ik).it(it_newton) = norm(rhs_gn)/norm(umeas(ik).data(:));
                    fprintf('Error GN=%d\n',norm(rhs_gn)/norm(umeas(ik).data(:)))
                    
                    flag_domain_step = 0;
                    fprintf('Newton step\n')

                    % domain
                    src_info   = src_info_gn;  
                    L = L_gn;	
                    n_bd = n_bd_gn;
                    h_bd = h_bd_gn;
                    t_bd = 0:h_bd:L-h_bd;

                    %field and operators            
                    uscat = uscat_gn;
                    u = u_gn;
                    uin = uin_gn;
                    dudn = dudn_gn;
                    dudnin = dudnin_gn;
                    Der = Der_gn;
                    inv_Fw = inv_Fw_gn;
                    D_tgt = D_tgt_gn;
                    S_tgt = S_tgt_gn;

                    %residual
                    rhs = rhs_gn(:);

                    step_flag(ik).iteration(it_newton) = 1;
                    filtering(ik).sigma(it_newton) = sigma_bd;
                    break;
                end
                
                if ratio_sd(ik).it(it_newton) < 10^-3
                    [uscat_sd,u_sd,uin_sd,dudn_sd,dudnin_sd,Der_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                        forward_transmission(kh,khi,rho,src_info_sd,L_sd,h_bd_sd,dir,tgt);        

                    rhs_sd = umeas(ik).data-uscat_sd;
                    rhs_sd = rhs_sd(:);
                    error_sd(ik).it(it_newton) = norm(rhs_sd)/norm(umeas(ik).data(:));
                    fprintf('Error SD=%d\n',norm(rhs_sd)/norm(umeas(ik).data(:)))
                    
                    flag_domain_step = 0;
                    fprintf('SD step\n')

                    % domain
                    src_info   = src_info_sd;  
                    L = L_sd;	
                    n_bd = n_bd_sd;
                    h_bd = h_bd_sd;
                    t_bd = 0:h_bd:L-h_bd;

                    %field and operators
                    uscat = uscat_sd;
                    u = u_sd;
                    uin = uin_sd;
                    dudn = dudn_sd;
                    dudnin = dudnin_sd;
                    Der = Der_sd;
                    inv_Fw = inv_Fw_sd;
                    D_tgt = D_tgt_sd;
                    S_tgt = S_tgt_sd;

                    %residual
                    rhs = rhs_sd(:);

                    step_flag(ik).iteration(it_newton) = 0;
                    filtering(ik).sigma(it_newton) = sigma_bd;
                    break;
                end
                
                if ((ratio_gn(ik).it(it_newton) >= 10^-3) && (ratio_sd(ik).it(it_newton) >= 10^-3))
                    %filtering should be used
                    sigma_bd = sigma_bd/10;
                    sol_aux   = 0;
                    hg        = 2/(2*N_var);
                    tg        = -1:hg:1;
                    gauss_val = exp(-tg.*tg/sigma_bd);
                    gauss_new = zeros(1,length(gauss_val));
                    gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                    gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                    delta_gn = (delta_gn'.*gauss_new)';
                    delta_sd = (delta_sd'.*gauss_new)';                    
                end
                
            end
        end
        
        if it_filtering == max_filtering
            % Newton method
            [uscat_gn,u_gn,uin_gn,dudn_gn,dudnin_gn,Der_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                forward_transmission(kh,khi,rho,src_info_gn,L_gn,h_bd_gn,dir,tgt);        

            rhs_gn = umeas(ik).data-uscat_gn;            
            rhs_gn = rhs_gn(:);
            error_gn(ik).it(it_newton) = norm(rhs_gn)/norm(umeas(ik).data(:));
            fprintf('Error GN=%d\n',norm(rhs_gn)/norm(umeas(ik).data(:)))

            %S-descent
            [uscat_sd,u_sd,uin_sd,dudn_sd,dudnin_sd,Der_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                    forward_transmission(kh,khi,rho,src_info_sd,L_sd,h_bd_sd,dir,tgt);        

            rhs_sd = umeas(ik).data-uscat_sd;
            rhs_sd = rhs_sd(:);
            error_sd(ik).it(it_newton) = norm(rhs_sd)/norm(umeas(ik).data(:));
            fprintf('Error SD=%d\n',norm(rhs_sd)/norm(umeas(ik).data(:)))

            if norm(rhs_sd)<norm(rhs_gn)
                flag_domain_step = 0;
                fprintf('SD step\n')

                % domain
                src_info   = src_info_sd;  
                L = L_sd;	
                n_bd = n_bd_sd;
                h_bd = h_bd_sd;
                t_bd = 0:h_bd:L-h_bd;

                %field and operators
                uscat = uscat_sd;
                u = u_sd;
                uin = uin_sd;
                dudn = dudn_sd;
                dudnin = dudnin_sd;
                Der = Der_sd;
                inv_Fw = inv_Fw_sd;
                D_tgt = D_tgt_sd;
                S_tgt = S_tgt_sd;

                %residual
                rhs = rhs_sd(:);

                step_flag(ik).iteration(it_newton) = 0;

            else
                flag_domain_step = 0;
                fprintf('Newton step\n')

                % domain
                src_info   = src_info_gn;  
                L = L_gn;	
                n_bd = n_bd_gn;
                h_bd = h_bd_gn;
                t_bd = 0:h_bd:L-h_bd;

                %field and operators            
                uscat = uscat_gn;
                u = u_gn;
                uin = uin_gn;
                dudn = dudn_gn;
                dudnin = dudnin_gn;
                Der = Der_gn;
                inv_Fw = inv_Fw_gn;
                D_tgt = D_tgt_gn;
                S_tgt = S_tgt_gn;

                %residual
                rhs = rhs_gn(:);

                step_flag(ik).iteration(it_newton) = 1;

            end
            filtering(ik).sigma(it_newton) = sigma_bd;
                        
        end
                
        iesc = 0;
        %stopping parameter
        if ~flag_domain_step

	        fprintf('Step domain=%d\n',norm(delta_var))
            bd_sols(ik,it_newton).xs  = src_info(1,:);
            bd_sols(ik,it_newton).ys  = src_info(2,:);
    
            %stopping parameter
            if (norm(delta_var) < eps_step)
              flag_newton = 0;
              iesc        = 1;
              disp("eps_step "+eps_step)
              disp("norm step "+norm(delta_var))
              fprintf('Step too small!\n')
            end
            
            if it_newton > max_it
              flag_newton = 0;
              iesc = 2;
              fprintf('Reached max iteration!\n')            
            end
        
            if norm(rhs)/norm(umeas(ik).data) < eps_res
               flag_newton = 0;
               iesc = 3;
               fprintf('RHS too small!\n')
            end
            if norm(rhs_old)<norm(rhs)                                
                src_info = src_info_old;                
		        L  = L_old;
                iesc = 4;
                fprintf('RHS increasing! %d -> %d\n',norm(rhs_old)/norm(umeas(ik).data(:)),norm(rhs)/norm(umeas(ik).data(:)))
                bd_sols(ik,it_newton).xs  = src_info(1,:);
                bd_sols(ik,it_newton).ys  = src_info(2,:);
                rhs = rhs_old;
                break;
            end
            rhs_old = rhs;
        
            fprintf('RHS =%d\n',norm(rhs)/norm(umeas(ik).data(:)))
        
            it_newton = it_newton + 1;
        else
            src_info = src_info_old;
            L        = L_info;
            iesc = 5;
            break;
        end
    end

    fprintf('After newton RHS =%d\n',norm(rhs)/norm(umeas(ik).data(:)))     
    bd_fsols(ik).xs  = src_info(1,:);
    bd_fsols(ik).ys  = src_info(2,:);
    bd_sols(ik).rnx = src_info(3,:);
    bd_sols(ik).rny = src_info(4,:);
    bd_sols(ik).ds  = src_info(5,:);
    iesc_flag(ik) = iesc;
    it_newtons(ik) = it_newton;
    rhs_mags(ik) = norm(rhs)/norm(umeas(ik).data);
    ik = ik + 1;
    
    if mod(ik,15) % ~mod(kh,10)
        disp("ik="+ik)
        save(res_file,...
        'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
        'error_gn','error_sd','filtering','iesc_flag','step_flag','ratio_gn',...
        'ratio_sd','xo','yo')
    end

end

save(res_file,...
'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
'error_gn','error_sd','filtering','iesc_flag','step_flag','ratio_gn',...
'ratio_sd','xo','yo')


