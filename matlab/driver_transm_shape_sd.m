%driver for inverse shape
clear

addpath('./src')

%set up bounday
N_bd          = 8;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1.;
coefs_bd(4)   = 0.2;
coefs_bd(5)   = 0.02;
coefs_bd(7)   = 0.1;
coefs_bd(9)   = 0.1;
n_bd  = 1000;
t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
h_bd  = 2*pi/n_bd;
c_t   = cos(t_bd);
s_t   = sin(t_bd);
p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
        sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
        bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
xs    = p_t .* c_t;
ys    = p_t .* s_t;
dxs   = dp_t .* c_t + p_t .* (-s_t);
dys   = dp_t .* s_t + p_t .* c_t;
ds    = sqrt(dxs.^2+dys.^2);
Lplane    = sum(ds)*h_bd;
        
%incident field frequency
k0    = 1.;
dk    = 0.25;
n_kh  = 101;
khv   = k0:dk:k0+(n_kh-1)*dk;

%domain jump
khi_fac   = 1.2; %interior
rho   = 1.1;

% incidence directions
n_dir = 8;
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

%choose to add noise
ifnoise   = 0;
noise_lvl = 0.02;

%generating data
generate = 0;

max_it      = 30;
bd_inv_refs = cell(max_it,n_kh);


if generate 

    for ik = 1 : length(khv)    
        %incident data
        kh = khv(ik);
        khi = khi_fac*kh;
%         wl = 2*pi/kh;
%         Nw = Lplane/wl;

        %generating the boundary
        n_bd  = ceil(60*kh);    
        if (n_bd < 600)
            n_bd = 600;
        end
        if mod(n_bd,2) 
            n_bd = n_bd+1;
        end
        t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
        h_bd     = 2*pi/n_bd;
        c_t   = cos(t_bd);
        s_t   = sin(t_bd);
        p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
            sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
        dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
            bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
        d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
            bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
        xs    = p_t .* c_t;
        ys    = p_t .* s_t;
        dxs   = dp_t .* c_t + p_t .* (-s_t);
        dys   = dp_t .* s_t + p_t .* c_t;
        d2xs  = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
        d2ys  = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);        
        ds    = sqrt(dxs.^2 + dys.^2);    
        H     = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
        L    = length(ds)*h_bd;    

        %reparametrization
        var_up    = zeros(1,2*N_bd+1);    
        src_info      = zeros(6,n_bd);
        src_info(1,:) = xs;
        src_info(2,:) = ys;
        src_info(3,:) = dys./ds;
        src_info(4,:) = -dxs./ds;
        src_info(5,:) = ds;    
        src_info(6,:) = H;
        [srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_bd,var_up',n_bd);   
        n_bd = size(srcout,2);
        h_bd = hout;
        L    = Lout;
        t_bd = 0:h_bd:L-h_bd;    

        zks = [complex(khi);complex(kh)];
        a = [complex(1.0); complex(1.0)];
        b = [complex(1.0); complex(rho)];
        q = (a(1)/b(1) + a(2)/b(2))*0.5;

        
        %setting up the bdry for the operators
        src_info = srcout;  
        
        %caulating field
        [uscat,~,~,~,~,~,inv_Fw,D_tgt,S_tgt] = forward_transmission(kh,khi,rho,src_info,L,h_bd,dir,tgt);        
        umeas(ik).data = uscat;        
        
        % Analytic solution test to make sure correct data is generated  
        %transmission_out
        src0 = [0.01;-0.07];
        uex_out = helm_c_p(zks(2),src0,tgt);
        uout = helm_c_p(zks(2),src0(:,1),src_info);
        dudnout = helm_c_gn(zks(2),src0(:,1),src_info);
        rhsa =[ (uout)/q; rho*dudnout ];
        sola = inv_Fw * rhsa;        
        utestout = (D_tgt *sola(1:n_bd) - S_tgt *  sola(n_bd+1:end))/rho;
        errs(ik) = norm(utestout-uex_out)/norm(uex_out);
        fprintf('For kh=%d with n_bd=%d, error=%d\n',kh,n_bd,errs(ik));
        if mod(ik,15) 
            save('./test-data/data_plane_k10.mat','umeas','errs','khv');
        end

    end
    save('./test-data/data_plane_k10.mat','umeas','errs','khv');
else
    load('./test-data/data_plane_k10.mat')
end

%add noise
if (ifnoise)
   for ik=1:length(khv)    
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                        * abs(umeas(ik).data).* noise./abs(noise);
   end
end

 
%Inverse problem
ik=1;

while ik <= length(khv)
    
    %incident data
    kh = khv(ik);
    khi = kh*khi_fac;
    
    fprintf('Wavenumber=%d\n',kh)
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    eps_step    = 1e-5;
    eps_res     = 1e-5;
    
        
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
    wl = 2*pi/kh;
    Nw = L/wl;

    n_bd  = ceil(60*Nw);
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

%     zks = [complex(khi);complex(kh)];
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
%         fprintf('Condition number for shape %d\n',cond(DFw_var))
        fprintf('Condition number for matrix %d\n',cond(Minv))
        %finding delta %check the sign
        delta_gn =  Minv \ [ real(rhs); imag(rhs) ];
        delta_sd =  Minv' * [ real(rhs); imag(rhs) ];
        
        Delta = 0.1 * pi/(kh);
        scale_sd = norm(delta_sd);
        scale_gn = norm(delta_gn);
        
        t = delta_sd' * delta_sd/norm(Minv*delta_sd)^2;               
        
        fprintf('Delta = %d, |gn|=%d, |t*sd|=%d\n',Delta,scale_gn,t*scale_sd)
                
        %old values for domain
        src_info_old   = src_info;  
        L_old = L;	
        n_bd_old = n_bd;
        h_bd_old = h_bd;
        t_bd_old = t_bd;
        
        % Gauss-Newton
        delta_var = delta_gn';
        wl = 2*pi/kh;
        Nw = L/wl;

        n_bd  = ceil(40*Nw);
        if (n_bd < 300)
            n_bd = 300;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        %resample
        [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,delta_var',n_bd);
        n_bd_gn = size(srcout,2);
        h_bd_gn = hout;
        L_gn    = Lout;
        t_bd_gn = 0:h_bd_gn:L_gn-h_bd_gn;
        src_info_gn = srcout;
        
        [uscat_gn,u_gn,uin_gn,dudn_gn,dudnin_gn,Der_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
            forward_transmission(kh,khi,rho,src_info_gn,L_gn,h_bd_gn,dir,tgt);        
               
        %right hand side
        
        rhs_gn = umeas(ik).data-uscat_gn;
        rhs_gn = rhs_gn(:);
        err_gn = norm(rhs_gn,'fro')/norm(umeas(ik).data(:),'fro');
        fprintf('Error GN=%d\n',err_gn)
        
        % SD
        src_info   = src_info_old;  
        L = L_old;	
        n_bd = n_bd_old;
        h_bd = h_bd_old;
        t_bd = t_bd_old;
        delta_var = Delta * delta_sd' / scale_sd;
        wl = 2*pi/kh;
        Nw = L/wl;

        n_bd  = ceil(40*Nw);
        if (n_bd < 300)
            n_bd = 300;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        %resample
        [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,delta_var',n_bd);
        n_bd_sd = size(srcout,2);
        h_bd_sd = hout;
        L_sd   = Lout;
        t_bd_sd = 0:h_bd_sd:L_sd-h_bd_sd;
        src_info_sd = srcout;
        
        [uscat_sd,u_sd,uin_sd,dudn_sd,dudnin_sd,Der_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
            forward_transmission(kh,khi,rho,src_info_sd,L_sd,h_bd_sd,dir,tgt);        
               
        %right hand side
        
        rhs_sd = umeas(ik).data-uscat_sd;
        rhs_sd = rhs_sd(:);
        err_sd = norm(rhs_sd,'fro')/norm(umeas(ik).data(:),'fro');
        fprintf('Error SD=%d\n',err_sd);
        istep_newton = 1;
        if(err_sd<err_gn) istep_newton = 0; end;
        
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

        bd_inv_ref0 = struct('L',L,'n_bd',n_bd,'h_bd',h_bd,'uscat',uscat,...
         'rhs_gn',rhs_gn,'rhs_sd',rhs_sd,'istep',istep_newton,...
         'delta_gn',delta_gn,'delta_sd',delta_sd,'fret_cond',cond(Minv),...
         'src_info_sd',src_info_sd,'src_info_gn',src_info_gn);
        bd_inv_refs{it_newton,ik} = bd_inv_ref0; 
%         bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
%         'dys',dys,'L',L,'H',H,'n_bd',n_bd,'rhs',rhs,'delta_var',delta_var,...
%         'delta_imp',delta_imp,'t_lambda',t_lambda,'lambda_imp',lambda_imp,...
%         'sigma_bd',sigma_bd,'it_filtering',it_filtering,'fret_imp_cond',...
%         cond(DFw_imp),'fret_dom_cond',cond(DFw_var),'fret_cond',cond(Minv),...
%         'rank_imp',rank(DFw_imp,1e-10),'rank_dom',rank(DFw_var,1e-10),...
%         'rank_fret',rank(Minv,1e-10));

                
        iesc = 0;
        %stopping parameter
        if ~flag_domain_step

	        fprintf('Step domain=%d\n',norm(delta_var))
                
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
    bd_sols(ik).xs  = src_info(1,:);
    bd_sols(ik).ys  = src_info(2,:);
    bd_sols(ik).rnx = src_info(3,:);
    bd_sols(ik).rny = src_info(4,:);
    bd_sols(ik).ds  = src_info(5,:);
    bd_sols(ik).H   = src_info(6,:);
    iesc_flag(ik) = iesc;
    it_newtons(ik) = it_newton;
    rhs_mags(ik) = norm(rhs,'fro')/norm(umeas(ik).data,'fro');
    ik = ik + 1;
    
     if mod(ik,15) % ~mod(kh,10)
       disp("ik="+ik)
       save('test-data/sol_plane.mat',...
       'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
       'bd_inv_refs')
     end
    
     

end

save('test-data/sol_plane.mat',...
       'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
       'bd_inv_refs')
     
% disp("ik="+ik)

% save('plane_2k .mat',...
%       'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
%       'bd_inv_refs')
% exit;
figure;
plot(src_info(1,:),src_info(2,:))

t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
c_t   = cos(t_bd);
s_t   = sin(t_bd);
p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
    sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
xo    = p_t .* c_t;
yo    = p_t .* s_t;

hold on
plot(xo,yo,'r')
