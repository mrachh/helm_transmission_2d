%driver for inverse shape
clear

addpath('./src')
addpath('..')

%incident field frequency
dk    = 0.25;
n_kh  = 360;
khv   = 1:dk:n_kh*dk;

% incidence directions
n_dir = 16;
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
load('data_t20_k120.mat')
%load('data_t70_k90.mat')
generate = 0;

if generate     
    %set up bounday
    load_image
    %genrating 
    for ik = 1 : length(khv)    
        %incident data
        kh = khv(ik);

        fprintf('Wavenumber=%d\n',kh)
        %generating the boundary
        wl = 2*pi/kh;
        Nw = Lplane/wl;

        %generating the boundary
        n_bd  = ceil(50*Nw);
        if (n_bd < 300)
            n_bd = 300;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        t_bd  = 0:Lplane/n_bd:Lplane-Lplane/n_bd;    
        xs    = fnval(t_bd,xplane);
        ys    = fnval(t_bd,yplane);
        dxs   = fnval(t_bd,dxplane);
        dys   = fnval(t_bd,dyplane);
        ds    = sqrt(dxs.^2 + dys.^2);
        H     = ones(1,n_bd);
        L     = Lplane;

        %reparametrization
        var_up    = zeros(1,2*N_bd+1);    
        src_info      = zeros(6,n_bd);
        src_info(1,:) = xs;
        src_info(2,:) = ys;
        src_info(3,:) = dys./ds;
        src_info(4,:) = -dxs./ds;
        src_info(5,:) = ds;    
        src_info(6,:) = H;
        src_old = src_info;
        [srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_bd,var_up',n_bd);   
        n_bd = size(srcout,2);
        h_bd = hout;
        L    = Lout;
        t_orig = t_bd;
        t_bd = 0:h_bd:L-h_bd;    
        xs   = srcout(1,:);
        ys   = srcout(2,:);
        ds   = srcout(5,:);
        dxs  = -srcout(4,:).*srcout(5,:);
        dys  = srcout(3,:).*srcout(5,:); 

        %setting up the bdry for the operators
        src = zeros(4,n_bd);
        src(1,:) = xs;
        src(2,:) = ys;
        src(3,:) = dxs;
        src(4,:) = dys;

        %setting up the bdry for the operators
        src_info = srcout;  

        %generating operators
        norder = 16;
        rsc = 2*pi/L;
        S  = slp_mat(kh,norder,h_bd,src_info);
        Sp = sprime_ext_mat(kh,norder,h_bd,src_info);
        D  = dlp_ext_mat(kh,norder,h_bd,src_info);
        Der = specdiffmat(n_bd,src_info)*rsc;    
        T = Der * S * Der  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',S),dys./ds) + ...
            bsxfun(@times,bsxfun(@times,(dxs./ds)',S),dxs./ds));

        %operators to target
        S_tgt = slmat_out(kh,h_bd,src,tgt);
        D_tgt = dlmat_out(kh,h_bd,src,tgt);       

        %generating the impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp_orig = lambda_imp_f(t_lambda);    

        %solving the system
        eta    = kh;    
        Fw_mat = (T + 1i* eta * Sp + ... %du part
            1i * kh * bsxfun(@times,lambda_imp_orig',D+1i*eta*S));%i lambda u part  
        inv_Fw = inv(Fw_mat);    

        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp_orig',1,n_dir) .* uinc);

        %calculating scattered field at target
        pot = inv_Fw * bd_data;    
        umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;

        %add noise
        if (ifnoise == 1)
            noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
            umeas(ik).data = umeas(ik).data + noise_lvl ...
                * abs(umeas(ik).data).* noise./abs(noise);
         end

    end
end

%filtering info
iffilter_bd = 1; 
sigma_bd = 0.1;
iffilter_imp = 1; 
sigma_imp = 0.01;

ik=1;

while ik <= 37
    
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    eps_step    = 1e-3;
    eps_res     = 1e-3;
    max_it      = 30;
    rhs_old     = 1e16;  
    
    %set up initial domain        
    if ik == 1
        if (n_bd < 280)
            n_bd = 280;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end        
        t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
        h_bd     = 2*pi/n_bd;
        c_t   = cos(t_bd+pi);
        s_t   = sin(t_bd+pi);
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

    N_var      = floor(3*kh);

    %%%%%%%%%%%Impedance%%%%%%%%%%%%%%
    if ik == 1
        %set initial guess for domain
        %impedance
        N_imp      = floor(0.5*kh);%kh
        var_imp    = zeros(1,2*N_imp+1);
        var_imp(1) = 0.5;
    else
        %impedance
        N_imp_old    = N_imp;
        var_imp_old  = var_imp;
        N_imp        = floor(0.5*kh);
        var_imp      = zeros(1,2*N_imp+1);
        var_imp(1)   = var_imp_old(1);
        var_imp( 2 : N_imp_old + 1 ) = var_imp_old( 2 : N_imp_old + 1 );
        var_imp( N_imp + 2 : N_imp + N_imp_old + 1 ) = var_imp_old( N_imp_old +2 : end );
    end 

    while flag_newton
        
        if it_newton == 1
            src_info_init = src_info;
            var_imp_init = var_imp;
        end

        fprintf('Iteration =%d\n',it_newton)
                
        %%%%%%%%%%%Shape%%%%%%%%%%%%%%
        %generating the boundary
        wl = 2*pi/kh;
        Nw = L/wl;    

        n_bd  = ceil(40*Nw);
        if (n_bd < 280)
            n_bd = 280;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end        
        %resample
        src_old = src_info;
        [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var,zeros(2*N_var+1,1),n_bd);   
        n_bd = size(srcout,2);
        h_bd = hout;
        L    = Lout;
        t_orig = t_bd;
        t_bd = 0:h_bd:L-h_bd;    
        xs   = srcout(1,:);
        ys   = srcout(2,:);
        ds   = srcout(5,:);
        dxs  = -srcout(4,:).*srcout(5,:);
        dys  = srcout(3,:).*srcout(5,:);    
        
        %setting up the bdry for the operators
        src = zeros(4,n_bd);
        src(1,:) = xs;
        src(2,:) = ys;
        src(3,:) = dxs;
        src(4,:) = dys;

        %setting up the bdry for the operators
        src_info = srcout;  

        %generating operators
        norder = 16;
        rsc = 2*pi/L;
        S  = slp_mat(kh,norder,h_bd,src_info);
        Sp = sprime_ext_mat(kh,norder,h_bd,src_info);
        D  = dlp_ext_mat(kh,norder,h_bd,src_info);
        Der = specdiffmat(n_bd,src_info)*rsc;    
        T = Der * S * Der  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',S),dys./ds) + ...
            bsxfun(@times,bsxfun(@times,(dxs./ds)',S),dxs./ds));

        %this is necessary - need to calculate the curvature
        H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );  

        %operators to target
        S_tgt = slmat_out(kh,h_bd,src,tgt);
        D_tgt = dlmat_out(kh,h_bd,src,tgt);          
                
        %impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*var_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp))*var_imp(N_imp+2:end)')';

        %fw for lambda
        eta = kh;
        Fw_mat = (T + 1i* eta * Sp + ... %du part
                1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S));%i lambda u part  
        inv_Fw = inv(Fw_mat);

        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
                
        %scattered field
        pot = inv_Fw * bd_data;        
        uscat = (D_tgt + 1i * eta * S_tgt) * pot;        
        ubd_aux  = (D + 1i * eta *S ) * pot;
        dubd_aux = (T + 1i * eta * Sp ) * pot;
        
        %right hand side
        rhs= umeas(ik).data-uscat;
        rhs = rhs(:);
        
        %constructing matrix for inverse problem
        % delta impedance        
        DFw_imp    = zeros(length(rhs),2*N_imp+1);        
        for iimp = 1 : (2*N_imp+1)
            %basis delta imp
%             t_lambda = t_bd*2*pi/L;
            delta_imp  = zeros(1,2*N_imp+1);    
            delta_imp(iimp) = 1;        
            delta_lambda_imp = (delta_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(2:N_imp+1)'+...
                                sin(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(N_imp+2:end)')';
                            
            %boundary data
            bd_data =  -1i * kh * repmat(delta_lambda_imp',1,n_dir) .* ...
                    ( uinc + ubd_aux );                
            
            %find potential
            pot = inv_Fw * bd_data;
        
            % measure delta
            DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
            DFw_imp(:,iimp) = DFw_col(:);
            
        end
        
        %delta shape
        DFw_var    = zeros(length(rhs),2*N_var+1);
        for ivar = 1 : (2*N_var+1)
            %basis delta shape            
            delta_var  = zeros(1,2*N_var+1);    
            delta_var(ivar) = 1;     
            t_h = t_bd*2*pi/L;
            h_t = (delta_var(1)+cos(bsxfun(@times,t_h',1:N_var))*delta_var(2:N_var+1)'+...
                    sin(bsxfun(@times,t_h',1:N_var))*delta_var(N_var+2:end)')';
            
            % right hand side - boundary data
            uinc    = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
            duincdn = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir)) .* uinc ./ repmat(ds',1,n_dir);            
            u      = ubd_aux + uinc;
            dudn   = dubd_aux + duincdn;
            du     = Der*u;
                                        
            h_nu     = repmat(h_t',1,n_dir); 

            bd_data_delta1 = kh^2 * h_nu .* u;

            bd_data_delta2 = Der*(h_nu.*du);

            bd_data_delta3 = -1i *kh * repmat(transpose(lambda_imp),1,n_dir) .* h_nu .* ...
                ( dudn  + 1* repmat(H',1,n_dir) .* u ) ;

            bd_data_delta = bd_data_delta1 + bd_data_delta2 + bd_data_delta3;
            
            %find potential
            pot = inv_Fw * bd_data_delta;
        
            % measure delta
            DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
            DFw_var(:,ivar) = DFw_col(:);
            
        end
        
        fprintf('Condition number for impedance %d\n',cond(DFw_imp))
        fprintf('Condition number for shape %d\n',cond(DFw_var))
        Minv = [real(DFw_var) real(DFw_imp);imag(DFw_var) imag(DFw_imp)];
        fprintf('Condition number for matrix %d\n',cond(Minv))
        %finding delta
        delta = Minv \ [ real(rhs); imag(rhs) ];
                
        delta = delta';
        
        delta_var = delta(1:2*N_var+1); 
        delta_imp = delta(2*N_var+2:end);
        
	    %impedance
        if (iffilter_imp == 1)
            sigma = sigma_imp;
            sol_aux   = 0;
            hg        = 2/(2*N_imp);
            tg        = -1:hg:1;
            gauss_val = exp(-tg.*tg/sigma);
            gauss_new = zeros(1,length(gauss_val));
            gauss_new(1:N_imp+1) = gauss_val(N_imp+1:end);
            gauss_new(N_imp+2:end) = gauss_val(N_imp:-1:1);
            delta_imp = delta_imp.*gauss_new;
        end

        %update impedance
        var_imp_old = var_imp;
        var_imp     = var_imp + delta_imp;
        
	%update domain
        flag_filtering = 1;
        it_filtering   = 1;
        N_it_filtering = 2;
        sigma_bd       = 1;
        src_info_old   = src_info;  
        L_info = L;	
        while it_filtering < N_it_filtering
        
            [srcout,~,Lout,~,~] = resample_curve(src_info,L_info,N_var,delta_var',n_bd);

            if (issimple(srcout(1,:),srcout(2,:)))
		        fprintf('Passed simple test with sigma=%d at iteration =%d\n',sigma_bd, it_filtering)                           
                L_old    = L;
                L        = Lout;
                xs  = srcout(1,:);
                ys  = srcout(2,:);
                ds  = srcout(5,:);
                dxs = -srcout(4,:).*srcout(5,:);
                dys = srcout(3,:).*srcout(5,:);
                rsc = 2*pi/L;
                Der = specdiffmat(n_bd,srcout)*rsc;
                H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );
                Hhat = fft(H);
                NH = N_var;
                Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];        
                filtering(ik).value(it_newton) = norm(Hhat_bounded)/norm(Hhat);
                filtering(ik).sigma(it_newton) = sigma_bd;
                ratio_k = 1 - norm(Hhat_bounded)/norm(Hhat);
		        fprintf('Ratio test=%d\n',ratio_k)
                if ratio_k<10^-2
                    fprintf('Passed ratio test with sigma=%d at iteration =%d with ratio=%d\n',sigma_bd, it_filtering,ratio_k)
                    src_info = srcout; 
		            L = Lout;
                    flag_domain_step = 0;
                    break;
                else
                    sigma_bd = sigma_bd/10;
                    sol_aux   = 0;
                    hg        = 2/(2*N_var);
                    tg        = -1:hg:1;
                    gauss_val = exp(-tg.*tg/sigma);
                    gauss_new = zeros(1,length(gauss_val));
                    gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                    gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                    delta_var = delta_var.*gauss_new;
                end
            else
                sigma_bd = sigma_bd/10;
                sol_aux   = 0;
                hg        = 2/(2*N_var);
                tg        = -1:hg:1;
                gauss_val = exp(-tg.*tg/sigma);
                gauss_new = zeros(1,length(gauss_val));
                gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                delta_var = delta_var.*gauss_new;                
            end
            flag_domain_step = 1;
            it_filtering     = it_filtering + 1;
        end

        %stopping parameter
        if ~flag_domain_step

           fprintf('Step impedance= %d\n',norm(delta_imp)/norm(var_imp))
	       fprintf('Step domain=%d\n',norm(delta_var))
        
            %stopping parameter
            if (norm(delta_imp)/norm(var_imp) < eps_step) && norm(delta_var) < eps_step
              flag_newton = 0;
              fprintf('Step too small!\n')  
            end
        
            if it_newton > max_it
              flag_newton = 0;
              fprintf('Reached max iteration!\n')            
            end
        
            if norm(rhs)/norm(umeas(ik).data) < eps_res
               flag_newton = 0;
               fprintf('RHS too small!\n')
            end
        
            if norm(rhs_old)<norm(rhs)                
                var_imp  = var_imp_old;
                src_info = src_info_old;
		        L = L_info;
                fprintf('RHS increasing!->%d\n',norm(rhs)/norm(umeas(ik).data))
                break;
            end
            rhs_old = rhs;
        
            fprintf('RHS =%d\n',norm(rhs)/norm(umeas(ik).data))
        
            it_newton = it_newton + 1;
        else            
            break;
        end
    end
   
    fprintf('After newton RHS =%d\n',norm(rhs)/norm(umeas(ik).data)) 
    lambda_vecs(ik).coefs = var_imp;
    bd_sols(ik).xs  = src_info(1,:);
    bd_sols(ik).ys  = src_info(2,:);
    bd_sols(ik).rnx = src_info(3,:);
    bd_sols(ik).rny = src_info(4,:);
    bd_sols(ik).ds  = src_info(5,:);
    bd_sols(ik).H   = src_info(6,:);
    bd_sols(ik).L = L;
    ik = ik + 1;
    
    %if ~mod(kh,10) 
    %        save('data_t70_testcurv.mat')
    %end

end

% save('data_t20_testcurv.mat')
% 
% for ii=1:length(khv)
%     n_lambda = 1000;
%     t_lambda = 0:2*pi/n_lambda:2*pi-2*pi/n_lambda;
%     lambda_imp_orig = lambda_imp_f(t_lambda);
%     kh     = khv(ii);
%     N_imp  = floor(kh);
%     if N_imp < 1
%         N_imp = 1;
%     end
%     lcoefs = lambda_vecs(ii).coefs;
%     lambda_imp = (lcoefs(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*lcoefs(2:N_imp+1)'+sin(bsxfun(@times,t_lambda',1:N_imp))*lcoefs(N_imp+2:end)')';
%     lambda_error(ii) = norm(lambda_imp_orig-lambda_imp);
%     fprintf('lambda error =%d\n',lambda_error(ii))
% end

