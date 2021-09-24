function [src_info,L,t_bd,h_bd,iesc,it_newton,rhs,step_flag, ...
          filtering,ratio_gn,ratio_sd,error_sd,error_gn,it_domain] = ...
          inverse_dirichlet_v1(data,kh,src_input,L_input,t_input,dir,tgt,...
               N_var,max_it,eps_step,eps_res,src0)
           
    %initializing variables 
    ratio_gn = 0;
    ratio_sd = 0;
    error_sd = 0;
    error_gn = 0;
    
    %setting variables
    src_info = src_input;
    L        = L_input;    
    t_bd     = t_input;
    n_dir    = size(dir,2);
    eta      = kh;
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    
    % reparametrize and resample the boundary
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
    src_old = src_info;
    [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var, ...
       zeros(2*N_var+1,1),n_bd);
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

    % setting up the bdry for the operators
    src = zeros(4,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;

    % setting up the bdry for the operators
    src_info = srcout;

    %first forward problem
    [uscat,ubd_aux,uinc,dubd_aux,duincdn,inv_Fw,D_tgt,S_tgt] = ...
        forward_dirichlet(kh,src_info,L,h_bd,dir,tgt,src0);
            
    %right hand side
    rhs= data-uscat;
    rhs = rhs(:);
    rhs_old = rhs;
        
    while flag_newton
        
        if it_newton == 1
            src_info_init = src_info;            
        end

        fprintf('Iteration =%d\n',it_newton)
                
        %constructing matrix for inverse problem               
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
            dudn   = dubd_aux + duincdn;
                                        
            h_nu     = repmat(h_t',1,n_dir); 

            bd_data_delta3 = - h_nu .* dudn;

            bd_data_delta = bd_data_delta3;
            
            %find potential
            pot = inv_Fw * bd_data_delta;
        
            % measure delta
            DFw_col         = (D_tgt + 1i * kh * S_tgt) * pot;            
            DFw_var(:,ivar) = DFw_col(:);
            
        end       

        Minv = [real(DFw_var); imag(DFw_var)];
        fprintf('Condition number for shape %d\n',cond(DFw_var))
        fprintf('Condition number for matrix %d\n',cond(Minv));        
        
        %finding delta
        delta_gn = Minv \ [ real(rhs); imag(rhs) ];
        delta_sd =  Minv' * [ real(rhs); imag(rhs) ];        

        scale_sd = norm(delta_sd);
        scale_gn = norm(delta_gn);        
        
        t = delta_sd' * delta_sd/norm(Minv*delta_sd)^2;
        delta_sd = t* delta_sd;
%         fprintf('|gn|=%d, |t*sd|=%d\n',scale_gn,t*scale_sd)

        %old values for domain
        src_info_old = src_info;  
        L_old = L;	
        n_bd_old = n_bd;
        h_bd_old = h_bd;
        t_bd_old = t_bd;

        %parameteres
        wl = 2*pi/kh;
        Nw = L/wl;
        n_bd  = ceil(40*Nw);
        if (n_bd < 300)
            n_bd = 300;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        
        sigma_bd  = 1;
        it_filtering  = 0;
        max_filtering = 10;

        delta_var_gn = delta_gn(1:2*N_var+1); 
        delta_var_sd = delta_sd(1:2*N_var+1); 
                
        while it_filtering < max_filtering
            
            it_filtering = it_filtering + 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%Checking domains%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % setting domain for gauss-newton
            delta_var = delta_var_gn';

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
            H_gn = H;
            Hhat = fft(H);
            NH = floor(N_var);
            Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];                
            ratio_gn(it_newton) = 1 - norm(Hhat_bounded)/norm(Hhat);
            fprintf('Gauss-Newton ratio of high frequency modes is %d\n',ratio_gn(it_newton))

            % setting domain for s-descent            
            delta_var = delta_var_sd';% is this the right one? That is the step size that they use to compare in doleg
            
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
            H_sd = H;
            Hhat = fft(H);
            NH = floor(N_var);
            Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];                
            ratio_sd(it_newton) = 1 - norm(Hhat_bounded)/norm(Hhat);
            fprintf('SD ratio of high frequency modes is %d\n',ratio_sd(it_newton))

            if ((ratio_gn(it_newton) < 10^-3) && (ratio_sd(it_newton) < 10^-3))

                % Newton method
                %calculating field for gn
                [uscat_gn,ubd_gn,uinc_gn,dubd_gn,duinc_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                    forward_dirichlet(kh,src_info_gn,L_gn,h_bd_gn,dir,tgt,src0);
                
                rhs_gn = data-uscat_gn;            
                rhs_gn = rhs_gn(:);
                error_gn(it_newton) = norm(rhs_gn)/norm(data(:));
                fprintf('Error GN=%d\n',norm(rhs_gn)/norm(data(:)))

                %S-descent
                %update the impedance for sd
                %calculating field for sd                
                [uscat_sd,ubd_sd,uinc_sd,dubd_sd,duinc_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                    forward_dirichlet(kh,src_info_sd,L_sd,h_bd_sd,dir,tgt,src0);
                
                rhs_sd = data-uscat_sd;
                rhs_sd = rhs_sd(:);
                error_sd(it_newton) = norm(rhs_sd)/norm(data(:));
                fprintf('Error SD=%d\n',norm(rhs_sd)/norm(data(:)))
                
                if norm(rhs_sd)<norm(rhs_gn)
                    flag_domain_step = 0;
                    fprintf('SD step - 1\n')

                    % domain
                    src_info   = src_info_sd;  
                    L = L_sd;	
                    n_bd = n_bd_sd;
                    h_bd = h_bd_sd;
                    t_bd = 0:h_bd_sd:L-h_bd_sd;     
                    H    = H_sd;
                    
                    %field and operators                    
                    uscat = uscat_sd;
                    ubd_aux = ubd_sd;
                    uinc = uinc_sd;
                    dubd_aux = dubd_sd;
                    duincdn = duinc_sd;
                    inv_Fw = inv_Fw_sd;
                    D_tgt = D_tgt_sd;
                    S_tgt = S_tgt_sd;

                    %residual
                    rhs = rhs_sd(:);

                    step_flag(it_newton) = 0;

                else
                    flag_domain_step = 0;
                    fprintf('Newton step - 1\n')

                    % domain
                    src_info   = src_info_gn;  
                    L = L_gn;	
                    n_bd = n_bd_gn;
                    h_bd = h_bd_gn;
                    t_bd = 0:h_bd_gn:L-h_bd_gn;
                    H    = H_gn;
                    
                    %field and operators            
                    uscat = uscat_gn;
                    ubd_aux = ubd_gn;
                    uinc = uinc_gn;
                    dubd_aux = dubd_gn;
                    duincdn = duinc_gn;                    
                    inv_Fw = inv_Fw_gn;
                    D_tgt = D_tgt_gn;
                    S_tgt = S_tgt_gn;

                    %residual
                    rhs = rhs_gn(:);

                    step_flag(it_newton) = 1;

                end
                filtering(it_newton) = sigma_bd;
                break;
                   
            else
                if ratio_gn(it_newton) < 10^-3

                    % Newton method
                    %calculating field for gn
                    [uscat_gn,ubd_gn,uinc_gn,dubd_gn,duinc_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                        forward_dirichlet(kh,src_info_gn,L_gn,h_bd_gn,dir,tgt,src0);
                    
                    rhs_gn = data-uscat_gn;            
                    rhs_gn = rhs_gn(:);
                    error_gn(it_newton) = norm(rhs_gn)/norm(data(:));
                    fprintf('Error GN=%d\n',norm(rhs_gn)/norm(data(:)))
                    
                    flag_domain_step = 0;
                    fprintf('Newton step - 2\n')

                    % domain
                    src_info   = src_info_gn;  
                    L = L_gn;	
                    n_bd = n_bd_gn;
                    h_bd = h_bd_gn;
                    t_bd = 0:h_bd_gn:L-h_bd_gn;
                    H    = H_gn;
                    
                    %field and operators  
                    uscat = uscat_gn;
                    ubd_aux = ubd_gn;
                    uinc = uinc_gn;
                    dubd_aux = dubd_gn;
                    duincdn = duinc_gn;
                    inv_Fw = inv_Fw_gn;
                    D_tgt = D_tgt_gn;
                    S_tgt = S_tgt_gn;

                    %residual
                    rhs = rhs_gn(:);

                    step_flag(it_newton) = 1;
                    filtering(it_newton) = sigma_bd;
                    break;
                end
                
                if ratio_sd(it_newton) < 10^-3

                    %calculating field for sd                                    
                    [uscat_sd,ubd_sd,uinc_sd,dubd_sd,duinc_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                        forward_dirichlet(kh,src_info_sd,L_sd,h_bd_sd,dir,tgt,src0);
                    
                    rhs_sd = data-uscat_sd;
                    rhs_sd = rhs_sd(:);
                    error_sd(it_newton) = norm(rhs_sd)/norm(data(:));
                    fprintf('Error SD=%d\n',norm(rhs_sd)/norm(data(:)))
                    
                    flag_domain_step = 0;
                    fprintf('SD step - 2\n')

                    % domain
                    src_info   = src_info_sd;  
                    L = L_sd;	
                    n_bd = n_bd_sd;
                    h_bd = h_bd_sd;
                    t_bd = 0:h_bd:L-h_bd;
                    H    = H_sd;
                    
                    %field and operators
                    uscat = uscat_sd;
                    ubd_aux = ubd_sd;
                    uinc = uinc_sd;
                    dubd_aux = dubd_sd;
                    duincdn = duinc_sd;                    
                    inv_Fw = inv_Fw_sd;
                    D_tgt = D_tgt_sd;
                    S_tgt = S_tgt_sd;                    
                    
                    %residual
                    rhs = rhs_sd(:);

                    step_flag(it_newton) = 0;
                    filtering(it_newton) = sigma_bd;
                    break;
                end
                
                if ((ratio_gn(it_newton) >= 10^-3) && (ratio_sd(it_newton) >= 10^-3))
                    %filtering should be used
                    sigma_bd = sigma_bd/10;
                    sol_aux   = 0;
                    hg        = 2/(2*N_var);
                    tg        = -1:hg:1;
                    gauss_val = exp(-tg.*tg/sigma_bd);
                    gauss_new = zeros(1,length(gauss_val));
                    gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                    gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                    delta_var_gn = (delta_var_gn'.*gauss_new)';
                    delta_var_sd = (delta_var_sd'.*gauss_new)';                    
                end
                
            end
        end
        
        if it_filtering == max_filtering
            % Newton method
            [uscat_gn,ubd_gn,uinc_gn,dubd_gn,duinc_gn,inv_Fw_gn,D_tgt_gn,S_tgt_gn] = ...
                forward_dirichlet(kh,src_info_gn,L_gn,h_bd_gn,dir,tgt,src0);
            
            rhs_gn = data-uscat_gn;            
            rhs_gn = rhs_gn(:);
            error_gn(it_newton) = norm(rhs_gn)/norm(data(:));
            fprintf('Error GN=%d\n',norm(rhs_gn)/norm(data(:)))

            %S-descent
            [uscat_sd,ubd_sd,uinc_sd,dubd_sd,duinc_sd,inv_Fw_sd,D_tgt_sd,S_tgt_sd] = ...
                forward_dirichlet(kh,src_info_sd,L_sd,h_bd_sd,dir,tgt,src0);            
            
            rhs_sd = data-uscat_sd;
            rhs_sd = rhs_sd(:);
            error_sd(it_newton) = norm(rhs_sd)/norm(data(:));
            fprintf('Error SD=%d\n',norm(rhs_sd)/norm(data(:)))

            if norm(rhs_sd)<norm(rhs_gn)
                flag_domain_step = 0;
                fprintf('SD step\n')

                % domain
                src_info   = src_info_sd;  
                L = L_sd;	
                n_bd = n_bd_sd;
                h_bd = h_bd_sd;
                t_bd = 0:h_bd:L-h_bd;
                H    = H_sd;

                %field and operators
                uscat = uscat_sd;
                ubd_aux = ubd_sd;
                uinc = uinc_sd;
                dubd_aux = dubd_sd;
                duincdn = duinc_sd;
                inv_Fw = inv_Fw_sd;
                D_tgt = D_tgt_sd;
                S_tgt = S_tgt_sd;   

                %residual
                rhs = rhs_sd(:);

                step_flag(it_newton) = 0;

            else
                flag_domain_step = 0;
                fprintf('Newton step\n')

                % domain
                src_info   = src_info_gn;  
                L = L_gn;	
                n_bd = n_bd_gn;
                h_bd = h_bd_gn;
                t_bd = 0:h_bd:L-h_bd;
                H    = H_gn;
                
                %field and operators 
                uscat = uscat_gn;
                ubd_aux = ubd_gn;
                uinc = uinc_gn;
                dubd_aux = dubd_gn;
                duincdn = duinc_gn;                
                inv_Fw = inv_Fw_gn;
                D_tgt = D_tgt_gn;
                S_tgt = S_tgt_gn;                

                %residual
                rhs = rhs_gn(:);

                step_flag(it_newton) = 1;

            end
            filtering(it_newton) = sigma_bd;

        end
        
        iesc = 0;
        %stopping parameter
        if ~flag_domain_step
                        
	        fprintf('Step domain=%d\n',norm(delta_var))
        
            %stopping parameter
            if (norm(delta_var) < eps_step)
                flag_newton = 0;
                iesc = 1;
                fprintf('esp_step=%d\n',eps_step)
                fprintf('norm shape step=%d\n',norm(delta_var))
                fprintf('Step shape too small!\n');
            end            
        
            if it_newton > max_it
              flag_newton = 0;
              iesc = 2;
              fprintf('Reached max iteration!\n')            
            end
        
            if norm(rhs(:))/norm(data(:)) < eps_res
               flag_newton = 0;
               iesc = 3;
               fprintf('RHS too small!\n')
            end
        
            if norm(rhs_old)<norm(rhs)                 
                src_info = src_info_old;
		        L  = L_old;
                iesc = 4;
                fprintf('RHS increasing! %d -> %d\n',norm(rhs_old)/norm(data(:)),norm(rhs)/norm(data(:)))
                it_domain(it_newton).src = src_info;
		break;
            end
            rhs_old = rhs;
        
            fprintf('RHS =%d\n',norm(rhs)/norm(data(:)))
        
            it_domain(it_newton).src = src_info;
            it_newton = it_newton + 1;
        else                     
            src_info = src_info_old;
            L        = L_old;
            iesc = 5;
            it_domain(it_newton).src = src_info;
            break;
        end
    end
    
    return
