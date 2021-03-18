function var_imp = imp_correction(kh,dir,umeas,tgt,t_bd,h_bd,src_info,L,var_imp0,N_imp)
        
	fprintf('Correction impedance!\n')
        %direction of incident wave
        n_dir = size(dir,2);
        x_dir = dir(1,:);
        y_dir = dir(2,:);

        %%%%%%%%%%%Shape%%%%%%%%%%%%%%        
        n_bd = size(src_info,2);
        xs   = src_info(1,:);
        ys   = src_info(2,:);
        ds   = src_info(5,:);
        dxs  = -src_info(4,:).*src_info(5,:);
        dys  = src_info(3,:).*src_info(5,:);      
        
        %setting up the bdry for the operators
        src = zeros(4,n_bd);
        src(1,:) = xs;
        src(2,:) = ys;
        src(3,:) = dxs;
        src(4,:) = dys;

        %N_imp 
        N_imp_tot = (length(var_imp0)-1)/2;
        var_imp = zeros(1,2*N_imp_tot+1);
        var_imp = var_imp0;
        t_lambda = t_bd*2*pi/L;

	%modes to search for
	%N_imp = floor(0.5*kh);
                
        %generating operators
        norder = 16;
        rsc = 2*pi/L;
        S  = slp_mat(kh,norder,h_bd,src_info);
        Sp = sprime_ext_mat(kh,norder,h_bd,src_info);
        D  = dlp_ext_mat(kh,norder,h_bd,src_info);
        Der = specdiffmat(n_bd,src_info)*rsc;    
        T = Der * S * Der  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',S),dys./ds) + ...
            bsxfun(@times,bsxfun(@times,(dxs./ds)',S),dxs./ds));

        Sik  = slp_mat(1i*kh,norder,h_bd,src_info);
        Spik = sprime_ext_mat(1i*kh,norder,h_bd,src_info);

        zpars2    = complex(zeros(2,1));
        zpars2(1) = kh;
        zpars2(2) = 1j*abs(kh);
        Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);

        %this is necessary - need to calculate the curvature
        H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );  

        %operators to target
        S_tgt = slmat_out(kh,h_bd,src,tgt);
        D_tgt = dlmat_out(kh,h_bd,src,tgt); 
        
        %incoming data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));

        %impedance
        lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
                sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';

        %fw for lambda
        eta = kh;
        Fw_mat = Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
                 1i * kh * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
        %Fw_mat = (T + 1i* eta * Sp + ... %du part
        %            1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S));%i lambda u part
        inv_Fw = inv(Fw_mat);

        bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);

        %scattered field
        pot = inv_Fw * bd_data;
        uscat =  (S_tgt + 1i * eta * D_tgt*Sik)*pot;
        ubd_aux  = (S + 1i * eta * D*Sik) *pot;	
        %uscat = (D_tgt + 1i * eta * S_tgt) * pot;
        %ubd_aux  = (D + 1i * eta *S ) * pot;

        %right hand side
        rhs= umeas-uscat;
        rhs = rhs(:);
        
        %newton variables
        flag_newton = 1;
        it_newton   = 1;
        eps_step    = 1e-2;
        eps_res     = 1e-2;
        max_it      = 3;
                
        while flag_newton
        
            fprintf('Iteration =%d\n',it_newton)
       
            %calculating the derivatve
            DFw_imp    = zeros(length(rhs),2*N_imp+1);
            for iimp = 1 : (2*N_imp+1)
                %basis delta imp
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
		DFw_col         = (S_tgt + 1i * eta * D_tgt*Sik)*pot;
                %DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
                DFw_imp(:,iimp) = DFw_col(:);

            end

            fprintf('Condition number for impedance %d\n',cond(DFw_imp))

            Minv = [real(DFw_imp);imag(DFw_imp)];
            fprintf('Condition number for matrix %d\n',cond(Minv))
            %finding delta
            delta = Minv \ [ real(rhs); imag(rhs) ];

            delta = delta';

            alpha_imp = 1;
            var_imp_old = var_imp;
            var_imp(1) = var_imp(1) + alpha_imp * delta(1);
            if N_imp>1
	       var_imp(2:N_imp+1) = var_imp(2:N_imp+1) + alpha_imp * delta(2:N_imp+1);
               var_imp(N_imp_tot+2:N_imp_tot+1+N_imp) = var_imp(N_imp_tot+2:N_imp_tot+1+N_imp) + alpha_imp * delta(N_imp+2:end);
            end
            %impedance           
            lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
                sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';

            %fw for lambda
            eta = kh;
            Fw_mat = Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
                    1i * kh * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
            %Fw_mat = (T + 1i* eta * Sp + ... %du part
            %        1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S));%i lambda u part  
            inv_Fw = inv(Fw_mat);

            bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
                
            %scattered field
            pot = inv_Fw * bd_data;        
            uscat =  (S_tgt + 1i * eta * D_tgt*Sik)*pot;
            ubd_aux  = (S + 1i * eta * D*Sik) *pot;
	    %uscat = (D_tgt + 1i * eta * S_tgt) * pot;        
            %ubd_aux  = (D + 1i * eta *S ) * pot;
        
            %right hand side
	    rhs_old = rhs;
            rhs= umeas-uscat;
            rhs = rhs(:);
            
            if (norm(delta_imp)/norm(var_imp) < eps_step)
                flag_newton = 0;
                fprintf('Step too small!\n')  
            end
        
            if it_newton > max_it
                flag_newton = 0;
                fprintf('Reached max iteration!\n')            
            end
        
            if norm(rhs)/norm(umeas(:)) < eps_res
                flag_newton = 0;
                fprintf('RHS too small!\n')
            end
        
            if norm(rhs_old)<norm(rhs)
                flag_newton = 0;
                var_imp = var_imp_old;                    
                fprintf('RHS increasing!\n')
                fprintf('RHS old = %d   --->  RHS new = %d\n', norm(rhs_old), norm(rhs))
            end           
            
	    fprintf('RHS old = %d   --->  RHS new = %d\n', norm(rhs_old), norm(rhs))
	    rhs_old = rhs;
            it_newton = it_newton + 1;
        end

	fprintf('Finished correction impedance!\n')
return

            
            
