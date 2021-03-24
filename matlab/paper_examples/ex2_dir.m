%driver for inverse shape
clear

addpath('../src')
addpath('../')

%incident field frequency
dk    = 0.25;
n_kh  = 1;
khv   = 1:dk:(1+n_kh*dk);

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


% Newton stopping criterion
eps_step    = 1e-3;
eps_res     = 1e-3;
max_it      = 200;

bd_inv_refs = cell(max_it,n_kh);

%choose to add noise
ifnoise   = 0;
noise_lvl = 0.02;

%generating data
generate = 0;

if generate 
    %genrating    
    N_bd          = 8;
    coefs_bd      = zeros(1,2*N_bd+1);
    coefs_bd(1)   = 1.;
    coefs_bd(4)   = 0.2;
    coefs_bd(5)   = 0.02;
    coefs_bd(7)   = 0.1;
    coefs_bd(9)   = 0.1;
    n_bd = 1000;
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
    xs    = p_t .* c_t;
    ys    = p_t .* s_t;
    dxs   = dp_t .* c_t + p_t .* (-s_t);
    dys   = dp_t .* s_t + p_t .* c_t;
    ds    = sqrt(dxs.^2 + dys.^2);
    Lplane    = length(ds)*h_bd;
    d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';	
    d2xs = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
    d2ys = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);
    ds   = sqrt(dxs.^2 + dys.^2);
	H    = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
    bd_ref = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
        'dys',dys,'ds','ds','Lplane',Lplane,'d2xs',d2xs,'d2ys',d2ys,'H',H,...
        'n_bd',n_bd);
    
    


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
        if (n_bd < 600)
            n_bd = 600;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        t_bd  = 0:Lplane/n_bd:Lplane-Lplane/n_bd;    
        h_bd  = 2*pi/n_bd;
        c_t   = cos(t_bd);
        s_t   = sin(t_bd);
        p_t   = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
                sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
        dp_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';
	    d2p_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(2:N_bd+1)' + ...
                bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd(N_bd+2:end)')';	
        xs   = p_t .* c_t;
        ys   = p_t .* s_t;
        dxs  = dp_t .* c_t + p_t .* (-s_t);
        dys  = dp_t .* s_t + p_t .* c_t;
	    d2xs = d2p_t .* c_t + 2 * dp_t .* (-s_t) + p_t .* (-c_t);
        d2ys = d2p_t .* s_t + 2 * dp_t .* c_t + p_t .* (-s_t);
        ds   = sqrt(dxs.^2 + dys.^2);
	    H    = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );
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

        Sik  = slp_mat(1i*kh,norder,h_bd,src_info);
        Spik = sprime_ext_mat(1i*kh,norder,h_bd,src_info);
    
        zpars2    = complex(zeros(2,1));
        zpars2(1) = kh;
        zpars2(2) = 1j*abs(kh);    
        Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);

        %operators to target
        S_tgt = slmat_out(kh,h_bd,src,tgt);
        D_tgt = dlmat_out(kh,h_bd,src,tgt);       

        %generating the impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp_orig = inf(size(t_lambda));    

        %solving the system
        eta    = kh;
        Fw_mat = D + 1i*kh*S;
        inv_Fw = inv(Fw_mat);    

        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(uinc);

        %calculating scattered field at target
        pot = inv_Fw * bd_data;    
        umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
        
        % Analytic solution test to make sure correct data is generated
        src0 = [0.01;-0.07];
        uex = helm_c_p(kh,src0,tgt);
        uin_a = helm_c_p(kh,src0,src_info);
        dudnin_a = helm_c_gn(kh,src0,src_info);
        rhs_a = uin_a;
        sol_a = inv_Fw * rhs_a;
        utest = (D_tgt + 1i * eta * S_tgt)*sol_a;
        errs(ik) = norm(utest-uex)/norm(uex);
        fprintf('error=%d\n',errs(ik));

        %umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
       if mod(ik,15)
           save('./dir-example-data/data_k20.mat','umeas','errs','t_lambda',...
           'lambda_imp_orig',...
           'N_bd','coefs_bd','khv','bd_ref');     
       end
    end
    save('./dir-example-data/data_k20.mat','umeas','errs','t_lambda',...
      'lambda_imp_orig',...
      'N_bd','coefs_bd','khv','bd_ref');
    return
else
   load('./dir-example-data/data_k20.mat')
end

if (ifnoise == 1)
   for ik=1:length(khv)
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                * abs(umeas(ik).data).* noise./abs(noise);
   end
end

%Inverse problem
ik=1;

while ik <= n_kh
    
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    
    
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
    N_var      = floor(3*kh);
    N_imp      = floor(0.5*kh);


    %%%%%%%%%%%Impedance%%%%%%%%%%%%%%
    if ik == 1 && it_newton == 1
	   N_imp_tot  = 50;
	   var_imp    = zeros(1,2*N_imp_tot+1);
       var_imp(1) = 0.5;
    end 

    %%%%%%%%%%%Shape%%%%%%%%%%%%%%
    %generating the boundary
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

    %impedance
    t_lambda = t_bd*2*pi/L;
    lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';

     %fw for lambda
     eta = kh;
     Fw_mat = Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
        1i * kh * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
     %Fw_mat = (T + 1i* eta * Sp + ... %du part
     %           1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S));%i lambda u part
     inv_Fw = inv(Fw_mat);

     %bd_data
     uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
     duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
     bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);

     %scattered field
     pot = inv_Fw * bd_data;
     uscat =  (S_tgt + 1i * eta * D_tgt*Sik)*pot;  
     ubd_aux  = (S + 1i * eta * D*Sik) *pot;
     dubd_aux = (Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4))*pot;
     %uscat = (D_tgt + 1i * eta * S_tgt) * pot;
     %ubd_aux  = (D + 1i * eta *S ) * pot;
     %dubd_aux = (T + 1i * eta * Sp ) * pot;

     %right hand side
     rhs      = umeas(ik).data-uscat;
     rhs      = rhs(:);
     rhs_old = rhs;
     
     while flag_newton
        
        if it_newton == 1
            src_info_init = src_info;
            var_imp_init = var_imp;
        end

        fprintf('Iteration =%d\n',it_newton)
                
        %constructing matrix for inverse problem
        % delta impedance        
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
            DFw_col         = (S_tgt + 1i * eta * D_tgt*Sik)*pot;
            %DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
            DFw_var(:,ivar) = DFw_col(:);
            
        end
        
        fprintf('Condition number for impedance %d\n',cond(DFw_imp))
        fprintf('Condition number for shape %d\n',cond(DFw_var))
        Minv = [real(DFw_var) real(DFw_imp);imag(DFw_var) imag(DFw_imp)];
        fprintf('Condition number for matrix %d\n',cond(Minv))
        %finding delta
        delta = Minv \ [ real(rhs); imag(rhs) ];
                
        delta = delta';
        
        rscale_delta = min(pi/abs(kh),1);
        disp("kh="+abs(kh));
        disp("rscale="+rscale_delta);
        delta_var = delta(1:2*N_var+1)*rscale_delta;
        delta_imp = delta(2*N_var+2:end)*rscale_delta;

        %update impedance
        var_imp_old = var_imp;
        var_imp(1) = var_imp(1)+delta_imp(1);
	    if N_imp>1
	      var_imp(2:N_imp+1) = var_imp(2:N_imp+1)+delta_imp(2:N_imp+1);
	      var_imp(N_imp_tot+2:N_imp_tot+1+N_imp) = var_imp(N_imp_tot+2:N_imp_tot+1+N_imp)+delta_imp(N_imp+2:end);
        end

	%update domain
        flag_filtering = 1;
        it_filtering   = 1;
        N_it_filtering = 10;
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
                NH = floor(N_var);
                Hhat_bounded = [ Hhat(1:NH+1) Hhat(end-NH+1:end)];        
                filtering(ik).value(it_newton).nh3(it_filtering) = norm(Hhat_bounded)/norm(Hhat);
                filtering(ik).sigma(it_newton) = sigma_bd;
                ratio_k = 1 - norm(Hhat_bounded)/norm(Hhat);
		        fprintf('Ratio test=%d\n',ratio_k)
                if ratio_k<1*10^-3
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
                    gauss_val = exp(-tg.*tg/sigma_bd);
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
                gauss_val = exp(-tg.*tg/sigma_bd);
                gauss_new = zeros(1,length(gauss_val));
                gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                delta_var = delta_var.*gauss_new;                
            end
            flag_domain_step = 1;
            it_filtering     = it_filtering + 1;
        end

	if it_filtering>N_it_filtering
		src_info = srcout;
		L = Lout;
	end

	%correction of the impedance
%	var_imp0 = var_imp;
%	var_imp = imp_correction(kh,dir,umeas(ik).data,tgt,t_bd,h_bd,src_info,L,var_imp0,N_imp);
%	impedance(ik).iteration(it_newton)=norm(var_imp0-var_imp);

        %%%%%%%%%%%Shape%%%%%%%%%%%%%%
        %generating the boundary
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

        %impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';

        %fw for lambda
        eta = kh;
        Fw_mat = Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
        1i * kh * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
	%Fw_mat = (T + 1i* eta * Sp + ... %du part
        %        1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S));%i lambda u part
        inv_Fw = inv(Fw_mat);

        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);

	%scattered field
        pot = inv_Fw * bd_data;
        uscat =  (S_tgt + 1i * eta * D_tgt*Sik)*pot;
        ubd_aux  = (S + 1i * eta * D*Sik) *pot;
        dubd_aux = (Sp + 1i* kh *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4))*pot;
        %uscat = (D_tgt + 1i * eta * S_tgt) * pot;
        %ubd_aux  = (D + 1i * eta *S ) * pot;
        %dubd_aux = (T + 1i * eta * Sp ) * pot;

        %right hand side
        rhs= umeas(ik).data-uscat;
        rhs = rhs(:);
        
        bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
        'dys',dys,'L',L,'H',H,'n_bd',n_bd,'rhs',rhs,'delta_var',delta_var,...
        'delta_imp',delta_imp,'t_lambda',t_lambda,'lambda_imp',lambda_imp,...
        'sigma_bd',sigma_bd,'it_filtering',it_filtering,'fret_imp_cond',...
        cond(DFw_imp),'fret_dom_cond',cond(DFw_var),'fret_cond',cond(Minv),...
        'rank_imp',rank(DFw_imp,1e-10),'rank_dom',rank(DFw_var,1e-10),...
        'rank_fret',rank(Minv,1e-10));
        bd_inv_refs{it_newton,ik} = bd_inv_ref0; 

        
        
        iesc = 0;
        %stopping parameter
        if ~flag_domain_step

            fprintf('Step impedance= %d\n',norm(var_imp-var_imp_old)/norm(var_imp))
	        fprintf('Step domain=%d\n',norm(delta_var))
        
            %stopping parameter
            if (norm(var_imp-var_imp_old)/norm(var_imp) < eps_step) %&& (norm(delta_var) < eps_step)
              flag_newton = 0;
              iesc = 1;
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
                var_imp  = var_imp_old;
                src_info = src_info_old;
		        L  = L_info;
                iesc = 4;
                fprintf('RHS increasing! %d -> %d\n',norm(rhs_old)/norm(umeas(ik).data),norm(rhs)/norm(umeas(ik).data))
                break;
            end
            rhs_old = rhs;
        
            fprintf('RHS =%d\n',norm(rhs)/norm(umeas(ik).data))
        
            it_newton = it_newton + 1;
        else
            var_imp  = var_imp_old;
            src_info = src_info_old;
            L        = L_info;
            iesc = 5;
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
    iesc_flag(ik) = iesc;
    it_newtons(ik) = it_newton;
    rhs_mags(ik) = norm(rhs)/norm(umeas(ik).data);
    ik = ik + 1;
    
    if mod(ik,15) % ~mod(kh,10)
      disp("ik="+ik)
      save('dir-example-data/sol_3kd05ki_impnocor_n000_damp_movie.mat',...
      'bd_sols','lambda_vecs','khv','iesc_flag','it_newtons','rhs_mags',...
      'bd_inv_refs')
    end

end
disp("ik="+ik)

save('dir-example-data/sol_3kd05ki_impnocor_n000_damp_movie.mat',...
      'bd_sols','lambda_vecs','khv','iesc_flag','it_newtons','rhs_mags',...
      'bd_inv_refs')
