%driver for inverse shape
clear

addpath('../kress')
addpath('./src')

%set up bounday
N_bd          = 5;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1;
coefs_bd(N_bd+1) = 0.3;

%incident field frequency
dk    = 0.5;
n_kh  = 40;
khv   = 1:dk:n_kh*dk;

% incidence directions
n_dir = 32;
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
ifnoise   = 1;
noise_lvl = 0.01;

%genrating 
for ik = 1 : length(khv)    
    %incident data
    kh = khv(ik);
                 
    %generating the boundary
    n_bd  = ceil(100*kh);
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 300)
        n_bd = 300;
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
    
    %fw for lambda
    eta = kh;
    Fw_mat = T + 1i* eta * Sp  + ... %du part
            1i * kh * bsxfun(@times,lambda_imp_orig',D+1i*eta*S);%i lambda u part  
    inv_Fw = inv(Fw_mat);    
    
    %bd_data
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
        exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    bd_data = -(duinc + 1i * kh * repmat(lambda_imp_orig',1,n_dir) .* uinc);
    
    %calculating the measure
    pot = inv_Fw * bd_data;
    umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
    
    %add noise
    if (ifnoise == 1)
        noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
        umeas(ik).data = umeas(ik).data + noise_lvl ...
            * abs(umeas(ik).data).* noise./abs(noise);
     end
    
end

%filter parameters
iffilter = 1; 
sigma = 0.2;%works for 5
% sigma = 0.3;

%RLA with Newton method
for ik = 1 : length(khv)
    
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
        
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    eps_step    = 1e-4;
    eps_res     = 1e-2;
    max_it      = 50;
    rhs_old     = 1e16;
%     alpha       = .5;
    
    while flag_newton
        
        fprintf('Iteration =%d\n',it_newton)
        
        %generating the boundary
        n_bd  = ceil(64*kh);        
        if (n_bd < 80)
            n_bd = 80;
        end
        if ik == 1
        %set up initial domain        
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
        
        %resample
        if kh<1
            N_var      = 1;
        else
            N_var      = floor(2*kh);
        end                
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
        lambda_imp = lambda_imp_f(t_lambda);

        %fw for lambda
        eta = kh;
        Fw_mat = T + 1i* eta * Sp  + ... %du part
                1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S);%i lambda u part  
        inv_Fw = inv(Fw_mat);    
       
        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
                
        %scattered field
        pot      = inv_Fw * bd_data;        
        uscat    = (D_tgt + 1i * eta * S_tgt) * pot;        
        ubd_aux  = (D + 1i * eta * S) * pot;
        dubd_aux = (T + 1i * eta * Sp) * pot;
        
        %right hand side
        rhs= umeas(ik).data-uscat;
        rhs = rhs(:);
        
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
        
%         fprintf('Condition number of inverse problem operator is %d\n',cond(DFw_imp))
        
        %finding delta
        delta = [real(DFw_var); imag(DFw_var)] \ [ real(rhs); imag(rhs) ];
                
        delta = delta';
        %filter        
        if (iffilter == 1)
            sol_aux   = 0;            
            hg        = 2/(2*N_var);
            tg        = -1:hg:1;
            gauss_val = exp(-tg.*tg/sigma);
            gauss_new = zeros(1,length(gauss_val)); 
            gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
            gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
            delta = delta.*gauss_new;
        end
                
        %update domain
        src_update = src_info;
        alpha      = 2;
        while 1
            fprintf('alpha=%d\n',alpha)                
            alpha = alpha/2;            
            [srcout,~,Lout,~,~] = resample_curve(src_info,L,N_var,alpha*delta',n_bd);
            if issimple(srcout(1,:),srcout(2,:))                
                fprintf('alpha=%d\n',alpha)                
                break;
            end
            if alpha <.1
                alpha=0;
                fprintf('alpha=%d\n',alpha)                
                break;
                srcout = src_info;
            end
        end                       
%         [srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_var,alpha*delta',n_bd);   
        src_info = srcout;
        L_old    = L;
        L        = Lout;
                
        %stopping parameter
        %for the step
        [srcold,hout,Lout,~,tt] = resample_curve(src_old,L,N_var,zeros(2*N_var+1,1),size(src_info,2));   
        fprintf('Step =%d\n',norm(srcold(1:2,:)-src_info(1:2,:))/norm(src_info(1:2,:)))
        if norm(srcold(1:2,:)-src_info(1:2,:))/norm(src_info(1:2,:)) < eps_step            
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
            flag_newton = 0;
            src_info = src_old;
            fprintf('RHS increasing!\n')
        end
        rhs_old = rhs;
                
        fprintf('RHS =%d\n',norm(rhs)/norm(umeas(ik).data))        
        fprintf('Delta size =%d\n',norm(delta))
        
        it_newton = it_newton + 1;
    end
    
    bd_sols(ik).xs = src_info(1,:);
    bd_sols(ik).ys = src_info(2,:);
    
end

figure; hold on;
h1 = plot(0,0);
for ik=1:length(khv)
    n_bd  = ceil(32*kh);
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 64)
        n_bd = 64;
    end
    t_bd   = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
    c_t    = cos(t_bd);
    s_t    = sin(t_bd);
    %old domain
    p_t    = (coefs_bd(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(2:N_bd+1)'+...
        sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd(N_bd+2:end)')';
    xs     = p_t .* c_t;
    ys     = p_t .* s_t;
        
    %new domain    
    xs1    = bd_sols(ik).xs;
    ys1    = bd_sols(ik).ys;
        
    h0 = plot(xs, ys,'b');
    delete(h1)        
    h1 = plot(xs1, ys1,'-.k');
    pause(1);
    
end

plot(xs1, ys1,'-.k')
