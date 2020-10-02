%driver for inverse impedance
clear

addpath('../kress')
addpath('./src')

%set up bounday
N_bd          = 5;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1;
coefs_bd(N_bd+1) = 0.3;

%set up impedance funcion
% N_imp          = 10;
% coefs_imp      = zeros(1,2*N_imp+1);
% coefs_imp(1)   = 1;
% coefs_imp(2)   = 0.5;
% coefs_imp(11)   = 0.3;

%incident field frequency
dk    = 0.25;
n_kh  = 20;
khv   = 1:dk:n_kh*dk;

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
ifnoise   = 1;
noise_lvl = 0.02;

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
    [srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_bd,var_up');   
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
        uinc;
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
sigma = 0.1;

%RLA with Newton method
for ik = 1 : length(khv)
    
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
    
    %target_points
%     r_tgt = 10;
%     n_tgt = 5;
%     t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
%     x_t   = r_tgt * cos(t_tgt);
%     y_t   = r_tgt * sin(t_tgt);    
%     tgt   = [ x_t; y_t];
    
    if ik == 1
        %set initial guess for domain       
        N_imp      = floor(2*kh);
        var_imp    = zeros(1,2*N_imp+1);
        var_imp(1) = 1;        
    else
        N_imp_old    = N_imp;
        var_imp_old  = var_imp;        
        N_imp        = floor(2*kh);
        var_imp      = zeros(1,2*N_imp+1);
        var_imp(1)   = var_imp_old(1);
        var_imp( 2 : N_imp_old + 1 ) = var_imp_old( 2 : N_imp_old + 1 );
        var_imp( N_imp + 2 : N_imp + N_imp_old + 1 ) = var_imp_old( N_imp_old +2 : end );
    end
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    eps_step    = 1e-3;
    eps_res     = 1e-2;
    max_it      = 20;
    rhs_old     = 1e16;
    
    while flag_newton
        
        %generating the boundary       
        n_bd  = ceil(80*kh);
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
        [srcout,hout,Lout,~,tt] = resample_curve(src_info,L,N_bd,var_up');   
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
                
        %impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*var_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp))*var_imp(N_imp+2:end)')';

        %fw for lambda
        eta = kh;
        Fw_mat = T + 1i* eta * Sp  + ... %du part
                1i * kh * bsxfun(@times,lambda_imp',D+1i*eta*S);%i lambda u part  
        inv_Fw = inv(Fw_mat);    
        
        %bd_data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            uinc;
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
        
        %scattered field
        pot = inv_Fw * bd_data;
        
        uscat = (D_tgt + 1i * eta * S_tgt) * pot;        
        ubd_aux  = (D + 1i* eta *S ) * pot;
        
        %right hand side
        rhs= umeas(ik).data-uscat;
        rhs = rhs(:);
        
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
            DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
            DFw_imp(:,iimp) = DFw_col(:);
            
        end
        
%         fprintf('Condition number of inverse problem operator is %d\n',cond(DFw_imp))
        
        %finding delta
        delta = [real(DFw_imp); imag(DFw_imp)] \ [ real(rhs); imag(rhs) ];
                
        delta = delta';
        %filter        
        if (iffilter == 1)
            sol_aux   = 0;            
            hg        = 2/(2*N_imp);
            tg        = -1:hg:1;
            gauss_val = exp(-tg.*tg/sigma);
            gauss_new = zeros(1,length(gauss_val)); 
            gauss_new(1:N_imp+1) = gauss_val(N_imp+1:end);
            gauss_new(N_imp+2:end) = gauss_val(N_imp:-1:1);
            delta = delta.*gauss_new;
%             var_imp = var_imp + delta;
%             var_imp = var_imp.*gauss_new;
        end
        
        
        %update domain
        var_imp = var_imp + delta;
                
        %stopping parameter
        if norm(delta)/norm(var_imp) < eps_step
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
            var_imp = var_imp - delta;
            fprintf('RHS increasing!\n')
        end
        rhs_old = rhs;
        
        fprintf('Iteration =%d\n',it_newton)
        fprintf('RHS =%d\n',norm(rhs)/norm(umeas(ik).data))
        fprintf('Step =%d\n',norm(rhs)/norm(umeas(ik).data))
        
        it_newton = it_newton + 1;
    end
    
    lambda_vecs(ik).coefs = var_imp;
    
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
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
    var_imp = lambda_vecs(ik).coefs;
    N_imp = (length(var_imp)-1)/2;
    lambda_imp = (var_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*var_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_bd',1:N_imp))*var_imp(N_imp+2:end)')';
    lambda_imp_orig = lambda_imp_f(t_bd);    
    h0 = plot(t_bd, lambda_imp_orig,'b');
    delete(h1)        
    h1 = plot(t_bd, lambda_imp,'r');
    pause(1);
    
end

