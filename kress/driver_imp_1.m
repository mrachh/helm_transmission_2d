%driver 
clear

%bounday
N_bd          = 5;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1;
coefs_bd(N_bd+1) = 0.3;

%impedance funcion
N_imp          = 2;
coefs_imp      = zeros(1,2*N_imp+1);
coefs_imp(1)   = 1;
coefs_imp(2) = 0.;

%incident data frequency
dk    = 0.5;
n_kh  = 6;
khv    = 1:dk:n_kh*dk;

for ik = 1 : n_kh-1
    %incident data
    kh = khv(ik);
    n_dir = 8;
    t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
    x_dir = cos(t_dir);
    y_dir = sin(t_dir);
    dir =[ x_dir; y_dir ];
    
    %target_points
    r_tgt = 10;
    n_tgt = 100;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
    x_t   = r_tgt * cos(t_tgt);
    y_t   = r_tgt * sin(t_tgt);    
    tgt   = [ x_t; y_t];
         
    %generating the boundary
    n_bd  = ceil(32*kh);
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 64)
        n_bd = 64;
    end
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
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
            
    %generating the impedance
    lambda_imp = (coefs_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(2:N_imp+1)'+...
        sin(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(N_imp+2:end)')';
         
    %setting up the bdry for the operators
    src = zeros(8,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;
    src(5,:) = d2xs;
    src(6,:) = d2ys;
    src(7,:) = d3xs;
    src(8,:) = d3ys;       
    
    %bd_data
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
        exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
    
    %calculating the measure
    [umeas(ik).data,ubd(ik).data,~,~] = forward_data_imp(kh,bd_data,src,tgt,t_bd,lambda_imp);
    
end

%derivative of lambda
%bounday

for ik = 1 : n_kh-1
    
    %incident data
    kh = khv(ik);
    n_dir = 8;
    t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
    x_dir = cos(t_dir);
    y_dir = sin(t_dir);
    dir =[ x_dir; y_dir ];
    
    fprintf('Wavenumber=%d\n',kh)
    
    %target_points
    r_tgt = 10;
    n_tgt = 100;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
    x_t   = r_tgt * cos(t_tgt);
    y_t   = r_tgt * sin(t_tgt);    
    tgt   = [ x_t; y_t];
         
    %generating the boundary
    n_bd  = ceil(32*kh);
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 64)
        n_bd = 64;
    end
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
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

    %setting up the bdry for the operators
    src = zeros(8,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;
    src(5,:) = d2xs;
    src(6,:) = d2ys;
    src(7,:) = d3xs;
    src(8,:) = d3ys;       
    
    %generating operators
    S  = slmat(kh,src,t_bd);
    D  = dlmat(kh,src,t_bd);
    Sp = sprimelmat(kh,src,t_bd);
    T  = dprimelmat(kh,src,t_bd);
    
    %generating the impedance
    lambda_imp = (coefs_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(2:N_imp+1)'+...
        sin(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(N_imp+2:end)')';
    
    %fw for lambda
    eta = kh;
    Fw_mat = (T + 1i* eta * (Sp - eye(n_bd)/2) + ... %du part
            1i * kh * bsxfun(@times,lambda_imp',D+eye(n_bd)/2+1i*eta*S));%i lambda u part  
    inv_Fw = inv(Fw_mat);
    
    %calculating scattered field at target
    S_tgt = slmat_out(kh,src,tgt);
    D_tgt = dlmat_out(kh,src,tgt);    
    
    %impedance funcion
    N_imp      = 2;        
    for iimp = 1: 2* N_imp+1
        %impedance data -> lambda+delta
        var_imp    = zeros(1,2*N_imp+1);
        delta = 1e-7;
        var_imp(1) = 1;    
        var_imp(iimp) = var_imp(iimp) + delta;        
        lambda_imp_1 = (var_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*var_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_bd',1:N_imp))*var_imp(N_imp+2:end)')';
        
        %forward operator
        eta    = kh;    
        Fw_mat1 = (T + 1i* eta * (Sp - eye(n_bd)/2) + ... %du part
            1i * kh * bsxfun(@times,lambda_imp_1',D+eye(n_bd)/2+1i*eta*S));%i lambda u part  
        inv_Fw1 = inv(Fw_mat1);

        %boundary data
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc + 1i * kh * repmat(lambda_imp_1',1,n_dir) .* uinc);
        
        %measured data
        pot = inv_Fw1 * bd_data;
        umeas1(ik,iimp).data = (D_tgt + 1i * eta * S_tgt)*pot;
        
        % delta impedance
        delta_imp  = zeros(1,2*N_imp+1);    
        delta_imp(iimp) = delta;        
        delta_lambda_imp = (delta_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*delta_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_bd',1:N_imp))*delta_imp(N_imp+2:end)')';
        
        % right hand side
        bd_data =  -1i * kh * repmat(delta_lambda_imp',1,n_dir) .* ...
            ( uinc + ubd(ik).data );                
        
        % finding potential
        pot = inv_Fw * bd_data;
        
        % measure delta
        umeas_delta(ik,iimp).data = (D_tgt + 1i * eta * S_tgt)*pot;
        
        fprintf('Difference for mode %d is %d\n', iimp, max(max(abs(umeas(ik).data+umeas_delta(ik,iimp).data-umeas1(ik,iimp).data))))
    end
        
end
