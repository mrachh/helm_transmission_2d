%driver 
clear

addpath('../kress')
addpath('./src')

%bounday
N_bd          = 4;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1.;
coefs_bd(N_bd+1) = 0.3;

%impedance funcion
N_imp          = 5;
coefs_imp      = zeros(1,2*N_imp+1);
coefs_imp(1)   = 0.7;
coefs_imp(N_imp+1) = 0.3;

%incident data frequency
khv = 3;

% incidence directions
n_dir = 5;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
x_dir = cos(t_dir);
y_dir = sin(t_dir);
dir =[ x_dir; y_dir ];

for ik = 1 : length(khv)
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
    
    %target_points
    r_tgt = 10;
    n_tgt = 5;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
    x_t   = r_tgt * cos(t_tgt);
    y_t   = r_tgt * sin(t_tgt);    
    tgt   = [ x_t; y_t];
         
    %generating the boundary
    n_bd  = ceil(32*kh);
    if (n_bd < 500)
        n_bd = 500;
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
    %it is not calculating the curvature have to diy
       
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
    lambda_imp = (coefs_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*coefs_imp(2:N_imp+1)'+...
        sin(bsxfun(@times,t_lambda',1:N_imp))*coefs_imp(N_imp+2:end)')';

    
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
    
    %calculating the measure
%     [umeas(ik).data,ubd(ik).data,~,~] = forward_data_imp(kh,bd_data,src,tgt,t_bd,lambda_imp);
    pot = inv_Fw * bd_data;
    umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
    ubd(ik).data   = (D + 1i * eta * S) * pot ;
    dubd(ik).data  = (T + 1i * eta * Sp) * pot;
    
    %impedance funcion     
    for iimp = 1: 2* N_imp+1
        %impedance data -> lambda+delta
        var_imp    = zeros(1,2*N_imp+1);
        delta = 1e-5;        
        var_imp(iimp) = delta;        
        coefs_imp1    = coefs_imp + var_imp;
        lambda_imp_1 = (coefs_imp1(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*coefs_imp1(2:N_imp+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp))*coefs_imp1(N_imp+2:end)')';
        
        %forward operator
        eta    = kh;    
        Fw_mat1 = T + 1i* eta * Sp  + ... %du part
                1i * kh * bsxfun(@times,lambda_imp_1',D+1i*eta*S);%i lambda u part  
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
        delta_lambda_imp = (delta_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(N_imp+2:end)')';
        
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
