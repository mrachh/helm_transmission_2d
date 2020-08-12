%driver 
clear

%bounday
N_bd          = 3;
coefs_bd      = zeros(1,2*N_bd+1);
coefs_bd(1)   = 1.;
coefs_bd(N_bd+1) = 0.3;

%impedance funcion
N_imp          = 3;
coefs_imp      = zeros(1,2*N_imp+1);
coefs_imp(1)   = 0;
% coefs_bd(N_bd+1) = 0.1;

%incident data frequency
% dk    = 0.5;
% n_kh  = 2;
khv    = 3;%:dk:n_kh*dk;

% incidence directions
n_dir = 1;
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
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 128)
        n_bd = 128;
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
    H     = ( dxs .* d2ys - dys .* d2xs )  ./ ( ds.^3 );

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
    
    %operators to target
    S_tgt = slmat_out(kh,src,tgt);
    D_tgt = dlmat_out(kh,src,tgt);    
    
    %generating the impedance
    lambda_imp = (coefs_imp(1)+cos(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(2:N_imp+1)'+...
        sin(bsxfun(@times,t_bd',1:N_imp))*coefs_imp(N_imp+2:end)')';
    
    %fw for lambda
    eta = kh;
    Fw_mat = (T + 1i* eta * (Sp - eye(n_bd)/2) + ... %du part
            1i * kh * bsxfun(@times,lambda_imp',D+eye(n_bd)/2+1i*eta*S));%i lambda u part  
    inv_Fw = inv(Fw_mat);    
    
    %bd_data
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
        uinc;
    bd_data = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
    
    %calculating the measure
    [umeas(ik).data,ubd(ik).data,dubd(ik).data,~] = forward_data_imp(kh,bd_data,src,tgt,t_bd,lambda_imp);
      
    % checking the domain update    
    for iup = 1: 2* N_bd+1
        %bd -> bd + delta
        var_up    = zeros(1,2*N_bd+1);
        delta = 1e-7;        
        var_up(iup) = delta;           
        coefs_bd1 = coefs_bd + var_up;
        
        
        %generating new domain
        p_t1   = (coefs_bd1(1)+cos(bsxfun(@times,t_bd',1:N_bd))*coefs_bd1(2:N_bd+1)'+...
            sin(bsxfun(@times,t_bd',1:N_bd))*coefs_bd1(N_bd+2:end)')';
        dp_t1  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(2:N_bd+1)' + ...
            bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(N_bd+2:end)')';
        d2p_t1 = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(2:N_bd+1)' + ...
            bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(N_bd+2:end)')';
        d3p_t1 = (bsxfun(@times,((1:N_bd).*(1:N_bd).*(1:N_bd))',sin(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(2:N_bd+1)' + ...
        bsxfun(@times,((1:N_bd).*(1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*coefs_bd1(N_bd+2:end)')';
        xs1    = p_t1 .* c_t;
        ys1    = p_t1 .* s_t;
        dxs1   = dp_t1 .* c_t + p_t1 .* (-s_t);
        dys1   = dp_t1 .* s_t + p_t1 .* c_t;
        d2xs1  = d2p_t1 .* c_t + 2 * dp_t1 .* (-s_t) + p_t1 .* (-c_t);
        d2ys1  = d2p_t1 .* s_t + 2 * dp_t1 .* c_t + p_t1 .* (-s_t);    
        d3xs1  = d3p_t1 .* c_t + d2p_t1 .*(-s_t) - 2 * (d2p_t1 .* s_t + dp_t1 .* c_t) - (dp_t1 .*c_t + p_t1 .*(-s_t));
        d3ys1  = d3p_t1 .* s_t + d2p_t1 .* c_t + 2 * (d2p_t1 .* c_t - dp_t1 .* s_t) - (dp_t1 .* s_t + p_t1 .* c_t);
        ds1    = sqrt(dxs1.^2 + dys1.^2);
        
        %setting up the bdry for the operators
        src1 = zeros(8,n_bd);
        src1(1,:) = xs1;
        src1(2,:) = ys1;
        src1(3,:) = dxs1;
        src1(4,:) = dys1;
        src1(5,:) = d2xs1;
        src1(6,:) = d2ys1;
        src1(7,:) = d3xs1;
        src1(8,:) = d3ys1;       
        
        %generating operators
        S1  = slmat(kh,src1,t_bd);
        D1  = dlmat(kh,src1,t_bd);
        Sp1 = sprimelmat(kh,src1,t_bd);
        T1  = dprimelmat(kh,src1,t_bd);
        
        %forward operator
        eta    = kh;    
        Fw_mat1 = (T1 + 1i* eta * (Sp1 - eye(n_bd)/2) + ... %du part
            1i * kh * bsxfun(@times,lambda_imp',D1+eye(n_bd)/2+1i*eta*S1));%i lambda u part  
        inv_Fw1 = inv(Fw_mat1);

        %boundary data for new domain
        uinc  = exp(1i *kh * (bsxfun(@times,xs1',x_dir)+bsxfun(@times,ys1',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys1',x_dir)-bsxfun(@times,dxs1',y_dir))./repmat(ds1',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs1',x_dir)+bsxfun(@times,ys1',y_dir)));
        bd_data1 = -(duinc + 1i * kh * repmat(lambda_imp',1,n_dir) .* uinc);
        
        %calculating scattered field at target
        S1_tgt = slmat_out(kh,src1,tgt);
        D1_tgt = dlmat_out(kh,src1,tgt);    
        
        %measured data
        pot = inv_Fw1 * bd_data1;
        umeas1(ik,iup).data = (D1_tgt + 1i * eta * S1_tgt)*pot;
        
        % delta impedance        
        h_t   = (var_up(1)+cos(bsxfun(@times,t_bd',1:N_bd))*var_up(2:N_bd+1)'+...
            sin(bsxfun(@times,t_bd',1:N_bd))*var_up(N_bd+2:end)')';
        dh_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(2:N_bd+1)' + ...
            bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(N_bd+2:end)')';
        dh2_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(2:N_bd+1)' + ...
            bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(N_bd+2:end)')';
        hxs    = h_t .* c_t;
        hys    = h_t .* s_t;
        dhxs   = dh_t .* c_t + h_t .* (-s_t);        
        dhys   = dh_t .* s_t + h_t .* c_t;
        

        % right hand side        
        uinc    = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duincdn = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir)) .* uinc ./ repmat(ds',1,n_dir);
        u      = ubd(ik).data + uinc;
        dudn   = dubd(ik).data + duincdn;
        du     = zeros(n_bd,n_dir);
        d2u    = zeros(n_bd,n_dir);
        
        for idir = 1 : n_dir 
            
            [ rduaux , rduaux2 ]  = spectral_derivatives(transpose(real(u(:,idir))));   
            [ iduaux , iduaux2 ]  = spectral_derivatives(transpose(imag(u(:,idir))));            
            du(:,idir)  = transpose(rduaux) + 1i * transpose(iduaux);
            d2u(:,idir) = transpose(rduaux2) + 1i * transpose(iduaux2);

        end
        
        h_nu = repmat((( hxs .* dys - hys .* dxs ) ./ ds)', 1, n_dir );
        
        der_ds = repmat((-( dxs .* d2xs + dys .* d2ys ) ./ ( ds.^3))',1,n_dir);
                
        der_h_nu = ( repmat((( dhxs.* dys + hxs .* d2ys - dhys .* dxs - hys.* d2xs ) ./ ds)', 1, n_dir) + ...
            repmat(( hxs .* dys - hys .* dxs )', 1, n_dir) .* der_ds );
                
        bd_data_delta1 = kh^2 * h_nu .* u;
        
        bd_data_delta2 = repmat((1./(ds))', 1, n_dir) .* der_h_nu .* repmat((1./(ds))', 1, n_dir) .* du + ...
            repmat((1./(ds))', 1, n_dir) .* h_nu .* der_ds .* du + ...
            repmat((1./(ds))',1,n_dir) .* h_nu .* repmat((1./(ds))', 1, n_dir).* d2u ;

        bd_data_delta3 = -1i *kh * repmat(transpose(lambda_imp),1,n_dir) .* h_nu .* ...
            ( dudn  + 1* repmat(H',1,n_dir) .* u ) ;
        
        bd_data_delta = bd_data_delta1 + bd_data_delta2 + bd_data_delta3;
        
        % finding potential
        pot = inv_Fw * bd_data_delta;
        
        % measure delta
        umeas_delta(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
        
        fprintf('Difference for mode %d is %d\n', iup, max(max(abs(umeas(ik).data+umeas_delta(ik).data-umeas1(ik,iup).data))))
        pause
    end
        
end