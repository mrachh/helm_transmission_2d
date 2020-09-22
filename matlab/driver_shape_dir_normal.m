%driver 
clear

%bounday
N_bd             = 3;
coefs_bd         = zeros(1,2*N_bd+1);
coefs_bd(1)      = 1.;
coefs_bd(N_bd+1) = 0.;

%incident data frequency
khv    = 3;

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
         
    %generating the object
    n_bd  = ceil(32*kh);
    if mod(n_bd,2) 
        n_bd = n_bd+1;
    end
    if (n_bd < 128)
        n_bd = 128;
    end
    t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
    h_bd  = 2*pi/n_bd;    
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
    src_old = src_info;
    [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_bd,var_up');   
    n_bd = size(srcout,2);
    h_bd = hout;
    L    = Lout;
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
    src_info      = srcout;    

    %generating operators
    norder = 16;
    S  = slp_mat(kh,norder,h_bd,src_info);
    Sp = sprime_ext_mat(kh,norder,h_bd,src_info);
    D  = dlp_ext_mat(kh,norder,h_bd,src_info);
    Der = specdiffmat(n_bd,src_info);
    T = Der * S * Der  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',S),dys./ds) + ...
        bsxfun(@times,bsxfun(@times,(dxs./ds)',S),dxs./ds));

    %operators to target
    S_tgt = slmat_out(kh,h_bd,src,tgt);
    D_tgt = dlmat_out(kh,h_bd,src,tgt);    
    
    %fw for lambda
    eta = kh;
    Fw_mat = D + 1i * eta * S;    
    inv_Fw = inv(Fw_mat);    
    
    %bd_data
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    bd_data = -uinc;
    
    %calculating the measure
    pot = inv_Fw * bd_data;
    umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
    ubd(ik).data   = (D + 1i * eta * S) * pot ;
    dubd(ik).data  = (T + 1i * eta * Sp) * pot;
      
    % checking the domain update    
    for iup = 1: 2* N_bd+1
        %bd -> bd + delta
        var_up    = zeros(1,2*N_bd+1);
        delta = 5e-1;        
        var_up(iup) = delta;           
        
        %generating new domain        
        srcinfo1 = src_info;
        [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_bd,var_up',n_bd);                
        tout = 0:hout:(Lout-hout);
        t_bd1 = tout;
        xs1  = srcout(1,:);
        ys1  = srcout(2,:);
        dxs1 = -srcout(4,:) .* srcout(5,:);
        dys1 = srcout(3,:) .* srcout(5,:);
        ds1  = srcout(5,:);  
        pause
        %setting up the bdry for the operators
        nout = size(srcout,2);        
        src1 = zeros(4,nout);
        src1(1,:) = xs1;
        src1(2,:) = ys1;
        src1(3,:) = dxs1;
        src1(4,:) = dys1;
                
        %setting up the bdry for the operators
        src_info1 = srcout;
    
        %generating operators
        norder = 16;
        S1  = slp_mat(kh,norder,hout,src_info1);
        Sp1 = sprime_ext_mat(kh,norder,hout,src_info1);
        D1  = dlp_ext_mat(kh,norder,hout,src_info1);
        Der1 = specdiffmat(nout,src_info1);
        T1 = Der1 * S1 * Der1  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys1./ds1)',S1),dys1./ds1) + ...
            bsxfun(@times,bsxfun(@times,(dxs1./ds1)',S1),dxs1./ds1));
        
        %forward operator
        eta    = kh;    
        Fw_mat1 = D1 + 1i * eta * S1;
        inv_Fw1 = inv(Fw_mat1);

        %boundary data for new domain
        uinc  = exp(1i *kh * (bsxfun(@times,xs1',x_dir)+bsxfun(@times,ys1',y_dir)));
        bd_data1 = -uinc;
        
        %calculating scattered field at target
        S1_tgt = slmat_out(kh,hout,src1,tgt);
        D1_tgt = dlmat_out(kh,hout,src1,tgt);    
        
        %measured data
        pot1 = inv_Fw1 * bd_data1;
        umeas1(ik,iup).data = (D1_tgt + 1i * eta * S1_tgt)*pot1;
        
        % delta shape        
        t_h  = t_bd;
        h_t  = (var_up(1)+cos(bsxfun(@times,t_h',1:N_bd))*var_up(2:N_bd+1)'+...
            sin(bsxfun(@times,t_h',1:N_bd))*var_up(N_bd+2:end)')';        
        
        % right hand side        
        uinc    = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duincdn = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir)) .* uinc ./ repmat(ds',1,n_dir);
        u      = ubd(ik).data + uinc;
        dudn   = dubd(ik).data + duincdn;        
        
        h_nu = repmat(h_t', 1, n_dir );
        
        bd_data_delta3 = -h_nu .* dudn  ;
        
        bd_data_delta = bd_data_delta3;
        
        % finding potential
        poth = inv_Fw * bd_data_delta;
        
        % measure delta
        umeas_delta(ik).data = (D_tgt + 1i * eta * S_tgt)*poth;
        
        fprintf('Difference for mode %d is %d\n', iup, max(max(abs(umeas(ik).data+umeas_delta(ik).data-umeas1(ik,iup).data))))
%         pause
    end
        
end