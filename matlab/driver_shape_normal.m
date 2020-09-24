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
N_imp          = 1;%for the test to work it has to be like this.
coefs_imp      = zeros(1,2*N_imp+1);
coefs_imp(1)   = 0.2;
% coefs_imp(N_bd+1) = 0.;

%incident data frequency
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

    %this is necessary - need to calculate the curvature
    H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );    
    
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
%     [umeas(ik).data,ubd(ik).data,dubd(ik).data,~] = forward_data_imp(kh,bd_data,src,tgt,t_bd,lambda_imp);
    pot = inv_Fw * bd_data;
    umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
    ubd(ik).data   = (D + 1i * eta * S) * pot ;
    dubd(ik).data  = (T + 1i * eta * Sp) * pot;
          
    % checking the domain update    
    for iup = 1: 2* N_bd+1
        %bd -> bd + delta
        var_up    = zeros(1,2*N_bd+1);
        delta = 1e-5;        
        var_up(iup) = delta;           
        
        %generating new domain
        [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_bd,var_up');                
        L1    = Lout;
        h_bd1 = hout;
        t_bd1 = 0:h_bd1:(L1-h_bd1);        
        xs1  = srcout(1,:);
        ys1  = srcout(2,:);
        dxs1 = -srcout(4,:) .* srcout(5,:);
        dys1 = srcout(3,:) .* srcout(5,:);
        ds1  = srcout(5,:);  
        H1   = srcout(6,:);
        
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
        rsc1 = 2*pi/L1;
        S1  = slp_mat(kh,norder,h_bd1,src_info1);
        Sp1 = sprime_ext_mat(kh,norder,h_bd1,src_info1);
        D1  = dlp_ext_mat(kh,norder,h_bd1,src_info1);
        Der1 = specdiffmat(nout,src_info1)*rsc1;        
        T1 = Der1 * S1 * Der1  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys1./ds1)',S1),dys1./ds1) + ...
            bsxfun(@times,bsxfun(@times,(dxs1./ds1)',S1),dxs1./ds1));
        
        %generating the impedance
        t_lambda1 = t_bd1*2*pi/L1;
        lambda_imp1 = (coefs_imp(1)+cos(bsxfun(@times,t_lambda1',1:N_imp))*coefs_imp(2:N_imp+1)'+...
            sin(bsxfun(@times,t_lambda1',1:N_imp))*coefs_imp(N_imp+2:end)')';

%         max(abs(lambda_imp))
        %forward operator
        eta    = kh;    
        Fw_mat1 = (T1 + 1i* eta * Sp1) + ... %du part
            1i * kh * bsxfun(@times,lambda_imp1',D1+1i*eta*S1);%i lambda u part  
        inv_Fw1 = inv(Fw_mat1);

        %boundary data for new domain
        uinc  = exp(1i *kh * (bsxfun(@times,xs1',x_dir)+bsxfun(@times,ys1',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys1',x_dir)-bsxfun(@times,dxs1',y_dir))./repmat(ds1',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs1',x_dir)+bsxfun(@times,ys1',y_dir)));
        bd_data1 = -(duinc + 1i * kh * repmat(lambda_imp1',1,n_dir) .* uinc);
        
        %calculating scattered field at target
        S1_tgt = slmat_out(kh,h_bd1,src1,tgt);
        D1_tgt = dlmat_out(kh,h_bd1,src1,tgt);    
        
        %measured data
        pot = inv_Fw1 * bd_data1;
        umeas1(ik,iup).data = (D1_tgt + 1i * eta * S1_tgt)*pot;
        
%         umeas1(ik,iup).data(1:5)
        % delta impedance  
        t_h = t_bd*2*pi/L;
        h_t   = (var_up(1)+cos(bsxfun(@times,t_h',1:N_bd))*var_up(2:N_bd+1)'+...
            sin(bsxfun(@times,t_h',1:N_bd))*var_up(N_bd+2:end)')';
%         dh_t  = (bsxfun(@times,(1:N_bd)',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(2:N_bd+1)' + ...
%             bsxfun(@times,(1:N_bd)',cos(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(N_bd+2:end)')';
%         dh2_t = (bsxfun(@times,((1:N_bd).*(1:N_bd))',-cos(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(2:N_bd+1)' + ...
%             bsxfun(@times,((1:N_bd).*(1:N_bd))',-sin(bsxfun(@times,(1:N_bd)',t_bd)))'*var_up(N_bd+2:end)')';
%         hxs    = h_t .* c_t;
%         hys    = h_t .* s_t;
%         dhxs   = dh_t .* c_t + h_t .* (-s_t);        
%         dhys   = dh_t .* s_t + h_t .* c_t;
        

        % right hand side        
        uinc    = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duincdn = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir)) .* uinc ./ repmat(ds',1,n_dir);
        u      = ubd(ik).data + uinc;
        dudn   = dubd(ik).data + duincdn;
        du     = Der*u;
                
        h_nu     = repmat(h_t',1,n_dir); 
                        
        bd_data_delta1 = kh^2 * h_nu .* u;
        
        bd_data_delta2 = Der*(h_nu.*du);

        bd_data_delta3 = -1i *kh * repmat(transpose(lambda_imp),1,n_dir) .* h_nu .* ...
            ( dudn  + 1* repmat(H',1,n_dir) .* u ) ;
        
        bd_data_delta = bd_data_delta1 + bd_data_delta2 + bd_data_delta3;
        
        % finding potential
        pot = inv_Fw * bd_data_delta;
        
        % measure delta
        umeas_delta(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
        
        fprintf('Difference for mode %d is %d\n', iup, max(max(abs(umeas(ik).data+umeas_delta(ik).data-umeas1(ik,iup).data))))
%         pause
    end
        
end