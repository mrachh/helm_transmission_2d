function [uscat,u,uin,dudn,dudnin,Der,inv_Fw,D_tgt,S_tgt] = forward_transmission(kh,khi,rho,src_info,L,h_bd,dir,tgt)

    %setting up x_dir and y_dir
    x_dir = dir(1,:);
    y_dir = dir(2,:);
    n_dir = length(x_dir);

    %%%%%%%%%%%Shape%%%%%%%%%%%%%%
    xs   = src_info(1,:);
    ys   = src_info(2,:);
    ds   = src_info(5,:);
    dxs  = -src_info(4,:).*src_info(5,:);
    dys  = src_info(3,:).*src_info(5,:);
    n_bd = length(xs);


    a = [complex(1.0); complex(1.0)];
    b = [complex(1.0); complex(rho)];
    q = (a(1)/b(1) + a(2)/b(2))*0.5;
    
    %setting up the bdry for the operators
    src = zeros(4,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;


    %generating operators
    norder = 16;
    rsc = 2*pi/L;    
    Se  = slp_mat(kh,norder,h_bd,src_info);    
    Spe = sprime_ext_mat(kh,norder,h_bd,src_info) + eye(n_bd)/2;    
    De  = dlp_ext_mat(kh,norder,h_bd,src_info) - eye(n_bd)/2;    
    Der = specdiffmat(n_bd,src_info)*rsc;    
    Te = Der * Se * Der  + kh^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',Se),dys./ds) + ...
            bsxfun(@times,bsxfun(@times,(dxs./ds)',Se),dxs./ds));    
    %generating interior operators
    Si = slp_mat(khi,norder,h_bd,src_info);    
    Spi = sprime_ext_mat(khi,norder,h_bd,src_info) + eye(n_bd)/2;    
    Di  = dlp_ext_mat(khi,norder,h_bd,src_info) - eye(n_bd)/2;     
    Ti = Der * Si * Der  + khi^2 * (bsxfun(@times,bsxfun(@times,(dys./ds)',Si),dys./ds) + ...
        bsxfun(@times,bsxfun(@times,(dxs./ds)',Si),dxs./ds));
    
    %difference operator
    zpars2    = complex(zeros(2,1));
    zpars2(1) = khi;
    zpars2(2) = kh;
    Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);    
    
    %scattered field u_i = D_i sigma - S_i mu, u_e = (1/rho)(D_e sigma - S_e mu)
    Tmat = [ eye(n_bd)+(De/rho-Di)/q (Si-Se/rho)/q; -Tdiff eye(n_bd)+Spi-Spe];

    %operators to target
    S_tgt = slmat_out(kh,h_bd,src,tgt);
    D_tgt = dlmat_out(kh,h_bd,src,tgt);
    %fw for transmission
    inv_Fw = inv(Tmat);   
    
    %bd_data
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    duincdn = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            uinc;
    bd_data = -[uinc/q; rho*duincdn ];

    %scattered field     
    soln = inv_Fw * bd_data;   
    
    % test solution at exterior target
    uscat = (D_tgt *soln(1:n_bd,:) - S_tgt *  soln(n_bd+1:end,:))/rho;

    uin    = (Di-eye(n_bd)/2)*soln(1:n_bd,:) - Si*soln(n_bd+1:end,:);
    uex    = ((De+eye(n_bd)/2)*soln(1:n_bd,:) - Se*soln(n_bd+1:end,:))/rho;
    u      = uex + uinc;

    fprintf('Error u=%d\n',max(max(abs(uin-u)))./max(max(abs(uin))))

    dudnin = Ti*soln(1:n_bd,:) - (Spi+eye(n_bd)/2)*soln(n_bd+1:end,:);
    dudnex = (Te*soln(1:n_bd,:) - (Spe-eye(n_bd)/2)*soln(n_bd+1:end,:))/rho;
    dudn   = dudnex + duincdn;

    fprintf('Error dudn=%d\n',max(max(abs(dudnin-rho*dudn)))./max(max(abs(dudnin))))

return
