function [umeas,uscat,uinc,duscatdn,duinc,inv_Fw_comb,D_tgt,S_tgt] = ...
  forward_dirichlet_fast(kh,src_info,L,h_bd,dir,tgt,src0,work,lw)

    fprintf('Wavenumber = %d\n',kh)

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
    
    %setting up the bdry for the operators
    src = zeros(4,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;
    
    % get the forward matrix
    norder = 16;
    [S,Sp,D] = get_dir_mats(kh,norder,h_bd,src_info,work,lw);
    D = D - eye(n_bd)/2;
    Sp = Sp + eye(n_bd)/2;
    
    
    Fw_mat_comb  = eye(n_bd)/2 + D + 1i* kh * S;
    Fw_mat_green = eye(n_bd)/2 + Sp - 1i * kh *S;

    %operators to target
    S_tgt = slmat_out(kh,h_bd,src,tgt);
    D_tgt = dlmat_out(kh,h_bd,src,tgt);               

    %solving the system        
    inv_Fw_comb = inv(Fw_mat_comb); 
    inv_Fw_green = inv(Fw_mat_green);
    
    %bd_data    
    uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    bd_data_green = duinc - 1i * kh *uinc;    
    bd_data_comb  = - uinc;

    %calculating scattered field at target
    pot_green = inv_Fw_green * bd_data_green;   
    pot_comb  = inv_Fw_comb * bd_data_comb;   
    umeas_green = (D_tgt + 1i * kh * S_tgt) * pot_comb;
    umeas_comb  = - S_tgt * pot_green;
    uscat  = - uinc;
    duscatdn = pot_green;
    
    fprintf('Error Green-Comb=%d\n',norm(umeas_green(:)-umeas_comb(:))/norm(umeas_green(:)))
    umeas = umeas_green;    
    
    % Analytic solution test to make sure correct data is generated        
%      src0 = [0.01;-0.07];%other domains
%    src0 = [0.01;-1.1];%horseshoe    
%     src0 = [-1.0;0.0];%horseshoe    
    uex = helm_c_p(kh,src0,tgt);        
    
    %data for checking    
    uin_a = helm_c_p(kh,src0,src_info(1:5,:));
    rhs_a = uin_a;
    sol_a = inv_Fw_comb * rhs_a;
    
    utest = (D_tgt + 1i * kh * S_tgt) * sol_a;    
%     zs = helm_c_p(kh,src_info(1:5,:),tgt);
%     zd = helm_d_p(kh,src_info(1:5,:),tgt);
%     z  = 1i*kh*zs + zd;
%     ucomp = (z.*src_info(5,:))*sol_a*h_bd;
%     norm(utest(:)-ucomp(:))/norm(utest(:))
%     uex
%     utest
%     ucomp
    
    errs = norm(utest(:)-uex(:))/norm(uex(:));
    fprintf('Checking Relative Error = %d\n',errs)
    fprintf('Checking Absolute Error = %d\n\n',norm(utest(:)-uex(:)))
    
return
