%driver for inverse shape
clear

addpath('../src')
addpath('../')

%incident field frequency
dk    = 0.25;
n_kh  = 77;
khv   = 1:dk:(1+(n_kh-1)*dk);

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
eps_step    = 1e-5;
eps_res     = 1e-5;
max_it      = 100;

bd_inv_refs = cell(max_it,n_kh);

%choose to add noise
ifnoise   = 0;
noise_lvl = 0.02;


% which inverse problem, 
% iinv_type = 1, recover shape, iinv_type = 2, recover shape  + impedance

iinv_type = 2;

%generating data
generate = 0;


% domain type, idom=1, star shaped plane
% idom=2, t20 plane, idom=3, t70 plane

idom = 1;

% bc_type, bc_type = 1, dirichet, bc_type = 2, neumann, bc type=3 impedance
% only relevant for data generation
bc_type = 3;

% iimp function, only relevant if bc_type = 3
% iimp func = 1, lambda = 1 + 0.1*cos(t) + 0.02*cos(9t);
% iimp func = 2, non smooth hat function between 0.5-0.6,
%  with discontinuity at t\pi
% iimp func = 3, smooth and scaled version of iimp func = 2, with
%     std dev of gaussian = 0.005 (impedance between 1-2);
% iimp func = 3, smooth and scaled version of iimp func = 2, with
%     std dev of gaussian = 0.05 (impedance between 1-2);

iimp_func = 1;
% file which contains the stored data if generate = 0,
% it also represents file in which data will be stored if generate=1;
gen_data_file = strcat('data/data_k',int2str(ceil(khv(n_kh))),...
    '_idom',int2str(idom),'_bctype',int2str(bc_type),'_iimp_func',...
     int2str(iimp_func),'.mat');
 
 
gen_data_file= 'data/data_k20_idom1_bctype3_iimp_func1.mat';
gen_data_file= 'data/data_k20_idom1_bctype3_iimp_func1_highres.mat';

 
ifmultifreq = 0; 
multifreqmax = 0;

cb = 2;
ci = 0.5;
scontrol = 0;
res_file = ...
strcat('data/sol_k',int2str(ceil(khv(n_kh))),...
    '_idom',int2str(idom),'_bctype',int2str(bc_type),'_iimp_func',...
     int2str(iimp_func),'_multifreqmax',int2str(ceil(multifreqmax)),...
    '_cb',num2str(cb),'_ci',num2str(ci),'_damp',num2str(scontrol),...
    '_iinv_type',int2str(iinv_type),'_additive_noise.mat');
res_file = 'data/sol_k20_idom1_bctype3_iimp_func1_highres_iiv_type2.mat';
ifexit = 0;

if generate 
    %genrating 
    
    if (idom == 1)
          load_starplane;
    end
    if (idom == 2)
         load_image20;
    end
    if(idom == 3)
        load_image70;
    end


    %genrating 
    for ik = 1 : length(khv)    
        %incident data
        kh = khv(ik);

        fprintf('Wavenumber=%d\n',kh)
        %generating the boundary
        wl = 2*pi/kh;
        Nw = Lplane/wl;

        %generating the boundary
        n_bd  = ceil(100*Nw);
        if (n_bd < 600)
            n_bd = 600;
        end
        if mod(n_bd,2)
            n_bd = n_bd+1;
        end
        t_bd  = 0:Lplane/n_bd:Lplane-Lplane/n_bd; 
        if (idom == 1)
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
        end
        
        if (idom == 2 || idom == 3)
            xs    = fnval(t_bd,xplane);
            ys    = fnval(t_bd,yplane);
            dxs   = fnval(t_bd,dxplane);
            dys   = fnval(t_bd,dyplane);
            ds    = sqrt(dxs.^2 + dys.^2);
            H     = ones(1,n_bd);
            L     = Lplane;
        end
        
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
        
        %generating the impedance
        t_lambda = t_bd*2*pi/L;
        lambda_imp_orig = lambda_imp_guru(t_lambda,iimp_func);  
        
        if(bc_type == 2)
            lambda_imp_orig = zeros(size(t_lambda));
        end
        
        % get the forward matrix
        S  = slp_mat(kh,norder,h_bd,src_info);
        Sp = sprime_ext_mat(kh,norder,h_bd,src_info);
        D  = dlp_ext_mat(kh,norder,h_bd,src_info);
        Sik  = slp_mat(1i*kh,norder,h_bd,src_info);
        Spik = sprime_ext_mat(1i*kh,norder,h_bd,src_info);
    
        zpars2    = complex(zeros(2,1));
        zpars2(1) = kh;
        zpars2(2) = 1j*abs(kh);    
        Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);
        eta    = kh;    
        if(bc_type == 2 || bc_type == 3)
           Fw_mat = Sp + 1i* kh *( Tdiff*Sik + ...
                 (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - ...
                 eye(n_bd)/4)+ ...
                 1i * kh * bsxfun(@times,lambda_imp_orig',S+1i*eta*D*Sik);
        end
       
        
        if(bc_type == 1)
            Fw_mat = D + 1i*eta*S;
        end
        

        %operators to target
        S_tgt = slmat_out(kh,h_bd,src,tgt);
        D_tgt = dlmat_out(kh,h_bd,src,tgt);       

        

        %solving the system
        
        inv_Fw = inv(Fw_mat);    

        %bd_data
        src0 = [0.01;-0.07];
        uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        uin_a = helm_c_p(kh,src0,src_info);
        dudnin_a = helm_c_gn(kh,src0,src_info);
        
        if(bc_type == 2 || bc_type == 3)
            bd_data = -(duinc + ...
              1i * kh * repmat(lambda_imp_orig',1,n_dir) .* uinc);
            rhs_a = dudnin_a + 1i*kh*lambda_imp_orig'.*uin_a;
        end
        
        if(bc_type == 1)
            bd_data = -uinc;
            rhs_a = uin_a;
        end

        %calculating scattered field at target
        pot = inv_Fw * bd_data;   
        sol_a = inv_Fw * rhs_a;
        
        
        % Analytic solution test to make sure correct data is generated
        
        uex = helm_c_p(kh,src0,tgt);
        
        
        if(bc_type == 2 || bc_type == 3)
            umeas(ik).data = (S_tgt + 1i * eta * D_tgt*Sik)*pot;
            utest = (S_tgt + 1i * eta * D_tgt*Sik)*sol_a;
        end
        
        if(bc_type == 1)
            umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;
            utest = (D_tgt + 1i * eta * S_tgt)*sol_a;
        end
        errs(ik) = norm(utest-uex)/norm(uex);
        fprintf('error=%d\n',errs(ik));
        if mod(ik,15) 
            save(gen_data_file,'umeas','errs','t_lambda',...
            'lambda_imp_orig','khv','bd_ref');
            
        end

        %umeas(ik).data = (D_tgt + 1i * eta * S_tgt)*pot;

    end
    save(gen_data_file,'umeas','errs','t_lambda',...
      'lambda_imp_orig','khv','bd_ref');
    if(ifexit) 
        exit;
    end
    return
else
   load(gen_data_file)
end



if (ifnoise == 1)
   for ik=1:length(khv)
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                * abs(umeas(ik).data).* noise./abs(noise);
   end
end

if (ifnoise == 2)
   for ik=1:length(khv)
       noise = randn(n_tgt,n_dir)+1i*randn(n_tgt,n_dir);
       umeas(ik).data = umeas(ik).data + noise_lvl ...
                .* noise./abs(noise);
   end
end


%Inverse problem
ik=1;

% Initialize domain to unit circle

n_bd = 300;
if mod(n_bd,2)
    n_bd = n_bd+1;
end        
t_bd  = 0:2*pi/n_bd:2*pi-2*pi/n_bd;
h_bd     = 2*pi/n_bd;
if(idom == 1)
    c_t   = cos(t_bd);
    s_t   = sin(t_bd);
else
    c_t = cos(t_bd+pi);
    s_t = sin(t_bd+pi);
end
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

% Setup initial impedance
if(iinv_type == 2)
    N_imp_tot  = 120;
    var_imp    = zeros(1,2*N_imp_tot+1);
    var_imp(1) = 0.5;
end



for ik=1:n_kh
    
    %incident data
    kh = khv(ik);
    
    fprintf('Wavenumber=%d\n',kh)
    
    %newton variables
    flag_newton = 1;
    it_newton   = 1;
    
    %modes to search for
    N_var      = floor(cb*kh);
    
    if(iinv_type == 2) 
        N_imp      = floor(ci*kh);
    else
        N_imp = 0;
    end


    %%%%%%%%%%%Shape%%%%%%%%%%%%%%
    % reparametrize and resample the boundary
    wl = 2*pi/kh;
    Nw = L/wl;

    n_bd  = ceil(60*Nw);
    if (n_bd < 300)
    n_bd = 300;
    end
    if mod(n_bd,2)
        n_bd = n_bd+1;
    end
    %resample
    src_old = src_info;
    [srcout,hout,Lout,~,~] = resample_curve(src_info,L,N_var, ...
       zeros(2*N_var+1,1),n_bd);
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

    % setting up the bdry for the operators
    src = zeros(4,n_bd);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;

    % setting up the bdry for the operators
    src_info = srcout;

    % generating operators
    norder = 16;
    rsc = 2*pi/L;
    Der = specdiffmat(n_bd,src_info)*rsc;
    %this is necessary - need to calculate the curvature
    H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );
    
    t_lambda = t_bd*2*pi/L;
    if(iinv_type == 2)
        lambda_imp = (var_imp(1)+...
            cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';
    else
        lambda_imp = lambda_imp_guru(t_lambda,iimp_func);
    end
    
    if( ifmultifreq && kh <= multifreqmax )
        nfreqs = ik;
        freqs = khv(1:ik);
        ifreqs = 1:ik;
    else
        nfreqs = 1;
        freqs = khv(ik);
        ifreqs = ik;
    end    
    
    fprintf('nfreqs=%d\n',nfreqs);
   
    inv_Fw_mats = complex(zeros(n_bd,n_bd,nfreqs));
    sol_to_receptor_mats = complex(zeros(n_tgt,n_bd,nfreqs));
    ubd_aux_all = complex(zeros(n_bd,n_dir,nfreqs));
    dubd_aux_all = complex(zeros(n_bd,n_dir,nfreqs));
    uinc_all = complex(zeros(n_bd,n_dir,nfreqs));
    duinc_all = complex(zeros(n_bd,n_dir,nfreqs));
    zpars2    = complex(zeros(2,1));
    rhs_all = complex(zeros(n_tgt,n_dir,nfreqs));
    rhs_recast_all = complex(zeros(n_tgt*n_dir*nfreqs,1));
    
    
    for jk=1:nfreqs
        khuse = freqs(jk);
        
        S  = slp_mat(khuse,norder,h_bd,src_info);
        Sp = sprime_ext_mat(khuse,norder,h_bd,src_info);
        D  = dlp_ext_mat(khuse,norder,h_bd,src_info);
        
        Sik  = slp_mat(1i*abs(khuse),norder,h_bd,src_info);
        Spik = sprime_ext_mat(1i*khuse,norder,h_bd,src_info);
        
        
        
        zpars2(1) = khuse;
        zpars2(2) = 1j*abs(khuse);
        Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);

    %operators to target
        S_tgt = slmat_out(khuse,h_bd,src,tgt);
        D_tgt = dlmat_out(khuse,h_bd,src,tgt);

     %fw for lambda
        eta = khuse;
        Fw_mat = Sp + 1i* khuse *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
                1i * khuse * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
        inv_Fw_mats(:,:,jk) = inv(Fw_mat);
     
     %bd_data
        uinc_all(:,:,jk)  = exp(1i *khuse * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        duinc_all(:,:,jk) = 1i* khuse * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *khuse * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
        bd_data = -(duinc_all(:,:,jk) + 1i * khuse * repmat(lambda_imp',1,n_dir) .* uinc_all(:,:,jk));
        
        sol_to_receptor_mats(:,:,jk) = S_tgt + 1i*eta*D_tgt*Sik;

     %scattered field
        pot = inv_Fw_mats(:,:,jk) * bd_data;
        uscat =  sol_to_receptor_mats(:,:,jk)*pot;  
        ubd_aux_all(:,:,jk)  = (S + 1i * eta * D*Sik) *pot;
        dubd_aux_all(:,:,jk) = (Sp + 1i* khuse *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4))*pot;
   

     %right hand side
        rhs_all(:,:,jk) = umeas(ifreqs(jk)).data-uscat;
    end
    rhs_recast_all(:) = rhs_all(:);
    rhs_old = rhs_recast_all;
     
     while flag_newton
        
        if it_newton == 1
            src_info_init = src_info;
            if(iinv_type == 2)
                var_imp_init = var_imp;
            end
        end

        fprintf('Iteration =%d\n',it_newton)
                
        %constructing matrix for inverse problem
        % delta impedance        
        if (iinv_type == 2)
            DFw_imp    = zeros(length(rhs_recast_all),2*N_imp+1); 
            for iimp = 1 : (2*N_imp+1)
                %basis delta imp
                delta_imp  = zeros(1,2*N_imp+1);    
                delta_imp(iimp) = 1;        
                delta_lambda_imp = (delta_imp(1)+cos(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(2:N_imp+1)'+...
                                sin(bsxfun(@times,t_lambda',1:N_imp))*delta_imp(N_imp+2:end)')';
                            
            %boundary data
                for jk=1:nfreqs
                    khuse = freqs(jk);
                    bd_data =  -1i * khuse * repmat(delta_lambda_imp',1,n_dir) .* ...
                        ( uinc_all(:,:,jk) + ubd_aux_all(:,:,jk) );                
            
               %find potential
            
                    pot = inv_Fw_mats(:,:,jk) * bd_data;
        
                % measure delta
                    DFw_col         = sol_to_receptor_mats(:,:,jk)*pot;
                    istart = (jk-1)*n_tgt*n_dir+1;
                    iend = jk*n_tgt*n_dir;
                
                
            %DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
                    DFw_imp(istart:iend,iimp) = DFw_col(:);
                end
            end
        end
        
        
        
        %delta shape
        DFw_var    = zeros(length(rhs_recast_all),2*N_var+1);
        for ivar = 1 : (2*N_var+1)
            %basis delta shape            
            delta_var  = zeros(1,2*N_var+1);    
            delta_var(ivar) = 1;     
            t_h = t_bd*2*pi/L;
            h_t = (delta_var(1)+cos(bsxfun(@times,t_h',1:N_var))*delta_var(2:N_var+1)'+...
                    sin(bsxfun(@times,t_h',1:N_var))*delta_var(N_var+2:end)')';
            
            % right hand side - boundary data
            for jk=1:nfreqs
                khuse = complex(freqs(jk));
                u      = ubd_aux_all(:,:,jk) + uinc_all(:,:,jk);
                dudn   = dubd_aux_all(:,:,jk) + duinc_all(:,:,jk);
                du     = Der*u;                        
                h_nu     = repmat(h_t',1,n_dir); 

                bd_data_delta1 = khuse^2 * h_nu .* u;

                bd_data_delta2 = Der*(h_nu.*du);

                bd_data_delta3 = -1i *khuse * ...
                     repmat(transpose(lambda_imp),1,n_dir) .* h_nu .* ...
                    ( dudn  + 1* repmat(H',1,n_dir) .* u ) ;

                bd_data_delta = bd_data_delta1 + bd_data_delta2 + bd_data_delta3;
            
            %find potential
                pot = inv_Fw_mats(:,:,jk) * bd_data_delta;
        
            % measure delta
                DFw_col         = sol_to_receptor_mats(:,:,jk)*pot;
                istart = (jk-1)*n_tgt*n_dir+1;
                iend = jk*n_tgt*n_dir;
            %DFw_col         = (D_tgt + 1i * eta * S_tgt)*pot;
                DFw_var(istart:iend,ivar) = DFw_col(:);
            end
        end
       
        if (iinv_type == 2)
            Minv = [real(DFw_var) real(DFw_imp);imag(DFw_var) imag(DFw_imp)];
        else
            Minv = [real(DFw_var); imag(DFw_var)];
        end
        if(nfreqs == 1)
            if (iinv_type == 2)
                fprintf('Condition number for impedance %d\n',cond(DFw_imp))
            end
            fprintf('Condition number for shape %d\n',cond(DFw_var))
            fprintf('Condition number for matrix %d\n',cond(Minv));
        end
        
        %finding delta
        delta = Minv \ [ real(rhs_recast_all); imag(rhs_recast_all) ];
        
                
        delta = delta';
        if( abs(scontrol)<0.001)
            rscale_delta = 1;
        else
           rscale_delta = min(scontrol*pi/abs(kh),1);
        end
        disp("kh="+abs(kh));
        disp("rscale="+rscale_delta);
        
        delta_var = delta(1:2*N_var+1)*rscale_delta;
        if(iinv_type == 2)
            delta_imp = delta(2*N_var+2:end)*rscale_delta;

        %update impedance
            var_imp_old = var_imp;
            var_imp(1) = var_imp(1)+delta_imp(1);
            if N_imp>1
                var_imp(2:N_imp+1) = var_imp(2:N_imp+1)+delta_imp(2:N_imp+1);
                var_imp(N_imp_tot+2:N_imp_tot+1+N_imp) = var_imp(N_imp_tot+2:N_imp_tot+1+N_imp)+delta_imp(N_imp+2:end);
            end
        end

	%update domain
        flag_filtering = 1;
        it_filtering   = 1;
        N_it_filtering = 10;
        sigma_bd       = 1;
        src_info_old   = src_info;  
        L_info = L;	
        
        delta_old = delta_var;
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
                    delta_var = delta_old.*gauss_new;
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
                delta_var = delta_old.*gauss_new;                
            end
            flag_domain_step = 1;
            it_filtering     = it_filtering + 1;
        end

        if it_filtering>N_it_filtering
            src_info = srcout;
            L = Lout;
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
        
        rsc = 2*pi/L;
        Der = specdiffmat(n_bd,src_info)*rsc;
        %this is necessary - need to calculate the curvature
        H = ( dxs .* (Der*dys')' - dys .* (Der*dxs')' )  ./ ( ds.^3 );
    
        t_lambda = t_bd*2*pi/L;
        if(iinv_type == 2)
            lambda_imp = (var_imp(1)+...
            cos(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(2:N_imp_tot+1)'+...
            sin(bsxfun(@times,t_lambda',1:N_imp_tot))*var_imp(N_imp_tot+2:end)')';
        else
            lambda_imp = lambda_imp_guru(t_lambda,iimp_func);
        end

        %generating operators
        norder = 16;
        
        inv_Fw_mats = complex(zeros(n_bd,n_bd,nfreqs));
        sol_to_receptor_mats = complex(zeros(n_tgt,n_bd,nfreqs));
        ubd_aux_all = complex(zeros(n_bd,n_dir,nfreqs));
        dubd_aux_all = complex(zeros(n_bd,n_dir,nfreqs));
        uinc_all = complex(zeros(n_bd,n_dir,nfreqs));
        duinc_all = complex(zeros(n_bd,n_dir,nfreqs));
    
        rhs_all = complex(zeros(n_tgt,n_dir,nfreqs));
        rhs_recast_all = complex(zeros(n_tgt*n_dir*nfreqs,1));
        
        
        for jk=1:nfreqs
            khuse = freqs(jk);
            S  = slp_mat(khuse,norder,h_bd,src_info);
            Sp = sprime_ext_mat(khuse,norder,h_bd,src_info);
            D  = dlp_ext_mat(khuse,norder,h_bd,src_info);
        
            Sik  = slp_mat(1i*abs(khuse),norder,h_bd,src_info);
            Spik = sprime_ext_mat(1i*khuse,norder,h_bd,src_info);
        
            zpars2(1) = khuse;
            zpars2(2) = 1j*abs(khuse);
            Tdiff     = ddiff_neu_mat(zpars2,norder,h_bd,src_info);

    %operators to target
            S_tgt = slmat_out(khuse,h_bd,src,tgt);
            D_tgt = dlmat_out(khuse,h_bd,src,tgt);

     %fw for lambda
            eta = khuse;
            Fw_mat = Sp + 1i* khuse *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4)+ ...
                1i * khuse * bsxfun(@times,lambda_imp',S+1i*eta*D*Sik);
            inv_Fw_mats(:,:,jk) = inv(Fw_mat);
     
     %bd_data
            uinc_all(:,:,jk)  = exp(1i *khuse * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
            duinc_all(:,:,jk) = 1i* khuse * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
                exp(1i *khuse * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
            bd_data = -(duinc_all(:,:,jk) + 1i * khuse * repmat(lambda_imp',1,n_dir) .* uinc_all(:,:,jk));
        
            sol_to_receptor_mats(:,:,jk) = S_tgt + 1i*eta*D_tgt*Sik;

     %scattered field
            pot = inv_Fw_mats(:,:,jk) * bd_data;
            uscat =  sol_to_receptor_mats(:,:,jk)*pot;  
            ubd_aux_all(:,:,jk)  = (S + 1i * eta * D*Sik) *pot;
            dubd_aux_all(:,:,jk) = (Sp + 1i* khuse *( Tdiff*Sik + (Spik + eye(n_bd)/2) * (Spik + eye(n_bd)/2) - eye(n_bd)/4))*pot;
   

     %right hand side
            rhs_all(:,:,jk) = umeas(ifreqs(jk)).data-uscat;
        end
        rhs_recast_all(:) = rhs_all(:);
        
    
        if(nfreqs == 1)
            if(iinv_type == 2)
                bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
                'dys',dys,'L',L,'H',H,'n_bd',n_bd,'rhs',rhs_recast_all,'delta_var',delta_var,...
                'delta_imp',delta_imp,'t_lambda',t_lambda,'lambda_imp',lambda_imp,...
                'sigma_bd',sigma_bd,'it_filtering',it_filtering,'fret_imp_cond',...
                cond(DFw_imp),'fret_dom_cond',cond(DFw_var),'fret_cond',cond(Minv),...
                'rank_imp',rank(DFw_imp,1e-10),'rank_dom',rank(DFw_var,1e-10),...
                'rank_fret',rank(Minv,1e-10));
            else
                bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
                'dys',dys,'L',L,'H',H,'n_bd',n_bd,'rhs',rhs_recast_all,'delta_var',delta_var,...
                't_lambda',t_lambda,'lambda_imp',lambda_imp,...
                'sigma_bd',sigma_bd,'it_filtering',it_filtering,...
                'fret_dom_cond',cond(DFw_var),'rank_dom',rank(DFw_var,1e-10));
            end
        else
            if(iinv_type == 2)
                bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
                'dys',dys,'L',L,'H',H,'n_bd',n_bd,'delta_var',delta_var,...
                'delta_imp',delta_imp,'t_lambda',t_lambda,'lambda_imp',lambda_imp,...
                'sigma_bd',sigma_bd,'it_filtering',it_filtering);
            else
                bd_inv_ref0 = struct('t_bd',t_bd,'h_bd',h_bd,'xs',xs,'ys',ys,'dxs',dxs,...
                'dys',dys,'L',L,'H',H,'n_bd',n_bd,'delta_var',delta_var,...
                't_lambda',t_lambda,'lambda_imp',lambda_imp,...
                'sigma_bd',sigma_bd,'it_filtering',it_filtering);
            end
        end
        bd_inv_refs{it_newton,ik} = bd_inv_ref0; 
 
        
        umeas_test = zeros(n_tgt,n_dir,nfreqs);
        for jk=1:nfreqs
            umeas_test(:,:,jk) = umeas(ifreqs(jk)).data;
        end
        
        
        iesc = 0;
        %stopping parameter
        if ~flag_domain_step

            if(iinv_type == 2) 
                fprintf('Step impedance= %d\n',norm(var_imp-var_imp_old)/norm(var_imp))
            end
	        fprintf('Step domain=%d\n',norm(delta_var))
        
            %stopping parameter
            if(iinv_type == 2)
                if (norm(var_imp-var_imp_old)/norm(var_imp) < eps_step*rscale_delta) %&& (norm(delta_var) < eps_step)
                    flag_newton = 0;
                    iesc = 1;
                    fprintf('Step too small!\n')  
                end
            else
                if (norm(delta_var) < eps_step*rscale_delta)
                    flag_newton = 0;
                    iesc = 1;
                    disp('esp_step='+eps_step);
                    disp('norm step='+norm(delta_var));
                    fprintf('Step too small!\n');
                end
            end
        
            if it_newton > max_it
              flag_newton = 0;
              iesc = 2;
              fprintf('Reached max iteration!\n')            
            end
        
            if norm(rhs_recast_all)/norm(umeas_test(:)) < eps_res
               flag_newton = 0;
               iesc = 3;
               fprintf('RHS too small!\n')
            end
        
            if norm(rhs_old)<norm(rhs_recast_all) 
                if(iinv_type == 2)
                    var_imp  = var_imp_old;
                end
                src_info = src_info_old;
		        L  = L_info;
                iesc = 4;
                fprintf('RHS increasing! %d -> %d\n',norm(rhs_old)/norm(umeas_test(:)),norm(rhs_recast_all)/norm(umeas_test(:)))
                break;
            end
            rhs_old = rhs_recast_all;
        
            fprintf('RHS =%d\n',norm(rhs_recast_all)/norm(umeas_test(:)))
        
            it_newton = it_newton + 1;
        else
            if(iinv_type == 2)
                var_imp  = var_imp_old;
            end
            src_info = src_info_old;
            L        = L_info;
            iesc = 5;
            break;
        end
    end

    fprintf('After newton RHS =%d\n',norm(rhs_recast_all)/norm(umeas(ik).data(:))) 
    if(iinv_type == 2)
        lambda_vecs(ik).coefs = var_imp;
    end
    bd_sols(ik).xs  = src_info(1,:);
    bd_sols(ik).ys  = src_info(2,:);
    bd_sols(ik).rnx = src_info(3,:);
    bd_sols(ik).rny = src_info(4,:);
    bd_sols(ik).ds  = src_info(5,:);
    bd_sols(ik).H   = src_info(6,:);
    bd_sols(ik).h_bd = h_bd;
    bd_sols(ik).L = L;
    iesc_flag(ik) = iesc;
    it_newtons(ik) = it_newton;
    rhs_mags(ik) = norm(rhs_recast_all)/norm(umeas_test(:));
    if(nfreqs>1) 
        rtmp = rhs_all(:,:,ik);
        rtmp = rtmp(:);
        umeas_tmp = umeas_test(:,:,ik);
        umeas_tmp = umeas_tmp(:);
        rhs_mags2(ik) = norm(rtmp)/norm(umeas_tmp);
    else
        rhs_mags2(ik) = norm(rhs_recast_all)/norm(umeas_test(:));
    end
    
    if mod(ik,15) % ~mod(kh,10)
      disp("ik="+ik)
      if(iinv_type == 2)
        save(res_file,...
        'bd_sols','lambda_vecs','khv','iesc_flag','it_newtons','rhs_mags',...
         'rhs_mags2','bd_inv_refs')
      else
          save(res_file,...
        'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
         'rhs_mags2','bd_inv_refs')
      end
    end

end
disp("ik="+ik)

if(iinv_type == 2)
    save(res_file,...
      'bd_sols','lambda_vecs','khv','iesc_flag','it_newtons','rhs_mags',...
      'rhs_mags2','bd_inv_refs')
else
    save(res_file,...
        'bd_sols','khv','iesc_flag','it_newtons','rhs_mags',...
         'rhs_mags2','bd_inv_refs')
end
if(ifexit)
   exit;
end
