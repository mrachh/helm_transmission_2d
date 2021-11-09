function [mats,varargout] = get_fw_mats(kh,src_info,bc,sensor_info,opts)
%
% get_fw_mats returns the relevant forward matrices for solving
% Dirichlet, Neumann/ Impedance boundary value problem
% along with an error estimate in the solution if requested.
%
% Input:
%   kh - Helmholtz wave number
%   src_info - source info struct;
%      src_info.xs = x coordinates;
%      src_info.ys = y coordinates;
%      src_info.dxs = dxdt;
%      src_info.dys = dydt;
%      src_info.ds = sqrt(dxdt^2 + dydt^2);
%      src_info.h = h in trapezoidal parametrization;
%      src_info.lambda - imepdance value at discretization nodes 
%           (optional if solving impedance boundary value problem);
%   bc - boundary condition struct;
%     bc.type = type of boundary condition;
%        'd' or 'Dirichlet' for dirichlet
%        'n' or 'Neumann' for Neumann
%        'i' or 'Impedance' for impedance
%   sensor_info - sensor_info struct
%       sensor_info.xtgt = x coordinates of targets;
%       sensor_info.ytgt = y coordinates of targets;
%
%   opts - options struct
%      opts.ifflam - whether to use FLAM to compress the matrices (false)
%      opts.verbose - flag for displaying detailed output (false)
%      opts.test_analytic - flag for whether to test for analytic solution
%         (false)
%      opts.src_in - interior point if testing analytic solution (0,0)
% Output
%
   if nargin < 5
       opts = [];
   end
   mats = [];
   ifflam = false;
   if isfield(opts,'ifflam')
       ifflam = opts.ifflam;
   end
   
   test_analytic = false;
   verbose = false;
   if(isfield(opts,'verbose'))
       verbose = opts.verbose;
   end
   if isfield(opts,'test_analytic')
       test_analytic = opts.test_analytic;
       if isfield(opts,'src_in')
          src_in = opts.src_in;
       else
          if(verbose)
            fprintf('Interior source points not set, assuming to be origin\n\n');
          end
          src_in = zeros(2,1);
       end
   end
   

   
   
   norder = 16;
   n_bd = length(src_info.xs);
   srctmp = zeros(6,n_bd);
   srctmp(1,:) = src_info.xs;
   srctmp(2,:) = src_info.ys;
   srctmp(3,:) = src_info.dys./src_info.ds;
   srctmp(4,:) = -src_info.dxs./src_info.ds;
   srctmp(5,:) = src_info.ds;
   if(isfield(src_info,'H'))
      srctmp(6,:) = src_info.H;
   else
      srctmp(6,:) = 0;
   end
   h_bd = src_info.h;
   
   
   src = zeros(4,n_bd);
   src(1,:) = src_info.xs;
   src(2,:) = src_info.ys;
   src(3,:) = src_info.dxs;
   src(4,:) = src_info.dys;
   
   tgt = [sensor_info.xtgt; sensor_info.ytgt];
   varargout{1} = 0;
   
   if ~ifflam
      bctype = bc.type;
      
      
      % matrices for dirichlet boundary value problem
      % Representation for solution:
      % u = D_{k} + ik S_{k}
      %
      % Integral equation:
      % u = (1/2 + D_{k} + ik S_{k})[\sigma]
      %
      % To get dirichlet data: 
      % u = Fw_dir_mat*sigma
      %
      % To get neumann data:
      % dudn = Fw_neu_mat*(dudn_inc - ik u_inc)
      if(strcmpi(bctype,'d') || strcmpi(bctype,'Dirichlet'))
          S  = slp_mat(kh,norder,h_bd,srctmp);
          Sp = sprime_ext_mat(kh,norder,h_bd,srctmp) + eye(n_bd)/2;
          D  = dlp_ext_mat(kh,norder,h_bd,srctmp) - eye(n_bd)/2;
          mats.Fw_mat = eye(n_bd)/2 + D + 1i* kh * S;
          mats.inv_Fw_mat = inv(mats.Fw_mat);
          mats.Fw_dir_mat = mats.Fw_mat;
          mats.Fw_neu_mat = inv(eye(n_bd)/2 + Sp - 1i * kh *S);
          
          
          
          %operators to target
          S_tgt = slmat_out(kh,h_bd,src,tgt);
          D_tgt = dlmat_out(kh,h_bd,src,tgt); 
          mats.sol_to_receptor = (D_tgt + 1i * kh * S_tgt);
          
          if (test_analytic)
              %data for checking  
            uin_a = helm_c_p(kh,src_in,srctmp);
            rhs_a = uin_a;
            sol_a = mats.inv_Fw_mat * rhs_a;
    
            utest = mats.sol_to_receptor * sol_a;  
            uex = helm_c_p(kh,src_in,tgt); 
            varargout{1} = norm(utest(:)-uex(:))/norm(uex(:));
          end
      end
      
      % Matrices for Neumann boudnary value problem
      % Representation for solution
      % u = (S_{k} + ik*D_{k} S_{ik})\sigma
      %
      % Integral equation:
      % dudn = -(2 + ik)/4 + S' + ik (D'_{k} - D'_{ik}) S_{ik} + S_{ik}'^2)
      %
      % To get dirichlet data: 
      % u = Fw_dir_mat*sigma
      %
      % To get neumann data:
      % dudn = Fw_neu_mat*sigma
      if(strcmpi(bctype,'n') || strcmpi(bctype,'Neumann'))
          S  = slp_mat(kh,norder,h_bd,srctmp);
          Sp = sprime_ext_mat(kh,norder,h_bd,srctmp) + eye(n_bd)/2;
          D  = dlp_ext_mat(kh,norder,h_bd,srctmp) - eye(n_bd)/2;
          Sik = slp_mat(1i*kh,norder,h_bd,srctmp);
          Spik = sprime_ext_mat(1i*kh,norder,h_bd,srctmp) + eye(n_bd)/2;
          zpars2 = complex(zeros(2,1));
          zpars2(1) = kh;
          zpars2(2) = 1i*abs(kh);
          ddiff = ddiff_neu_mat(zpars2,norder,h_bd,srctmp);
          mats.Fw_mat = -eye(n_bd)/2 + Sp + 1i* kh *( ddiff*Sik + ...
                 Spik * Spik - ...
                 eye(n_bd)/4);
          mats.inv_Fw_mat = inv(mats.Fw_mat);
          mats.Fw_dir_mat = S+1i*kh*(D+eye(n_bd)/2)*Sik;
          mats.Fw_neu_mat = mats.Fw_mat;
          
          %operators to target
          S_tgt = slmat_out(kh,h_bd,src,tgt);
          D_tgt = dlmat_out(kh,h_bd,src,tgt); 
          mats.sol_to_receptor = (S_tgt + 1i * kh * D_tgt*Sik);
          if (test_analytic)
            %data for checking    
            
            dudnin_a = helm_c_gn(kh,src_in,srctmp);
            rhs_a = dudnin_a;
            sol_a = mats.inv_Fw_mat * rhs_a;
    
            utest = mats.sol_to_receptor * sol_a;  
            uex = helm_c_p(kh,src_in,tgt); 
            varargout{1} = norm(utest(:)-uex(:))/norm(uex(:));
          end 
      end
      
      
      % Matrices for Impedance boudnary value problem
      % Representation for solution
      % u = S_{k} + ik*D_{k} S_{ik}
      %
      % Integral equation:
      % dudn + ik\lambda u = -(2 + ik)/4 + S' + ik (D'_{k} - D'_{ik}) S_{ik} + S_{ik}'^2)
      %                     +ik \lambda(ik S_{ik}/2 + S_{k} + ik D_{k}
      %                     S_{ik})
      %
      % To get dirichlet data: 
      % u = Fw_dir_mat*sigma
      %
      % To get neumann data:
      % dudn = Fw_neu_mat*sigma
      if(strcmpi(bctype,'i') || strcmpi(bctype,'Impedance'))
          S  = slp_mat(kh,norder,h_bd,srctmp);
          Sp = sprime_ext_mat(kh,norder,h_bd,srctmp) + eye(n_bd)/2;
          D  = dlp_ext_mat(kh,norder,h_bd,srctmp) - eye(n_bd)/2;
          Sik = slp_mat(1i*kh,norder,h_bd,srctmp);
          Spik = sprime_ext_mat(1i*kh,norder,h_bd,srctmp) + eye(n_bd)/2;
          zpars2 = complex(zeros(2,1));
          zpars2(1) = kh;
          zpars2(2) = 1i*abs(kh);
          ddiff = ddiff_neu_mat(zpars2,norder,h_bd,srctmp);
          mats.Fw_mat = -eye(n_bd)/2 + Sp + 1i* kh *( ddiff*Sik + ...
                 Spik * Spik - ...
                 eye(n_bd)/4) + ...
                 1i*kh*bsxfun(@times,src_info.lambda,S+1i*kh*(D+eye(n_bd)/2)*Sik);
          mats.inv_Fw_mat = inv(mats.Fw_mat);
          mats.Fw_dir_mat = S+1i*kh*(D+eye(n_bd)/2)*Sik;
          mats.Fw_neu_mat = -eye(n_bd)/2 + Sp + 1i* kh *( ddiff*Sik + ...
                 Spik * Spik - ...
                 eye(n_bd)/4);
          
          %operators to target
          S_tgt = slmat_out(kh,h_bd,src,tgt);
          D_tgt = dlmat_out(kh,h_bd,src,tgt); 
          mats.sol_to_receptor = (S_tgt + 1i * kh * D_tgt*Sik);
          if (test_analytic)
            %data for checking    
            uin_a = helm_c_p(kh,src_in,srctmp);
            dudnin_a = helm_c_gn(kh,src_in,srctmp);
            rhs_a = dudnin_a+1i*kh*src_info.lambda.*uin_a;
            sol_a = mats.inv_Fw_mat * rhs_a;
    
            utest = mats.sol_to_receptor * sol_a;  
            uex = helm_c_p(kh,src_in,tgt); 
            varargout{1} = norm(utest(:)-uex(:))/norm(uex(:));
          end 
      end
   end
end