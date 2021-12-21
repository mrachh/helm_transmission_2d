function fields = compute_fields(kh,src_info,mats,sensor_info,bc,opts)
   fields = [];
   
   if(nargin < 6)
       opts = [];
   end
   
   ifflam = false;
   if(isfield(opts,'ifflam'))
       ifflam = opts.ifflam;
   end
   
   [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
   [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
   nt_uni = length(tgt_uni(:,1));
   induse = itgt + (idir-1)*nt_uni;
   
   
   t_dir_uni = t_dir_uni(:);
   x_dir = cos(t_dir_uni)';
   y_dir = sin(t_dir_uni)';
   n_dir = length(x_dir);
   xs = src_info.xs(:)';
   ys = src_info.ys(:)';
   ds = src_info.ds(:)';
   dxs = src_info.dxs(:)';
   dys = src_info.dys(:)';
   
   fields.uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   fields.dudninc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   
   if(strcmpi(bc.type,'d') || strcmpi(bc.type,'Dirichlet'))
      bd_data = -fields.uinc;
   end

   if(strcmpi(bc.type,'n') || strcmpi(bc.type,'Neumann'))
      bd_data = -fields.dudninc;
   end

   if(strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
     lambda_rep = repmat(src_info.lambda,1,n_dir);
     bd_data = -(fields.dudninc + 1i*kh*lambda_rep.*fields.uinc);
   end


   if(~ifflam)
     uscat_tgt_all = mats.bdrydata_to_receptor*bd_data;
     uscat_tgt_all = uscat_tgt_all(:);
     fields.uscat_tgt = uscat_tgt_all(induse);
     sigma = mats.inv_Fw_mat*bd_data;
     fields.uscat = mats.Fw_dir_mat*sigma;
     if(strcmpi(bc.type,'d') || strcmpi(bc.type,'Dirichlet'))
         bd_data2 = fields.dudninc - 1i * kh *fields.uinc;
         fields.dudnscat = mats.Fw_neu_mat*bd_data2 - fields.dudninc;
     end
     
     if(strcmpi(bc.type,'n') || strcmpi(bc.type,'Neumann') || strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
        fields.dudnscat = mats.Fw_neu_mat*sigma; 
     end
   end
   
end