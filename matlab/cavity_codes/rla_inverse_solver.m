function [inv_data_all,src_info_out] = rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,rla_path_opts)
   if(nargin<6)
       rla_path_opts = [];
   end
   if(nargin<5)
       src_init = [];
   end
   if(nargin<4)
       opts = [];
   end
   if(nargin<3)
       optim_opts = [];
   end
   
   if(isfield(rla_path_opts,'ik_list'))
       ik_list = rla_path_opts.ik_list;
   else
       ik_list = 1:1:length(u_meas);
   end
   npath = length(ik_list);
   
   nppw = 10;
   if(isfield(opts,'nppw'))
       nppw = opts.nppw;
   end
   
   if(isempty(src_init))
       kh0 = u_meas{ik_list(1)}.kh;
       n = max(300,ceil(abs(kh0)*nppw));
       src_init = geometries.ellipse(1,1,n);
   end
   inv_data_all = cell(1,npath);
   src_info = src_init;
   for i=1:npath
       kh = u_meas{ik_list(i)}.kh;
       [inv_data_all{i},src_out] = inverse_solver(kh,src_info,bc, ...
          u_meas{ik_list(i)},optim_opts,opts);
       src_info = src_out;
   end
   src_info_out = src_info;
end