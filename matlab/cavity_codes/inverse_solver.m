function [inverse_sol_data,src_info_out] = inverse_solver(kh,src_info,bc, ...
     u_meas,optim_opts,opts)
    inverse_sol_data = [];
    
    if(nargin<5)
      opts = [];
    end

    if(nargin<4)
      optim_opts = [];
    end
    
    verbose = false;
    if(isfield(opts,'verbose'))
        verbose = opts.verbose;
    end
    
    opts_use = opts;
    if(isfield(opts_use,'ncoeff_boundary_mult'))
        opts_use.ncoeff_boundary = floor(kh*ncoeff_boundary_mult);
    end
    
    if(isfield(opts_use,'ncoeff_impedance_mult'))
        opts_use.ncoeff_impedance = floor(kh*ncoeff_impedance_mult);
    end
    
    maxit = 100;
    if(isfield(optim_opts,'maxit'))
        maxit = optim_opts.maxit;
    end
    
   
    optim_type = 'gn';
    if(isfield(optim_opts,'optim_type'))
        optim_type = optim_opts.optim_type;
    end
    
    sd_iter = 0;
    
    if(strcmpi(optim_type,'sd-gn') || strcmpi(optim_type,'sd-min(sd,gn)') || strcmpi(optim_type,'sd-min(gn,sd)'))
        sd_iter = 50;
        if(isfield(optim_opts,'sd_iter'))
           sd_iter = optim_opts.sd_iter;
        end
    end
    
    eps_res = 1e-5;
    if(isfield(optim_opts,'eps_res'))
        eps_res = optim_opts.eps_res;
    end
    
    eps_upd = 1e-5;
    if(isfield(optim_opts,'eps_upd'))
        eps_upd = optim_opts.eps_upd;
    end
    
    src_info_all = cell(maxit,1);
    res_all = zeros(maxit,1);
    exit_criterion = 0;
    fields_all = cell(maxit,1);
    deltas = cell(maxit,1);
    ier = zeros(maxit,1);
    store_fields = false;
    if(isfield(opts,'store_fields'))
        store_fields = opts.store_fields;
    end
    
    store_src_info = false;
    if(isfield(opts,'store_src_info'))
        store_src_info = opts.store_src_info;
    end
    
%
%  Initialize matrices and fields for iteration 0
%
    src_use = src_info;
    mats = get_fw_mats(kh,src_use,bc,u_meas,opts);
    fields = compute_fields(kh,src_use,mats,u_meas,bc,opts);
    size(fields);
    
    if(verbose)
      fprintf('-------------------------------------\n')
      fprintf('Starting inverse interation for frequency kh: %d\n',kh);
    end
    
    optim_opts_use0 = optim_opts;
    optim_opts_use0.optim_type = optim_type;
    optim_opts_use0.eps_res = eps_res;
    optim_opts_use0.eps_upd = eps_upd;
    optim_opts_use0.sd_iter = sd_iter;
    optim_opts_use0.maxit = maxit;
    ict = 1;
    for iter=1:maxit
        optim_opts_use = optim_opts_use0;
        if(strcmpi(optim_type,'sd-gn') || strcmpi(optim_type,'sd-min(sd,gn)') || strcmpi(optim_type,'sd-min(gn,sd)'))
            if(iter<=sd_iter)
                optim_opts_use.optim_type = 'sd';
            else
                optim_opts_use.optim_type=optim_type(4:end);
            end
        end
        
        [deltas{iter},src_info_all{iter},mats_out,fields_all{iter},res_all(iter),ier(iter)] = ...
           update_inverse_iterate(kh,src_use,mats,fields,u_meas,bc,optim_opts,opts);
        if(verbose)
            fprintf('iter number: %d \t optim_type: %s \t residue: %d \t ier: %d\n',iter,optim_opts_use.optim_type,res_all(iter),ier(iter));
        end
        src_use = src_info_all{iter};
        mats = mats_out;
        fields = fields_all{iter};
        
        if(ier(iter) ~=0) 
            exit_criterion = -1;
            break;
        end
        
        if(res_all(iter) <= eps_res)
            exit_criterion = 1;
            break;
        end
        
        res_upd_norm = 0;
        if(isfield(deltas{iter},'delta_bdry'))
            res_upd_norm = res_upd_norm + norm(deltas{iter}.delta_bdry(:));
        end
        
        if(isfield(deltas{iter},'delta_impedance'))
            res_upd_norm = res_upd_norm + norm(deltas{iter}.delta_impedance(:));
        end
        
        if(res_upd_norm <= eps_upd)
            exit_criterion = 2;
            break;
        end
        ict = ict + 1;
    end
    if(ict == (maxit +1)) 
        exit_criterion = 3;
    end
    if(verbose)
        fprintf('Exit criterion: %d \t iteration count: %d\n',exit_criterion,iter);
    end
    src_info_out = src_use;
    inverse_sol_data.res_all = res_all(1:iter);
    inverse_sol_data.iter_count = iter;
    inverse_sol_data.res_opt = res_all(iter);
    if(store_src_info)
        inverse_sol_data.src_info_all = src_info_all(1:iter);
    end
    inverse_sol_data.src_info_opt = src_use;
    if(store_fields)
        inverse_sol_data.fields_all = fields_all(1:iter); 
    end
    inverse_sol_data.deltas = deltas(1:iter);
    inverse_sol_data.fields_opt = fields_all{iter};
    inverse_sol_data.exit_criterion = exit_criterion;
    inverse_sol_data.optim_opts = optim_opts;
    if(verbose)
      fprintf('Completing inverse interation for frequency kh: %d\n',kh);
      fprintf('-------------------------------------\n');
    end
    
end