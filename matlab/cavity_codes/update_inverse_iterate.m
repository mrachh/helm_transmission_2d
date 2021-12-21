function [deltas,src_out,mats_out,fields_out,res,ier] = ...
   update_inverse_iterate(kh,src_info,mats,fields,u_meas,bc,optim_opts,opts)
%
%  This subroutine updates the boundary and/or the impedance function
%  of the curve, based on one optimization step of 
%  minimizing the objective function of the scattered field
%  at a collection of target locations, and incident directions
%  for a single frequency kh.
%
%  Let F(k,\bx_{j}, d_{j}, \Gamma, \lambda) denote the scattered
%  field at frequency k, at target location \bx_{j}, due to incident 
%  field given by plane wave with direction $d_{j}$, where
%  the boundary of the curve is given by $\Gamma$ and the impedance
%  function given by $\lambda$ (note that the impedance function)
%  is optional, and let u_meas denote the corresponding measured data.
%
%  We support 3 optimization methods for updating the boundary/impedance
%   sd - Steepest descent with step size set to the Cauchy point
%   gn - Gauss Newton
%   min(gn,sd) or min(sd,gn) - Use the better of steepest descent or Gauss-newton
%
%  The code sequentially updates the boundary of the obstacle first if requested,
%  followed by the impedance if requested.
%
%  Suppose J_b is the Frechet derivative of F with respect to the boundary
%  where the update is assumed to be a scalar function \delta_b in the
%  normal direction, and the scalar function is discretized in a fourier
%  basis.
%  Then if optimization method = sd, then
%    \delta_{b} = t J_{b}^{T} (u_meas - F); 
%      with t = (u_meas - F)^{T} J_{b} J_{b}^{T} (u_meas - F)/|J_{b}
%      J_{b}^{T} (u_meas-F)|^2;
%  
%   If the optimization method = gn, then
%     \delta_{b} = J_{b}^{+} (u_meas -F); 
%
%   where J_{b}^{+} is the pseudo inverse of J_{b}
%
%   If the updated curve is self-intersecting or if the curvature
%   is requested to have a nearly band limited tail then,
%   the above updates are modified by either scaling them by 2^i 
%   for i<=maxit_filtering, so that the conditions are met, or
%   by damping the high frequency components of the update using 
%   a Gaussian until the conditions are met.
%
%   If the conditions aren't met with the maximum filtering iterations, 
%   then the input curve is returned with an error code
%
%   If the optimization scheme is min(sd,gn) or min(gn,sd), then
%   both the gn and the sd updates are computed, and the updated curve
%   which has a smaller residue |u_meas - F(.,.,\Gamma_new,.)|
%   is returned
%
%   For the impedance, a similar procedure is applied after updating the
%   boundary. However, there are no filtering involved in the case
%   of the impedance update.
%-------
%   NOTE: IMPEDANCE currently only runs in gauss-newton mode. This needs
%    to be fixed
%------
%   
%
%
%
%
%
%
%
%
%
%
%
%
%
%

    deltas = [];
    if(nargin<8)
        opts = [];
    end

    if(nargin<7)
        optim_opts = [];
    end
   
    opts_use = opts;
    if(~isfield(opts,'ncoeff_boundary'))
        ncoeff_boundary = ceil(2*abs(kh));
        opts_use.ncoeff_boundary = ceil(2*abs(kh));
    else
        ncoeff_boundary = opts.ncoeff_boundary;
    end
    if(~isfield(opts,'ncoeff_impedance'))
        opts_use.ncoeff_impedance = ceil(0.5*abs(kh));
        ncoeff_impedance = ceil(0.5*abs(kh));
    else
        ncoeff_impedance = opts.ncoeff_impedance;
    end
    
    
    nppw = 10;
    rlam = Inf;
    if(isfield(opts,'nppw'))
        nppw = opts.nppw;
        rlam = 2*pi/real(kh);
    end
    
    eps_curv = Inf;
    if(isfield(optim_opts,'eps_curv'))
        eps_curv = optim_opts.eps_curv;
        
    end
    
    
    optim_type = 'gn';
    if(isfield(optim_opts,'optim_type'))
        optim_type = optim_opts.optim_type;
    end
    
    filter_type = 'gauss-conv';
    if(isfield(optim_opts,'filter_type'))
        filter_type = optim_opts.filter_type;
    end
    
    maxit_filter = 10;
    if(isfield(optim_opts,'maxit_filter'))
        maxit_filter = optim_opts.maxit_filter;
    end
    
    verbose = false;
    if(isfield(opts,'verbose'))
        verbose = opts.verbose;
    end
    
    rhs = u_meas.uscat(:) - fields.uscat_tgt(:);
    frechet_mats = get_frechet_ders(kh,mats,src_info,u_meas,fields,bc,opts_use);
    
    
    opts_update_geom = [];
    opts_update_geom.eps_curv = eps_curv;
    if(isfield(optim_opts,'n_curv'))
       opts_update_geom.n_curv = optim_opts.n_curv;
    else
        opts_update_geom.n_curv = max(30,ncoeff_boundary);
    end
 
    opts_update_geom.nppw = nppw;
    opts_update_geom.rlam = rlam;
    
    if(strcmpi(bc.invtype,'o') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
        deltas.nmodes_bdry = opts_use.ncoeff_boundary;
        Minv = [real(frechet_mats.bdry); imag(frechet_mats.bdry)];
        if(strcmpi(optim_type,'gn') || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_bdry_gn = Minv \ [real(rhs(:)); imag(rhs(:))];
            delta_bdry_gn0 = delta_bdry_gn;
            ier_gn = 10;
            iter_filter_bdry_gn = -1;
            for iter=1:maxit_filter
                [src_out_gn,ier_gn] = update_geom(src_info,ncoeff_boundary,delta_bdry_gn,opts_update_geom);
                if(ier_gn == 0) 
                    iter_filter_bdry_gn = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_boundary);
                        N_var = ncoeff_boundary;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_bdry_gn = (delta_bdry_gn0'.*gauss_new)';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_bdry_gn = delta_bdry_gn0/(2^(iter));
                    end
                end
            end
            
            if(iter_filter_bdry_gn == -1)
                iter_filter_bdry_gn = maxit_filter;
            end
            
            if(verbose)
                fprintf('Post gn filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_bdry_gn,ier_gn);
            end
            if(ier_gn ~=0)
                mats_out_gn = mats;
                fields_out_gn = fields;
            else
                [mats_out_gn] = get_fw_mats(kh,src_out_gn,bc,u_meas,opts);
                fields_out_gn = compute_fields(kh,src_out_gn,mats_out_gn,u_meas,bc,opts);
            end
        end
        if(strcmpi(optim_type,'sd')  || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_bdry_sd = Minv'*[real(rhs(:)); imag(rhs(:))];
            t = delta_bdry_sd'*delta_bdry_sd/norm(Minv*delta_bdry_sd)^2;
            delta_bdry_sd = t*delta_bdry_sd;
            delta_bdry_sd0 = delta_bdry_sd;
            %ier_sd = 10;
            iter_filter_bdry_sd = -1;
            for iter=1:maxit_filter
                ier_sd = 10;
                [src_out_sd,ier_sd] = update_geom(src_info,ncoeff_boundary,delta_bdry_sd,opts_update_geom);
                if(ier_sd == 0) 
                    iter_filter_bdry_sd = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_boundary);
                        N_var = ncoeff_boundary;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_bdry_sd = (delta_bdry_sd0'.*gauss_new)';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_bdry_sd = delta_bdry_sd0/(2^(iter));
                    end
                end
            end
            if(iter_filter_bdry_sd == -1)
                iter_filter_bdry_sd = maxit_filter;
            end
               
            
            if(verbose)
                fprintf('Post sd filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_bdry_sd,ier_sd);
            end
            
            if(ier_sd ~=0)
                mats_out_sd = mats;
                fields_out_sd = fields;
            else
                [mats_out_sd] = get_fw_mats(kh,src_out_sd,bc,u_meas,opts);
                fields_out_sd = compute_fields(kh,src_out_sd,mats_out_sd,u_meas,bc,opts);
            end
        end
        
        
        if(strcmpi(optim_type,'gn'))
            deltas.iter_type = 'gn';
        elseif(strcmpi(optim_type,'sd'))
            deltas.iter_type = 'sd';
        elseif(strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            rhs_gn = u_meas.uscat(:) - fields_out_gn.uscat_tgt(:);
            rhs_sd = u_meas.uscat(:) - fields_out_sd.uscat_tgt(:);
            if(norm(rhs_gn(:)) < norm(rhs_sd(:)))
                deltas.iter_type = 'gn';
            else
                deltas.iter_type = 'sd';
            end
            
            fprintf('optimization method used = %s\n',deltas.iter_type);
        end
        
        if(strcmpi(deltas.iter_type,'gn'))
            rhs_gn = u_meas.uscat(:) - fields_out_gn.uscat_tgt(:);
            deltas.delta_bdry = delta_bdry_gn;
            deltas.iter_filter_bdry = iter_filter_bdry_gn;
            mats_out = mats_out_gn;
            fields_out = fields_out_gn;
            ier = ier_gn;
            src_out = src_out_gn;
            res = norm(rhs_gn(:))/norm(u_meas.uscat(:));
        elseif(strcmpi(deltas.iter_type,'sd'))
            rhs_sd = u_meas.uscat(:) - fields_out_sd.uscat_tgt(:);
            deltas.delta_bdry = delta_bdry_sd;
            deltas.iter_filter_bdry = iter_filter_bdry_sd;
            
            mats_out = mats_out_sd;
            fields_out = fields_out_sd;
            ier = ier_sd;
            src_out = src_out_sd;
            res = norm(rhs_sd(:))/norm(u_meas.uscat(:));
        end      
    end
    if(strcmpi(bc.invtype,'i') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
        deltas.nmodes_imp = ncoeff_impedance;
        Minv = [real(frechet_mats.impedance); imag(frechet_mats.impedance)];
        deltas.delta_impedance = Minv \ [real(rhs(:)); imag(rhs(:))];
        nh = ncoeff_impedance;
        hcoefs_use = deltas.delta_impedance(:);
        n = length(src_out.xs);
        t = 0:2*pi/n:2*pi*(1.0-1.0/n); 
        h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
        src_out.lambda = src_out.lambda + h_upd';  
        
        
%       update matrices and fields after updated impedance        
        mats_out = get_fw_mats(kh,src_out,bc,u_meas,opts);
        fields_out = compute_fields(kh,src_out,mats_out,u_meas,bc,opts);
        rhs = u_meas.uscat(:) - fields_out.uscat_tgt(:);
        res = norm(rhs(:))/norm(u_meas.uscat(:));
    end 
end