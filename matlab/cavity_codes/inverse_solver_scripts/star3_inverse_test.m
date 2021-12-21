load('../data/star3_ik1_nk9_tensor_data_Dirichlet.mat');

optim_opts = [];
opts = [];
opts.verbose=true;
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';
optim_opts.optim_type = 'min(gn,sd)';
optim_opts.eps_curv = 1e-3;
optim_opts.filter_type = 'step_length';
opts.store_src_info = true;
[inv_data_all,src_info_out] = rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts);