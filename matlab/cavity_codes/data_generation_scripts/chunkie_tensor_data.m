% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction


n  = 300;
addpath('../');
nc = 9;
coefs = zeros(2*nc+1,1);
coefs(1) = 1;
coefs(nc+1) = 0.3;
src_info.xs = [1,2,2,1,-1,-2,-2,-1,-1,-1.5,-1.5,1.5,1.5,1];
src_info.ys = [-1,0,1,2,2,1,0,-1,0,0,1,1,0,0];
plot(src_info.xs,src_info.ys);



nk = 9;


% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';

fname = ['../data/knight9_ik1_nk' int2str(nk) '_tensor_data_' bc.type '.mat'];
save(fname,'src_info');


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions
n_dir = 100;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

[t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
t_tgt_grid = t_tgt_grid(:);
t_dir_grid = t_dir_grid(:);
xtgt = r_tgt*cos(t_tgt_grid);
ytgt = r_tgt*sin(t_tgt_grid);
tgt   = [ xtgt'; ytgt'];


sensor_info = [];
sensor_info.tgt = tgt;
sensor_info.t_dir = t_dir_grid;



%src_info = geometries.starn(coefs,nc,n);


dk = 0.25;
kh = 1:dk:(1+(nk-1)*dk);

u_meas = cell(nk,1);



nppw = 20;

for ik=1:nk
    ik
   %n = ceil(nppw*L*abs(kh(ik))/2/pi);
   %n = max(n,300);
   %src_info = geometries.starn(coefs,nc,n);
   
   %[mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
   %fields = rla.compute_fields_chunkie(kh(ik),src_info,sensor_info,bc,opts);
   [fields] = helm_dirichlet_solver(kh(ik),src_info,sensor_info);
   u_meas0 = [];
   u_meas0.kh = kh(ik);
   u_meas0.uscat_tgt = fields.uscat_tgt;
   u_meas0.tgt = sensor_info.tgt;
   u_meas0.t_dir = sensor_info.t_dir;
   %u_meas0.err_est = erra;
   u_meas{ik} = u_meas0;
end


save(fname,'u_meas','-append');
