% Testing script for working with more complicated sensor + direction
% info


% set target locations
%receptors
r_tgt = 10;
n_tgt = 3;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions
n_dir = 3;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

[t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
t_tgt_grid = t_tgt_grid(:);
t_dir_grid = t_dir_grid(:);
xtgt = r_tgt*cos(t_tgt_grid);
ytgt = r_tgt*sin(t_tgt_grid);
tgt   = [xtgt'; ytgt'];

sensor_info = [];
sensor_info.tgt = tgt;
sensor_info.t_dir = t_dir_grid;