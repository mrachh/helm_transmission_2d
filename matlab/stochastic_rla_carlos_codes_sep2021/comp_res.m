data = load('larry_01_pi3_nk097_nsource100_ntarget100.mat');
load('combined_sol_k25_cb2_ns1_stoc1-45.mat');

igood = [2;21;28;35;38;41;44];
ibad = [6;9;12;18;19;24;25;30;31];


res_good = zeros(97,1000,length(igood));
ikh_good = zeros(length(igood),1);
jkh_good = zeros(1000,length(igood));

res_bad = zeros(97,1000,length(ibad));
ikh_bad = zeros(length(ibad));
jkh_bad = zeros(1000,length(ibad));


% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];

% Set interior point
src0 = [0.01;-1.1];%horseshoe


% incidence directions
n_dir = 100;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
x_dir = cos(t_dir);
y_dir = sin(t_dir);
dir =[ x_dir; y_dir ];

[work,lw] = hank106_workarray();

for igeom=1:1
    i = igood(igeom);
    nkh = length(S{i}.kh);
    ikh_good(igeom) = nkh;
    tic,
    for ikh = nkh:nkh
        % Define geometry here
        n = length(S{i}.bd(ikh).xs);
        src_info = zeros(6,n);
        src_info(1,:) = S{i}.bd(ikh).xs;
        src_info(2,:) = S{i}.bd(ikh).ys;
        src_info(5,:) = S{i}.bd(ikh).ds;
        src_info(3,:) = S{i}.bd(ikh).rnx;
        src_info(4,:) = S{i}.bd(ikh).rny;
        h_bd  = S{i}.bd(ikh).h_bd;
        L     = S{i}.bd(ikh).L;
        kh_cur = S{i}.kh(ikh);
        jkhmax = find(data.khv == kh_cur);
        jkh_good(ikh,igeom) = jkhmax;
        for jkh=1:jkhmax
            khuse = data.khv(jkh);
            [ucomp,~,~,~,~,~,~,~] = forward_dirichlet_fast(khuse,...
               src_info,L,h_bd,dir,tgt,src0,work,lw);
           uex = data.umeas(jkh).data;
           rhs = ucomp - uex;
           res_good(jkh,ikh,igeom) = norm(rhs(:))/norm(uex);
        end
    end
    t = toc;
    fprintf('\n\n\n\n Total time taken for geometry %d\n\n\n',t);
end