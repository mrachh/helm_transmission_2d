x2 = load('test-data/plane_2k_contrast1.2.mat');
x3 = load('test-data/sol_k20_idom1_bctype3_iimp_func1_highres.mat');
x4 = load('test-data/sol_k20_idom1_bctype3_iimp_func1_highres_iiv_type2.mat');


figure(1), 
clf()
semilogy(x2.khv,x2.rhs_mags,'kx'); hold on;
semilogy(x3.khv,x3.rhs_mags,'rx');
semilogy(x4.khv,x4.rhs_mags,'bx');


figure(2), 
clf()
plot(x2.khv,x2.it_newtons,'kx'); hold on;
plot(x3.khv,x3.it_newtons,'rx');
plot(x4.khv,x4.it_newtons,'bx');


figure(3),
clf()
plot(x2.khv,x2.iesc_flag,'kx'); hold on;
plot(x3.khv,x3.iesc_flag,'rx');
plot(x4.khv,x4.iesc_flag,'bx');


figure(4)
clf()

ikp = 17;
subplot(2,2,1)
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.');
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.');
plot(x2.xo,x2.yo,'k--')
axis equal


subplot(2,2,2)
ikp = 37;
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.')
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.');
plot(x2.xo,x2.yo,'k--')

axis equal


subplot(2,2,3)
ikp = 57;
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.')
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.');
plot(x2.xo,x2.yo,'k--')

axis equal


subplot(2,2,4)
ikp = 77;
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.')
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.');
plot(x2.xo,x2.yo,'k--')

axis equal