x2 = load('test-data/plane_2k.mat');
x3 = load('test-data/plane_3k.mat');
x4 = load('test-data/plane_4k.mat');

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


max_it = 100;
x2step = zeros(max_it,length(x2.khv));
x3step = zeros(max_it,length(x2.khv));
x4step = zeros(max_it,length(x2.khv));

for i=1:length(x2.khv)
    x2step(1:length(x2.step_flag(i).iteration),i) = x2.step_flag(i).iteration;
    x3step(1:length(x3.step_flag(i).iteration),i) = x3.step_flag(i).iteration;
    x4step(1:length(x4.step_flag(i).iteration),i) = x4.step_flag(i).iteration;
end

figure(4)
clf()
subplot(3,1,1)
spy(x2step')

subplot(3,1,2)
spy(x3step')

subplot(3,1,3)
spy(x4step')


figure(5)
clf()

ikp = 17;
subplot(2,1,1)
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.')
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.')
axis equal


subplot(2,1,2)
ikp = 37;
plot(x2.bd_sols(ikp).xs,x2.bd_sols(ikp).ys,'k.'), hold on;
plot(x3.bd_sols(ikp).xs,x3.bd_sols(ikp).ys,'r.')
plot(x4.bd_sols(ikp).xs,x4.bd_sols(ikp).ys,'b.')
axis equal
