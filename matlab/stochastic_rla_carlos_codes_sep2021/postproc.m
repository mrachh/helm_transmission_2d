
load('combined_sol_k25_cb2_ns1_stoc1-45.mat')


a = 0.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/3;  % controls approx opening angle in radians (keep small for resonant)

kh = 25;

% Estimate length of curve
n = 300;
nhalf = ceil(n/2);
s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
z = z*1.2;  % vert stretch! makes ellipse cavity
Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve



figure(1)
clf
plot(real(Z),imag(Z),'k.'), hold on;
figure(2)
clf
hold on;
for ii=1:length(S)
    i = ii;
    fprintf('icur=%d\n',i);
    kk =length(S{i}.kh);
    figure(1)
    plot(S{i}.bd(kk).xs,S{i}.bd(kk).ys);
    figure(2)
    plot(S{i}.kh)
    %pause
end



igood = [2;21;28;35;38;41;44];

figure(3)
clf
plot(real(Z),imag(Z),'k.'), hold on;
figure(4)
clf
hold on;
for ii=1:length(igood)
    i = igood(ii);
    fprintf('icur=%d\n',i);
    kk =length(S{i}.kh);
    figure(3)
    plot(S{i}.bd(kk).xs,S{i}.bd(kk).ys);
    figure(4)
    plot(S{i}.kh)
    %pause
end


ibad = [6;9;12;18;19;24;25;30;31];

figure(5)
clf
plot(real(Z),imag(Z),'k.'), hold on;
figure(6)
clf
hold on;
for ii=1:length(ibad)
    i = ibad(ii);
    fprintf('icur=%d\n',i);
    kk =length(S{i}.kh);
    figure(5)
    plot(S{i}.bd(kk).xs,S{i}.bd(kk).ys);
    figure(6)
    plot(S{i}.kh)
    %pause
end


