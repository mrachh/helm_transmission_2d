S1 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc1-1.mat');
S2 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc2-8.mat');
S3 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc9-15.mat');
S4 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc16-22.mat');
S5 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc23-23.mat');
S6 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc24-24.mat');
S7 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc25-25.mat');
S8 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc26-30.mat');
S9 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc31-35.mat');
S10 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc41-45.mat');
S11 = load('multi_larry_01_pi3_k25_cb2_ns1_nostoc46-50.mat');

S = cell(45,1);
S{1} = S1.stoc;
istart = 1
for i =1:length(S2.stoc)
    S{istart+i} = S2.stoc(i);
end
istart = istart + length(S2.stoc);
for i =1:length(S3.stoc)
    S{istart+i} = S3.stoc(i);
end
istart = istart + length(S3.stoc);
for i =1:length(S4.stoc)
    S{istart+i} = S4.stoc(i);
end
istart = istart + length(S4.stoc);
for i =1:length(S5.stoc)
    S{istart+i} = S5.stoc(i);
end
istart = istart + length(S5.stoc);

for i =1:length(S6.stoc)
    S{istart+i} = S6.stoc(i);
end
istart = istart + length(S6.stoc);

for i =1:length(S7.stoc)
    S{istart+i} = S7.stoc(i);
end
istart = istart + length(S7.stoc);


for i =1:length(S8.stoc)
    S{istart+i} = S8.stoc(i);
end
istart = istart + length(S8.stoc);

for i =1:length(S9.stoc)
    S{istart+i} = S9.stoc(i);
end
istart = istart + length(S9.stoc);

for i =1:length(S10.stoc)
    S{istart+i} = S10.stoc(i);
end
istart = istart + length(S10.stoc);

for i =1:length(S11.stoc)
    S{istart+i} = S11.stoc(i);
end
istart = istart + length(S11.stoc);

xs0 = S1.xs_orig;
ys0 = S1.ys_orig;
save('combined_sol_k25_cb2_ns1_stoc1-45.mat','S','xs0','ys0');