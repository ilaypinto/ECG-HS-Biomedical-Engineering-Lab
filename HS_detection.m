% HS detection %
%Based on code of another student inn an earlier year.
%Thank you!
% import data
clear all; clc;
data = fopen('Recording3.txt','rt');
C = textscan(data, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',1);
fclose(data); s_r = 500; dt = 1/500;
time_vec = [0:length(C{1})-1]/s_r; fc = 8; fc2 = 50; fs = length(C{1})/30;
b1 = fir1(48,[fc/(fs/2) fc2/(fs/2)],'bandpass'); a=1; %Define filter for ECG
fc = 50; fc2 = 220; fs = length(C{1})/30;
b2 = fir1(48,[fc/(fs/2) fc2/(fs/2)],'bandpass'); a=1; %Define filter for HS
%%
% graphs for intuition
figure(1);
subplot(4,1,1);plot(time_vec,C{1});
subplot(4,1,2);plot(time_vec,filter(b2,a,C{1}));
subplot(4,1,3);plot(time_vec,C{2});
subplot(4,1,4);plot(time_vec,filter(b1,a,C{2}));

%%
% filtering and slicing the data to rest and active periods
hs_filt = filter(b2,a,C{1}); ecg_filt = filter(b1,a,C{2});
hs_rest = hs_filt(1:20*s_r);
hs_active = hs_filt(20*s_r+1:end);
ECG_rest = ecg_filt(1:20*s_r);
ECG_active = ecg_filt(20*s_r+1:end);
%%
% HS processing
% normalization relatively to biggest value
hs_rest_norm = hs_rest./max(hs_rest);
hs_active_norm = hs_active./max(hs_active);
%%
% Signal Shannon Entropy
for i = 1:length(hs_rest_norm)
    Ehs_rest(i) = -((hs_rest_norm(i))^2)*(log((hs_rest_norm(i))^2));
end

for i = 1:length(hs_active_norm)
    Ehs_active(i) = -((hs_active_norm(i))^2)*(log(hs_active_norm(i))^2);
end
%%
% ECG QRS detection - using AF2 algorithm from 'QRS_detection'

% rest
thresh1 = 0.4 * max(ECG_rest);
Y0_1 = abs(ECG_rest); Y1_1 = zeros(length(ECG_rest),1);
for i = 1:length(ECG_rest)
    if ECG_rest(i) >= thresh1  
        Y1_1(i) = Y0_1(i);
    else
        Y1_1(i) = thresh1;
    end
end

Y2_1 = diff(Y1_1);
for j = 1:length(Y2_1)
    if Y2_1(j) < 0.3 * max(Y2_1)  % QRS candidates are set
        Y2_1(j) = 0;
    end
end
Y2_1 = islocalmax(Y2_1); m_rest=[];
for j = 1:length(Y2_1)
    if Y2_1(j) == 1
       m_rest(end+1) = j;
    end
end

% active
thresh1 = 0.4*max(ECG_active);
Y0_1 = abs(ECG_active); Y1_1 = zeros(length(ECG_active),1);
for i = 1:length(ECG_active)
    if ECG_active(i) >= thresh1  
        Y1_1(i) = Y0_1(i);
    else
        Y1_1(i) = thresh1;
    end
end
Y2_1 = diff(Y1_1);
for j = 1:length(Y2_1)
    if Y2_1(j) < 0.3*max(Y2_1)  % QRS candidates are set
        Y2_1(j) = 0;
    end
end
Y2_1 = islocalmax(Y2_1); m_active=[];
for j = 1:length(Y2_1)
    if Y2_1(j) == 1
       m_active(end+1) = j;
    end
end
%%
% S1,S2 detection
peaks_rest = zeros(size(Ehs_rest));
peaks_active = zeros(size(Ehs_active));

%Threshold that isolates high energy parts from low energy parts
for i = 1:length(Ehs_rest)
    if Ehs_rest(i) >= 0.2 * max(Ehs_rest)
        peaks_rest(i) = 1;
    else
        peaks_rest(i) = 0;
    end
end

for i = 1:length(Ehs_active)
    if Ehs_active(i) >= 0.2*max(Ehs_active)
        peaks_active(i) = 1;
    else
        peaks_active(i) = 0;
    end
end

%%
%Looking for high energy part in the physiological time window, appropriate
%to each QRS detected, respectively.

% rest
s1_idx_rest = zeros(1, length(hs_rest));
s1_val_rest = zeros(1, length(hs_rest));
s2_idx_rest = zeros(1, length(hs_rest));
s2_val_rest = zeros(1, length(hs_rest));
R_peak_of_S1_rest = zeros(1, length(m_rest));
R_peak_of_S2_rest = zeros(1, length(m_rest));

for i = 1:length(m_rest)-1
    vec_Ehs_s1 = Ehs_rest(m_rest(i):m_rest(i)+100);
    [val_s1, idx_s1] = max(vec_Ehs_s1);
    s1_idx_rest(m_rest(i)+idx_s1) = 1;
    s1_val_rest(m_rest(i)+idx_s1) = val_s1;
    R_peak_of_S1_rest(i) = idx_s1 * dt;
end

for i = 1:length(m_rest)-1
    vec_Ehs_s2 = Ehs_rest(m_rest(i) + 100:m_rest(i)+250);
    [val_s2, idx_s2] = max(vec_Ehs_s2);
    s2_idx_rest(m_rest(i) + 100 + idx_s2) = 1;
    s2_val_rest(m_rest(i) + 100 + idx_s2) = val_s2;
    R_peak_of_S2_rest(i) = (idx_s2 + 100) * dt;
end

% active

s1_idx_active = zeros(1, length(hs_active));
s1_val_active = zeros(1, length(hs_active));
s2_idx_active = zeros(1, length(hs_active));
s2_val_active = zeros(1, length(hs_active));
R_peak_of_S1_active = zeros(1, length(m_active));
R_peak_of_S2_active = zeros(1, length(m_active));

for i = 1:length(m_active)
    vec_Ehs_s1 = Ehs_active(m_active(i):m_active(i)+80);
    [val_s1, idx_s1] = max(vec_Ehs_s1);
    s1_idx_active(m_active(i)+idx_s1) = 1;
    s1_val_active(m_active(i)+idx_s1) = val_s1;
    R_peak_of_S1_active(i) = idx_s1 * dt;
end

for i = 1:length(m_active)
    vec_Ehs_s2 = Ehs_active(m_active(i) + 80:m_active(i) + 200);
    [val_s2, idx_s2] = max(vec_Ehs_s2);
    s2_idx_active(m_active(i) + 80 + idx_s2) = 1;
    s2_val_active(m_active(i) + 80 + idx_s2) = val_s2;
    R_peak_of_S2_active(i) = (idx_s2 + 100) * dt;
end

s1_rest = max(hs_rest) + 0.3 * max(hs_rest) * s1_idx_rest;
s1_rest(s1_idx_rest == 0) = nan;

s2_rest = max(hs_rest) + 0.3 * max(hs_rest) * s2_idx_rest;
s2_rest(s2_idx_rest == 0) = nan;

s1_active = max(hs_active) + 0.3 * max(hs_active) * s1_idx_active;
s1_active(s1_idx_active == 0) = nan;

s2_active = max(hs_active) + 0.3 * max(hs_active) * s2_idx_active;
s2_active(s2_idx_active == 0) = nan;

%%
t_rest = time_vec(1:length(hs_rest_norm));
t_active = time_vec(length(hs_rest_norm)+1:end);

% rest plot
figure ('Name', 'S1 & S2 Detection - Rest');
subplot(2,1,1);
plot(t_rest, hs_rest);
xlabel('Time [sec]'); ylabel('Voltage [mV]'); title('S1 & S2 Detection - Rest'); xlim([0 22]);

hold on
scatter(t_rest, s1_rest, 'm*');
scatter(t_rest, s2_rest, 'c*');
legend('Stethoscope', 'S1', 'S2', 'Location', 'eastoutside');
hold off

subplot(2,1,2)
plot(t_rest, ECG_rest);
xlabel('Time [sec]'); ylabel('Voltage [mV]'); title('ECG Signal - Rest'); xlim([0 22]);

legend('Lead I'); 

% active plot
figure ('Name', 'S1 & S2 Detection - Active');
subplot(2,1,1);
plot(t_active, hs_active); 
xlabel('Time [sec]'); ylabel('Voltage [mV]'); title('S1 & S2 Detection - Active'); xlim([22.4 33.6]);

hold on
scatter(t_active, s1_active, 'm*');
scatter(t_active, s2_active, 'c*');
legend('Stethoscope', 'S1', 'S2', 'Location', 'eastoutside');
hold off

subplot(2,1,2)
plot(t_active, ECG_active);
xlabel('Time [sec]'); ylabel('Voltage [mV]'); title('ECG Signal - Active'); xlim([22.4 33.6]);

legend('Lead I'); 

%% Table
%See table in the attached text file.
%rest
BPM_rest = time_vec(m_rest); 
BPM_rest = diff(BPM_rest);
BPM_rest = 60./BPM_rest;
BPM_std_rest = std(BPM_rest);
BPM_mean_rest = mean(BPM_rest);

dt_R_to_S1_rest = time_vec(find(s1_idx_rest==1)) - time_vec(m_rest(1:end-1));
dt_R_to_S2_rest = time_vec(find(s2_idx_rest==1)) - time_vec(m_rest(1:end-1));
dt_S1_to_S2_rest = time_vec(find(s2_idx_rest==1)) - time_vec(find(s1_idx_rest==1));
dt_S2_to_next_S1_rest = diff(time_vec(m_rest(2:end))) - dt_R_to_S2_rest(1:end-1) + dt_R_to_S1_rest(2:end);

% active
BPM_active = time_vec(m_active); 
BPM_active = diff(BPM_active);
BPM_active = 60./BPM_active;
BPM_std_active = std(BPM_active);
BPM_mean_active = mean(BPM_active);

dt_R_to_S1_active = time_vec(find(s1_idx_active==1)) - time_vec(m_active);
dt_R_to_S2_active = time_vec(find(s2_idx_active==1)) - time_vec(m_active);
dt_S1_to_S2_active = time_vec(find(s2_idx_active==1)) - time_vec(find(s1_idx_active==1));
dt_S2_to_next_S1_active = diff(time_vec(m_active)) - dt_R_to_S2_active(1:end-1) + dt_R_to_S1_active(2:end);

% Delta t
mean_dt_rest_S1 = mean(dt_R_to_S1_rest(1:end-1));
STD_dt_rest_S1 = std(dt_R_to_S1_rest(1:end-1));
mean_dt_rest_S2 = mean(dt_R_to_S2_rest(1:end-1));
STD_dt_rest_S2 = std(dt_R_to_S2_rest(1:end-1));

S1_to_S2_rest = dt_R_to_S2_rest-dt_R_to_S1_rest;
S1_to_S2_mean_rest = mean(S1_to_S2_rest);
S1_to_S2_std_rest = std(S1_to_S2_rest);

mean_dt_active_S1 = mean(dt_R_to_S1_active(1:end-1));
STD_dt_active_S1 = std(dt_R_to_S1_active(1:end-1));
mean_dt_active_S2 = mean(dt_R_to_S2_active(1:end-1));
STD_dt_active_S2 = std(dt_R_to_S2_active(1:end-1));

S1_to_S2_active = dt_R_to_S2_active-dt_R_to_S1_active;
S1_to_S2_mean_active = mean(S1_to_S2_active(1:end-1));
S1_to_S2_std_active = std(S1_to_S2_active(1:end-1));

mean_dt_S2_to_S1_rest = mean(dt_S2_to_next_S1_rest);
STD_dt_S2_to_S1_rest = std(dt_S2_to_next_S1_rest);
mean_dt_S2_to_S1_active = mean(dt_S2_to_next_S1_active);
STD_dt_S2_to_S1_active = std(dt_S2_to_next_S1_active);

% table - Peak to Peak from stethoscope

%rest
S1_vec_maxval_rest_t = zeros(1,length(hs_rest));
S1_vec_minval_rest_t = zeros(1,length(hs_rest));
S2_vec_maxval_rest_t = zeros(1,length(hs_rest));
S2_vec_minval_rest_t = zeros(1,length(hs_rest));

for i = 1:length(m_rest)-1
    vec_s = hs_rest(m_rest(i):m_rest(i)+100);
    [max_s1,ind_max_s1] = max(vec_s);
    [min_s1,ind_min_s1] = min(vec_s);
    S1_vec_maxval_rest_t(m_rest(i) + ind_max_s1) = max_s1;
    S1_vec_minval_rest_t(m_rest(i) + ind_min_s1) = min_s1;
end

for i=1:length(m_rest)-1
    vec_s = hs_rest(m_rest(i)+100:m_rest(i)+250);
    [max_s2,ind_max_s2] = max(vec_s);
    [min_s2,ind_min_s2] = min(vec_s);
    S2_vec_maxval_rest_t(m_rest(i)+100+ind_max_s2) = max_s2;
    S2_vec_minval_rest_t(m_rest(i)+100+ind_min_s2) = min_s2;
end

f_max1 = find(S1_vec_maxval_rest_t);
f_min1 = find(S1_vec_minval_rest_t);
peak_to_peak_s1_rest_mean = mean(S1_vec_maxval_rest_t(f_max1)-S1_vec_minval_rest_t(f_min1));
peak_to_peak_s1_rest_std = std(S1_vec_maxval_rest_t(f_max1)-S1_vec_minval_rest_t(f_min1));
f_max2 = find(S2_vec_maxval_rest_t);
f_min2 = find(S2_vec_minval_rest_t);
peak_to_peak_s2_rest_mean = mean(S2_vec_maxval_rest_t(f_max2)-S2_vec_minval_rest_t(f_min2));
peak_to_peak_s2_rest_std = std(S2_vec_maxval_rest_t(f_max2)-S2_vec_minval_rest_t(f_min2));

%active
s1_maxval_active_t=zeros(1,length(t_active));
s1_minval_active_t=zeros(1,length(t_active));
s2_maxval_active_t=zeros(1,length(t_active));
s2_minval_active_t=zeros(1,length(t_active));

for i=1:length(m_active)-1
    vec_s=t_active(m_active(i):m_active(i)+80);
    [max_s1,ind_max_s1]=max(vec_s);
    [min_s1,ind_min_s1]=min(vec_s);
    s1_maxval_active_t(m_active(i)+ind_max_s1)=max_s1;
    s1_minval_active_t(m_active(i)+ind_min_s1)=min_s1;
end

for i=1:length(m_active)-1
    vec_s=t_active(m_active(i)+80:m_active(i)+200);
    [max_s2,ind_max_s2]=max(vec_s);
    [min_s2,ind_min_s2]=min(vec_s);
    s2_maxval_active_t(m_active(i)+80+ind_max_s2)=max_s2;
    s2_minval_active_t(m_active(i)+80+ind_min_s2)=min_s2;
end

f_max1=find(s1_maxval_active_t);
f_min1=find(s1_minval_active_t);
peak_to_peak_s1_active_mean=mean(s1_maxval_active_t(f_max1)-s1_minval_active_t(f_min1));
peak_to_peak_s1_active_std=std(s1_maxval_active_t(f_max1)-s1_minval_active_t(f_min1));
f_max2=find(s2_maxval_active_t);
f_min2=find(s2_minval_active_t);
peak_to_peak_s2_active_mean=mean(s2_maxval_active_t(f_max2)-s2_minval_active_t(f_min2));
peak_to_peak_s2_active_std=std(s2_maxval_active_t(f_max2)-s2_minval_active_t(f_min2));





