%HRV&SNR%
%import raw data
clear all; clc;
data = fopen('Hila&Ilay_part_b.txt','rt');
B = textscan(data, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',1);
fclose(data);
time_vec=[0:length(B{1})-1]/500;fc = 8;fc2=50;fs=length(B{1})/60.244;
b=fir1(48,[fc/(fs/2) fc2/(fs/2)],'bandpass');a=1;%define filter
%%
%finding QRS for later analysis- Lead I
%using algorithm from "QRS detection"-AF2
filt_lead_1=filter(b,a,B{1});thresh1=0.4*max(filt_lead_1);
Y0_1=abs(filt_lead_1);Y1_1=zeros(length(filt_lead_1),1);
for i=1:length(filt_lead_1)
    if filt_lead_1(i)>=thresh1  
        Y1_1(i)=Y0_1(i);
    else
        Y1_1(i)=thresh1;
    end
end
Y2_1=diff(Y1_1);
for j=1:length(Y2_1)
    if Y2_1(j)<0.3*max(Y2_1)  %QRS candidates are set
        Y2_1(j)=0;
    end
end
Y2_1=islocalmax(Y2_1);m=[];
for j=1:length(Y2_1)
    if Y2_1(j)==1
       m(end+1)=j;
    end
end
RR=time_vec(m);RR=diff(RR);HR=60./RR;
HR_mean=mean(HR);HR_std=std(HR);
figure(1);
plot(1:15,HR(1:15),'o');yline([HR_mean,HR_mean-HR_std,HR_mean+HR_std],'-', ...
    {'HR mean','HR mean-std','HR mean+std'});
xlabel('Beat number');ylabel('Heart Rate(Bpm)');title('HR/beats');
%%
%HRV
HRV_cell={};
for i=m
    HRV_cell{end+1}=B{1}(i-140:i+210);
end
mean_sig=mean(cell2mat(HRV_cell),2);
figure(2);
subplot(3,1,1);plot(time_vec(1:length(mean_sig)),mean_sig);xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Mean Heart beat cycle')
subplot(3,1,2);plot(time_vec(1:length(mean_sig)),HRV_cell{1});xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('First Heart beat cycle')
subplot(3,1,3);plot(time_vec(1:length(mean_sig)),HRV_cell{64});xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Last Heart beat cycle')

%%
%mean voltage and STD - mean signal
mean_volt_mean = mean(mean_sig(185:225));
std_volt_mean = std(mean_sig(185:225));

%mean voltage and STD - first signal
mean_volt_first = mean(HRV_cell{1}(185:225));
std_volt_first = std(HRV_cell{1}(185:225));

%mean voltage and STD - last signal
mean_volt_last = mean(HRV_cell{64}(185:225));
std_volt_last = std(HRV_cell{64}(185:225));

%%
%SNR
% mean signal vs. first signal
mean_max_std = max(mean_sig)/std_volt_mean;
first_max_std = max(HRV_cell{1})/std_volt_first;

new_mean_sig = mean_sig./std_volt_mean;
figure(3);
subplot(3,1,1);plot(time_vec(1:length(mean_sig)),mean_sig);xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Heart beat cycle');
hold on;
plot(time_vec(1:length(HRV_cell{1})),HRV_cell{1});legend('Average cycle','First cycle');
subplot(3,1,2);plot(time_vec(1:length(mean_sig)),new_mean_sig);xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Heart beat cycle / Isoelectric STD')
hold on; 
plot(time_vec(1:length(HRV_cell{1})),HRV_cell{1}./std_volt_first);legend('Average cycle','First cycle');
subplot(3,1,3);plot(time_vec(1:length(mean_sig)),new_mean_sig-mean_sig);xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('SNR Improvment')
hold on;
plot(time_vec(1:length(HRV_cell{1})),HRV_cell{1}./std_volt_first-HRV_cell{1});legend('Average cycle','First cycle');
