%QRS detection%
%import raw data
clear all; clc;
data = fopen('Hila&Ilay.txt','rt');
A = textscan(data, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',1);
time_vec=[0:length(A{1})-1]/500;
fclose(data);
fc = 8;fc2=50;fs=length(A{1})/49.962; 
%for IIR filter
%[b,a] = butter(2,[fc/(fs/2) fc2/(fs/2)],'bandpass');%freqz(b,a); %to see the filter's graph
%for FIR filter
b=fir1(48,[fc/(fs/2) fc2/(fs/2)],'bandpass');a=1;

%% 
%ECG before manipulation
figure(1);
subplot(6,1,1);plot(time_vec,A{1});xlabel('Time(sec)');ylabel('Amplitude(mV)');title('Lead I');
subplot(6,1,2);plot(time_vec(1:5000),A{1}(1:5000));xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead I-Closeup');
subplot(6,1,3);plot(time_vec,A{2});xlabel('Time(sec)');ylabel('Amplitude(mV)');title('Lead III');
subplot(6,1,4);plot(time_vec(1:5000),A{2}(1:5000));xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead III-Closeup');
subplot(6,1,5);plot(time_vec,A{3});xlabel('Time(sec)');ylabel('Amplitude(mV)');title('Lead II');
subplot(6,1,6);plot(time_vec(1:5000),A{3}(1:5000));xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead II-Closeup');
%% 
%filtering
figure(2);
subplot(3,1,1);filt_lead_1=filter(b,a,A{1});plot(time_vec,filt_lead_1);hold on;
plot(time_vec,A{1});legend('filtered','unfiltered');xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead I');
subplot(3,1,2);filt_lead_3=filter(b,a,A{2});plot(time_vec,filt_lead_3);hold on;
plot(time_vec,A{2});legend('filtered','unfiltered');xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead III');
subplot(3,1,3);filt_lead_2=filter(b,a,A{3});plot(time_vec,filt_lead_2);hold on;
plot(time_vec,A{3});legend('filtered','unfiltered');xlabel('Time(sec)');
ylabel('Amplitude(mV)');title('Lead II');
%%
%ECG fourier tranform
L=length(A{1});
lead_1_fourier=fft(A{1});filt_lead_1_fourier=fft(filt_lead_1);
lead_3_fourier=fft(A{2});filt_lead_3_fourier=fft(filt_lead_3);
lead_2_fourier=fft(A{3});filt_lead_2_fourier=fft(filt_lead_2);
P2org1 = abs(lead_1_fourier/L);P2org2 = abs(filt_lead_1_fourier/L);P2org3 = abs(lead_3_fourier/L);
P2org4 = abs(filt_lead_3_fourier/L);P2org5 = abs(lead_2_fourier/L);P2org6 = abs(filt_lead_2_fourier/L);
P1org1 = P2org1(1:L/2+1);P1org2 = P2org2(1:L/2+1);P1org3 = P2org3(1:L/2+1);
P1org4 = P2org4(1:L/2+1);P1org5 = P2org5(1:L/2+1);P1org6 = P2org6(1:L/2+1);
P1org1(2:end-1) = 2*P1org1(2:end-1);P1org2(2:end-1) = 2*P1org2(2:end-1);
P1org3(2:end-1) = 2*P1org3(2:end-1);P1org4(2:end-1) = 2*P1org4(2:end-1);
P1org5(2:end-1) = 2*P1org5(2:end-1);P1org6(2:end-1) = 2*P1org6(2:end-1);
freq_vec = fs*(0:(L/2))/L; % Define the frequency domain f

figure(3);
subplot(3,1,1);plot(freq_vec,P1org1);hold on;
plot(freq_vec,P1org2);legend('unfiltered','filtered');xlabel('Frequency(Hz)');
ylabel('Amplitude(mV)');title('Lead I');
subplot(3,1,2);plot(freq_vec,P1org3);hold on;
plot(freq_vec,P1org4);legend('unfiltered','filtered');xlabel('Frequency(Hz)');
ylabel('Amplitude(mV)');title('Lead III');
subplot(3,1,3);plot(freq_vec,P1org5);hold on;
plot(freq_vec,P1org6);legend('unfiltered','filtered');xlabel('Frequency(Hz)');
ylabel('Amplitude(mV)');title('Lead II');

%%
%QRS recognition- AF2
%only one lead is required, Lead II was chosen
thresh1=0.4*max(filt_lead_1);thresh3=0.4*max(filt_lead_3);thresh2=0.4*max(filt_lead_2);%setting thresh
Y0_1=abs(filt_lead_1);Y0_2=abs(filt_lead_2);Y0_3=abs(filt_lead_3);%taking the absolute val of signal
Y1_1=zeros(length(filt_lead_2),1);Y1_2=zeros(length(filt_lead_2),1);Y1_3=zeros(length(filt_lead_2),1);
%condition 1 of AF2 algorithm
for i=1:length(filt_lead_1)
    if filt_lead_1(i)>=thresh1  
        Y1_1(i)=Y0_1(i);
    else
        Y1_1(i)=thresh1;
    end
end
for i=1:length(filt_lead_2)
    if filt_lead_2(i)>=thresh2
        Y1_2(i)=Y0_2(i);
    else
        Y1_2(i)=thresh2;
    end
end
for i=1:length(filt_lead_3)
    if filt_lead_3(i)>=thresh3
        Y1_3(i)=Y0_3(i);
    else
        Y1_3(i)=thresh3;
    end
end
%derivative of the signal after condition 1
Y2_1=diff(Y1_1);Y2_2=diff(Y1_2);Y2_3=diff(Y1_3);
%condition 2 of AF2 algorithm, different threshold was set then the
%original algorithm, appropriate to our signal
figure(4);
for j=1:length(Y2_2)
    if Y2_2(j)<0.3*max(Y2_2)  %QRS candidates are set
        Y2_2(j)=0;
    end
end
%taking the local maximum from candidates gives us the extract high peaks of
%candidates
Y2_2=islocalmax(Y2_2);m=[];
for j=1:length(Y2_2)
    if Y2_2(j)==1
       m(end+1)=j;
    end
end
plot(time_vec,A{3},m/500,ones(length(m),1)./2.1,'o');
xline([30 40],'-',{'Standing','Standing & Deep Breaths'});
xlabel('Time(sec)');ylabel('Amplitude(mV)');title('Lead II with QRS');
%%
%confusion matrix
pred=zeros(length(A{1}),1);true=zeros(length(A{1}),1);
for i=1:length(m)
    pred(m(i))=1;
    true(m(i))=1;
end
true([17414 21292])=1;
figure(6);
confusionchart(true,pred);
title('Confusion Matrix: AF2 QRS detection algorithm')
set(gca,'FontSize',30);
%%
%HR
m(end+1)=17414;m(end+1)=21292;
m=sort(m);%QRS who werent detected added manually to compute HR
g1=Y2_2(1:15000);g2=Y2_2(1:20000);g3=Y2_2;
g1=sum(g1);g2=sum(g2);g3=sum(g3);
m1=m(1:g1);m2=m(g1+1:g2);m3=m(g2+1:g3);
HR_sitting=time_vec(m1);HR_standing=time_vec(m2);HR_standing_b=time_vec(m3);
HR_sitting=diff(HR_sitting);HR_standing=diff(HR_standing);HR_standing_b=diff(HR_standing_b);
HR_sitting=60./HR_sitting;HR_standing=60./HR_standing;HR_standing_b=60./HR_standing_b;
std1=std(HR_sitting);std2=std(HR_standing);std3=std(HR_standing_b);
HR_sitting=mean(HR_sitting,2);HR_standing=mean(HR_standing,2);HR_standing_b=mean(HR_standing_b,2);
%%
%Einthoven's law confirmation
Einth=A{1}+A{2}-A{3};
figure(5);
subplot(3,1,1);plot(time_vec,A{1}+A{2});xlabel('time(sec)');ylabel('Amplitude(mV)');
title('Lead I+Lead III');
subplot(3,1,2);plot(time_vec,A{3});xlabel('time(sec)');ylabel('Amplitude(mV)');
title('Lead II');
subplot(3,1,3);plot(time_vec,Einth);xlabel('time(sec)');ylabel('Amplitude(mV)');
title('Lead I+Lead III-Lead II');
%%
%T-test: standard breathing vs heavy breathing while standing
HR_sitting=time_vec(m1);HR_standing=time_vec(m2);HR_standing_b=time_vec(m3);
HR_sitting=diff(HR_sitting);HR_standing=diff(HR_standing);HR_standing_b=diff(HR_standing_b);
HR_sitting=60./HR_sitting;HR_standing=60./HR_standing;HR_standing_b=60./HR_standing_b;
[h,p]=ttest2(HR_standing,HR_standing_b(1:9))
