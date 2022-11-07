% Kalman filter based on autoregressive model

clc
clear
close all;

load('SampleECG1.mat');
% load('dataEEG.txt'); data = dataEEG(:,11)';
% data = load('eeg.txt')'; %data = data(5001:20000)';

load LPFilterCoefs50Hzfs250Hz h
fs = 250;

% x = filtfilt(h,1,data);
x = data;


% x = data + (5.5+sin(2*pi*n*.1/fs)).*sin(2*pi*n*f0/fs)/std(data);


% [b,a] = ellip(20,1,50,20/125);
% load tempcoefs a
% a = [1 -.9];
% a = [1 -2*(.99) (.99)^2];
% a = [.1 2 3 4 3 2 .1];
% b = 1;
% beta = .5;

mdl = arx(x', 12);
a = get(mdl, 'a');
b = 1;
beta = get(mdl, 'NoiseVariance');
gamma = 0.1;%0.9;
wlen = round(fs*.250);

[y1f,y1s,Pbar,Phat,PSmoothed,Kgain] = KalmanARFilter(x, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');
% [y2f,y2s,Pbar,Phat,PSmoothed,Kgain] = KalmanARFilter(x, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'cntrl');

t = (0:length(data)-1)/fs;

figure;
hold on;
lbl = {};
plot(t, x, 'b'); lbl = cat(2, lbl, 'Noisy signal');
plot(t, x - y1s, 'r'); lbl = cat(2, lbl, 'Kalman filter denoised');
plot(t, x - y1f, 'm');  lbl = cat(2, lbl, 'Kalman smoother denoised');
grid;
xlabel('time (sec.)');
legend(lbl);

figure;
hold on;
lbl = {};
pwelch(x, length(x)/2, length(x)/8, 1024, fs); lbl = cat(2, lbl, 'Noisy signal');
pwelch(x - y1f, length(x)/2, length(x)/8, 1024, fs); lbl = cat(2, lbl, 'Kalman filter denoised');
pwelch(x - y1s, length(x)/2, length(x)/8, 1024, fs); lbl = cat(2, lbl, 'Kalman smoother denoised');
legend(lbl);

figure;
freqz(b,a,256,fs);


% % % psd(y2,256,fs);

% % % fid = fopen('KalmanFilter.txt','wt');
% % % fprintf(fid,'%f\n',y1);
% % % fclose(fid);
% % % fid = fopen('KalmanSmoother.txt','wt');
% % % fprintf(fid,'%f\n',y2);
% % % fclose(fid);


% % % psd(y1,256,fs);

% % % figure;
% % % plot(Kgain);
% % % grid;

% % % title('noisy spectrum');
% % % 
% % % % % % figure;
% % % % % % psd(y1,1000,fs);
% % % % % % title('Kalman filter output spectrum');
% % % 
% % % figure;
% % % psd(y2,1000,fs);
% % % title('Kalman smoother output spectrum');
