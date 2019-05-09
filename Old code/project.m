clear all;
close all;

[y,Fs] = audioread('./train/A01_X01.wav');
load('outputAngle.mat');
y = y(:,1:6);

d = 4.75e-2;
freq = 600;
c = 343;
lambda = c/freq;
passfreq = [freq-0.001 freq+0.001];
x_loc=4; y_loc=4;

zmag=[];
zang=[];

for idx = 1:6 
    yF(idx,:) = fft(y(:,idx), Fs);
    zmag = [zmag abs(yF(idx,freq))];
    zang = [zang angle(yF(idx,freq))];
end

%figure;
%plot(1:size(y,1),y(:,1));
%hold on;
%plot(1:size(fy,1),fy(:,1));
%yF = fftshift(fft(fy(:,1),Fs));
%figure;
%plot(-Fs/2:(Fs/2)-1,yF);]]

t=0:2*freq;
freq_y = [];

for idx = 1:6
    fy = zmag(idx) * exp((1i*2*pi*freq*t) + zang(idx));
    freq_y = [freq_y fy'];
end

[first_mic, mic_delays] = first_mic_solver(y);

calc_ref_mic = first_mic(1);
ref_mic = 5;
ang = [120 60 0 300 240 180];
mic_loc = [x_loc+d*cosd(ang') y_loc+ d*sind(ang')];

mic_dis = [];
for i=1:6
    dis_vec = [mic_loc(i,1)-mic_loc(ref_mic,1) mic_loc(i,2)-mic_loc(ref_mic,2)];
    mic_dis = [mic_dis; dis_vec];
end

direc_vecs = [];
for tau = 0:359
    tau_vec = [cosd(tau) sind(tau)];
    direc_vecs = [direc_vecs; tau_vec/norm(tau_vec)];
end

mic_dis_projections = mic_dis * direc_vecs';
phase_diff = mic_dis_projections*2*pi./lambda;

M = exp(1i*phase_diff);

spectrum = freq_y * M;
spec_result= [];
for idx =1:size(spectrum,2)
    spec_result = [spec_result sum(real(spectrum(:,idx)))];
end

figure;
plot(0:359,spec_result);
[pks,locs] = findpeaks(spec_result);
disp([locs, calc_ref_mic, ang(calc_ref_mic)]);
disp(outputAngle(1,:))


function [first_mic, mic_delays] = first_mic_solver(y)
    mic_delays = [];
    for i=1:6
        [res, lags] = xcorr(y(:,1), y(:,i));
        [max_v, max_i] = max(res);
        mic_delays = [mic_delays lags(max_i)];
    end
    [least_lag] = max(mic_delays);
    first_mic = find(mic_delays == least_lag);
end

