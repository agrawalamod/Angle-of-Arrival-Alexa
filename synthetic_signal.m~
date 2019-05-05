close all;
clear all;

freq=1000;
c = 343;
lambda = c/freq;
TxD = 20*lambda;
t = 1;

x_loc = 2;
y_loc = 4;
d = 4.5e-2;
Fs = 44100;
ref_mic = 5;
ang = [120 60 0 300 240 180];
mic_loc = [x_loc+d*cosd(ang') y_loc+ d*sind(ang')];

theta=50;
source_loc = [x_loc+TxD*cosd(theta) y_loc+TxD*sind(theta)];

figure;
scatter(mic_loc(:,1), mic_loc(:,2));
hold on;
scatter(source_loc(:,1),source_loc(:,2));

mic_dist = mic_loc - source_loc;
dist = sqrt(sum(mic_dist.^2,2));
time_diff = dist/c;
phase_d = dist*2*pi/lambda;

fy = phase_d;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

final_result = aoa(fy, freq, c, Fs, mic_dis_projections);
[pks, locs] = findpeaks(final_result, 'MinPeakDistance', 20);
disp(locs);
figure;
plot(0:359, final_result);
hold on;
scatter(locs,pks, 'v');


%[first_mic, mic_delays] = first_mic_solver(y);

function [result] = aoa(fy, freq, c, Fs, mic_dis_projections)
    lambda = c/freq;
  
    phase_diff = mic_dis_projections*2*pi/lambda;
    M = exp(1i*-1*phase_diff');
    spectrum = (M * fy)'; % 360x6 * 6x1
    spec_result= [];
    for idx =1:size(spectrum,2)
        spec_result = [spec_result real(spectrum(idx))];
    end
    
    result = spec_result;
    
    %figure;
    %plot(0:359,spec_result);
    %[pks,locs] = findpeaks(spec_result);
    %result = locs;
    %disp([locs, calc_ref_mic, ang(calc_ref_mic)]);
    %disp(outputAngle(2,:));
end


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