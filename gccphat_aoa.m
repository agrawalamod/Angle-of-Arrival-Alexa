close all;
clear all;

[raw_y,Fs] = audioread('./train/A03_X02.wav');
ans_X=2;

x_loc = 5;
y_loc = 5;
d = 4.75e-2;
c = 343;
ref_mic = 5;
ang = [120 60 0 300 240 180];
mic_loc = [x_loc+d*cosd(ang') y_loc+ d*sind(ang')];

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
save('mic_dis_projections.mat', 'mic_dis_projections');

time_diff = mic_dis_projections/c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('outputAngle.mat');
raw_y = raw_y(:,1:6);

[first_mic, mic_delays] = first_mic_solver(raw_y);

ref_sig = raw_y(:,first_mic(1));
[max_val, max_indx] = max(ref_sig);
timeband = [max_indx-500 max_indx+25000];

%chops the signal
%fy = raw_y(timeband(1):timeband(2),:);
fy = raw_y;

refsig = fy(:,ref_mic);

fyt = [];
fyt2 = [];
for idx=1:6
    tau = gccphat(refsig, fy(:,idx), Fs);
    tau2 = interpolated_gccphat(refsig, fy(:,idx), Fs);
    fyt= [fyt;tau];
    fyt2 = [fyt2; tau2];
end

result = [];
result2 = [];
for alpha = 1:360
    result = [result sum(abs(time_diff(:,alpha)-fyt))];
    result2 = [result2 sum(abs(time_diff(:,alpha)-fyt2))];
end

figure;
plot(0:359, result);
hold on;
plot(0:359, result2);
[min_val, min_indx] = min(result);
[min_val2, min_indx2] = min(result2);
disp([min_indx min_indx2]);
disp(outputAngle(ans_X,:));
dist = 343*100*fyt;
dist2 = 343*100*fyt2;
disp(dist);
disp(dist2);

function [tau] = interpolated_gccphat(y1, y2, Fs)
    esp = 0.000001;
    N = length(y1);
    
    fft_y1 = real(fft(y1,2*N));
    fft_y2= real(fft(y2,2*N));
    fft_prod = fft_y1 .* conj(fft_y2);
    fft_prod_mag = sum(abs(fft_prod));
    Gphat = fft_prod./(fft_prod_mag+esp);
    Rphat = ifft(Gphat, 2*N)';
    
    x = horzcat(Rphat(N+1:end),Rphat(1:N));
    lag = -N:N;
    
    
   %interpolation
   %[corr_val, corr_index] = max(x(N-8:N+8));
   %corr_index = corr_index + (N-8);
   [corr_val, corr_index] = max(x);
   tau_approx = lag(corr_index)/Fs;
   
   dk = 0.5*(x(corr_index-1) - x(corr_index+1)) / (x(corr_index+1) + x(corr_index-1) - 2* x(corr_index));
   if(abs(dk)>=1)
       dk=0;
   end
   t_offset = dk/Fs;
   tau = tau_approx + t_offset;
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
