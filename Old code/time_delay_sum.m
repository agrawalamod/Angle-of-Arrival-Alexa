clear all;
close all;
[y,Fs] = audioread('./train/A01_X01.wav');
y = y(10000:13500,1:6);


d = 4.75e-2;
freq = 500;
c = 343;
lambda = c/freq;
passfreq = 500;
x_loc=4; y_loc=4;

for idx = 1:6 
    y(:,idx) = lowpass(y(:,idx), passfreq, Fs);
end

[first_mic, mic_delays] = first_mic_solver(y);

ref_mic = first_mic(1);
ang = [120 60 0 300 240 180];
mic_loc = [x_loc+d*cosd(ang') y_loc+d*sind(ang')];

mic_dis = [];
for i=1:6
    dis_vec = [mic_loc(i,1)-mic_loc(ref_mic,1) mic_loc(i,2)-mic_loc(ref_mic,2)];
    %mic_dis = [mic_dis; dis_vec/norm(dis_vec)];
    mic_dis = [mic_dis; dis_vec];
end

direc_vecs = [];
for tau = 0:359
    tau_vec = [cosd(tau) sind(tau)];
    direc_vecs = [direc_vecs; tau_vec/norm(tau_vec)];
end

mic_dis_projections = mic_dis * direc_vecs';
time_diff = mic_dis_projections/c;
phase_diff = mic_dis_projections*2*pi./lambda;

ref_signal = y(:,ref_mic(1))';

spec_result= [];

for idx=1:360
    delays=time_diff(:,idx);
    delayed_signal=[];
    for idx2=1:6
        delayed_signal = [delayed_signal; delayseq(ref_signal,delays(idx2),Fs)];
    end
    signal_sum = norm(sum(delayed_signal,1));
    spec_result = [spec_result signal_sum];
    disp(idx);
end


%{
M = exp(1i*phase_diff);

spectrum = y * M;
for idx =1:size(spectrum,2)
    spec_result = [spec_result sum(real(spectrum(:,idx)))];
end
%}
figure;
plot(0:359,spec_result);


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

