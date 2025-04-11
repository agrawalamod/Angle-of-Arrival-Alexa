% aoa_pipeline.m
clear; close all; clc;

% Load ground truth
load('outputAngle.mat');  % matrix of shape [10 x 3]

% Parameters
num_signals = 10;
num_arrays = 3;
mic_count = 6;
mic_spacing = 0.0475; % meters
sound_speed = 343;    % m/s
results = zeros(num_signals, num_arrays, 3);  % [signal, array, algorithm]

for signal_index = 1:num_signals
    for array_index = 1:num_arrays
        file_name = sprintf('./train/A%02d_X%02d.wav', array_index, signal_index);
        fprintf('Processing %s\n', file_name);
        [raw_y, Fs] = audioread(file_name);
        raw_y = raw_y(:, 1:mic_count);  % keep first 6 mics

        [deg_gcc, deg_interp] = estimate_gccphat_full(raw_y, Fs, mic_spacing, sound_speed, 0, '');
        %deg_music = estimate_music(raw_y, Fs, mic_spacing, sound_speed);
        deg_music = estimate_music_wideband(raw_y, Fs, mic_spacing, sound_speed);


        true_deg = outputAngle(signal_index, array_index);
        err_gcc = angular_diff(deg_gcc, true_deg);
        err_interp = angular_diff(deg_interp, true_deg);
        err_music = angular_diff(deg_music, true_deg);

        results(signal_index, array_index, :) = [err_gcc, err_interp, err_music];
    end
end
% Plot CDF
errors = reshape(results, [], 3);
figure;
hold on;
for i = 1:3
    [f, x] = ecdf(errors(:,i));
    plot(x, f, 'LineWidth', 2);
end
legend({'GCC-PHAT', 'GCC-PHAT Interp', 'MUSIC'});
xlabel('AoA Error (degrees)');
ylabel('CDF');
title('AoA Estimation Error CDF');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Local functions start below this line ==============
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau = gccphat_generated(y1, y2, Fs)
    N = length(y1);
    Y1 = fft(y1, 2*N);
    Y2 = fft(y2, 2*N);
    G = Y1 .* conj(Y2);
    G = G ./ (abs(G) + eps);
    R = ifft(G);
    [~, idx] = max(abs(R));
    lag = idx - 1;
    if lag > N
        lag = lag - 2*N;
    end
    tau = lag / Fs;
end

function diff = angular_diff(est, true_angle)
    diff = abs(mod(est - true_angle + 180, 360) - 180);
end

function [deg1, deg2] = estimate_gccphat_full(y, Fs, d, c, plotting, figname)
    ang = [120 60 0 300 240 180];
    mic_pos = d * [cosd(ang); sind(ang)]';
    ref_mic = 1;
    mic_dis = [];
    for i = 1:size(mic_pos, 1)
        if i ~= ref_mic
            delta = mic_pos(i,:) - mic_pos(ref_mic,:);
            mic_dis = [mic_dis; delta, 0];  % make it 3D by padding z = 0
        end
    end
    res = 1;
    azi_grid = 0:res:359;
    ele_grid = linspace(0, 90, 360);
    [azi_2d, ele_2d] = meshgrid(azi_grid, ele_grid);

    [deg1, deg2] = gcc_phat_aoa(y(:,1:6), Fs, c, res, azi_2d, ele_2d, mic_dis, plotting, figname);
end





function angle = estimate_music(y, Fs, d, c)
    ang = [120 60 0 300 240 180];
    mic_pos = d * [cosd(ang); sind(ang)]';
    snapshot = y(1:2048, :);
    snapshot = snapshot - mean(snapshot);
    R = snapshot' * snapshot / size(snapshot, 1);
    [E, D] = eig(R);
    [~, idx] = sort(diag(D), 'descend');
    En = E(:, 2:end);  % assume 1 source

    theta_scan = 0:1:359;
    f0 = 1000;
    lambda = c / f0;
    P = zeros(size(theta_scan));
    for i = 1:length(theta_scan)
        theta = theta_scan(i) * pi / 180;
        d_vec = [cos(theta); sin(theta)];
        delays = (mic_pos * d_vec) / c;
        a = exp(-1j * 2 * pi * f0 * delays);
        P(i) = 1 / (a' * (En * En') * a);
    end
    [~, max_idx] = max(P);
    angle = theta_scan(max_idx);
end

function angle = estimate_music_wideband(y, Fs, d, c)
    % Estimate AoA using wideband MUSIC for a hexagonal mic array

    % === Mic geometry ===
    ang = [120 60 0 300 240 180];  % angles in degrees
    mic_pos = d * [cosd(ang); sind(ang)]';  % 6 x 2 positions in meters

    % === Frame & FFT settings ===
    frame_len = 1024;        % can be 512 or 2048 depending on SNR
    fft_len = 2048;          % zero-padded FFT
    hop_size = frame_len;    % non-overlapping
    num_mics = size(y, 2);

    % === Sliding window FFT ===
    num_frames = floor((size(y, 1) - frame_len) / hop_size) + 1;
    Y = zeros(fft_len/2, num_mics, num_frames);  % freq x mic x time

    for i = 1:num_frames
        idx_start = (i-1)*hop_size + 1;
        frame = y(idx_start:idx_start+frame_len-1, :);
        frame = frame - mean(frame, 1);  % zero mean per channel
        fft_frame = fft(frame, fft_len);
        Y(:, :, i) = fft_frame(1:fft_len/2, :);
    end

    % === Frequency energy to select top bins ===
    energy_per_bin = squeeze(sum(sum(abs(Y).^2, 2), 3));  % freq x 1
    [~, sorted_idx] = sort(energy_per_bin, 'descend');
    top_bins = sorted_idx(1:20);  % pick top 20 frequency bins

    % === MUSIC Spectrum over angles ===
    theta_scan = 0:1:359;
    P_avg = zeros(size(theta_scan));

    for bin = top_bins'
        f = (bin-1) * Fs / fft_len;
        lambda = c / f;
        k = 2 * pi / lambda;

        % === Average spatial covariance R over all frames ===
        R = zeros(num_mics, num_mics);
        for t = 1:num_frames
            snap = squeeze(Y(bin, :, t)).';  % mic x 1
            R = R + (snap * snap');
        end
        R = R / num_frames;

        % === Eigen-decomposition ===
        [E, D] = eig(R);
        [~, idx] = sort(diag(D), 'descend');
        En = E(:, 2:end);  % noise subspace

        % === MUSIC spectrum scan ===
        for i = 1:length(theta_scan)
            theta = theta_scan(i) * pi / 180;
            direction = [cos(theta); sin(theta)];
            delays = (mic_pos * direction) / c;  % 6 x 1
            a = exp(-1j * 2 * pi * f * delays);
            a = a / norm(a);  % normalize steering vector
            P_avg(i) = P_avg(i) + 1 / abs(a' * (En * En') * a);
        end
    end

    % === Normalize and pick peak ===
    P_avg = real(P_avg);
    P_avg = P_avg / max(P_avg);
    [~, max_idx] = max(P_avg);
    angle = theta_scan(max_idx);
end

function angle = estimate_music_wideband_voting(y, Fs, d, c)
    % Wideband MUSIC with per-bin AoA voting for hex mic array

    % === Hexagonal mic geometry ===
    ang = [120 60 0 300 240 180];
    mic_pos = d * [cosd(ang); sind(ang)]';  % 6 x 2 in meters
    num_mics = size(mic_pos, 1);

    % === STFT parameters ===
    frame_len = 1024;
    fft_len = 2048;
    hop = frame_len;
    num_frames = floor((size(y,1) - frame_len) / hop) + 1;

    % === FFT buffer ===
    Y = zeros(fft_len/2, num_mics, num_frames);
    for t = 1:num_frames
        idx = (t-1)*hop + 1;
        frame = y(idx:idx+frame_len-1, :);
        frame = frame - mean(frame, 1);
        frame_fft = fft(frame, fft_len);
        Y(:,:,t) = frame_fft(1:fft_len/2, :);
    end

    % === Bin energy calculation ===
    bin_energy = squeeze(sum(sum(abs(Y).^2, 2), 3));  % freq x 1
    [~, sorted_bins] = sort(bin_energy, 'descend');
    top_bins = sorted_bins(1:20);  % use top 20 bins only

    % === Per-bin AoA voting ===
    theta_scan = 0:1:359;
    angle_votes = zeros(size(theta_scan));

    for bin = top_bins'
        f = (bin - 1) * Fs / fft_len;
        lambda = c / f;

        % === Build R by averaging across frames with energy threshold ===
        R = zeros(num_mics, num_mics);
        used_frames = 0;
        for t = 1:num_frames
            snap = squeeze(Y(bin,:,t)).';
            if norm(snap)^2 < 1e-6  % frame energy too low
                continue;
            end
            R = R + (snap * snap');
            used_frames = used_frames + 1;
        end
        if used_frames < 2
            continue;  % skip bin if not enough good frames
        end
        R = R / used_frames;

        % === Eigen-decomposition ===
        [E, D] = eig(R);
        [~, idx] = sort(diag(D), 'descend');
        En = E(:, 2:end);  % assume 1 source

        % === MUSIC spectrum scan ===
        P = zeros(size(theta_scan));
        for i = 1:length(theta_scan)
            theta = theta_scan(i) * pi/180;
            dir = [cos(theta); sin(theta)];
            delays = (mic_pos * dir) / c;
            a = exp(-1j * 2 * pi * f * delays);
            a = a / norm(a);  % normalize
            P(i) = 1 / abs(a' * (En * En') * a);
        end

        % === Vote for the angle with max response (weighted) ===
        [~, peak_idx] = max(real(P));
        angle_votes(peak_idx) = angle_votes(peak_idx) + bin_energy(bin);
    end

    % === Final AoA decision ===
    [~, best_idx] = max(angle_votes);
    angle = theta_scan(best_idx);
end



function [deg1, deg2] = gcc_phat_aoa(fy, Fs, c, res, azi_2d, ele_2d, mic_dis, plotting, figname)
    ref_mic = 1;
    refsig = fy(:, ref_mic);
    fy_tau1 = [];
    fy_tau2 = [];
    
    for idx = 1:6
        if ref_mic == idx
            tau1 = 0;
            tau2 = 0;
        else
            %tau2 = interpolated_gccphat(refsig, fy(:,idx), Fs);
            tau2 = interpolated_gccphat_generated(refsig, fy(:,idx), Fs);
            tau1 = gccphat_generated(refsig, fy(:,idx), Fs);
        end
        fy_tau1 = [fy_tau1; tau1];
        fy_tau2 = [fy_tau2; tau2];
    end
    
    fy_dis1 = fy_tau1(2:end) * c;
    fy_dis2 = fy_tau2(2:end) * c;
    scores1 = zeros(size(azi_2d));
    scores2 = zeros(size(azi_2d));
    for A=1:size(azi_2d,2)
        for E=1:size(azi_2d,1)
            azi = azi_2d(E,A)*pi/180;
            ele = ele_2d(E,A)*pi/180;
            vec = [cos(azi)*cos(ele) sin(azi)*cos(ele) sin(ele)];
            scores1(E,A)= 1/norm(mic_dis*vec' - fy_dis1);
            scores2(E,A)= 1/norm(mic_dis*vec' - fy_dis2);            
        end
    end

    if plotting == 1
        fig = figure('visible','off');
        imagesc(scores2);
        xlabel('Azimuth (0-360)');
        ylabel('Elevation (0-90)');
        title(['AoA space search: 1/distance - ' figname]);
        saveas(fig, ['./results/images_updated/' figname '.png']);
    end

    result1 = sum(scores1,1);
    result2 = sum(scores2,1);

    [~, max_indx1] = max(result1);
    [~, max_indx2] = max(result2);

    deg1 = res * max_indx1;
    deg2 = res * max_indx2;
end
function tau = interpolated_gccphat_generated(y1, y2, Fs)
    % Robust interpolated GCC-PHAT using 3-point quadratic peak fitting

    N = length(y1);
    eps_val = 1e-6;

    % === FFTs and GCC-PHAT computation ===
    Y1 = fft(y1, 2*N);
    Y2 = fft(y2, 2*N);
    G = Y1 .* conj(Y2);
    G = G ./ (abs(G) + eps_val);
    R = real(ifft(G));
    R = fftshift(R);  % center the zero delay
    lags = (-N):(N-1);  % integer lags

    % === Find peak index ===
    [~, peak_idx] = max(R);

    % === Check if we can interpolate safely ===
    if peak_idx <= 1 || peak_idx >= length(R)
        % Out of bounds — fallback to sample resolution
        tau = lags(peak_idx) / Fs;
        return;
    end

    % === 3-point parabolic interpolation ===
    y_m1 = R(peak_idx - 1);
    y_0  = R(peak_idx);
    y_p1 = R(peak_idx + 1);

    denom = y_m1 - 2*y_0 + y_p1;
    if abs(denom) < eps_val
        dk = 0;  % flat parabola, fallback to peak
    else
        dk = 0.5 * (y_m1 - y_p1) / denom;
    end

    % === Clip interpolation offset if extreme ===
    if abs(dk) > 1
        dk = 0;
    end

    % Final interpolated delay (in seconds)
    tau = (lags(peak_idx) + dk) / Fs;
end

function tau = interpolated_gccphat(y1, y2, Fs)
    esp = 1e-6;
    N = length(y1);
    ffta_y1 = fft(y1, 2*N);
    ffta_y2 = fft(y2, 2*N);
    fft_y1 = ffta_y1(1:N+1);
    fft_y2 = ffta_y2(1:N+1);
    fft_prod = fft_y1 .* conj(fft_y2);
    fft_prod_mag = abs(fft_prod);
    Gphat = fft_prod ./ (fft_prod_mag + esp);
    Rphat = irfft(Gphat, 1, 2*N);
    x = [Rphat(N+1:end), Rphat(1:N)];
    lag = -N:N;
    [~, corr_index] = max(x);
    tau_approx = lag(corr_index)/Fs;
    if corr_index > 1 && corr_index < length(x)
        dk = 0.5 * (x(corr_index-1) - x(corr_index+1)) / ...
            (x(corr_index+1) + x(corr_index-1) - 2 * x(corr_index));
    else
        dk = 0;
    end
    dk = min(max(dk, -1), 1); % clip extreme values
    t_offset = dk / Fs;
    tau = tau_approx + t_offset;
end

function y = irfft(x, even, N)
    if even
        n = 2 * (length(x) - 1);
        s = length(x) - 1;
    else
        n = 2 * (length(x) - 1) + 1;
        s = length(x);
    end
    xn = zeros(1, n);
    xn(1:s+1) = x;
    xn(s+2:end) = conj(x(s:-1:2));
    y = ifft(xn, N);
end