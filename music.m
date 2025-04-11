function music(file_path, d)
    % file_path : path to multichannel .wav file (6 channels)
    % d         : distance from center to each mic (hexagon size) in meters (e.g., 0.0475)
    
    % Constants
    c = 343;  % Speed of sound (m/s)
    [audio_data, Fs] = audioread(file_path);
    %disp(audio_data)
    audio_data = audio_data(:, 1:6);
    [num_samples, num_mics] = size(audio_data);
    
    if num_mics ~= 6
        error('Expected 6-channel audio (hexagonal array)');
    end

    % --- Mic Geometry (Hexagon around origin) ---
    angles_deg = 0:60:300;
    mic_positions = d * [cosd(angles_deg); sind(angles_deg)];  % 2 x 6

    % --- Preprocess Signal (use a snapshot) ---
    snapshot = audio_data(1:2048, :);
    snapshot = snapshot - mean(snapshot);  % zero-mean

    % --- Spatial Covariance Matrix ---
    R = snapshot' * snapshot / size(snapshot,1);

    % --- Eigen-decomposition ---
    [E, D] = eig(R);
    [~, idx] = sort(diag(D), 'descend');
    E = E(:, idx);
    
    % Assume 1 source → noise subspace has 5 vectors
    En = E(:, 2:end);  % Noise subspace (6-1)

    % --- MUSIC Spectrum over 360 degrees ---
    theta_scan = 0:1:359;
    Pmusic = zeros(size(theta_scan));

    % Choose a nominal frequency (e.g., 1 kHz) for narrowband approximation
    f0 = 1000;  % Hz
    lambda = c / f0;
    k = 2 * pi / lambda;

    for i = 1:length(theta_scan)
        theta = theta_scan(i) * pi / 180;
        doa_vec = [cos(theta); sin(theta)];  % unit vector direction
        delays = (mic_positions' * doa_vec) / c;  % delay per mic in seconds
        steering_vec = exp(-1j * 2 * pi * f0 * delays);  % 6 x 1

        Pmusic(i) = 1 / (steering_vec' * (En * En') * steering_vec);
    end

    % Normalize and Plot
    Pmusic_dB = 10*log10(abs(Pmusic) / max(abs(Pmusic)));
    polarplot(deg2rad(theta_scan), Pmusic_dB);
    title('MUSIC AoA Spectrum (Hexagonal Array)');
    rlim([-20 0]);  % dynamic range in dB
    
    % Estimate AoA
    [~, max_idx] = max(Pmusic);
    estimated_angle = theta_scan(max_idx);
    fprintf('Estimated AoA: %d degrees\n', estimated_angle);
end
