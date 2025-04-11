function run_triangulation_pipeline()
    % === Constants ===
    c = 343;                  % speed of sound (m/s)
    d = 0.0475;               % mic spacing (hex array side)
    res = 1;                  % AoA search resolution
    azi_grid = 0:res:359;
    ele_grid = 0;
    [azi_2d, ele_2d] = meshgrid(azi_grid, ele_grid);

    % === Array locations in meters ===
    arrays = [2.38 4.90; 1.24 3.38; 2.93 1.30];  % A01, A02, A03

    % === Load ground truth locations ===
    load('outputLocation.mat');  % should load variable outputLocation [10 x 2]

    % === Storage ===
    est_angles = zeros(10, 3);       % signal x array (AoAs)
    est_locs = zeros(10, 2);         % signal x (x, y)
    errors = zeros(10, 1);           % localization errors

    for signal_idx = 1:10
        angles = zeros(1, 3);
        for array_idx = 1:3
            % Load file
            file_name = sprintf('./train/A%02d_X%02d.wav', array_idx, signal_idx);
            [raw_y, Fs] = audioread(file_name);
            fprintf('Processing %s\n', file_name);

            % Extract 6-mic signals
            fy = raw_y(:,1:6);

            % Estimate AoA using interpolated GCC-PHAT
            ref_mic = 1;
            refsig = fy(:, ref_mic);
            fy_tau2 = [];
            for idx = 1:6
                if idx ~= ref_mic
                    tau = interpolated_gccphat_stable(refsig, fy(:,idx), Fs);
                    fy_tau2 = [fy_tau2; tau * c];
                end
            end

            % Get mic positions (relative)
            ang = [120 60 0 300 240 180];
            mic_pos = d * [cosd(ang); sind(ang)]';
            mic_dis = [];
            for i = 2:6
                mic_dis = [mic_dis; mic_pos(i,:) - mic_pos(1,:)];
            end

            % Search AoA
            scores = zeros(size(azi_2d));
            for A=1:length(azi_grid)
                azi = azi_grid(A)*pi/180;
                ele = 0;
                vec = [cos(azi)*cos(ele) sin(azi)*cos(ele) sin(ele)];
                scores(1,A) = 1 / norm(mic_dis * vec' - fy_tau2);
            end

            [~, max_idx] = max(scores);
            angles(array_idx) = azi_grid(max_idx);
        end

        % Store AoAs
        est_angles(signal_idx,:) = angles;

        % === Triangulate location ===
        plotting = 1;
        sig_name = sprintf('X%02d', signal_idx);
        pointCoord = findLocation(angles, arrays, d, plotting, sig_name);
        est_locs(signal_idx,:) = pointCoord';

        % === Compute localization error ===
        gt = outputLocation(signal_idx, :);
        errors(signal_idx) = norm(gt - pointCoord');
    end

    % === Save Results ===
    save('./results/estimated_locations.mat', 'est_angles', 'est_locs', 'errors');

    % === Plot Estimated vs Ground Truth ===
    figure;
    scatter(outputLocation(:,1), outputLocation(:,2), 50, 'g', 'filled');
    hold on;
    scatter(est_locs(:,1), est_locs(:,2), 50, 'r', 'filled');
    legend('Ground Truth', 'Estimated');
    title('Localization Results');
    xlabel('X (m)'); ylabel('Y (m)'); grid on;

    % === Plot Error CDF ===
    figure;
    [f, x] = ecdf(errors);
    plot(x, f, 'LineWidth', 2);
    xlabel('Localization Error (m)'); ylabel('CDF');
    title('CDF of Localization Error'); grid on;
end