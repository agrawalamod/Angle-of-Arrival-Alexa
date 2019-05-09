close all;
clear all;

%% Contants
c = 343;
d = 4.5e-2;
res = 0.25;

%% Run once
[mic_dis, azi_2d, ele_2d] = system_setup(d, res);

%{
%% Training Data
path = './train/';
mic_arrays=[2.38 4.90; 1.27 3.38; 2.93 1.30];
train_result = [];
load('outputAngle.mat');
load('outputLocation.mat');

for train_signal = 1:10
    angles1 = [];
    angles2 = [];
    if(train_signal<10)
        sig_name = ['X0' num2str(train_signal)];
    else
        sig_name = ['X' num2str(train_signal)];
    end
    for array_index=1:3
        filename = ['A0' num2str(array_index) '_' sig_name '.wav'];
        filepath = [path filename];
        disp(filepath);
        [y,Fs] = audioread(filepath);
        
        [calc_angle1, calc_angle2] = gcc_phat_aoa(y, Fs, c, res, azi_2d, ele_2d, mic_dis);
        angles1 = [angles1 calc_angle1];
        angles2 = [angles2 calc_angle2];
    end    
    location = findLocation(angles2, mic_arrays, d, 1, sig_name)';
    angle_error = abs(angles2-outputAngle(train_signal,:));
    location_error = norm(location - outputLocation(train_signal,:));
    result = [train_signal angles2 angle_error location location_error];
    train_result = [train_result; result];
end

train_location_result = [train_result(:,1) train_result(:,8:10)];
train_angle_result = [train_result(:,1) train_result(:,2:7)];
dlmwrite('./results/train/location.csv', train_location_result, '\t');
dlmwrite('./results/train/angles.csv', train_angle_result, '\t');
%}


%% Testing Data 1
path = './test/';
mic_arrays=[2.38 4.90; 1.27 3.38; 2.93 1.30];
test_result = [];

for test_signal = 11:15
    angles1 = [];
    angles2 = [];
    if(test_signal<10)
        sig_name = ['X0' num2str(test_signal)];
    else
        sig_name = ['X' num2str(test_signal)];
    end
    for array_index=1:3
        filename = ['A0' num2str(array_index) '_' sig_name '.wav'];
        filepath = [path filename];
        disp(filepath);
        [y,Fs] = audioread(filepath);
        [calc_angle1, calc_angle2] = gcc_phat_aoa(y, Fs, c, res, azi_2d, ele_2d, mic_dis);
        angles1 = [angles1 calc_angle1];
        angles2 = [angles2 calc_angle2];
    end    
    location = findLocation(angles2, mic_arrays, d, 1, sig_name)';
    result = [test_signal angles2 location];
    test_result = [test_result; result];
end

test_location_result = [test_result(:,1) test_result(:,5:6)];
test_angle_result = [test_result(:,1) test_result(:,2:4)];

dlmwrite('./results/test/location.csv', test_location_result, '\t');
dlmwrite('./results/test/angles.csv', test_angle_result, '\t');

%{
%% Testing Data 3
path = './test_exam/';
mic_arrays=[2.38 4.90; 1.27 3.38; 2.93 1.30];
testex_result = [];

for test_signal = 16:20
    angles1 = [];
    angles2 = [];
    if(test_signal<10)
        sig_name = ['X0' num2str(test_signal)];
    else
        sig_name = ['X' num2str(test_signal)];
    end
    for array_index=1:3
        filename = ['A0' num2str(array_index) '_' sig_name '.wav'];
        filepath = [path filename];
        disp(filepath);
        [y,Fs] = audioread(filepath);
        [calc_angle1, calc_angle2] = gcc_phat_aoa(y, Fs, c, res, azi_2d, ele_2d, mic_dis);
        angles1 = [angles1 calc_angle1];
        angles2 = [angles2 calc_angle2];
    end    
    location = findLocation(angles2, mic_arrays, d, 0, sig_name)';
    result = [test_signal angles2 location];
    testex_result = [testex_result; result];
end

testex_location_result = [testex_result(:,1) testex_result(:,5:6)];
testex_angle_result = [testex_result(:,1) testex_result(:,2:4)];

dlmwrite('./results/test_exam/location.csv', testex_location_result, '\t');
dlmwrite('./results/test_exam/angles.csv', testex_angle_result, '\t');
%}

