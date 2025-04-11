close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343;
d = 4.5e-2;
res = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mic_dis, azi_2d, ele_2d] = system_setup(d, res);
load('outputAngle.mat');
load('outputLocation.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run AoA on three signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path = './test_exam/';
%signal_num = '20';

path = './train/';
type = 'train-';Li
signal_num = '01';


angles1 = [];
angles2 = [];
for array_index=1:3
    filename = ['A0' num2str(array_index) '_X' signal_num '.wav'];
    filepath = [path filename];
    [y,Fs] = audioread(filepath);
    [calc_angle1, calc_angle2] = gcc_phat_aoa(y, Fs, c, res, azi_2d, ele_2d, mic_dis, 1, ['TEST-' filename]);
    angles1 = [angles1 calc_angle1];
    angles2 = [angles2 calc_angle2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(outputAngle(str2num(signal_num),:));
%disp(angles1);
disp(angles2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and plot location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mic_arrays=[2.38 4.90; 1.27 3.38; 2.93 1.30];
[location] = findLocation(angles2,mic_arrays, d, 1, [type signal_num]);
%disp(outputLocation(str2num(signal_num),:));
disp(location');


