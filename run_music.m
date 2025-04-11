close all;
clear all;

d = 0.0475
angles = 0:60:300;  % in degrees
mic_positions = d * [cosd(angles); sind(angles)];  % 2 x 6 matrix
disp([mic_positions])
music('./train/A01_X01.wav', d);