clear all;
close all;
% given parameters
f = 100e3; % frequency
 
% other parameters
Rs = 1e4;
Ns = 1e3;
c = 3e8;
lambda = c/f;
 
% simulation
N_ant = 5;
a_axis = 0:N_ant-1;
%THETA = rand()*180;
THETA = 60;
r_o = lambda*20;
x_o = r_o*cos(THETA/180*pi);
y_o = r_o*sin(THETA/180*pi);
 
% time axis
t_axis = [0:(Ns-1)]/Rs;
 
%% Part a
DIS = [lambda, lambda/2, lambda/4];
t_est = [0:5:180];
t_spec = zeros(length(DIS),length(t_est));
 
for kD = 1:3
    Xant = [0:N_ant-1].'*DIS(kD);
    ant_tau = sqrt((x_o-Xant).^2+y_o^2)/c;
    ant_sig = exp(1j*2*pi*f*(t_axis-ant_tau));
    
    for kt = 1:length(t_est)
        ang_comp = exp(-1j*2*pi/lambda*DIS(kD)*cos(t_est(kt)/180*pi).*a_axis); ang_comp = ang_comp.';
    t_spec(kD,kt) = abs(sum(repmat(ang_comp,[1,Ns]).*ant_sig,[1,2]));
    end
 
end
 
%figure;
%plot(t_est,t_spec(1,:)/1000); hold on;
%plot(t_est,t_spec(2,:)/1000); hold on;
%plot(t_est,t_spec(3,:)/1000); hold on;
%legend('full','half','quarter');
 
%% Part b
t_spec = zeros(length(DIS),length(t_est));
for kD = 2
    Xant = [0:N_ant-1].'*DIS(kD);
    ant_tau = sqrt((x_o-Xant).^2+y_o^2)/c;
    ant_sig = exp(1j*2*pi*f*(t_axis-ant_tau))+rand()*0.5;
 
    for kt = 1:length(t_est)
        ang_comp = exp(-1j*2*pi/lambda*DIS(kD)*cos(t_est(kt)/180*pi).*a_axis); ang_comp = ang_comp.';
        t_spec(kD,kt) = abs(sum(repmat(ang_comp,[1,Ns]).*ant_sig,[1,2]));
    end
end
 
%figure;
%plot(t_est,t_spec(2,:)/1000); hold on;
%title('Theta = 60 degrees, d = lambda/2');
%xlabel('Degrees');
%ylabel('Result of steering matrix')

%}
 
%% Part c
% THETA = [rand()*180;rand()*180;rand()*180;rand()*180;];
THETA = linspace(20,160,4);
No = length(THETA);
x_o = r_o*cos(THETA/180*pi);
y_o = r_o*sin(THETA/180*pi);
 
t_spec = zeros(length(DIS),length(t_est));
for kD = 2
    Xant = [0:N_ant-1].'*DIS(kD);
    ant_sig = zeros(N_ant,Ns);
    for ko = 1:No
        ant_tau = sqrt((x_o(ko)-Xant).^2+y_o(ko)^2)/c;
        ant_sig = ant_sig+exp(1j*2*pi*f*(t_axis-ant_tau))+rand();
    end
 
    for kt = 1:length(t_est)
        ang_comp = exp(-1j*2*pi/lambda*DIS(kD)*cos(t_est(kt)/180*pi).*a_axis); ang_comp = ang_comp.';
        t_spec(kD,kt) = abs(sum(repmat(ang_comp,[1,Ns]).*ant_sig,[1,2]));
    end
end
 
figure;
plot(t_est,t_spec(2,:)/1000); hold on;
title('Theta = 60 degrees, d = lambda/2');
xlabel('Degrees');
ylabel('Result of steering matrix')
