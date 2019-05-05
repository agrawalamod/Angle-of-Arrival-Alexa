clear all;
close all;

%Keep fixed
freq = 100e3;
c = 3e8;
lambda = c/freq;

%Knobs
d = 40*lambda;
dis = lambda/2;
theta = 60;

%Variables
N = 5;
ant_center = [0 0];

%Setting up system
ant_location=[];
for i=ceil(-N/2):1:floor(N/2)
    ant_location=[ant_location; [i*dis 0]];
end

dist_gt = [];
tx_loc = [d*cosd(theta) d*sind(theta)];
y = [];

for ant_index=1:1:N    
    dist_ant = sqrt((tx_loc(1)-ant_location(ant_index,1)).^2 + (tx_loc(2)-ant_location(ant_index,2)).^2);
    time_ant = dist_ant/c;
    dist_gt = [dist_gt; dist_ant];
    
    y = [y; exp(1i*2*pi*freq*time_ant)+rand()*0.5];
 
end


M = [];
for alpha = 0:1:180
    phi = dis*cosd(alpha)*2*pi/lambda;
    temp_M=[];
    for j=floor(N/2):-1:ceil(-N/2)
        temp_M=[temp_M; exp(1i*phi*j)];
    end
    M = horzcat(M, temp_M); 
end

result = [];
for k = 1:size(M,2)
    result = [result; real(dot(M(:,k),y))];
end

plot(0:1:180,result)
hold on;
title('Theta = 30 degrees and Alpha = 80, d = lambda/2');
xlabel('Degrees');
ylabel('Result of steering matrix')


%{

alpha = 150;
tx_loc =[d*cosd(theta) d*sind(theta)];
tx_loc2 =[d*cosd(alpha) d*sind(alpha)];
y = [];
for ant_index=N:-1:1
    disp(ant_index)
    dist_ant = sqrt((tx_loc(1)-ant_locs(ant_index)).^2 + tx_loc(2).^2);
    time_ant = dist_ant/c;
    
    dist_ant2 = sqrt((tx_loc2(1)-ant_locs(ant_index)).^2 + tx_loc2(2).^2);
    time_ant2 = dist_ant2/c;
    
    y = [y; exp(1i*2*pi*freq*time_ant) + exp(1i*2*pi*freq*time_ant2)];
    
end


%}















