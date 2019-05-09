
%mic2: 0
%mic1: dcos(t)
%mic6: H * sin(120-t) = root3 * d * sin(120-t)
%mic5: 2d sin(150-t)
%mic4: H * sin(t) = 2d sin60 * sin(t) 
%mic3: d sin(t - 30)

d = 4.75e-2;
freq = 44.1e3;
c = 343;
lambda = c/freq;
M = [];
for t = 0:360
    d_array = [d*cosd(t) 0 d*sind(t-30) 2*d*sind(60)*sind(t) 2*d*sind(150-t) sqrt(3)*d*sind(120-t))];
    phi = d_array*2*pi/lambda;
    sig = exp(1i*phi');
    M = horzcat(M, sig);
end
    
