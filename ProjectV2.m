d = 4.75;
h = sqrt(3)*d/2;
 xlim = [100 500];

shift = [-d/2 -sqrt(3)*d/2; d/2 -sqrt(3)*d/2; -d 1; d 1; -d/2 sqrt(3)*d/2; d/2 sqrt(3)*d/2];
mic1 = [238 490; 238 490; 238 490; 238 490; 238 490; 238 490;];
%mic2 = [238 490; 238 490; 238 490; 238 490; 238 490; 238 490;];
%mic3 = [238 490; 238 490; 238 490; 238 490; 238 490; 238 490;];
mic2 = [124 338; 124 338; 124 338; 124 338; 124 338; 124 338;];
mic3 = [293 130; 293 130; 293 130; 293 130; 293 130; 293 130;];

mic1 = mic1+shift;
mic2 = mic2+shift;
mic3 = mic3+shift;


line1point1 = [238 490];
%line2point1 = [238 490];
%line3point1 = [238 490]; 
line2point1 = [124,338];
line3point1 = [293,130];

scatter(mic1(:,1),mic1(:,2));
hold on
scatter(mic2(:,1),mic2(:,2));
scatter(mic3(:,1),mic3(:,2));
scatter(line1point1(1),line1point1(2));
scatter(line2point1(1),line2point1(2));
scatter(line3point1(1),line3point1(2));


alpha1=270;
alpha2=315;
alpha3=293;

%%%%% second points for lines
line1p2 = [line1point1(1)+10*cosd(alpha1) line1point1(2)+10*sind(alpha1)];
line2p2 = [line2point1(1)+10*cosd(alpha2) line2point1(2)+10*sind(alpha2)];
line3p2 = [line3point1(1)+10*cosd(alpha3) line3point1(2)+10*sind(alpha3)];

vector1=line1p2-line1point1;
vector2=line2p2-line2point1;
vector3=line3p2-line3point1;


t1=100;
t2=-100;
%two points  for t1 and t2
%---line 1---
line1r1 = line1point1 + t1*vector1;
line1r2 = line1point1 + t2*vector1;

%---line 2---
line2r1 = line2point1 + t1*vector2;
line2r2 = line2point1 + t2*vector2;

%---line 3---
line3r1 = line3point1 + t1*vector3;
line3r2 = line3point1 + t2*vector3;


hold on
line([line1r1(1),line1r2(1)], [line1r1(2),line1r2(2)])
hold on
line([line2r1(1),line2r2(1)], [line2r1(2),line2r2(2)])
line([line3r1(1),line3r2(1)], [line3r1(2),line3r2(2)])

hold off

