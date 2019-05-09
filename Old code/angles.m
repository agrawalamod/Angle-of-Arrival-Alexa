close all;
clear all;
x=[1.47 4.34; 0.45 4.11; 0.94 2.14; 1.11 1.07; 2 0.86; 3.27 0.98; 3.17 1.99; 3.12 2.85; 2.22 1.57; 2.79 4.39];

X01=[1.47 4.34];
X02=[0.45 4.11];
X03=[0.94 2.14];
X04=[1.11 1.07];
X05=[2 0.86];
X06=[3.27 0.98];
X07=[3.17 1.99];
X08=[3.12 2.85];
X09=[2.22 1.57];
X10=[2.79 4.39];

A1=[2.38 4.90];
A2=[1.24 3.38];
A3=[2.93 1.30];

u = [1,0];
i=1;
outputAngle = [0 0 0];
while i<11
   v=[x(i,1),x(i,2)];
   
   vec1 = v-A1;
   vec2 = v-A2;
   vec3 = v-A3;
   
   angle1 = acosd(dot(vec1,u)/(norm(vec1)*norm(u)));
   angle2 = acosd(dot(vec2,u)/(norm(vec2)*norm(u)));
   angle3 = acosd(dot(vec3,u)/(norm(vec3)*norm(u)));
   
   if(vec1(2)<0)
    angle1 = 360-angle1;
   end
   
    if(vec2(2)<0)
    angle2 = 360-angle2;
    end
   if(vec3(2)<0)
    angle3 = 360-angle3;
   end
   
   outputAngle(i,1)=angle1;
   outputAngle(i,2)=angle2;
   outputAngle(i,3)=angle3;
   i=i+1;
end



% %v = X01-A1;
% v=[-1,-1];
% %angle:
% angle = acosd(dot(u,v)/(norm(u)*norm(v)));
% %find direction angle
% if(v(2)<0)
%    angle=360-angle
% end
% 
% 
% i=1;
% 
% plot([A1(1) X01(1)], [A1(2) X01(2)]);
% scatter(2.38,4.90)
% hold on
% plot([0,1], [0,0])
% plot([0,v(1)], [0,v(2)])