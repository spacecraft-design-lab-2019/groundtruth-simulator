%code for showing whether there is an eclipse between satellite, earth, and
%sun. If earth is between sun and satellite, there is an eclipse.
%Purpose of this code is to generate conditions for pytests

%All distances in kilometers
close all

earth = [0,0,0];
Re = 6371
r_sat = [5900,-6500,-6500]
r_Earth2Sun = [1e7,  1.50147817e8,  1.4e8]

[x,y,z] = sphere(100);
figure
surf(Re*x,Re*y,Re*z)
line = [r_Earth2Sun;r_sat];
hold on
plot3(r_sat(1),r_sat(2),r_sat(3),'o')
plot3(line(:,1),line(:,2),line(:,3))
limit = 9000;
axis([-limit,limit,-limit,limit,-limit,limit])