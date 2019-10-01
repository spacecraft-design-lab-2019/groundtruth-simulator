clc
clear all
close all

X = load('constants.mat');

figure

dist = 1.9658e+08;
x = linspace(8.42, 250, 100);

fsrp = (4.57e-6*1*0.02^2/200).*ones(length(x), 1);



for i = 1:length(x)
    fmain(i) = norm(accelPointMass([dist + x(i);0;0], [dist;0;0], X.GM_ast));
    fsun(i) = X.GM_Sun/((dist + x(i))^2);
    
end
semilogy(x, fmain)
hold on
semilogy(x, fsun)
semilogy(x, fsrp)


planeta = [149598023, 227939186, 778298361, 1429394133];
for planetid = 1:4
    for i = 1:length(x)
        f2nd(i) = norm(accelPointMass([dist + x(i);0;0], [planeta(planetid);0;0], X.GM_Planets(planetid)));
    end
    semilogy(x, f2nd)
end



xlabel('Distance from Asteroid Center [km]')
ylabel('Accleration [km/s^2]')
legend('GM', 'Sun', 'SRP', 'Earth', 'Mars', 'Jupiter', 'Saturn')
hold off
grid on
title('Accelerations about 433 Eros Asteroid')

