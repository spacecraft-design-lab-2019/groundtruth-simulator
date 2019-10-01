%test R_eci2mci
d_moon  = 384000;  %km
r_earth = 6371;    %km
r_moon  = 1737.1;  %km
r_moon_eci = d_moon * [0; cosd(23.4+5.14); sind(23.4+5.14)];
%x_eci = d_moon * [0; 1; 0];
%x_eci = r_moon_eci / 2;
x_eci = 

obliq_e = 23.4;
obliq_m = 6.68;
ang = (pi/180) * (obliq_e - obliq_m);
R_eci2mci = [1,0,0;
             0,cos(ang),-1*sin(ang);
             0,sin(ang),   cos(ang)];
x_mci = R_eci2mci * (x_eci - r_moon_eci);


n=720;
angle = linspace(0, 2*pi, n); 
earth_pts = r_earth * [cos(angle); sin(angle)];
moon_pts  = r_moon_eci(2:3) + r_moon * [cos(angle); sin(angle)];
figure, grid on, hold on
plot(earth_pts(1,:), earth_pts(2,:))
plot( moon_pts(1,:),  moon_pts(2,:))
plot(x_eci(2), x_eci(3), 'ro')
xlabel('ECI_Y')
ylabel('ECI_Z')
axis('equal')