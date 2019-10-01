%calculate inertial conditions (r_eci, v_eci) for equatorial orbit around
%the moon

%50 x 125 km?
R_eci2mci = [1,0,0;0,0.957722128749237,-0.287694845459613;0,0.287694845459613,0.957722128749237];
R_mci2eci = R_eci2mci';
moon_eci = [61893860.5957580;348960520.146146;133767644.928460];
v_moon_eci = [-1027.49637842178;95.1436650753021;130.763601660728];
%r_mci = [1737150; 0;0 ];
r_mci = [1790220; 0;0 ];
%v_mci = [0; 1680.4; 0];
v_mci = [0; 1745; 0];

r_eci = R_mci2eci * (r_mci);
r_eci = r_eci + moon_eci;

v_eci = R_mci2eci * v_mci;
v_eci = v_eci + v_moon_eci;

%want initial attitude to be with the +Z axis pointed normal to the orbit
%normal
q_i = rotm2quat(R_mci2eci);
e_i = rotm2eul( R_mci2eci);
w_i = (pi/180)*(360/60)*[0;0;10];