%--------------------------------------------------------------------------
% SCRIPT: mass_properties_script.m
% AUTHOR: Adam Zufall
% DATE:   April 9, 2019

% NOTES:
% v1.0 for PS1
% Script to calculate mass, MoI, moment arm to each face
%--------------------------------------------------------------------------

% clc; clear all; close all;

%% s/c is a cylindrical octogon with 4 solar panels
l_body    =  24; %m      length   of main body
d_body    =   6; %m      diameter of main body
rho_body  = 180; %kg/m3  density  of main body
l_panel   =   6; %m      length   of solar panels
w_panel   =   4; %m      width    of solar panels
rho_panel = 1.4; %kg/m2  density  of solar panels
z_offset  =   8; %m      distance from face solar panels are mounted down
r_offset  =   2; %m      distance from face solar panels are mounted out

m_panel   = rho_panel * w_panel * l_panel; %kg  mass of one solar panel
sl        = d_body*tand(22.5);             %m   side length of octogonal face
a_body    = 2*(1+sqrt(2))*sl^2;            %m2         area of octogonal face
v_body    = a_body * l_body;               %m3  volume of main body
m_body    = rho_body * v_body;             %kg    mass of main body
m_total   = m_body + 4*m_panel;            %kg    mass of entire s/c 


%% define initial body coordinate system as centered on closest face
cg_body   = [0,0,l_body/2]; 
r_pan = d_body/2 + r_offset+ l_panel/2;
cg_panels = [   r_pan,    0, z_offset;
                0,    r_pan, z_offset;
             -1*r_pan,    0, z_offset;
                0, -1*r_pan, z_offset];
cg_total = m_body * cg_body + sum(m_panel*cg_panels);
cg_total = cg_total / m_total;

%% calculate principle moments of inertia
sns = sin(pi/8)*sin(pi/8);
r_body = d_body / (2*cosd(22.5));
Iz_bod = 0.5*m_body*r_body^2*(1-(2/3)*sns);    %kg*m2, through center of mass
Iz_rec = (1/12)*m_panel*(l_panel^2+w_panel^2); %kg*m2, through center of mass
Iz_pan = Iz_rec + m_panel*(r_pan)^2;           %kg*m2, parallel axis theorem
Izz = Iz_bod + 4*Iz_pan;                       %around long axis of s/c

Ix_pan = m_panel*(cg_total(3) - z_offset)^2;     %kg*m2, through center of rotation
Ix_rec = (1/12)*m_body*(l_body^2+d_body^2);    %kg*m2, through center of mass
Ix_bod = Ix_rec + m_body*(cg_body(3)-cg_total(3))^2; %parallel axis theorem
Ixx = Ix_bod + 4*Ix_pan;
Iyy = Ixx;

%% add arbitrary CG offset (to account for manufacturing / design wiggle room)
% cg_offset = [0.5, -0.1, 1.0];
% I = [Ixx,0,0; 0,Iyy,0; 0,0,Izz];
% R_CG = norm(cg_offset);
% R2 = R_CG^2;
% J = zeros(3,3);
% for i=1:3
%     for j=1:3
%         if(i == j)
%             d = 0;
%         else
%             d = 1;
%         end
%         J(i,j) = I(i,j) + m_total*(R2*d-cg_offset(i)*cg_offset(j));
%     end
% end
% cg_total = cg_total + cg_offset;

%% calculate normal vectors for each surface and moment arm to CG
c45 = cosd(45);
s45 = sind(45);
r_b = d_body/2;
r_b2 = r_b*c45;
face_centers = [     0,       0,        0;  %bod 1
                     0,       0, l_body  ;  %bod 2
                  r_b ,       0, l_body/2;  %oct 1
                  r_b2,    r_b2, l_body/2;  %oct 2
                     0,    r_b , l_body/2;  %oct 3
               -1*r_b2,    r_b2, l_body/2;  %oct 4
               -1*r_b ,       0, l_body/2;  %oct 5
               -1*r_b2, -1*r_b2, l_body/2;  %oct 6
                     0, -1*r_b , l_body/2;  %oct 7
                  r_b2, -1*r_b2, l_body/2;  %oct 8
                 r_pan,       0, z_offset;  %panel 1+
                 r_pan,       0, z_offset;  %panel 1-
                     0,   r_pan, z_offset;  %panel 2+
                     0,   r_pan, z_offset;  %panel 2-
              -1*r_pan,       0, z_offset;  %panel 3+
              -1*r_pan,       0, z_offset;  %panel 3-
                     0,-1*r_pan, z_offset;  %panel 4+
                     0,-1*r_pan, z_offset]; %panel 4-
face_norm = [  0,      0, -1;  %bod 1
               0,      0,  1;  %bod 2
               1,      0,  0;  %oct 1
             c45,    s45,  0;  %oct 2
               0,      1,  0;  %oct 3
          -1*c45,    s45,  0;  %oct 4
              -1,      0,  0;  %oct 5
          -1*c45, -1*s45,  0;  %oct 6
               0,     -1,  0;  %oct 7
             c45, -1*s45,  0;  %oct 8
               0,      0,  1;  %panel 1+
               0,      0, -1;  %panel 1-
               0,      0,  1;  %panel 2+
               0,      0, -1;  %panel 2-
               0,      0,  1;  %panel 3+
               0,      0, -1;  %panel 3-
               0,      0,  1;  %panel 4+
               0,      0, -1]; %panel 4-
 mom_arm = zeros(18,3);
 for i=1:18
     rad = face_centers(i,:) - cg_total;
     mom_arm(i,:) = cross(rad, face_norm(i,:)); %check sign when we define solar radiation pressure!
 end
                