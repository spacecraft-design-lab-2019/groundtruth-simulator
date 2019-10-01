function [Ip_dry, R_dry, cg_dry] = moi_fxn()
%--------------------------------------------------------------------------
% SCRIPT: mass_properties_script.m
% AUTHOR: Adam Zufall & Hridu Jain
% DATE:   April 23, 2019

% NOTES:
% v2.0 for PS3
% Script to calculate mass, MoI, moment arm to each face
%--------------------------------------------------------------------------

clc; clear all; close all;

%% s/c is a rectangular prism with four cylinders and antennae
r_off  = [20,15,10]/100; %m, arbitrary mass offset
h_bus = 80/100;  %m     length of avionics body
w_bus = 20/100;  %m     width  of avionics body
r_tank = 10/100; %m     radius of fuel tanks
h_tank = 60/100; %m     height of fuel tanks
l_wire =  3;     %m     length of wire antennae
r_wire = 0.2/100;%m     radius of wire antennae
w_sc   =  0.9;   %m      width of the s/c
h_sc   =  0.8;   %m     height of the s/c
m_bus  = 56;     %kg       dry mass of s/c
m_full = 48;     %kg    filled mass of fuel tanks
m_empt =  1;     %kg     empty mass of fuel tanks
m_ball =  1;     %kg     point mass to cause offset
m_wire =  0.1;   %kg           mass of wire antennae

m_wet = m_bus + 4*m_full + 4*m_wire + m_ball; %wet mass of s/c
m_dry = m_bus + 4*m_empt + 4*m_wire + m_ball; %dry mass of s/c

cg_wet = r_off / m_wet;
cg_dry = r_off / m_dry;

%% define initial body coordinate system as geometric center
cg_bus = [0,0,0];
cg_tank = [   r_tank+w_bus, 0, 0;
           0,    r_tank+w_bus, 0;
           -1*r_tank-w_bus, 0, 0;
           0, -1*r_tank-w_bus, 0];
cg_wire = [   l_wire/2+w_bus, 0, 0;
           0,    l_wire/2+w_bus, 0;
           -1*l_wire/2-w_bus, 0, 0;
           0, -1*l_wire/2-w_bus, 0];

%% calculate principle moments of inertia       
I_bus = (m_bus/12)* [w_bus^2 + h_bus^2, 0, 0;
                     0, w_bus^2 + h_bus^2, 0;
                     0, 0, 2*w_bus^2];
I_tank = 4*[(1/12)*(3*r_tank^2+h_tank^2) + 0.5*(r_tank+w_bus/2)^2, 0 ,0;
            0, (1/12)*(3*r_tank^2+h_tank^2) + 0.5*(r_tank+w_bus/2)^2, 0;
            0, 0, (1/2)*(r_tank^2)+(r_tank+w_bus/2)^2];
I_full = m_full*I_tank;
I_empt = m_empt*I_tank;
I_wire = m_wire*[2*((1/3)*l_wire^2+(w_sc/2)^2),0,0;
                 0,2*((1/3)*l_wire^2+(w_sc/2)^2),0;
                 0,0,4*((1/3)*l_wire^2+(w_sc/2)^2)];
I_ball = m_ball*[r_off(2)^2+r_off(3)^2,-1*r_off(1)*r_off(2),-1*r_off(1)*r_off(3);
                 -1*r_off(1)*r_off(2),r_off(1)^2+r_off(3)^2,-1*r_off(2)*r_off(3);
                 -1*r_off(1)*r_off(3),-1*r_off(2)*r_off(3),r_off(1)^2+r_off(2)^2];
I_wet = I_bus + I_full          + I_ball;
I_dry = I_bus + I_empt + I_wire + I_ball;

%% calculate rotation matrix between principal and body axis
[~,Ip_wet] = eig(I_wet);
[~,Ip_dry] = eig(I_dry);
R_wet = I_wet*inv(Ip_wet);
R_dry = I_dry*inv(Ip_dry);
%% calculate moment arm to each face
sc_w = 0.6; sc_h = 0.8;
face_centers = [   sc_w/2,0,0; %+x
                -1*sc_w/2,0,0; %-x
                0,   sc_w/2,0; %+y
                0,-1*sc_w/2,0; %-y
                0,0,   sc_h/2; %+z
                0,0,-1*sc_h/2];%-z
face_norm = [ 1, 0, 0;         %+x
             -1, 0, 0;         %-x  
              0, 1, 0;         %+y
              0,-1, 0;         %-y 
              0, 0, 1;         %+z
              0, 0,-1];        %-z
mom_arm_wet = zeros(6,3);
mom_arm_dry = zeros(6,3);
 for i=1:6
     rad = face_centers(i,:) - cg_wet;
     mom_arm_wet(i,:) = cross(rad, face_norm(i,:)); %check sign when we define solar radiation pressure!
     rad = face_centers(i,:) - cg_dry;
     mom_arm_dry(i,:) = cross(rad, face_norm(i,:)); %check sign when we define solar radiation pressure!
 end
     
end