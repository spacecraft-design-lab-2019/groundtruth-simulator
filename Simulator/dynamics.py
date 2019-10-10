# -*- coding: utf-8 -*-
import numpy as np


#-------------------------Torques--------------------------------




#-------------------------Forces---------------------------------

def accelPointMass(r_sat, r_body, GM):
    #--------currently Earth gravity only assuming r_sat is in ECI
    #----------will eventually change to any body given r_body
    accel = -GM * r_sat / (np.linalg.norm(r_sat)**3)
    return accel


def accelHarmonic(state, R, n, m):

    x = 0
    y = 0
    z = 0
    V = np.zeros(50)
    W = np.zeros(50)

    r_vec = state[1:3]
    x_pos = r_vec[1]
    y_pos = r_vec[2]
    z_pos = r_vec[3]
    r = np.linalg.norm(r_vec)

    for j in range(0,m):
        for i in range(j,n):
            if i == 0:
                V[0,0] = R/r
                W[0,0] = 0
            elif i == j:
                V[i+1,j+1] = (2*j - 1)*((x_pos*R/r^2)*V[j,j] - (y_pos*R/r^2)*W[j,j])
                W[i+1,j+1] = (2*j - 1)*((x_pos*R/r^2)*W[j,j] - (y_pos*R/r^2)*V[j,j])
            elif i == j + 1:
                V[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V[i,j+1]
                W[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W[i,j+1]

            else:
                V[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V[i,j+1] - ((i+j-1)/(i-j))*(R^2/r^2)*V[i-1,j+1]
                W[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W[i,j+1] - ((i+j-1)/(i-j))*(R^2/r^2)*W[i-1,j+1]


    for j in range(0,m):
        for i in range(j,n):
            if j == 0:
                x = x + (GM/R^2)*(-C[i+1,1]*V[i+2,2]);
                y = y + (GM/R^2)*(-C[i+1,1]*W[i+2,2]);
            else:
                x = x + (GM/R^2)*0.5*((-C[i+1,j+1]*V[i+2,j+2] - S[i+1,j+1]*W[i+2,j+2]) + (np.math.factorial(i-j+2)/np.math.factorial(i-j))*(C[i+1,j+1]*V[i+2,j] + S[i+1,j+1]*W[i+2,j]))
                y = y + (GM/R^2)*0.5*((-C[i+1,j+1]*W[i+2,j+2] - S[i+1,j+1]*V[i+2,j+2]) + (np.math.factorial(i-j+2)/np.math.factorial(i-j))*(-C[i+1,j+1]*W[i+2,j] + S[i+1,j+1]*V[i+2,j]))

            z = z + (GM/R^2)*((i-j+1)*(-C[i+1,j+1]*V[i+2,j+1] - S[i+1,j+1]*W[i+2,j+1]))

    return np.array([x,y,z])

"""
function f_sphergrav = accelHarmonic2(state, GM, R, n, m)

X = load('constants.mat');
C = X.C;
S = X.S;

x = 0;
y = 0;
z = 0;
V = zeros(50);
W = zeros(50);

r_vec = state(1:3);
x_pos = r_vec(1);
y_pos = r_vec(2);
z_pos = r_vec(3);
r = norm(r_vec);

for j = 0:m
    for i = j:n
        if i == 0
            V(1,1) = R/r;
            W(1,1) = 0;

        elseif i == j
            V(i+1,j+1) = (2*j - 1)*((x_pos*R/r^2)*V(j,j) - (y_pos*R/r^2)*W(j,j));
            W(i+1,j+1) = (2*j - 1)*((x_pos*R/r^2)*W(j,j) - (y_pos*R/r^2)*V(j,j));

        elseif i == j + 1
            V(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V(i,j+1);
            W(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W(i,j+1);

        else
            V(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V(i,j+1) - ((i+j-1)/(i-j))*(R^2/r^2)*V(i-1,j+1);
            W(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W(i,j+1) - ((i+j-1)/(i-j))*(R^2/r^2)*W(i-1,j+1);
        end
    end
end

for j = 0:m
    for i = j:n
        if j == 0
            x = x + (GM/R^2)*(-C(i+1,1)*V(i+2,2));
            y = y + (GM/R^2)*(-C(i+1,1)*W(i+2,2));
        else
            x = x + (GM/R^2)*0.5*((-C(i+1,j+1)*V(i+2,j+2) - S(i+1,j+1)*W(i+2,j+2)) +...
                (factorial(i-j+2)/factorial(i-j))*(C(i+1,j+1)*V(i+2,j) + S(i+1,j+1)*W(i+2,j)));
            y = y + (GM/R^2)*0.5*((-C(i+1,j+1)*W(i+2,j+2) - S(i+1,j+1)*V(i+2,j+2)) +...
                (factorial(i-j+2)/factorial(i-j))*(-C(i+1,j+1)*W(i+2,j) + S(i+1,j+1)*V(i+2,j)));
        end

        z = z + (GM/R^2)*((i-j+1)*(-C(i+1,j+1)*V(i+2,j+1) - S(i+1,j+1)*W(i+2,j+1)));
    end
end

f_sphergrav = [x; y; z];
end

"""
