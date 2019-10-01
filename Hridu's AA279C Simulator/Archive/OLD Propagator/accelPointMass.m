function f_grav = accelPointMass(state, state_body, GM_body)
r_sat = state(1:3);
r_body = state_body(1:3);
r_sat_body = r_body - r_sat;

f_grav = GM_body * (r_sat_body/(norm(r_sat_body)^3) - r_body/(norm(r_body)^3));
end