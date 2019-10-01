function f_srp = accelSRP(state, state_sun, state_2nd, r_shade, Csrp, Msrp, Asrp, Psrp)

f_srp = zeros(3, 1);

r_sat = state(1:3);
r_sun = state_sun(1:3);
r_2nd = state_2nd(1:3);

r_sat_sun = r_sat - r_sun;
r_sat_2nd = r_sat - r_2nd;

if r_shade == 0
    f_srp = Psrp*Csrp*Asrp/Msrp.*r_sat_sun./norm(r_sat_sun);
else
    beta = acosd(dot(r_sat_2nd, r_sat_sun)/(norm(r_sat_2nd)*norm(r_sat_sun)));
    if beta > 90
        if norm(r_sat_2nd)*sind(beta) > r_shade
            f_srp = Psrp*Csrp*Asrp/Msrp.*r_sat_sun./norm(r_sat_sun);
        end
    end
end
end