function [ sun_eci ] = calc_sun_pos(JD_tt, moon_eci, sol_emb_eph_coeff, sol_sun_eph_coeff, sol_sun_jd_lookup)
%calculate the position of the sun in the eci frame based on the lookup
%coefficients and table provided.
persistent coeff2;
persistent fract_v2;

if(JD_tt > sol_sun_jd_lookup(end))
    disp('Error, JD is past final sun table value.')
end

%index and fraction of JD
k = -1;
f = -1;
for(i=1:length(sol_sun_jd_lookup-1))
    if(JD_tt >= sol_sun_jd_lookup(i) && JD_tt <= sol_sun_jd_lookup(i+1))
        k = i;
        f = (JD_tt-sol_sun_jd_lookup(i))/(sol_sun_jd_lookup(i+1)-sol_sun_jd_lookup(i)); %f in range 0,1
        break
    end
end
f = 2*f-1; %fraction scaled in valid range -1,1
if(k<0)    %past end of table, k never reset
    k = length(mop_moon_jd_lookup)-1;
    f = 1;
end
%cheby
order = 12;
y = zeros(order+1,1);
if(isempty(coeff2))
    coeff2 = zeros(order+1, order+1);
    coeff2(1,1) = 1;
    coeff2(2,2) = 1;
    for(idx = 3:order+1)
        coeff2(idx,:) = 2*[0 coeff2(idx-1,1:end-1)]-coeff2(idx-2,:);
    end
    fract_v2 = zeros(order+1,1);
    fract_v2(1) = 1;
end

y(1) = 1;
for(ind = 2:order+1)
    fract_v2(ind) = fract_v2(ind-1)*f;
    for(sub_ind = ind:-2:1)
        y(ind) = y(ind)+fract_v2(sub_ind)*coeff2(ind,sub_ind);
    end
end
y11 = y(1:11);

coeff_x1 = sol_emb_eph_coeff(k,  1:13);
coeff_y1 = sol_emb_eph_coeff(k, 14:26);
coeff_z1 = sol_emb_eph_coeff(k, 27:39);

coeff_x2 = sol_sun_eph_coeff(k,  1:11);
coeff_y2 = sol_sun_eph_coeff(k, 12:22);
coeff_z2 = sol_sun_eph_coeff(k, 23:33);

moon_pos = 0.012150585609624*moon_eci;

sun_eci(1) = 1000*dot(coeff_x1, y);
sun_eci(2) = 1000*dot(coeff_y1, y);
sun_eci(3) = 1000*dot(coeff_z1, y);

sun_eci = sun_eci-moon_pos;

sun_eci(1) = 1000*dot(coeff_x2, y11)-sun_eci(1);
sun_eci(2) = 1000*dot(coeff_y2, y11)-sun_eci(2);
sun_eci(3) = 1000*dot(coeff_z2, y11)-sun_eci(3);

end