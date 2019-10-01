function [ moon_pos ] = calc_moon_pos(JD_tt, mop_moon_eph_coeff, mop_moon_jd_lookup)
%calculate the position of the moon in the eci frame based on the lookup
%coefficients and table provided.
persistent coeff;
persistent fract_v;

if(JD_tt > mop_moon_jd_lookup(end))
    disp('Error, JD is past final moon table value.')
end

%index and fraction of JD
k = -1;
f = -1;
for(i=1:length(mop_moon_jd_lookup-1))
    if(JD_tt >= mop_moon_jd_lookup(i) && JD_tt <= mop_moon_jd_lookup(i+1))
        k = i;
        f = (JD_tt-mop_moon_jd_lookup(i))/(mop_moon_jd_lookup(i+1)-mop_moon_jd_lookup(i)); %f in range 0,1
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
if(isempty(coeff))
    coeff = zeros(order+1, order+1);
    coeff(1,1) = 1;
    coeff(2,2) = 1;
    for(idx = 3:order+1)
        coeff(idx,:) = 2*[0 coeff(idx-1,1:end-1)]-coeff(idx-2,:);
    end
    fract_v = zeros(order+1,1);
    fract_v(1) = 1;
end

y(1) = 1;
for(ind = 2:order+1)
    fract_v(ind) = fract_v(ind-1)*f;
    for(sub_ind = ind:-2:1)
        y(ind) = y(ind)+fract_v(sub_ind)*coeff(ind,sub_ind);
    end
end

coeff_x = mop_moon_eph_coeff(k,  1:13);
coeff_y = mop_moon_eph_coeff(k, 14:26);
coeff_z = mop_moon_eph_coeff(k, 27:39);

moon_pos(1) = 1000*dot(coeff_x, y);
moon_pos(2) = 1000*dot(coeff_y, y);
moon_pos(3) = 1000*dot(coeff_z, y);
end