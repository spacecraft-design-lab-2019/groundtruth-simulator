function [ J ] = calc_I( CG )
%Calculate I, given a CG offset
%s/c is 1 meter cube, approximately 200 kg
m = 200;
Ixx = 17; %kg*m2
Iyy = 18;
Izz = 22;
%off diagonal terms are zero when center of mass is through geometric center
I = [Ixx,0,0; 0,Iyy,0; 0,0,Izz];
R = norm(CG);
R2 = R^2;
J = zeros(3,3);
for i=1:3
    for j=1:3
        if(i == j)
            d = 0;
        else
            d = 1;
        end
        J(i,j) = I(i,j) + m*(R2*d-CG(i)*CG(j));
    end
end

end

