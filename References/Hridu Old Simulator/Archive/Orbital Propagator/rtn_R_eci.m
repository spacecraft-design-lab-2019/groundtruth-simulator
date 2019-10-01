function R = rtn_R_eci(state, GM)
% Provides the rotation matrix (not the translation) to go from ECI to RTN
% given the current state of the s/c.

oe = eci2oe(state, GM);
Omega = oe(4); i = oe(3); omega = oe(5); nu = oe(6);

if isnan(Omega), Omega = 0; end
if isnan(omega), omega = 0; end
if isnan(nu), nu = ang; end

ROmega = [cos(-Omega), sin(-Omega), 0;
         -sin(-Omega), cos(-Omega), 0;
          0,            0,            1];
Ri = [1,  0,          0;
      0,  cos(-i),   sin(-i);
      0,  -sin(-i),  cos(-i)];
Romega = [cos(-omega), sin(-omega), 0;
         -sin(-omega), cos(-omega), 0;
          0,            0,            1];
rotation = ROmega*Ri*Romega;

Rnu = [cos(nu), sin(nu), 0;
      -sin(nu), cos(nu), 0;
       0      , 0      , 1];

R = Rnu*rotation';

end
