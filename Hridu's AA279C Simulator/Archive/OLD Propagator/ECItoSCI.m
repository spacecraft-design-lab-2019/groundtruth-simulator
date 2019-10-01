function state_SCI = ECItoSCI(state_ECI, inc)
rotation = [1,  0,       0;
      0,  cos(inc), sin(inc);
      0, -sin(inc), cos(inc)];
rSCI = rotation*state_ECI(1:3);
vSCI = rotation*state_ECI(4:6);
state_SCI = [rSCI; vSCI];
end