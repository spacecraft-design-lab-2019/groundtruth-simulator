function roe = oe2roe(oe_sat, oe_ast)
da = (oe_sat(1) - oe_ast(1))/oe_ast(1);
u_sat = oe_sat(5) + nu2M(oe_sat(6), oe_sat(2));
u_ast = oe_ast(5) + nu2M(oe_ast(6), oe_ast(2));
dlambda = (u_sat - u_ast) + (oe_sat(4)- oe_ast(4))*cos(oe_ast(3));
dex = oe_sat(2)*cos(oe_sat(5)) - oe_ast(2)*cos(oe_ast(5));
dey = oe_sat(2)*sin(oe_sat(5)) - oe_ast(2)*sin(oe_ast(5));
dix = oe_sat(3) - oe_ast(3);
diy = (oe_sat(4) - oe_ast(4))*sin(oe_ast(3));

roe = zeros(6, 1);
roe(1) = da;
roe(2) = dlambda;
roe(3) = dex;
roe(4) = dey;
roe(5) = dix;
roe(6) = diy;
end

function M = nu2M(nu, e)
E = 2*atan(sqrt((1-e)/(1+e))*tan(nu/2));
M = E - e*sin(E);
end
