function statedot = f(state, GM)
statedot = zeros(6,1);
statedot(1:3) = state(4:6);
statedot(4:6) = -GM*state(1:3)./norm(state(1:3))^3;
end