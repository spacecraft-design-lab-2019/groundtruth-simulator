function Lt = get_momentum(Is,Iw,Ws,Ww,ax)
I1 = Iw(1);
I2 = Iw(2);
w1 = Ww(1);
w2 = Ww(2);
A = getA(ax);

Lw = [I1*w1;I2*w2];
Lt = Is*Ws + A*Lw;
end