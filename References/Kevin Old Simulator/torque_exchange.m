function Ma = torque_exchange(I, w, wdot, A, Lw, Lwdot)
Ma = I*wdot + A*Lwdot + cross(w, I*w) + cross(w, A*Lw);
end


