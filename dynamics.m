function dxdt = dynamics(x, u, p)

	V = x(1); beta = x(2); dpsidt = x(3); psi = x(4);

    delta = u(1); omega_F = u(2); omega_R = u(3);

    [f_Fx, f_Rx, f_Fy, f_Ry] = tires(V, beta, dpsidt, delta, omega_F, omega_R, p);

    dVdt = (f_Fx*cos(delta - beta) - f_Fy*sin(delta - beta) + f_Rx*cos(beta) + f_Ry*sin(beta))/p.m;
    dbetadt = (f_Fx*sin(delta - beta) + f_Fy*cos(delta - beta) - f_Rx*sin(beta) + f_Ry*cos(beta))/(p.m*V) - dpsidt;
    dpsiddt = (p.l_F*(f_Fx*sin(delta) + f_Fy*cos(delta)) - p.l_R*f_Ry)/p.I_z;
    dXdt = V*cos(beta + psi);
    dYdt = V*sin(beta + psi);

    dxdt = [dVdt; dbetadt; dpsiddt; dpsidt; dXdt; dYdt];

end
