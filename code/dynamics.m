function dxdt = dynamics(x, u, p)

	V = x(1); beta = x(2); dpsidt = x(3); psi = x(4);

    delta = u(1); omega_F = u(2); omega_R = u(3);

    V_Fx = V*cos(beta - delta) + dpsidt*p.l_F*sin(delta);
    V_Rx = V*cos(beta);

    s_Fx = V_Fx/(p.R_w*omega_F) - 1;
    s_Rx = V_Rx/(p.R_w*omega_R) - 1;

    alpha_F = atan((V*sin(beta) + dpsidt*p.l_F)/(V*cos(beta))) - delta;
    alpha_R = atan((V*sin(beta) - dpsidt*p.l_R)/(V*cos(beta)));

    s_Fy = (1 + s_Fx)*tan(alpha_F);
    s_Ry = (1 + s_Rx)*tan(alpha_R);

    [mu_Fx, mu_Fy] = pacejka(s_Fx, s_Fy, p);
    [mu_Rx, mu_Ry] = pacejka(s_Rx, s_Ry, p);

    f_Fz = (p.l_R - p.h*mu_Rx)*p.m*p.g/(p.l_F + p.l_R + p.h*(mu_Fx*cos(delta) - mu_Fy*sin(delta) - mu_Rx));
    f_Rz = p.m*p.g - f_Fz;

    f_Fx = mu_Fx*f_Fz;
    f_Rx = mu_Rx*f_Rz;

    f_Fy = mu_Fy*f_Fz;
    f_Ry = mu_Ry*f_Rz;

    dVdt = (f_Fx*cos(delta - beta) - f_Fy*sin(delta - beta) + f_Rx*cos(beta) + f_Ry*sin(beta))/p.m;
    dbetadt = (f_Fx*sin(delta - beta) + f_Fy*cos(delta - beta) - f_Rx*sin(beta) + f_Ry*cos(beta))/(p.m*V) - dpsidt;
    dpsiddt = (p.l_F*(f_Fx*sin(delta) + f_Fy*cos(delta)) - p.l_R*f_Ry)/p.I_z;
    dXdt = V*cos(beta + psi);
    dYdt = V*sin(beta + psi);

    dxdt = [dVdt; dbetadt; dpsiddt; dpsidt; dXdt; dYdt];

end
