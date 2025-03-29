function y = measurements(x, u, p)

    V = x(1); beta = x(2); dpsidt = x(3); X = x(5); Y = x(6);

    delta = u(1); omega_F = u(2); omega_R = u(3);

    [f_Fx, f_Rx, f_Fy, f_Ry] = tires(V, beta, dpsidt, delta, omega_F, omega_R, p);

    a_x = (f_Fx*cos(delta) - f_Fy*sin(delta) + f_Rx)/p.m + V*sin(beta)*dpsidt;

    a_y = (f_Fx*sin(delta) + f_Fy*cos(delta) + f_Ry)/p.m - V*cos(beta)*dpsidt;

    y = [a_x; a_y; dpsidt; X; Y];

end