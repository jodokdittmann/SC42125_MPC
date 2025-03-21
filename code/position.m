function dzetadt = position(zeta, x)

    V = x(1);
    beta = x(2);
    psi = zeta(1);

    dpsidt = x(3);
    dXdt = V*cos(beta + psi);
    dYdt = V*sin(beta + psi);

    dzetadt = [dpsidt; dXdt; dYdt];

end
