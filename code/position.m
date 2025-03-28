function detadt = position(eta, zeta)

    V = zeta(1);
    beta = zeta(2);
    psi = eta(1);

    dpsidt = zeta(3);
    dXdt = V*cos(beta + psi);
    dYdt = V*sin(beta + psi);

    detadt = [dpsidt; dXdt; dYdt];

end
