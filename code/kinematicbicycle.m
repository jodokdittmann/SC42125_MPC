function dzetadt = kinematicbicycle(t, zeta, u, p)

	X = zeta(1); Y = zeta(2); Psi = zeta(3); V = zeta(4);
	deltaf = u(1); a = u(2);

	Beta = atan(p.lr/(p.lf + p.lr)*atan(deltaf));

	dXdt = V*cos(Psi + Beta);
	dYdt = V*sin(Psi + Beta);
	dPsidt = V*sin(Beta)/p.lr;
	dVdt = a;

	dzetadt = [dXdt; dYdt; dPsidt; dVdt];

end
