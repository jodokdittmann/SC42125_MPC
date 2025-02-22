function dzetadt = dynamicbicycle(t, zeta, u, p)

	X = zeta(1); Y = zeta(2); Psi = zeta(3); dxdt = zeta(4); dydt = zeta(5); dPsidt = zeta(6);
	deltaf = u(1); ax = u(2);

	%Fxf = pacejka(p.Bx, p.Cx, p.Dx, p.Ex, kappaf)
	%Fxr = pacejka(p.Bx, p.Cx, p.Dx, p.Ex, kappar)

	alphaf = atan((dydt + p.lf*dPsidt)/dxdt) - deltaf;
	alphar = atan((dydt + p.lr*dPsidt)/dxdt);

	Fyf = pacejka(p.By, p.Cy, p.Dy, p.Ey, alphaf);
	Fyr = pacejka(p.By, p.Cy, p.Dy, p.Ey, alphar);

	dXdt = dxdt*cos(Psi) - dydt*sin(Psi);
	dYdt = dxdt*sin(Psi) - dydt*cos(Psi);
	dxddt = dydt*dPsidt + ax;
	dyddt = - dxdt*dPsidt + 2/p.m*(Fyf*cos(deltaf) + Fyr);
	dPsiddt = 2/p.Iz*(p.lf*Fyf - p.lr*Fyr);

	dzetadt = [dXdt; dYdt; dPsidt; dxddt; dyddt; dPsiddt];

end
