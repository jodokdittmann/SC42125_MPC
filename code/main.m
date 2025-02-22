h = 1e-3;
tf = 5;

m = 2;
n = 6; 

zeta0 = zeros(n, 1);
zeta0(4, 1) = 0.01;

U = nan(tf/h + 1, m);
Zeta = nan(tf/h + 1, n);
Zeta(1, :) = zeta0;

for k = 1:tf/h
	U(k, :) = [1 2];
	[t, zeta] = ode45(@(t, zeta) dynamicbicycle(t, zeta, U(k, :), p), [(k - 1)*h k*h], Zeta(k, :));
	Zeta(k + 1, :) = zeta(end, :);
end

figure(1)
plot(Zeta(:,1),Zeta(:,2))



