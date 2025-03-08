clear
init 

h = 1e-3;
tf = 10;

m = 3;
n = 6; 

zeta0 = zeros(n, 1);
zeta0(1) = 0.1;

U = nan(tf/h + 1, m);
Zeta = nan(tf/h + 1, n);
Zeta(1, :) = zeta0;

for k = 1:tf/h
	U(k, :) = [0.3 100 100];
	[t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(k, :), p), [(k - 1)*h k*h], Zeta(k, :));
	Zeta(k + 1, :) = zeta(end, :);
end

figure(1)
plot(Zeta(:,5),Zeta(:,6))

figure(2)
plot(Zeta(:,2), Zeta(:,4))


