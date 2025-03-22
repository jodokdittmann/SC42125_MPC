function trim = positionMPC(Zeta, x_ss, p)

R = 10;  
n = 50;
l = 30;

X_ref = [linspace(0, l, n), R*sin(linspace(0, pi, n)) + l, l - linspace(0, l, n), - R*sin(linspace(0, pi, n))];
Y_ref = [zeros(1, n) + R, R*cos(linspace(0, pi, n)), zeros(1, n) - R, - R*cos(linspace(0, pi, n))];

progress = zeros(1, 15);
objective = zeros(1, 15);
zeta = cell(1, 15);

h = 1.5;

zeta0 = Zeta;

for k = 1:p.n_equi
        [t, Zeta] = ode45(@(t, zeta0) position(zeta0, x_ss{k}), [0, h], zeta0);
        zeta{k} = Zeta(end, :)';

        distX = X_ref - zeta{k}(2)*ones(1, 200);
        distY = Y_ref - zeta{k}(3)*ones(1, 200);
        dist = distX.^2 + distY.^2;
        [deviation, progress(k)] = min(dist);
        objective(k) = deviation - progress(k);
end

[~, trim] = min(objective);

end