clear
init

V_ss = 7;
R_ss = 7;
beta_ss = -10*pi/180;
dpsidt_ss = V_ss/R_ss;

[delta_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, R_ss, beta_ss, p);

syms zeta [6 1]
syms u [3 1]

dzetadt = dynamics(zeta, u, p);

A_c = jacobian(dzetadt([1, 2, 4]), zeta([1, 2, 4]));
B_c = jacobian(dzetadt([1, 2, 4]), u);

A_c = double(subs(A_c, {zeta1, zeta2, zeta4, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, omega_F_ss, omega_R_ss}));
B_c = double(subs(B_c, {zeta1, zeta2, zeta4, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, omega_F_ss, omega_R_ss}));

h = 1e-3;
tf = 50;

[A_d, B_d] = c2d(A_c, B_c, h);

Q = eye(3);
R = eye(3);

[K, S, P] = dlqr(A_d, B_d, Q, R);

m = 3;
n = 6;

zeta0 = zeros(n, 1);  
zeta0(1) = 7;

U = nan(m, tf/h + 1); 
Zeta = nan(n, tf/h + 1); 
Zeta(:, 1) = zeta0;

for k = 1:tf/h
    U(:, k) = [delta_ss; omega_F_ss; omega_R_ss] + K*(Zeta([1, 2, 4], k) - [V_ss; beta_ss; dpsidt_ss]);
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(:, k), p), [(k - 1)*h k*h], Zeta(:, k));
    Zeta(:, k + 1) = zeta(end, :);  
end

figure(1)
plot(Zeta(5, :), Zeta(6, :))

figure(3)
plot(Zeta(4, :))

figure(4)
plot(atan(Zeta(1, :)./Zeta(2, :)))

figures(Zeta(5, :)', Zeta(6, :)', Zeta(3, :)', p, tf, h)  