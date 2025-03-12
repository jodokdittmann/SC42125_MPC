clear
init

V_ss = 7;
R_ss = 7;
beta_ss = -51*pi/180;
dpsidt_ss = V_ss/R_ss;

x0 = [V_ss, beta_ss, 0, dpsidt_ss, 0, 0, -40.7*pi/180, 0, 20.44, 58.33];

[delta_ss, ddeltadt_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p, x0);

m = 3;
n = 7;

syms zeta [n 1]
syms u [m 1]

dzetadt = dynamics(zeta, u, p);

A_c = jacobian(dzetadt([1, 2, 4, 5]), zeta([1, 2, 4, 5]));
B_c = jacobian(dzetadt([1, 2, 4, 5]), u);

A_c = double(subs(A_c, {zeta1, zeta2, zeta4, zeta5, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, ddeltadt_ss, omega_F_ss, omega_R_ss}));
B_c = double(subs(B_c, {zeta1, zeta2, zeta4, zeta5, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, ddeltadt_ss, omega_F_ss, omega_R_ss}));

h = 1e-3;
tf = 30;

[A_d, B_d] = c2d(A_c, B_c, h);

Q = eye(4);
R = eye(3);

[K, S, P] = dlqr(A_d, B_d, Q, R);

zeta0 = zeros(n, 1);  
zeta0(1) = 7;

U = nan(m, tf/h + 1); 
Zeta = nan(n, tf/h + 1); 

Zeta(:, 1) = [V_ss; beta_ss; 0; dpsidt_ss; delta_ss; 0; 0] + rand(n, 1);

for k = 1:tf/h
    U(:, k) = [ddeltadt_ss; omega_F_ss; omega_R_ss] - K*(Zeta([1, 2, 4, 5], k) - [V_ss; beta_ss; dpsidt_ss; delta_ss]);
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(:, k), p), [(k - 1)*h k*h], Zeta(:, k));
    Zeta(:, k + 1) = zeta(end, :);  
end

figure(1)
plot(Zeta(6, :), Zeta(7, :))

figure(2)
plot(Zeta(4, :))

figure(3)
plot(U(1, :))

figures(Zeta(6, :)', Zeta(7, :)', Zeta(3, :)', p, tf, h)  