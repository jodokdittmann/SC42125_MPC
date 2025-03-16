clear
init

V_ss = 7;
R_ss = 7;
beta_ss = -51*pi/180;
dpsidt_ss = V_ss/R_ss;

% Equilibria

x0 = [V_ss, beta_ss, 0, dpsidt_ss, 0, -40.7*pi/180, 0, 20.44, 58.33];
[delta_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p, x0);

m = 3;
n = 6;

% Linearisation

syms zeta [n 1]
syms u [m 1]

dzetadt = dynamics(zeta, u, p);

A_c = jacobian(dzetadt([1, 2, 4]), zeta([1, 2, 4]));
B_c = jacobian(dzetadt([1, 2, 4]), u);

A_c = double(subs(A_c, {zeta1, zeta2, zeta4, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, omega_F_ss, omega_R_ss}));
B_c = double(subs(B_c, {zeta1, zeta2, zeta4, u1, u2, u3}, {V_ss, beta_ss, dpsidt_ss, delta_ss, omega_F_ss, omega_R_ss}));

% Discretisation

ts = 1e-3;
tf = 30;

[A_d, B_d] = c2d(A_c, B_c, ts);

% LQR

l = 3;

Q = eye(l);
R = eye(m);

K = dlqr(A_d, B_d, Q, R);

% MPC

N = 10;

% Initialisation

U = nan(m, tf/ts + 1); 
Zeta = nan(n, tf/ts + 1); 
Zeta(:, 1) = [V_ss; beta_ss; 0; dpsidt_ss; 0; 0] + rand(n, 1);

% Reference

u_r = [delta_ss; omega_F_ss; omega_R_ss];
zeta_r = [V_ss; beta_ss; dpsidt_ss];

% Simulation

for k = 1:tf/ts
    % U(:, k) = u_r - K*(Zeta([1, 2, 4], k) - zeta_r);
    U(:, k) = rhc(A_d, B_d, Q, R, zeta_r, u_r, Zeta([1, 2, 4], k));
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(:, k), p), [(k - 1)*ts k*ts], Zeta(:, k));
    Zeta(:, k + 1) = zeta(end, :);  
end

% Figures

figure(1)
plot(Zeta(5, :), Zeta(6, :))

figure(2)
plot(Zeta(4, :))

figure(3)
plot(U(1, :))

figures(Zeta(5, :)', Zeta(6, :)', Zeta(3, :)', p, tf, ts)  