clear
init

V_ss = 7;
R_ss = 7;
beta_ss = -51*pi/180;
dpsidt_ss = V_ss/R_ss;

% Equilibria

zeta_0 = [V_ss, beta_ss, 0, dpsidt_ss, 0, -40.7*pi/180, 0, 20.44, 58.33];
[delta_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p, zeta_0);

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
tf = 10;

[A_d, B_d] = c2d(A_c, B_c, ts);

% Reference

u_r = [delta_ss; omega_F_ss; omega_R_ss];
zeta_r = [V_ss; beta_ss; dpsidt_ss];

% MPC

l = 3;

Q = 2*eye(l);
R = eye(m);

N = 15;

u_min = [-Inf; -Inf; -Inf];
u_max = [Inf; Inf; Inf];
zeta_min = [-Inf; -Inf; -Inf];
zeta_max = [Inf; Inf; Inf];

P = 10*Q;

zeta_0 = zeta_r + 2*rand(3, 1);
   
H = blkdiag( kron(speye(N), Q), P, kron(speye(N), R) );
f = [repmat(-Q*zeta_r, N, 1); -P*zeta_r; repmat(-R*u_r, N, 1)];
   
G_eq = [kron(speye(N+1), -speye(l)) + kron(sparse(diag(ones(N, 1), -1)), A_d), kron([sparse(1, N); speye(N)], B_d)];
g_eq = [-zeta_0; zeros(N*l, 1)];

u_eq = g_eq;
    
G_ineq = speye((N+1)*l + N*m);  
g_ineq = [repmat(zeta_min, N+1, 1); repmat(u_min, N, 1)];

u_ineq = [repmat(zeta_max, N+1, 1); repmat(u_max, N, 1)];
    
G_c = [G_eq; G_ineq];
g_c = [g_eq; g_ineq];

u_c = [u_eq; u_ineq];

prob = osqp;

prob.setup(H, f, G_c, g_c, u_c, 'verbose', 0);
 
% Initialisation

U = nan(m, tf/ts + 1); 
Zeta = nan(n, tf/ts + 1); 
Zeta(:, 1) = [zeta_0(1); zeta_0(2); rand; zeta_0(3); rand; rand];

% Simulation

for k = 1:tf/ts
    sol = prob.solve();
    U(:, k) = sol.x((N + 1)*l + 1:(N + 1)*l + m);
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(:, k), p), [(k - 1)*ts k*ts], Zeta(:, k));
    Zeta(:, k + 1) = zeta(end, :);
    g_c(1:l) = -zeta(end, [1, 2, 4]);
    u_c(1:l) = -zeta(end, [1, 2, 4]);
    prob.update('l', g_c, 'u', u_c);
end

% Figures

figures(Zeta, U, p, tf, ts)  