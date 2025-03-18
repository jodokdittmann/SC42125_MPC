clear
init

V_ss = 7;
R_ss = 7;
beta_ss = -51*pi/180;
dpsidt_ss = V_ss/R_ss;

% Equilibria

x_ss = [V_ss; beta_ss; dpsidt_ss];
u_init = [-40.7*pi/180; 20.44; 58.33];

u_ss = equilibria(x_ss, u_init, p);

delta_ss = u_ss(1);
omega_F_ss = u_ss(2);
omega_R_ss = u_ss(3);

% Linearisation

syms x [p.nx 1]
syms u [p.nu 1]

dxdt = dynamics(x, u, p);

A_c = jacobian(dxdt, x);
B_c = jacobian(dxdt, u);

A_c = double(subs(A_c, [x; u], [x_ss; u_ss]));
B_c = double(subs(B_c, [x, u], [x_ss, u_ss]));

% Discretisation

[A_d, B_d] = c2d(A_c, B_c, p.ts);

% Reference

x_r = x_ss;
u_r = u_ss;

% MPC

Q = 2*eye(p.nx);
R = eye(p.nu);

u_min = [-Inf; -Inf; -Inf];
u_max = [Inf; Inf; Inf];
x_min = [-Inf; -Inf; -Inf];
x_max = [Inf; Inf; Inf];

P = 10*Q;

x_0 = x_r;

H = blkdiag(kron(speye(p.N), Q), P, kron(speye(p.N), R));
f = [repmat(-Q*x_r, p.N, 1); -P*x_r; repmat(-R*u_r, p.N, 1)];

G_eq = [kron(speye(p.N + 1), -speye(p.nx)) + kron(sparse(diag(ones(p.N, 1), -1)), A_d), kron([sparse(1, p.N); speye(p.N)], B_d)];
g_eq = [-x_0; zeros(p.N*p.nx, 1)];

u_eq = g_eq;

G_ineq = speye((p.N + 1)*p.nx + p.N*p.nu);  
g_ineq = [repmat(x_min, p.N + 1, 1); repmat(u_min, p.N, 1)];

u_ineq = [repmat(x_max, p.N + 1, 1); repmat(u_max, p.N, 1)];

G_c = [G_eq; G_ineq];
g_c = [g_eq; g_ineq];

u_c = [u_eq; u_ineq];

prob = osqp;

prob.setup(H, f, G_c, g_c, u_c, 'verbose', 0);

% Initialisation

X = nan(p.nx, p.tf/p.ts + 1); 
X(:, 1) = x_0 + 0.6*rand(p.nx, 1);

Zeta = nan(p.nzeta, p.tf/p.ts + 1);
Zeta(:, 1) = zeros(p.nzeta, 1);

U = nan(p.nu, p.tf/p.ts + 1); 

% Simulation

for k = 1:p.tf/p.ts
    sol = prob.solve();
    U(:, k) = sol.x((p.N + 1)*p.nx + 1:(p.N + 1)*p.nx + p.nu);
    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X(:, k));
    [t, zeta] = ode45(@(t, zeta) position(zeta, x(t == t, :)), t, Zeta(:, k));
    X(:, k + 1) = x(end, :);
    Zeta(:, k + 1) = zeta(end, :);
    g_c(1:p.nx) = -x(end, :);
    u_c(1:p.nu) = -x(end, :);
    prob.update('l', g_c, 'u', u_c);
end

% Figures

figures(Zeta, X, U, p)  