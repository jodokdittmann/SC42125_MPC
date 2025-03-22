clear
init

tic

% Equilibria        
V_ss = 7;
beta_ss = linspace(-40*pi/180, 4*pi/180, 6); 
R_ss = linspace(7, 14, 4); 
x_ss = cell(length(beta_ss), length(R_ss)); 
u_ss = cell(length(beta_ss), length(R_ss)); 
u_ss{1, 1} = [-40.7*pi/180; 20.44; 58.33]; 

% Linearisation        
syms x [p.nx 1]
syms u [p.nu 1]
dxdt = dynamics(x, u, p);
dfdx = matlabFunction(jacobian(dxdt, x), 'Vars', {x, u});
dfdu = matlabFunction(jacobian(dxdt, u), 'Vars', {x, u});

% Initialisation
A = cell(length(beta_ss), length(R_ss)); 
B = cell(length(beta_ss), length(R_ss)); 

for i = 1:length(beta_ss)
    for j = 1:length(R_ss)
        
        % Equilibria
        x_ss{i, j} = [V_ss; beta_ss(i); V_ss/R_ss(j)];
        u_ss{i, j} = equilibria(x_ss{i, j}, u_ss{max(i-1,1), max(j-1,1)}, p);
       
        % Linearisation
        A{i, j} = dfdx(x_ss{i, j}, u_ss{i, j});
        B{i, j} = dfdu(x_ss{i, j}, u_ss{i, j});
        
        % Discretisation
        [A{i, j}, B{i, j}] = c2d(A{i, j}, B{i, j}, p.ts);

    end
end

toc

% Reference

x_r = x_ss{1, 1};
u_r = u_ss{1, 1};
A = A{1, 1};
B = B{1, 1};

% MPC

Q = 10*eye(p.nx);
R = eye(p.nu);

u_min = [-Inf; -Inf; -Inf];
u_max = [Inf; Inf; Inf];
x_min = [-Inf; -Inf; -Inf];
x_max = [Inf; Inf; Inf];

P = 10*Q;

x_0 = [1; 0; 0];

H = blkdiag(kron(speye(p.N), Q), P, kron(speye(p.N), R));
f = [repmat(-Q*x_r, p.N, 1); -P*x_r; repmat(-R*u_r, p.N, 1)];

G_eq = [kron(speye(p.N + 1), -speye(p.nx)) + kron(sparse(diag(ones(p.N, 1), -1)), A), kron([sparse(1, p.N); speye(p.N)], B)];
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
X(:, 1) = x_0; % + 0.6*rand(p.nx, 1);

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