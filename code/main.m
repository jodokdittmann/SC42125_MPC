clear
init
        
% Equilibria 

V_ss = 7;
beta_ss = linspace(-40*pi/180, 40*pi/180, p.n_equi); 
dpsidt_ss = linspace(1, -1, p.n_equi); 

[zeta_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p);

% Motion Primitives

[A, B] = matrices(zeta_ss, u_ss, p);   

Delta = trajectories(zeta_ss, p);

% Initialisation

U = nan(p.n_u, p.tf/p.ts + 1); 

X = nan(p.n_x, p.tf/p.ts + 1); 
X(:, 1) = [1; 0; 0; 0; 0; 10]; 

% EKF

X_est = nan(p.n_x, p.tf/p.ts + 1); 
X_est(:, 1) = [0.1; 0; 0; 0; 0; 10]; 

H = [1, zeros(1, 5); zeros(1, 2), 1, zeros(1, 3); zeros(1, 4), 1, 0; zeros(1, 5), 1];
Y = H*X;

P = eye(p.n_x);

% Simulation

for k = 1:p.tf/p.ts
    trim = pathplanning(X_est(p.n_zeta + 1:p.n_x, k), Delta, p);
    U(:, k) = MPC(A{trim}, B{trim}, zeta_ss{trim}, u_ss{trim}, X(1:p.n_zeta, k), p);
    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X(:, k));
    X(:, k + 1) = x(end, :);
    Y(:, k + 1) = H*X(:, k + 1) + [normrnd(0, 0.01); normrnd(0, 0.01); normrnd(0, 0.1); normrnd(0, 0.1)];
    [X_est(:, k + 1), P] = EKF(U(:, k), X_est(:, k), H, Y(:, k + 1), P, p);
end

% Figures

% figures(X(1:p.n_zeta, :), X(p.n_zeta + 1:p.n_x, :), U, p)

figures(X_est(1:p.n_zeta, :), X_est(p.n_zeta + 1:p.n_x, :), U, p)