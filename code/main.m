clear
init
        
V_ss = 7;
beta_ss = linspace(-40*pi/180, 40*pi/180, p.n_equi); 
dpsidt_ss = linspace(1, -1, p.n_equi); 

[zeta_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p);

[A, B] = matrices(zeta_ss, u_ss, p);   

Delta = trajectories(zeta_ss, p);

% Initialisation

X_est = nan(p.n_x, p.tf/p.ts + 1); 
X_est(:, 1) = [0.1; 0; 0; 0; 0; 10]; 

P = eye(p.n_x);

Zeta = nan(p.n_zeta, p.tf/p.ts + 1); 
Zeta(:, 1) = [1; 0; 0]; 

Eta = nan(p.n_eta, p.tf/p.ts + 1);
Eta(:, 1) = [0; 0; 10];

Eta_meas = nan(p.n_eta, p.tf/p.ts + 1);
Eta_meas(:, 1) = Eta(:, 1) + normrnd(0, 0.1, [p.n_eta, 1]);

U = nan(p.n_u, p.tf/p.ts + 1); 

% Simulation

for k = 1:p.tf/p.ts
    trim = positionMPC(Eta_meas(:, k), Delta, p);
    U(:, k) = dynamicsMPC(A{trim}, B{trim}, zeta_ss{trim}, u_ss{trim}, X_est(1:p.n_zeta, k), p);
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(:, k), p), [(k - 1)*p.ts k*p.ts], Zeta(:, k));
    [t, eta] = ode45(@(t, eta) position(eta, zeta(t == t, :)), t, Eta(:, k));
    Zeta(:, k + 1) = zeta(end, :);
    Eta(:, k + 1) = eta(end, :);
    Eta_meas(:, k + 1) = eta(end, :) + normrnd(0, 0.1, [1, p.n_eta]);
    [X_est(:, k + 1), P] = EKF(X_est(:, k), P, U(:, k), Eta_meas(:, k), p);
end

% Figures

figures(Zeta, Eta, U, p)

% figures(X_est(1:p.n_zeta, :), X_est(p.n_zeta + 1:p.n_x, :), U, p)