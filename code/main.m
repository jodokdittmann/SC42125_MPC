clear
clc
init
        
%% Equilibria 

V_ss = 7;
beta_ss = linspace(-40*pi/180, 40*pi/180, p.n_equi); 
dpsidt_ss = linspace(1, -1, p.n_equi); 

[zeta_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p);

%% Motion Primitives

[A, B, C, D, P_MPC] = matrices(zeta_ss, u_ss, p);   

Delta_eta = trajectories(zeta_ss, p);

%% Simulation

% Choose a steady state trajectory (traj ∈ [1, p.n_equi]) and state (1)
% or output feedback (2).

traj = 1;
feedback = 2; 

%% Steady State

U = nan(p.n_u, p.tf/p.ts + 1);
U(:, 1) = u_ss{traj};

Zeta = nan(p.n_zeta, p.tf/p.ts + 1); 
Zeta(:, 1) = zeta_ss{traj};

Zeta_est = nan(p.n_zeta, p.tf/p.ts + 1); 
Zeta_est(:, 1) = zeta_ss{traj};

Z = nan(p.n_z, p.tf/p.ts + 1);

P_KF = eye(p.n_zeta);

u_min = [-pi/4; 0; 0];
u_max = [pi/4; 50; 50];
u_Delta = [Inf; Inf; Inf];

sys = LTISystem('A', A{traj}, 'B', B{traj});
sys.u.min = u_min;
sys.u.max = u_max;
opt.maxIterations = 4000;
invset = sys.invariantSet(opt);
figure('Name', 'Invariant Set', 'NumberTitle', 'off')
invset.plot()

for k = 1:p.tf/p.ts
    switch feedback
        case 1 
            U(:, k) = MPC(A{traj}, B{traj}, P_MPC{traj}, zeta_ss{traj}, u_ss{traj}, Zeta(:, k), U(:, max(k - 1, 1)), u_min, u_max, u_Delta, p);
            Zeta(:, k + 1) = zeta_ss{traj} + A{traj}*(Zeta(:, k) - zeta_ss{traj}) + B{traj}*(U(:, k) - u_ss{traj});
        case 2 
            Z(:, k) = C{traj}*Zeta(:, k) + D{traj}*U(:, max(k - 1, 1));   
            [Zeta_est(:, k), P_KF] = KF(A{traj}, B{traj}, C{traj}, D{traj}, u_ss{traj}, zeta_ss{traj}, U(:, max(k - 1, 1)), Zeta_est(:, max(k - 1, 1)), Z(:, k), P_KF, p);   
            U(:, k) = MPC(A{traj}, B{traj}, P_MPC{traj}, zeta_ss{traj}, u_ss{traj}, Zeta_est(:, k), U(:, max(k - 1, 1)), u_min, u_max, u_Delta, p);
            Zeta(:, k + 1) = zeta_ss{traj} + A{traj}*(Zeta(:, k) - zeta_ss{traj}) + B{traj}*(U(:, k) - u_ss{traj});
    end    
end

%% Figures

figures(U, [Zeta; nan(p.n_eta, p.tf/p.ts + 1)], [Zeta_est; nan(p.n_eta, p.tf/p.ts + 1)], p)

%% Oval Track

U(:, 1) = [0; 10; 10];

X = nan(p.n_x, p.tf/p.ts + 1); 
X(:, 1) = [1; 0; 0; 0; 0; 5]; 

X_est = nan(p.n_x, p.tf/p.ts + 1); 
X_est(:, 1) = [2; pi/4; 0.1; 0; 0; 0]; 

Y = nan(p.n_y, p.tf/p.ts + 1);

P_EKF = eye(p.n_x);

u_Delta = [0.05; 0.1; 0.1];

for k = 1:p.tf/p.ts
    switch feedback
        case 1            
            traj = pathplanning(X(p.n_zeta + 1:p.n_x, k), Delta_eta, p);
            U(:, k) = MPC(A{traj}, B{traj}, P_MPC{traj}, zeta_ss{traj}, u_ss{traj}, X(1:p.n_zeta, k), U(:, max(k - 1, 1)), u_min, u_max, u_Delta, p);
        case 2
            Y(:, k) = measurements(X(:, k), U(:, max(k - 1, 1)), p) + [normrnd(0, 0.01); normrnd(0, 0.01); normrnd(0, 0.01); normrnd(0, 0.1); normrnd(0, 0.1)];
            [X_est(:, k), P_EKF] = EKF(U(:, max(k - 1, 1)), X_est(:, max(k - 1, 1)), Y(:, k), P_EKF, p);
            traj = pathplanning(X_est(p.n_zeta + 1:p.n_x, k), Delta_eta, p);
            U(:, k) = MPC(A{traj}, B{traj}, P_MPC{traj}, zeta_ss{traj}, u_ss{traj}, X_est(1:p.n_zeta, k), U(:, max(k - 1, 1)), u_min, u_max, u_Delta, p);   
    end
    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X(:, k));
    X(:, k + 1) = x(end, :);
end

%% Figures

figures(U, X, X_est, p)