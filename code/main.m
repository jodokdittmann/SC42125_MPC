clear
init
        
V_ss = 7;
beta_ss = linspace(-40*pi/180, 40*pi/180, 15); 
dpsidt_ss = linspace(1, -1, 15); 

[x_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p);

[A, B] = matrices(x_ss, u_ss, p);   

% Initialisation

X = nan(p.n_x, p.tf/p.ts + 1); 
X(:, 1) = [1; 0; 0]; 
X_real = X; 

X_estimated = X;
X_time_update = X;
P_estimated = zeros(3,3,3001);
P_estimated(:,:, 1) = eye(3,3);
P_time_update = P_estimated;
Fk = eye(3);  % I think it should be like this because its linearized at ts

Zeta = nan(p.n_zeta, p.tf/p.ts + 1);
Zeta(:, 1) = [0; 0; 10];
Zeta_real = Zeta;
Zeta_noised = Zeta;
Zeta_EKF = Zeta;

U = nan(p.n_u, p.tf/p.ts + 1); 
U_real = U;

R = eye(3);
Q = eye(3);

% Simulation

for k = 1:p.tf/p.ts

    zeta_noise = normrnd(0,1,[3,1]);  % additive noises
    x_noise = normrnd(0,1,[3,1]);

    % Real states like original code, this is left for reference
    trim = positionMPC(Zeta_real(:, k), x_ss, p); 
    U_real(:, k) = dynamicsMPC(A{trim}, B{trim}, x_ss{trim}, u_ss{trim}, X_real(:, k), p);
    [t, x_real] = ode45(@(t, x_real) dynamics(x_real, U_real(:, k), p), [(k - 1)*p.ts k*p.ts], X_real(:, k));
    [t, zeta_real] = ode45(@(t, zeta_real) position(zeta_real, x_real(t == t, :)), t, Zeta_real(:, k));
    Zeta_real(:, k + 1) = zeta_real(end, :);
    X_real(:, k + 1) = x_real(end, :); 

    %EKF states
    U(:, k) = dynamicsMPC(A{trim}, B{trim}, x_ss{trim}, u_ss{trim}, X_estimated(:, k), p); % Calculates U using X_hat_k|k
    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X(:, k));
    X(:, k + 1) = x(end, :); % This is the IRL state
    [t, zeta] = ode45(@(t, zeta) position(zeta, x(t == t, :)), t, Zeta(:, k));
    Zeta(:, k + 1) = zeta(end, :); % This is the IRL position

    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X_time_update(:, k));
    X_time_update(:, k + 1) = x(end, :); % This is the dynamics for EKF
    [t, zeta] = ode45(@(t, zeta) position(zeta, x(t == t, :)), t, Zeta_EKF(:, k));
    Zeta_EKF(:, k + 1) = zeta(end, :); % This is the EKF position

    Zeta_noised(:, k + 1) = Zeta(:, k + 1) + zeta_noise; % yk measurement used in EKF
        k
     [X_estimated(:, k+1), P_estimated(:,:, k+1), P_time_update(:,:,k+1) ]= ekf(Zeta_noised(:, k + 1) , Zeta_EKF(:,k+1), X_time_update(:, k), U(:,k), p, P_time_update(:,:,k), R, Q, Fk);
end

% Figures

figures(Zeta_real, X_real, U_real, p)  
figures(Zeta_noised, X_estimated, U, p)
