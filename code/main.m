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

Zeta = nan(p.n_zeta, p.tf/p.ts + 1);
Zeta(:, 1) = [0; 0; 10];

U = nan(p.n_u, p.tf/p.ts + 1); 

% Simulation

for k = 1:p.tf/p.ts
    trim = positionMPC(Zeta(:, k), x_ss, p);
    U(:, k) = dynamicsMPC(A{trim}, B{trim}, x_ss{trim}, u_ss{trim}, X(:, k), p);
    [t, x] = ode45(@(t, x) dynamics(x, U(:, k), p), [(k - 1)*p.ts k*p.ts], X(:, k));
    [t, zeta] = ode45(@(t, zeta) position(zeta, x(t == t, :)), t, Zeta(:, k));
    X(:, k + 1) = x(end, :);
    Zeta(:, k + 1) = zeta(end, :);
end

% Figures

figures(Zeta, X, U, p)  