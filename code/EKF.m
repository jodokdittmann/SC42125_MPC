function [x_est, P] = EKF(x_est, P, u, eta_meas, p)

f = @(~, x) [dynamics(x(1:p.n_zeta), u, p); position(x(p.n_zeta + 1:p.n_x), x(1:p.n_zeta))];
F = eye(p.n_x) + p.ts*numjac(f, [], x_est, f([], x_est), 1e-6, [], 0);

H = [zeros(p.n_zeta, p.n_zeta), eye(p.n_eta)];

% Measurement Update

K = P*H'/(H*P*H' + p.R_EKF);

x_est = x_est + K*(eta_meas - H*x_est);
P = P - K*H*P;

% Time Update

x_est = x_est + p.ts*f([], x_est);
P = F*P*F' + p.Q_EKF;

end