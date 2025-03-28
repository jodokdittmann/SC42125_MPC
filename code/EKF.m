function [x_est, P] = EKF(u, x_est, H, y, P, p)

f = @(~, x) dynamics(x, u, p);
F = eye(p.n_x) + p.ts*numjac(f, [], x_est, f([], x_est), 1e-6, [], 0);

% Measurement Update

K = P*H'/(H*P*H' + p.R_EKF);

x_est = x_est + K*(y - H*x_est);
P = P - K*H*P;

% Time Update

x_est = x_est + p.ts*f([], x_est);
P = F*P*F' + p.Q_EKF;

end