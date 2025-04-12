function [zeta_est, P] = KF(A, B, C, D, u_ss, zeta_ss, u, zeta_est, z_tilde, P, p)

    u_tilde = u - u_ss;
    zeta_tilde = zeta_est - zeta_ss;

    % Measurement Update

    K = P*C'/(C*P*C' + p.R_KF);

    zeta_tilde = zeta_tilde + K*(z_tilde - C*zeta_tilde - D*u_tilde);
    P = P - K*C*P;

    % Time Update

    zeta_tilde = A*zeta_tilde + B*u_tilde;
    P = A*P*A' + p.Q_KF;

    zeta_est = zeta_ss + zeta_tilde;

end