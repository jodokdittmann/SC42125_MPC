function [X_estimated_new, P_estimated_new, P_time_new] = ekf(Zeta_noised, Zeta_EKF, X_time, U, p, P_time_update, R, Q, Fk) % Takes them as 3d vectors

    dzetadt = position(Zeta_noised, X_time);
    dxdt = dynamics(X_time, U, p);
    Hk = dzetadt / dxdt; % this is Hk that we get from chain rule

    Kk = P_time_update * Hk' * 1/(Hk * P_time_update * Hk' + R);

    X_estimated_new = X_time + Kk*(Zeta_noised - Zeta_EKF);
    P_estimated_new = P_time_update - P_time_update * Hk' * 1/(Hk*P_time_update*Hk' + R)* Hk * P_time_update;

    % X time new in main
    
    P_time_new = Fk*P_estimated_new*Fk' + Q;


end