function [x_est, sqrtPx] = UKF(u, x_est, y, sqrtPx, p)  

    f = @(x) x + p.ts*dynamics(x, u, p);
    h = @(x) measurements(x, u, p);
   
    % Time Update

    sigma_X = [x_est, (repmat(x_est, 1, p.n_x) + p.gamma*sqrtPx), (repmat(x_est, 1, p.n_x) - p.gamma*sqrtPx)];
  
    for i = 1:1 + 2*p.n_x
        sigma_X(:,i) = f(sigma_X(:,i));
    end
 
    x_est = p.Wm*sigma_X(:,1) + p.W*sum(sigma_X(:,2:end), 2);
  
    [~, R] = qr([sqrt(p.W)*(sigma_X(:,2:end) - repmat(x_est, 1, p.n_x*2)), p.sqrtQ_UKF]', 0);
    sqrtPx = cholupdate(R, sqrt(p.Wc)*(sigma_X(:,1) - x_est), '-')';

    % Measurement Update
    
    sigma_X = [x_est, (repmat(x_est, 1, p.n_x) + p.gamma*sqrtPx), (repmat(x_est, 1, p.n_x) - p.gamma*sqrtPx)];
    
    sigma_Y = nan(p.n_y, 1 + 2*p.n_x);
    for i = 1:1 + 2*p.n_x
        sigma_Y(:,i) = h(sigma_X(:,i));
    end
    
    y_est = p.Wm*sigma_Y(:,1) + p.W*sum(sigma_Y(:,2:end), 2);

    [~, R] = qr([sqrt(p.W)*(sigma_Y(:,2:end) - repmat(y_est, 1, p.n_x*2)), p.sqrtR_UKF]', 0);
    sqrtPy = cholupdate(R, sqrt(p.Wc)*(sigma_Y(:,1) - y_est), '-')';
    
    delta_X = sigma_X - repmat(x_est, 1, p.n_x*2 + 1);
    delta_Z = sigma_Y - repmat(y_est, 1, p.n_x*2 + 1);
  
    Pxy = p.Wc*delta_X(:,1)*delta_Z(:,1)' + p.W*delta_X(:,2:end)*delta_Z(:,2:end)';
    
    K = (Pxy/sqrtPy')/sqrtPy;
    
    x_est = x_est + K*(y - y_est);

    U = K*sqrtPy;
  
    R = sqrtPx';
    for i = 1:size(U, 2)
        R = cholupdate(R, U(:,i), '-');
    end
    sqrtPx = R';
  
end

