function traj = pathplanning(eta, Delta_eta, p)

    progress = zeros(1, p.n_equi);
    objective = zeros(1, p.n_equi);

    for k = 1:p.n_equi
        X = eta(2) + cos(eta(1))*Delta_eta{k}(2) - sin(eta(1))*Delta_eta{k}(3);
        Y = eta(3) + sin(eta(1))*Delta_eta{k}(2) + cos(eta(1))*Delta_eta{k}(3);
        distX = p.X_ref - X;
        distY = p.Y_ref - Y;
        dist = distX.^2 + distY.^2;
        [deviation, progress(k)] = min(dist);
        objective(k) = 10*deviation - progress(k);
    end

    [~, traj] = min(objective);

end