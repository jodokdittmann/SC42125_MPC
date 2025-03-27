function trim = positionMPC(Zeta, Delta, p)

progress = zeros(1, p.n_equi);
objective = zeros(1, p.n_equi);

for k = 1:p.n_equi
       
        X = Zeta(2) + cos(Zeta(1))*Delta{k}(2) - sin(Zeta(1))*Delta{k}(3);
        Y = Zeta(3) + sin(Zeta(1))*Delta{k}(2) + cos(Zeta(1))*Delta{k}(3);

        distX = p.X_ref - X*ones(1, 200);
        distY = p.Y_ref - Y*ones(1, 200);
        dist = distX.^2 + distY.^2;
        [deviation, progress(k)] = min(dist);
        objective(k) = deviation - progress(k);
end

[~, trim] = min(objective);

end