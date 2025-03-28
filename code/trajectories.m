function Delta = trajectories(x_ss, p)

Delta = cell(p.n_equi, 1);

for k = 1:p.n_equi
    [~, eta] = ode45(@(t, eta) position(eta, x_ss{k}), [0, p.th], zeros(3,1));
    Delta{k} = eta(end, :)';
end

end