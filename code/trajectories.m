function Delta = trajectories(zeta_ss, p)

Delta = cell(p.n_equi, 1);

for k = 1:p.n_equi
    [~, eta] = ode45(@(t, eta) f(zeta_ss{k}, eta, p), [0, p.th], zeros(p.n_eta, 1));
    Delta{k} = eta(end, :);
end

function detadt = f(zeta, eta, p)
    dxdt = dynamics([zeta; eta], zeros(p.n_u, 1), p);
    detadt = dxdt(p.n_zeta + 1:p.n_x);
end

end