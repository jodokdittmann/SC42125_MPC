function u = MPC(A, B, zeta_ref, u_ref, zeta_0, p)

    u_min = [-Inf; -Inf; -Inf];
    u_max = [Inf; Inf; Inf];
    zeta_min = [-Inf; -Inf; -Inf];
    zeta_max = [Inf; Inf; Inf];

    H = blkdiag(kron(speye(p.N), p.Q_MPC), p.P_MPC, kron(speye(p.N), p.R_MPC));
    f = [repmat(-p.Q_MPC*zeta_ref, p.N, 1); -p.P_MPC*zeta_ref; repmat(-p.R_MPC*u_ref, p.N, 1)];

    G_eq = [kron(speye(p.N + 1), -speye(p.n_zeta)) + kron(sparse(diag(ones(p.N, 1), -1)), A), kron([sparse(1, p.N); speye(p.N)], B)];
    g_eq = [-zeta_0; zeros(p.N*p.n_zeta, 1)];

    u_eq = g_eq;

    G_ineq = speye((p.N + 1)*p.n_zeta + p.N*p.n_u);  
    g_ineq = [repmat(zeta_min, p.N + 1, 1); repmat(u_min, p.N, 1)];

    u_ineq = [repmat(zeta_max, p.N + 1, 1); repmat(u_max, p.N, 1)];

    G_c = [G_eq; G_ineq];
    g_c = [g_eq; g_ineq];

    u_c = [u_eq; u_ineq];

    prob = osqp;

    prob.setup(H, f, G_c, g_c, u_c, 'verbose', 0);

    sol = prob.solve();

    u = sol.x((p.N + 1)*p.n_zeta + 1:(p.N + 1)*p.n_zeta + p.n_u);

end