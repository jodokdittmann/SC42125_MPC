function u = MPC(A, B, P_MPC, zeta_ref, u_ref, zeta_0, u_0, u_min, u_max, u_Delta, p)
    
    H = blkdiag(kron(speye(p.N), p.Q_MPC), P_MPC, kron(speye(p.N), p.R_MPC));
    f = [repmat(-p.Q_MPC*zeta_ref, p.N, 1); -P_MPC*zeta_ref; repmat(-p.R_MPC*u_ref, p.N, 1)];

    G_dyn = [kron(speye(p.N + 1), -speye(p.n_zeta)) + kron(sparse(diag(ones(p.N, 1), -1)), A), kron([sparse(1, p.N); speye(p.N)], B)];
    g_dyn = [-zeta_0; zeros(p.N*p.n_zeta, 1)];

    G_u = [zeros(p.N*p.n_zeta, (p.N + 1)*p.n_zeta), speye(p.N*p.n_u)];  
    g_umin = repmat(u_min, p.N, 1);
    g_umax = repmat(u_max, p.N, 1);

    G_Delta = [zeros(p.N*p.n_zeta, (p.N + 1)*p.n_zeta), kron(speye(p.N), speye(p.n_u)) - kron(sparse(diag(ones(p.N - 1, 1), -1)), speye(p.n_u))];

    g_Deltamin = [u_0 - u_Delta; repmat(-u_Delta, p.N - 1, 1)];
    g_Deltamax = [u_0 + u_Delta; repmat(u_Delta, p.N - 1, 1)];

    prob = osqp;

    prob.setup(H, f, [G_dyn; G_u; G_Delta], [g_dyn; g_umin; g_Deltamin], [g_dyn; g_umax; g_Deltamax], 'verbose', 0, 'polish', 1);

    sol = prob.solve();

    u = sol.x((p.N + 1)*p.n_zeta + 1:(p.N + 1)*p.n_zeta + p.n_u);

end