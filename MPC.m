function u = MPC(A, B, P_LQR, zeta_ss, u_ss, zeta_0, u_0, zeta_min, zeta_max, u_min, u_max, G_f, g_fmin, g_fmax, u_Delta, p)

    H = blkdiag(kron(speye(p.N), p.Q_MPC), P_LQR, kron(speye(p.N), p.R_MPC));
    
    G_dyn = [kron(speye(p.N + 1), -speye(p.n_zeta)) + kron(sparse(diag(ones(p.N, 1), -1)), A), kron([sparse(1, p.N); speye(p.N)], B)];
    g_dyn = [zeta_ss - zeta_0; zeros(p.N*p.n_zeta, 1)];

    G_zeta = [speye(p.N*p.n_zeta), zeros(p.N*p.n_zeta, p.N*p.n_u + p.n_zeta); zeros(size(G_f, 1), p.N*p.n_zeta), G_f, zeros(size(G_f, 1), p.N*p.n_u)];  
    g_zetamin = [kron(ones(p.N, 1), zeta_min - zeta_ss); g_fmin];
    g_zetamax = [kron(ones(p.N, 1), zeta_max - zeta_ss); g_fmax];

    G_u = [zeros(p.N*p.n_u, (p.N + 1)*p.n_zeta), speye(p.N*p.n_u)];  
    g_umin = kron(ones(p.N, 1), u_min - u_ss);
    g_umax = kron(ones(p.N, 1), u_max - u_ss);

    G_Delta = [zeros(p.N*p.n_zeta, (p.N + 1)*p.n_zeta), kron(speye(p.N), speye(p.n_u)) - kron(sparse(diag(ones(p.N - 1, 1), -1)), speye(p.n_u))];

    g_Deltamin = [u_0 - u_ss - u_Delta; repmat(-u_Delta, p.N - 1, 1)];
    g_Deltamax = [u_0 - u_ss + u_Delta; repmat(u_Delta, p.N - 1, 1)];

    prob = osqp;

    prob.setup(H, [], [G_dyn; G_zeta; G_u; G_Delta], [g_dyn; g_zetamin; g_umin; g_Deltamin], [g_dyn; g_zetamax; g_umax; g_Deltamax], 'verbose', 0, 'polish', 0);

    sol = prob.solve();

    u = u_ss + sol.x((p.N + 1)*p.n_zeta + 1:(p.N + 1)*p.n_zeta + p.n_u);

end