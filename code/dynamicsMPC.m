function u = dynamicsMPC(A, B, x_ref, u_ref, x_0, p)

u_min = [-Inf; -Inf; -Inf];
u_max = [Inf; Inf; Inf];
x_min = [-Inf; -Inf; -Inf];
x_max = [Inf; Inf; Inf];

H = blkdiag(kron(speye(p.N), p.Q), p.P, kron(speye(p.N), p.R));
f = [repmat(-p.Q*x_ref, p.N, 1); -p.P*x_ref; repmat(-p.R*u_ref, p.N, 1)];

G_eq = [kron(speye(p.N + 1), -speye(p.n_x)) + kron(sparse(diag(ones(p.N, 1), -1)), A), kron([sparse(1, p.N); speye(p.N)], B)];
g_eq = [-x_0; zeros(p.N*p.n_x, 1)];

u_eq = g_eq;

G_ineq = speye((p.N + 1)*p.n_x + p.N*p.n_u);  
g_ineq = [repmat(x_min, p.N + 1, 1); repmat(u_min, p.N, 1)];

u_ineq = [repmat(x_max, p.N + 1, 1); repmat(u_max, p.N, 1)];

G_c = [G_eq; G_ineq];
g_c = [g_eq; g_ineq];

u_c = [u_eq; u_ineq];

prob = osqp;

prob.setup(H, f, G_c, g_c, u_c, 'verbose', 0);

sol = prob.solve();

u = sol.x((p.N + 1)*p.n_x + 1:(p.N + 1)*p.n_x + p.n_u);

end