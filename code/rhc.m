function ctrl = rhc(A_d, B_d, Q, R, zeta_r, u_r, zeta)

    [n, m] = size(B_d);

    % Constraints
    umin = [-Inf; -Inf; -Inf];
    umax = [Inf; Inf; Inf];
    xmin = [-Inf; -Inf; -Inf];
    xmax = [Inf; Inf; Inf];

    % Objective function
    QN = 10*Q;

    % Initial and reference states
    x0 = zeta;
   
    % Prediction horizon
    N = 10;

    % Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
    % - quadratic objective
    H = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
    % - linear objective
    f = [repmat(-Q*zeta_r, N, 1); -QN*zeta_r; repmat(-R*u_r, N, 1)];
    % - linear dynamics
    Ax = kron(speye(N+1), -speye(n)) + kron(sparse(diag(ones(N, 1), -1)), A_d);
    Bu = kron([sparse(1, N); speye(N)], B_d);
    Aeq = [Ax, Bu];
    leq = [-x0; zeros(N*n, 1)];
    ueq = leq;
    % - input and state constraints
    Aineq = speye((N+1)*n + N*m);
    lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
    uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
    % - OSQP constraints
    A = [Aeq; Aineq];
    l = [leq; lineq];
    u = [ueq; uineq];

    % Create an OSQP object
    prob = osqp;

    % Setup workspace
    prob.setup(H, f, A, l, u);

    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end

    % Apply first control input to the plant
    ctrl = res.x((N+1)*n+1:(N+1)*n+m);
    
    % % Update initial state
    % l(1:nx) = -x0;
    % u(1:nx) = -x0;
    % prob.update('l', l, 'u', u);

end