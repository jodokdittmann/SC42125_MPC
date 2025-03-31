function [A, B, C, D, P_MPC] = matrices(zeta_ss, u_ss, p)

    syms zeta [p.n_zeta 1]
    syms u [p.n_u 1]
    dxdt = dynamics([zeta; zeros(p.n_eta, 1)], u, p);
    dzetadt = dxdt(1:p.n_zeta);
    dfdzeta = matlabFunction(jacobian(dzetadt, zeta), 'Vars', {zeta, u});
    dfdu = matlabFunction(jacobian(dzetadt, u), 'Vars', {zeta, u});
    y = measurements([zeta; zeros(p.n_eta, 1)], u, p);
    z = y(1:p.n_z);
    dhdzeta = matlabFunction(jacobian(z, zeta), 'Vars', {zeta, u});
    dhdu = matlabFunction(jacobian(z, u), 'Vars', {zeta, u});

    A = cell(p.n_equi, 1); 
    B = cell(p.n_equi, 1);
    C = cell(p.n_equi, 1);
    D = cell(p.n_equi, 1);
    P_MPC = cell(p.n_equi, 1);

    for i = 1:p.n_equi
        A{i} = dfdzeta(zeta_ss{i}, u_ss{i});
        B{i} = dfdu(zeta_ss{i}, u_ss{i});
        [A{i}, B{i}] = c2d(A{i}, B{i}, p.ts);
        C{i} = dhdzeta(zeta_ss{i}, u_ss{i});
        D{i} = dhdu(zeta_ss{i}, u_ss{i});
        [~, P_MPC{i}, ~] = dlqr(A{i}, B{i}, p.Q_MPC, p.R_MPC);
    end

end