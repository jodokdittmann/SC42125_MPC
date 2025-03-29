function [A, B] = matrices(zeta_ss, u_ss, p)

    syms zeta [p.n_zeta 1]
    syms u [p.n_u 1]
    dxdt = dynamics([zeta; zeros(p.n_eta, 1)], u, p);
    dzetadt = dxdt(1:p.n_zeta);
    dfdzeta = matlabFunction(jacobian(dzetadt, zeta), 'Vars', {zeta, u});
    dfdu = matlabFunction(jacobian(dzetadt, u), 'Vars', {zeta, u});

    A = cell(p.n_equi, 1); 
    B = cell(p.n_equi, 1); 

    for i = 1:p.n_equi
        A{i} = dfdzeta(zeta_ss{i}, u_ss{i});
        B{i} = dfdu(zeta_ss{i}, u_ss{i});
        [A{i}, B{i}] = c2d(A{i}, B{i}, p.ts);
    end

end