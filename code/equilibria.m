function [delta_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p, x0)

xsol = fsolve(@(x) f(x, V_ss, beta_ss, dpsidt_ss, p), x0);

delta_ss = xsol(7);
omega_F_ss = xsol(8);
omega_R_ss = xsol(9);

function y = f(x, V_ss, beta_ss, dpsidt_ss, p)
    zeta(1) = x(1); % V
    zeta(2) = x(2); % beta
    zeta(3) = x(3); % psi
    zeta(4) = x(4); % dpsidt
    zeta(5) = x(5); % X
    zeta(6) = x(6); % Y
    u(1) = x(7); % ddeltadt
    u(2) = x(8); % omega_F
    u(3) = x(9); % omega_R
    dzetadt = dynamics(zeta, u, p);
    y(1) = x(1) - V_ss;
    y(2) = x(2) - beta_ss;
    y(3) = x(4) - dpsidt_ss;
    y(4) = dzetadt(1); % dVdt
    y(5) = dzetadt(2); % dbetadt
    y(6) = dzetadt(4); % dpsiddt
end

end

