function u_ss = equilibria(x_ss, u_init, p)

gamma_ss = fsolve(@(gamma) f(gamma, x_ss, p), [x_ss; u_init]);

u_ss = gamma_ss(p.nx + 1:p.nx + p.nu);

function kappa = f(gamma, x_ss, p)
    x = gamma(1:p.nx);
    u = gamma(p.nx + 1:p.nx + p.nu);
    dxdt = dynamics(x, u, p);
    kappa(1:p.nx) = gamma(1:p.nx) - x_ss;
    kappa(p.nx + 1:p.nx + p.nu) = dxdt;
end

end

