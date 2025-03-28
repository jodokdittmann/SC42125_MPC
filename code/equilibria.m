function [zeta_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p)

zeta_ss = cell(p.n_equi, 1); 
u_ss = cell(p.n_equi, 1); 
u_ss{1} = [-40.7*pi/180; 20.44; 58.33]; 
  
for i = 1:p.n_equi
    zeta_ss{i} = [V_ss; beta_ss(i); dpsidt_ss(i)];
    gamma_ss = fsolve(@(gamma) f(gamma, zeta_ss{i}, p), [zeta_ss{i}; u_ss{max(i-1, 1)}], optimset('Display', 'off'));
    u_ss{i} = gamma_ss(p.n_zeta + 1:p.n_zeta + p.n_u);
end

function kappa = f(gamma, zeta_ss, p)
    zeta = gamma(1:p.n_zeta);
    u = gamma(p.n_zeta + 1:p.n_zeta + p.n_u);
    dxdt = dynamics([zeta; zeros(p.n_eta, 1)], u, p);
    dzetadt = dxdt(1:p.n_zeta);
    kappa(1:p.n_zeta) = gamma(1:p.n_zeta) - zeta_ss;
    kappa(p.n_zeta + 1:p.n_zeta + p.n_u) = dzetadt;
end

end

