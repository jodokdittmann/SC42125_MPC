function [x_ss, u_ss] = equilibria(V_ss, beta_ss, dpsidt_ss, p)

x_ss = cell(p.n_equi, 1); 
u_ss = cell(p.n_equi, 1); 
u_ss{1} = [-40.7*pi/180; 20.44; 58.33]; 
  
for i = 1:p.n_equi
    x_ss{i} = [V_ss; beta_ss(i); dpsidt_ss(i)];
    gamma_ss = fsolve(@(gamma) f(gamma, x_ss{i}, p), [x_ss{i}; u_ss{max(i-1, 1)}], optimset('Display', 'off'));
    u_ss{i} = gamma_ss(p.n_x + 1:p.n_x + p.n_u);
end

function kappa = f(gamma, x_ss, p)
    x = gamma(1:p.n_x);
    u = gamma(p.n_x + 1:p.n_x + p.n_u);
    dxdt = dynamics(x, u, p);
    kappa(1:p.n_x) = gamma(1:p.n_x) - x_ss;
    kappa(p.n_x + 1:p.n_x + p.n_u) = dxdt;
end

end

