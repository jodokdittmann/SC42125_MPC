function [G_f, g_fmin, g_fmax] = terminalset(A, B, K_LQR, zeta_min, zeta_max, u_min, u_max, p)

G = [eye(p.n_zeta); K_LQR];
G_f = [];
g_fmin = [];
g_fmax = [];
opt = optimoptions('linprog','Display','iter');

for k = 0:Inf
    G_f = [G_f; G*(A - B*K_LQR)^k];
    g_fmin = [g_fmin; [zeta_min; u_min]];
    g_fmax = [g_fmax; [zeta_max; u_max]];
    min = nan(2*p.n_zeta, 1);
    max = nan(2*p.n_zeta, 1);
    for i = 1:2*p.n_zeta
        [~, min(i)] = linprog(G(i, :)*(A - B*K_LQR)^(k+1), [-G_f; G_f], [-g_fmin; g_fmax], [], [], [], [], opt);
        [~, max(i)] = linprog(-G(i, :)*(A - B*K_LQR)^(k+1), [-G_f; G_f], [-g_fmin; g_fmax], [], [], [], [], opt);
    end
    if all([-min; -max] <= [-zeta_min; -u_min; zeta_max; u_max])
        break
    end
end
    
