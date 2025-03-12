function [delta_ss, omega_F_ss, omega_R_ss] = equilibria(V_ss, R_ss, beta_ss, p)

    f_Ry_ss = p.m*V_ss^2/R_ss*cos(beta_ss)*p.l_F/(p.l_F + p.l_R);

    f_Fz_ss = (p.m*p.g*p.l_R + p.m*p.h*V_ss^2*sin(beta_ss)/R_ss)/(p.l_F + p.l_R);
    f_Rz_ss = (p.m*p.g*p.l_F - p.m*p.h*V_ss^2*sin(beta_ss)/R_ss)/(p.l_F + p.l_R);

    mu_Ry_ss = f_Ry_ss/f_Rz_ss;

    alpha_R_ss = atan((V_ss*sin(beta_ss) - V_ss*p.l_R/R_ss)/(V_ss*cos(beta_ss)));

    sol_R = [];
    
    for n = 0:9
        for m = 0:9
            init_R = [n*0.2 - 0.9, m*0.2 - 0.9,  sqrt((n*0.2 - 0.9)^2 + (m*0.2 - 0.9)^2)];

            res_R = fsolve(@(var_R) [tan(alpha_R_ss) - (var_R(2)/(1 + var_R(1))); ...
                                     var_R(3) - sqrt(var_R(1)^2 + var_R(2)^2); ...
                                     mu_Ry_ss + var_R(2) / var_R(3)*p.D*sin(p.C*atan(p.B*var_R(3)))], ...
                                     init_R);
            if isreal(res_R) && (isempty(sol_R) || min(vecnorm(sol_R - res_R, 2, 2)) > 1e-4) && norm(res_R) < 1e2
                sol_R = [sol_R; res_R(:)']; 
            end
        end
    end

    s_Rx_ss = res_R(1);
    s_R_ss = res_R(3);

    mu_R_ss = p.D*sin(p.C*atan(p.B*s_R_ss));
    mu_Rx_ss = -s_Rx_ss/s_R_ss*mu_R_ss;

    f_Rx_ss = mu_Rx_ss*f_Rz_ss;
    f_F_ss = sqrt(p.m^2*V_ss^4/R_ss^2 + f_Rx_ss^2 + f_Ry_ss^2 + 2*p.m*V_ss^2/R_ss*(f_Rx_ss*sin(beta_ss) - f_Ry_ss*cos(beta_ss)));

    mu_F_ss = f_F_ss/f_Fz_ss;

    s_F_ss = tan(asin(mu_F_ss/p.D)/p.C)/p.B;

    init_F = [0.0244, -0.1, 3.2*pi/180];

    sol_F = fsolve(@(var_F) [p.m*V_ss^2/R_ss*sin(beta_ss) - f_F_ss/s_F_ss*(var_F(1)*cos(var_F(3)) - var_F(2)*sin(var_F(3))) + f_Rx_ss; ...
                             var_F(2)/(1 + var_F(1)) - (V_ss*sin(beta_ss - var_F(3)) + V_ss*p.l_F*cos(var_F(3))/R_ss)/(V_ss*cos(beta_ss - var_F(3)) + V_ss*p.l_F*sin(var_F(3))/R_ss); ...
                             s_F_ss - sqrt(var_F(1)^2 + var_F(2)^2)], ...
                             init_F);

    s_Fx_ss = sol_F(1);
    delta_ss = sol_F(3);

    omega_F_ss = (V_ss*cos(beta_ss - delta_ss) + V_ss*p.l_F*sin(delta_ss)/R_ss)/((1 + s_Fx_ss)*p.R_w);
    omega_R_ss = V_ss*cos(beta_ss)/((1 + s_Rx_ss)*p.R_w);

end
