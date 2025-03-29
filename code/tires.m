function [f_Fx, f_Rx, f_Fy, f_Ry] = tires(V, beta, dpsidt, delta, omega_F, omega_R, p)

    V_Fx = V*cos(beta - delta) + dpsidt*p.l_F*sin(delta);
    V_Rx = V*cos(beta);

    s_Fx = V_Fx/(p.R_w*omega_F) - 1;
    s_Rx = V_Rx/(p.R_w*omega_R) - 1;

    alpha_F = atan((V*sin(beta) + dpsidt*p.l_F)/(V*cos(beta))) - delta;
    alpha_R = atan((V*sin(beta) - dpsidt*p.l_R)/(V*cos(beta)));

    s_Fy = (1 + s_Fx)*tan(alpha_F);
    s_Ry = (1 + s_Rx)*tan(alpha_R);

    [mu_Fx, mu_Fy] = pacejka(s_Fx, s_Fy, p);
    [mu_Rx, mu_Ry] = pacejka(s_Rx, s_Ry, p);

    f_Fz = (p.l_R - p.h*mu_Rx)*p.m*p.g/(p.l_F + p.l_R + p.h*(mu_Fx*cos(delta) - mu_Fy*sin(delta) - mu_Rx));
    f_Rz = p.m*p.g - f_Fz;

    f_Fx = mu_Fx*f_Fz;
    f_Rx = mu_Rx*f_Rz;

    f_Fy = mu_Fy*f_Fz;
    f_Ry = mu_Ry*f_Rz;

    function [mu_x, mu_y] = pacejka(s_x, s_y, p)

        s = sqrt(s_x^2 + s_y^2);
       
	    mu = p.D*sin(p.C*atan(p.B*s));

        mu_x = - s_x/s*mu;
        mu_y = - s_y/s*mu;

    end

end
