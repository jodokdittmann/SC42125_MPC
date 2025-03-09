function [mu_x, mu_y] = pacejka(s_x, s_y, p)

    s = sqrt(s_x^2 + s_y^2);

	mu = p.D*sin(p.C*atan(p.B*s));

    mu_x = - s_x/s*mu;
    mu_y = - s_y/s*mu;

end
