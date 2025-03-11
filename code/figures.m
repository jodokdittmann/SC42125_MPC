function figures(x, y, psi, p, tf, h)
    
    figure(5); clf; hold on;
    
    for i = 1:1000:tf/h + 1
       
        car = [- p.l_F*cos(psi(i)) + p.w/2*sin(psi(i)) + x(i), p.l_R*cos(psi(i)) + p.w/2*sin(psi(i)) + x(i), p.l_R*cos(psi(i)) - p.w/2*sin(psi(i)) + x(i), - p.l_F*cos(psi(i)) - p.w/2*sin(psi(i)) + x(i);
              - p.l_F*sin(psi(i)) - p.w/2*cos(psi(i)) + y(i), p.l_R*sin(psi(i)) - p.w/2*cos(psi(i)) + y(i), p.l_R*sin(psi(i)) + p.w/2*cos(psi(i)) + y(i), - p.l_F*sin(psi(i)) + p.w/2*cos(psi(i)) + y(i)];
        
        plot([car(1, :), car(1, 1)], [car(2, :), car(2, 1)], 'b-', 'LineWidth', 2);

        quiver(x(i), y(i), cos(psi(i)), sin(psi(i)), 1.5, 'r', 'LineWidth', 2);
    end
    
    axis equal;
    xlim([min(x)-5, max(x)+5]);
    ylim([min(y)-5, max(y)+5]);
    grid on;
end
