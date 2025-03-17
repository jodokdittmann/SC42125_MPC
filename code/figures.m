function figures(Zeta, U, p, tf, h)

    clf;
    
    figure(1)
    plot(Zeta(5, :), Zeta(6, :))

    figure(2)
    plot(Zeta(4, :))

    figure(3)
    plot(U(1, :))

    figure(4)
    for i = 1:500:tf/h + 1

        plot(Zeta(5, :), Zeta(6, :), 'g--');
        hold on; 
  
        rear = [Zeta(5, i) - p.l_R*cos(Zeta(3, i)), Zeta(6, i) - p.l_R*sin(Zeta(3, i))];
        front = [Zeta(5, i) + p.l_F*cos(Zeta(3, i)), Zeta(6, i) + p.l_F*sin(Zeta(3, i))];
        line([rear(1) front(1)], [rear(2) front(2)], 'Color', 'b', 'LineWidth', 2);
  
        d = 0.6; 
        reartire = [rear(1) - d/2*cos(Zeta(3, i)), rear(1) + d/2*cos(Zeta(3, i)), rear(2) - d/2*sin(Zeta(3, i)), rear(2) + d/2*sin(Zeta(3, i))];
        line([reartire(1) reartire(2)], [reartire(3) reartire(4)], 'Color', 'k', 'LineWidth', 8);
  
        fronttire = [front(1) - d/2*cos(Zeta(3, i) + U(1, i)), front(1) + d/2*cos(Zeta(3, i) + U(1, i)), front(2) - d/2*sin(Zeta(3, i) + U(1, i)), front(2) + d/2*sin(Zeta(3, i) + U(1, i))];
        line([fronttire(1) fronttire(2)],[fronttire(3) fronttire(4)], 'Color','k','LineWidth',8);
  
    end

    axis equal;
    xlim([min(Zeta(5, :)) - 5, max(Zeta(5, :)) + 5]);
    ylim([min(Zeta(6, :)) - 5, max(Zeta(6, :)) + 5]);
    grid on;

end

