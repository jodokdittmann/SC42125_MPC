function figures(Zeta, X, U, p)

    clf;
    
    figure(1)
    plot(Zeta(2, :), Zeta(3, :))

    figure(2)
    plot(Zeta(1, :))

    figure(3)
    plot(U(1, :))

    figure(4)
    for i = 1:500:p.tf/p.ts + 1

        plot(Zeta(2, :), Zeta(3, :), 'k--');
        hold on; 
  
        rear = [Zeta(2, i) - p.l_R*cos(Zeta(1, i)), Zeta(3, i) - p.l_R*sin(Zeta(1, i))];
        front = [Zeta(2, i) + p.l_F*cos(Zeta(1, i)), Zeta(3, i) + p.l_F*sin(Zeta(1, i))];
        line([rear(1) front(1)], [rear(2) front(2)], 'Color', 'k', 'LineWidth', 2);
  
        d = 0.6; 
        reartire = [rear(1) - d/2*cos(Zeta(1, i)), rear(1) + d/2*cos(Zeta(1, i)), rear(2) - d/2*sin(Zeta(1, i)), rear(2) + d/2*sin(Zeta(1, i))];
        line([reartire(1) reartire(2)], [reartire(3) reartire(4)], 'Color', 'k', 'LineWidth', 8);
  
        fronttire = [front(1) - d/2*cos(Zeta(1, i) + U(1, i)), front(1) + d/2*cos(Zeta(1, i) + U(1, i)), front(2) - d/2*sin(Zeta(1, i) + U(1, i)), front(2) + d/2*sin(Zeta(1, i) + U(1, i))];
        line([fronttire(1) fronttire(2)],[fronttire(3) fronttire(4)], 'Color', 'k','LineWidth', 8);
  
    end

    axis equal;
    xlim([min(Zeta(2, :)) - 5, max(Zeta(2, :)) + 5]);
    ylim([min(Zeta(3, :)) - 5, max(Zeta(3, :)) + 5]);
    grid on;

end

