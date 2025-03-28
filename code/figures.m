function figures(Zeta, Eta, U, p)

    clf;
    
    figure('Name', 'Velocity', 'NumberTitle', 'off');
    plot(Zeta(1, :))

    figure('Name', 'Slip Angle', 'NumberTitle', 'off');
    plot(Zeta(2, :))

    figure('Name', 'Yaw Rate', 'NumberTitle', 'off');
    plot(Zeta(3, :))

    figure('Name', 'Position', 'NumberTitle', 'off');
    plot(Eta(2, :), Eta(3, :))

    figure('Name', 'Orientation', 'NumberTitle', 'off');
    plot(Eta(1, :))

    figure('Name', 'Steering Angle', 'NumberTitle', 'off');
    plot(U(1, :))

    figure('Name', 'Trajectory', 'NumberTitle', 'off');
    for i = 1:500:p.tf/p.ts + 1

        plot(Eta(2, :), Eta(3, :), 'k--');
        hold on; 
  
        rear = [Eta(2, i) - p.l_R*cos(Eta(1, i)), Eta(3, i) - p.l_R*sin(Eta(1, i))];
        front = [Eta(2, i) + p.l_F*cos(Eta(1, i)), Eta(3, i) + p.l_F*sin(Eta(1, i))];
        line([rear(1) front(1)], [rear(2) front(2)], 'Color', 'k', 'LineWidth', 2);
  
        d = 0.6; 
        reartire = [rear(1) - d/2*cos(Eta(1, i)), rear(1) + d/2*cos(Eta(1, i)), rear(2) - d/2*sin(Eta(1, i)), rear(2) + d/2*sin(Eta(1, i))];
        line([reartire(1) reartire(2)], [reartire(3) reartire(4)], 'Color', 'k', 'LineWidth', 8);
  
        fronttire = [front(1) - d/2*cos(Eta(1, i) + U(1, i)), front(1) + d/2*cos(Eta(1, i) + U(1, i)), front(2) - d/2*sin(Eta(1, i) + U(1, i)), front(2) + d/2*sin(Eta(1, i) + U(1, i))];
        line([fronttire(1) fronttire(2)],[fronttire(3) fronttire(4)], 'Color', 'k','LineWidth', 8);
  
    end

    axis equal;
    xlim([min(Eta(2, :)) - 5, max(Eta(2, :)) + 5]);
    ylim([min(Eta(3, :)) - 5, max(Eta(3, :)) + 5]);
    grid on;
    hold off

end

