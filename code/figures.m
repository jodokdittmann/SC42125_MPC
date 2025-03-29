function figures(U, X, X_est, Y, p)

    clf;
    
    figure('Name', 'Velocity', 'NumberTitle', 'off')
    plot(X(1, :))
    hold on
    plot(X_est(1, :))
    hold off

    figure('Name', 'Slip Angle', 'NumberTitle', 'off')
    plot(X(2, :))
    hold on
    plot(X_est(2, :))
    hold off

    figure('Name', 'Yaw Rate', 'NumberTitle', 'off')
    plot(X(3, :))
    hold on
    plot(X_est(3, :))
    hold off

    figure('Name', 'Position', 'NumberTitle', 'off')
    plot(X(5, :), X(6, :))
    hold on
    plot(X_est(5, :), X_est(6, :))
    hold off

    figure('Name', 'Orientation', 'NumberTitle', 'off')
    plot(X(4, :))
    hold on
    plot(X_est(4, :))
    hold off

    figure('Name', 'Steering Angle', 'NumberTitle', 'off')
    plot(U(1, :))

    figure('Name', 'Front Wheel Speed', 'NumberTitle', 'off')
    plot(U(2, :))

    figure('Name', 'Rear Wheel Speed', 'NumberTitle', 'off')
    plot(U(3, :))

    figure('Name', 'Trajectory', 'NumberTitle', 'off')
    for i = 1:500:p.tf/p.ts + 1

        plot(X(5, :), X(6, :), 'k--')
        hold on
  
        rear = [X(5, i) - p.l_R*cos(X(4, i)), X(6, i) - p.l_R*sin(X(4, i))];
        front = [X(5, i) + p.l_F*cos(X(4, i)), X(6, i) + p.l_F*sin(X(4, i))];
        line([rear(1) front(1)], [rear(2) front(2)], 'Color', 'k', 'LineWidth', 2);
  
        d = 0.6; 
        reartire = [rear(1) - d/2*cos(X(4, i)), rear(1) + d/2*cos(X(4, i)), rear(2) - d/2*sin(X(4, i)), rear(2) + d/2*sin(X(4, i))];
        line([reartire(1) reartire(2)], [reartire(3) reartire(4)], 'Color', 'k', 'LineWidth', 8);
  
        fronttire = [front(1) - d/2*cos(X(4, i) + U(1, i)), front(1) + d/2*cos(X(4, i) + U(1, i)), front(2) - d/2*sin(X(4, i) + U(1, i)), front(2) + d/2*sin(X(4, i) + U(1, i))];
        line([fronttire(1) fronttire(2)],[fronttire(3) fronttire(4)], 'Color', 'k','LineWidth', 8);
  
    end

    axis equal;
    xlim([min(X(5, :)) - 5, max(X(5, :)) + 5]);
    ylim([min(X(6, :)) - 5, max(X(6, :)) + 5]);
    grid on;
    hold off

end

