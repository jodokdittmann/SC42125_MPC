clear
init 

h = 1e-3;
tf = 50;

m = 3;
n = 6; 

zeta0 = zeros(n, 1);
zeta0(1) = 7;

U = nan(tf/h + 1, m);
Zeta = nan(tf/h + 1, n);
Zeta(1, :) = zeta0;

for k = 1:tf/h
	U(k, :) = [0.0761 24.9092 16.5773];
    [t, zeta] = ode45(@(t, zeta) dynamics(zeta, U(k, :), p), [(k - 1)*h k*h], Zeta(k, :));
	Zeta(k + 1, :) = zeta(end, :);
end

figure(1)
plot(Zeta(:,5),Zeta(:,6))

figure(3)
plot(Zeta(:,4))

figure(4)
plot(atan(Zeta(:,1)./Zeta(:,2)))

visualize_vehicle([Zeta(:,5), Zeta(:,6)], Zeta(:,3), p)

function visualize_vehicle(position, heading, p)
  
    car_length = p.l_F + p.l_R; 
    car_width = 2.0; 
    
    corners = [-car_length/2, -car_width/2;
               car_length/2, -car_width/2;
               car_length/2,  car_width/2;
              -car_length/2,  car_width/2]';
    
    figure(5); clf; hold on;
    
    for i = 1:1000:size(position, 1) 
        pos = position(i, :); 
        theta = heading(i);    

        R = [cos(theta), -sin(theta);
             sin(theta),  cos(theta)];
        
        rotated_corners = R * corners;
        car_corners = rotated_corners + pos'; 
        
        car_corners = [car_corners, car_corners(:,1)];
        
        plot(car_corners(1,:), car_corners(2,:), 'b-', 'LineWidth', 2);
        
        quiver(pos(1), pos(2), cos(theta), sin(theta), 1.5, 'r', 'LineWidth', 2);
      
    end
    
   axis equal;
    xlim([min(position(:,1))-5, max(position(:,1))+5]);
    ylim([min(position(:,2))-5, max(position(:,2))+5]);
    grid on;
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Single Track Vehicle Visualization');
end
