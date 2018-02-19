function [  ] = plot_multiple_trajectories(num, States_multiple, part, params)

nx = params.nx;
nu = params.nu;
lane = params.lane;
car.d = params.car.d ;
x0 = params.x0 ;
truck.d = params.truck.d ;
truck.start = params.truck.start ;
xmin = params.xmin ;
xmax = params.xmax ;
middle = params.middle ;
P = params.P;


Colors =   1*[ 0    0   0 ;...
               1    0   0 ;...
               0    1   0 ;...
               0    0   1 ;...
               0.5  0   0.5;...
               0.5  0.5  0;...
               0    0.5  0.5;...
               0.5  0.5  0.5;...
               0.3  0.8  0.3];
%%
figure();
patch([xmin(1) xmin(1) xmax(1) xmax(1)],[xmin(2), xmax(2) ,xmax(2) , xmin(2)] , 9*[0.1 0.1 0.1])
    hold on
    plot(0:1:38, zeros(39,1), '--', 'Color','black', 'Linewidth', 2)
    plot(0:1:38, -lane*ones(39,1), '--', 'Color','black', 'Linewidth', 2)
    
x1_min_Truck = truck.start;
x1_max_Truck = truck.start + truck.d(1);
x2_min_Truck = ( middle - truck.d(2)/2 );
x2_max_Truck = ( middle + truck.d(2)/2 );
patch([x1_min_Truck x1_min_Truck  x1_max_Truck x1_max_Truck],[x2_min_Truck, x2_max_Truck ,x2_max_Truck , x2_min_Truck] , 6*[0.1 0.1 0.1])

x1_min_Car = x0(1) - car.d(1)/2;
x1_max_Car = x0(1) + car.d(1)/2;
x2_min_Car = x0(2) - car.d(2)/2;
x2_max_Car = x0(2) + car.d(2)/2;
patch([x1_min_Car x1_min_Car  x1_max_Car x1_max_Car],[x2_min_Car, x2_max_Car ,x2_max_Car , x2_min_Car] , 'r')
for i = 1:1000
    plot(States_multiple{i}(1,:),States_multiple{i}(2,:), 'Color', Colors(part{i},:) ,'Linewidth', 0.2)
end
axis equal

end

