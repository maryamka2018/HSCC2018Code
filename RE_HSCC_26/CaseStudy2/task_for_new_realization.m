function [car_states, truck_states  ] = task_for_new_realization( H_star, h_star, params, w_realiz )

N = params.N;
nx = params.nx;
nu = params.nu;
A = params.A;
B = params.B;
lane = params.lane;
car.d = params.car.d ;
car.v = params.car.v ;
x0 = params.x0 ;
truck.d = params.truck.d ;
truck.start = params.truck.start ;
xmin = params.xmin ;
xmax = params.xmax ;
vmin  =  params.vmin ;
vmax =  params.vmax ;
umin  =  params.umin ;
umax =  params.umax ;
middle = params.middle ;
c1 = params.c1;
c2 = params.c2;
c3 = params.c3;
t1 = params.t1;
t2 = params.t2;
t3 = params.t3;
P = params.P;
W = params.W;
v = params.v;
Matrix = params.Matrix;

%% compute States and Inputs, given w 
Abar = eye(N * nx) - [zeros(nx, N * nx); kron(eye(N-1),A) zeros((N-1)*nx,nx)];
Bbar = kron(eye(N),B);
Gamma = [zeros(nx, N*nu); Abar \ Bbar]; % Abar*Gamma = Bbar except for first rows

States_free_init  = [];
for k = 1:N+1
    States_free_init = [ States_free_init ; (A^(k-1))*x0 ];
end
for p = 1:P
    States_free{p} = States_free_init +  Gamma*h_star{p};
end

for p = 1:P
    Gamma_H{p} = Gamma*H_star{p} ;
end

for p = 1:P
    if all(W{p}*w_realiz <= v{p})
        States = States_free{p} + Gamma_H{p}*w_realiz;
        States =  reshape(States, nx,[]);
        Inputs =  reshape( h_star{p} + H_star{p}*w_realiz, nu,[]);
        break
    end
end

car_states  = States;
truck_states = [repmat(truck.start + w_realiz(1), 1, N+1); [middle +  Matrix*w_realiz(:)]'];


% %% plots
% x1_min_Truck = truck.start + w_realiz(1) ;
% x1_max_Truck = truck.start + w_realiz(1)+ truck.d(1);
% x2_min_Truck = ( middle - truck.d(2)/2 ) + Matrix*w_realiz;
% x2_max_Truck = ( middle + truck.d(2)/2 )+ Matrix*w_realiz;
% 
% figure
% fig = 1;
% for k = [1 5 6 10] 
%     subplot(2,2,fig)
%     title(['time k =',num2str(k-1)])
%     patch([xmin(1) xmin(1) xmax(1) xmax(1)],[xmin(2), xmax(2) ,xmax(2) , xmin(2)] , 9*[0.1 0.1 0.1])
%     hold on
%     plot(0:1:38, zeros(39,1), '--', 'Color','black', 'Linewidth', 2)
%     plot(0:1:38, -lane*ones(39,1), '--', 'Color','black', 'Linewidth', 2)
%     
%     
%     %% car box
%     x1_min_Car = States(1,k) - car.d(1)/2;
%     x1_max_Car = States(1,k) + car.d(1)/2;
%     x2_min_Car = States(2,k) - car.d(2)/2;
%     x2_max_Car = States(2,k) + car.d(2)/2;
%     
%     patch([x1_min_Car x1_min_Car  x1_max_Car x1_max_Car],[x2_min_Car, x2_max_Car ,x2_max_Car , x2_min_Car] , 'r')
%     
%     %% truck box
%     patch([x1_min_Truck x1_min_Truck  x1_max_Truck x1_max_Truck],[x2_min_Truck(k), x2_max_Truck(k) ,x2_max_Truck(k) , x2_min_Truck(k)] , 6*[0.1 0.1 0.1])
%     plot(States(1,1:k),States(2,1:k),'-.r','Linewidth', 1.5)
%     axis([0, 34, -8.5, 5 ])
%     axis equal
%     fig = fig + 1;
% end
% 
% figure
% subplot(2,1,1)
% plot( States(3,:), States(4,:),'*')
% hold on
% x = vmin(1):0.1:vmax(1);
% plot( x, c1*x - c2, 'Linewidth',1)
% plot( x, -c1*x + c2, 'Linewidth',1)
% plot( x, c1*x - c3, 'Linewidth',1)
% plot( x, -c1*x + c3, 'Linewidth',1)
% % axis([vmin(1) , vmax(1), vmin(2), vmax(2)])
% title('States x4 vs. x3 ')
% 
% % inputs
% subplot(2,1,2)
% plot( Inputs(1,:),Inputs(2,:),'*')
% hold on
% u = umin(1):0.1:umax(1);
% plot( u, t1*u - t2, 'Linewidth',1)
% plot( u, -t1*u + t2, 'Linewidth',1)
% plot( u, t1*u - t3, 'Linewidth',1)
% plot( u, -t1*u + t3, 'Linewidth',1)
% % axis([umin(1) , umax(1), umin(2), umax(2)])
% title('Inputs u2 vs. u1 ')
% Last_horizontal_position = States(1,end)


end

