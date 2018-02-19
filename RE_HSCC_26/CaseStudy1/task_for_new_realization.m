function [States, truck_poses ] = task_for_new_realization( H_star, h_star, params,w )

N = params.N;
truck.theta0 = params.truck.theta0 ;
truck.start = params.truck.start ;
Ts = params.Ts;
%%
% compute positions and orientations of the Truck 
pose_init = [truck.start; truck.theta0];
[tout, pose_out] = ode45( @(tout, pose_out) truck_dyn(tout, pose_out, w.v, w.theta/Ts, Ts, N) , [0:Ts:N*Ts], pose_init);
truck_poses = pose_out';

nw = 3;
w_vec = zeros(nw*(N+1),1);
w_vec(1:nw:end) = truck_poses(1,:);
w_vec(2:nw:end) = truck_poses(2,:);
w_vec(3:nw:end) = truck_poses(3,:);  %nw*(N+1) dimensional uncertainty

% compute the applied inputs and positions of the Car
nx = params.nx;
nu = params.nu;
Gamma = params.Gamma;
States_free_init = params.States_free_init;

Inputs =  h_star + H_star*w_vec;
States = States_free_init + Gamma*Inputs;
States =  reshape(States, nx,[]);
Inputs =  reshape(Inputs, nu,[]);

% plot Car and Truck for different time steps 
lane = params.lane;
car.d = params.car.d ;
car.v = params.car.v ;
x0 = params.x0 ;
truck.d = params.truck.d ;
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


end

