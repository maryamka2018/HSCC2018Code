function [ H_star, h_star, DIAGNOSTIC ] = Solve_Scenario_Program(params, w , selected_solver)
%%  Load needed parameters
lane = params.lane ;
car.d = params.car.d ;
car.v = params.car.v;
x0 = params.x0;
truck.d = params.truck.d ;
truck.theta0 = params.truck.theta0;
truck.start = params.truck.start ;
xmin = params.xmin;
xmax = params.xmax;
xmax_bold = params.xmax_bold;
xmin_bold = params.xmin_bold;
vmin = params.vmin ;
vmax = params.vmax ;
umin = params.umin ;
umax = params.umax ;
N = params.N;
middle = params.middle;
nx = params.nx;
nu = params.nu;
nw = params.nw;
A = params.A;
B = params.B;
Gamma = params.Gamma ;
States_free_init = params.States_free_init;
c1 = params.c1;
c2 = params.c2;
c3 = params.c3;
t1 = params.t1;
t2 = params.t2;
t3 = params.t3;
diag = params.diag;
Ts = params.Ts;
K_tot = params.K_tot;
K = params.K;

%% Optimization variables

h = sdpvar(nu*N,1,'full');
H = sdpvar(nu*(N+1) , nw*(N+1) , 'full');
H(logical(kron(triu(ones(N+1)), ones(nu,nw)))) = 0;
H = H(nu+1:end,:);
H(7:20,1:3)= 0;
H(9:20,4:6)= 0;
H(11:20,7:9)= 0;
H(13:20,10:12)= 0;
H(15:20,13:15)= 0;
H(17:20,16:18)= 0;
H(19:20,19:21)= 0;

States_free = States_free_init +  Gamma*h;
Gamma_H = Gamma*H ;
% so that...      All_states = States_free + Gamma_H*w     nx*(N+1)

% binary variables for obstacle
delta1_1 = binvar(N,1);
delta1_2 = binvar(N,1);
delta1_3 = binvar(N,1);
delta1_4 = binvar(N,1);


big_M = 200; % conservative upper-bound

% sampling
w_vec = zeros(nw*(N+1),K);

for k = 1:K
    % compute truck poses, for the given w{k}
    pose_init = [truck.start; truck.theta0];
    [tout, pose_out] = ode45( @(tout, pose_out) truck_dyn(tout, pose_out, w{k}.v, w{k}.theta/Ts, Ts, N) , [0:Ts:N*Ts], pose_init);
    truck_poses = pose_out';
    
    w_vec(1:nw:end,k) = truck_poses(1,:)';
    w_vec(2:nw:end,k) = truck_poses(2,:)';
    w_vec(3:nw:end,k) = truck_poses(3,:)';  %nw*(N+1) dimensional uncertainty
end
%% Constraints
constraints = [];
constraints_obs = [];

x_bold = sdpvar(nx*N,K,'full');

% state box constraints
A_cell = repmat({A},1,N);
B_cell = repmat({B},1,N);

constraints = [constraints,...
    x_bold == blkdiag(A_cell{:})* [repmat(x0,1,K); x_bold(1:end-nx,:)] + repmat(blkdiag(B_cell{:})*h,1,K) +  blkdiag(B_cell{:})*H*w_vec(:,:) ];

x1_bold = x_bold(1:nx:end,:);
x2_bold = x_bold(2:nx:end,:);
x3_bold = x_bold(3:nx:end,:);
x4_bold = x_bold(4:nx:end,:);

constraints = [constraints , x1_bold(:,1:K_tot(1))  <=  repmat(xmax_bold(nx+1:nx:end),1,K_tot(1)) , ...
    - x1_bold(:,1:K_tot(1)) <= - repmat(xmin_bold(nx+1:nx:end),1,K_tot(1)) , ...
    x2_bold(:,1:K_tot(1))  <=  repmat(xmax_bold(nx+2:nx:end),1,K_tot(1)) , ...
    - x2_bold(:,1:K_tot(1)) <= - repmat(xmin_bold(nx+2:nx:end),1,K_tot(1)) ];

%%%%%%% 3rd and 4th states have COUPLED CONSTRAINTS: ...
%       x4 - c1*x3 <= - c3
%       x4 - c1*x3 >= - c2
%       x4 + c1*x3 <= c2
%       x4 + c1*x3 >= c3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constraints = [constraints , x4_bold(:,1:K_tot(1)) - c1*x3_bold(:,1:K_tot(1))  <=  - repmat(c3, N, K_tot(1)) , ...
    - x4_bold(:,1:K_tot(1)) + c1*x3_bold(:,1:K_tot(1))  <=   repmat(c2, N, K_tot(1)),...
    x4_bold(:,1:K_tot(1)) + c1*x3_bold(:,1:K_tot(1))  <=   repmat(c2, N, K_tot(1))
    - x4_bold(:,1:K_tot(1)) - c1*x3_bold(:,1:K_tot(1))  <=  - repmat(c3, N, K_tot(1)) ];

%% Inputs have COUPLED CONSTRAINTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       u2 - t1*u1 <= - t3
%       u2 - t1*u1 >= - t2
%       u2 + t1*u1 <= t2
%       u2 + t1*u1 >= t3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1_bold = repmat(h(1:nu:end),1,K_tot(1)) + H(1:nu:end,:)*w_vec(:,1:K_tot(1));
u2_bold = repmat(h(2:nu:end),1,K_tot(1)) + H(2:nu:end,:)*w_vec(:,1:K_tot(1));

constraints = [constraints , u2_bold - t1*u1_bold  <=  - repmat(t3, N, K_tot(1)) , ...
    - u2_bold + t1*u1_bold  <=   repmat(t2, N, K_tot(1)),...
    u2_bold + t1*u1_bold  <=   repmat(t2, N, K_tot(1))
    - u2_bold - t1*u1_bold  <=  - repmat(t3, N, K_tot(1)) ];

for k = K_tot(1)+1:K
    
    % truck positions and orientations, given the sample w, from
    % the second time step
    
    truck_y1 = w_vec(nw+1:nw:end,k);
    truck_y2 = w_vec(nw+2:nw:end,k);
    truck_theta = w_vec(nw+3:nw:end,k);
    
    constraints_obs = [constraints_obs , ...
        (  x2_bold(:,k) - truck_y2 ).*cos(truck_theta) + diag <= (sin(truck_theta).*(x1_bold(:,k) - truck_y1  ) - truck.d(2)/2 ) + big_M*(1-delta1_1) , ...
        (- x2_bold(:,k) + truck_y2 ).*cos(truck_theta) + diag <= ( - sin(truck_theta).*(x1_bold(:,k) - truck_y1 ) - truck.d(2)/2) + big_M*(1-delta1_2) , ...
        ( x2_bold(:,k) - truck_y2 ).*sin(truck_theta) + diag <= ( - cos(truck_theta).*(x1_bold(:,k) - truck_y1 ) - truck.d(1)/2) + big_M*(1-delta1_3) , ...
        (- x2_bold(:,k) + truck_y2 ).*sin(truck_theta) + diag <= ( cos(truck_theta).*(x1_bold(:,k) - truck_y1 ) - truck.d(1)/2) + big_M*(1-delta1_4) , ...
        
        
        delta1_1 + delta1_2 + delta1_3+ delta1_4 == 1];
    
end

Constraints = [constraints, constraints_obs];


%% objective

objective = - ( States_free(N*nx+1)+ Gamma_H(N*nx+1,:)*mean(w_vec,2));


%% Solve the problem
ops = sdpsettings;
ops.solver = selected_solver;
ops.verbose = 3;
ops.debug = 1;
DIAGNOSTIC = optimize(Constraints, objective, ops)


%% Save results

h_star = value(h);
H_star = value(H);
H_star(isnan(H_star)) = 0;
h_star(isnan(h_star)) = 0;
end

