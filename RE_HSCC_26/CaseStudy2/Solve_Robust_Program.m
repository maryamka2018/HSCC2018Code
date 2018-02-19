function [ Optimizers , DIAGNOSTIC, Solvertime ] = Solve_Robust_Program(selected_solver, P, W, v, weight, samples, params)
%%  Load needed parameters
lane = params.lane ;
car.d = params.car.d ;
car.v = params.car.v;
x0 = params.x0;
truck.d = params.truck.d ;
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
c1 = params.c1;
c2 = params.c2;
c3 = params.c3;
t1 = params.t1;
t2 = params.t2;
t3 = params.t3;
Ts = params.Ts;
big_M = params.big_M;
Matrix = params.Matrix;

%% Optimization variables
for p = 1:P
    h{p} = sdpvar(nu*N,1,'full');
    H{p} = sdpvar(nu*(N+1) , nw*(N+1) , 'full');
    H{p}(logical(kron(triu(ones(N+1)), ones(nu,nw)))) = 0;
    H{p} = H{p}(nu+1:end,:);
end

Abar = eye(N * nx) - [zeros(nx, N * nx); kron(eye(N-1),A) zeros((N-1)*nx,nx)];
Bbar = kron(eye(N),B);
Gamma = [zeros(nx, N*nu); Abar \ Bbar]; % Abar*Gamma = Bbar except for first rows

States_free_init  = [];
for k = 1:N+1
    States_free_init = [ States_free_init ; (A^(k-1))*x0 ];
end
for p = 1:P
    States_free{p} = States_free_init +  Gamma*h{p};
end

for p = 1:P
    Gamma_H{p} = Gamma*H{p} ;
end
% so that...      All_states{p} = States_free{p} + Gamma_H{p}*w    nx*(N+1)

for p = 1:P
    %car
    delta_1{p} = binvar(N,1);
    delta_2{p} = binvar(N,1);
    delta_3{p} = binvar(N,1);
    delta_4{p} = binvar(N,1);

end

constraints_state_input = [];
constraints_obs = [];


for p = 1:P
   
%% state constraints
    % 1st state 
    for j = 1:N
   
    
    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');

        
    constraints_state_input = [constraints_state_input , States_free{p}(j*nx+1) + v{p}'*z1   <=  xmax_bold(j*nx+1) , ...
            W{p}'*z1 == Gamma_H{p}(j*nx+1,:)' ,...
            z1 >= 0   ,...
            
            v{p}'*z2 - States_free{p}(j*nx+1) <= - xmin_bold(j*nx+1) , ...
            W{p}'*z2 ==  - Gamma_H{p}(j*nx+1,:)' ,...
            z2 >= 0    ];
    end
     % 2nd state 
    for j = 1:N
   
    
    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');
    
        
    constraints_state_input = [constraints_state_input , States_free{p}(j*nx+2) + v{p}'*z1    <=  xmax_bold(j*nx+2) , ...
            W{p}'*z1 == Gamma_H{p}(j*nx+2,:)' ,...
            z1 >= 0  ,...
            
            v{p}'*z2 - States_free{p}(j*nx+2) <= - xmin_bold(j*nx+2) , ...
            W{p}'*z2 ==  -Gamma_H{p}(j*nx+2,:)' ,...
            z2 >= 0             ];
    end

%%%%%%% 3rd and 4th states have COUPLED CONSTRAINTS: ...
%       x4 - c1*x3 <= - c3
%       x4 - c1*x3 >= - c2
%       x4 + c1*x3 <= c2
%       x4 + c1*x3 >= c3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:N
   
    
    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');
    z3 = sdpvar(length(v{p}), 1, 'full');
    z4 = sdpvar(length(v{p}), 1, 'full');
    
    constraints_state_input = [constraints_state_input , States_free{p}(j*nx+4) - c1*States_free{p}(j*nx+3) + v{p}'*z1   <=  - c3 , ...
            W{p}'*z1 == (Gamma_H{p}(j*nx+4,:) - c1*Gamma_H{p}(j*nx+3,:))' ,...
            z1 >= 0   , ...
            
            v{p}'*z2 - States_free{p}(j*nx+4) + c1*States_free{p}(j*nx+3)<= c2 , ...
            W{p}'*z2 ==  -( Gamma_H{p}(j*nx+4,:) - c1*Gamma_H{p}(j*nx+3,:))' ,...
            z2 >= 0             ,...
           
            States_free{p}(j*nx+4) + c1*States_free{p}(j*nx+3) + v{p}'*z3   <=  c2 , ...
            W{p}'*z3 == (Gamma_H{p}(j*nx+4,:) + c1*Gamma_H{p}(j*nx+3,:))' ,...
            z3 >= 0  ,...
            
            v{p}'*z4 - States_free{p}(j*nx+4) - c1*States_free{p}(j*nx+3)<= - c3 , ...
            W{p}'*z4 ==  -( Gamma_H{p}(j*nx+4,:) + c1*Gamma_H{p}(j*nx+3,:))' ,...
            z4 >= 0             ];
    end
     
    
% Inputs have COUPLED CONSTRAINTS:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       u2 - t1*u1 <= - t3
        %       u2 - t1*u1 >= - t2
        %       u2 + t1*u1 <= t2
        %       u2 + t1*u1 >= t3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for j = 1:N
    
    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');
    z3 = sdpvar(length(v{p}), 1, 'full');
    z4 = sdpvar(length(v{p}), 1, 'full');
    
        constraints_state_input = [constraints_state_input , h{p}(nu*(j-1)+2) - t1*h{p}(nu*(j-1)+1)   + v{p}'*z1   <=  -t3 , ...
             W{p}'*z1 == H{p}(nu*(j-1)+2,:)' - t1*H{p}(nu*(j-1)+1,:)' ,...
            z1 >= 0  ,...
            
            -h{p}(nu*(j-1)+2) + t1*h{p}(nu*(j-1)+1)   + v{p}'*z2   <=  t2 , ...
             W{p}'*z2 == - H{p}(nu*(j-1)+2,:)' + t1*H{p}(nu*(j-1)+1,:)' ,...
            z2 >= 0   ,... 
            
            h{p}(nu*(j-1)+2) + t1*h{p}(nu*(j-1)+1)   + v{p}'*z3   <=  t2 , ...
             W{p}'*z3 == H{p}(nu*(j-1)+2,:)' + t1*H{p}(nu*(j-1)+1,:)' ,...
            z3 >= 0  ,...
            
             -h{p}(nu*(j-1)+2) - t1*h{p}(nu*(j-1)+1)   + v{p}'*z4   <=  -t3 , ...
            W{p}'*z4 == - H{p}(nu*(j-1)+2,:)' - t1*H{p}(nu*(j-1)+1,:)' ,...
            z4 >= 0   ];
    end   
    
     
%% obstacle 
    % on the first state.. 
   
    for j = 1:N

    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');

        constraints_obs = [constraints_obs , States_free{p}(j*nx+1) + v{p}'*z1  <=  truck.start - car.d(1)/2 + big_M*(1-delta_1{p}(j)) , ...
            W{p}'*z1 == Gamma_H{p}(j*nx+1,:)' - [1, zeros(1,nw*(N+1)-1)]' ,...
            z1 >= 0   ,...
            
            v{p}'*z2 - States_free{p}(j*nx+1) <= - truck.start - truck.d(1) - car.d(1)/2 + big_M*(1-delta_2{p}(j)) , ...
            W{p}'*z2 ==  -Gamma_H{p}(j*nx+1,:)' + [1, zeros(1,nw*(N+1)-1)]' ,...
            z2 >= 0     ];
    end    
    
    
    % on the second state.. 
    
    for j = 1:N
        
    z1 = sdpvar(length(v{p}), 1, 'full');
    z2 = sdpvar(length(v{p}), 1, 'full');
        
        constraints_obs = [constraints_obs , States_free{p}(j*nx+2) + v{p}'*z1   <=  middle - truck.d(2)/2 - car.d(2)/2 + big_M*(1-delta_3{p}(j)) , ...
            W{p}'*z1 == (Gamma_H{p}(j*nx+2,:) - Matrix(j+1,:))' ,...
            z1 >= 0   ,...
            
            v{p}'*z2 - States_free{p}(j*nx+2) <= - (middle + truck.d(2)/2 + car.d(2)/2) + big_M*(1- delta_4{p}(j)) , ...
            W{p}'*z2 ==  - (Gamma_H{p}(j*nx+2,:) - Matrix(j+1,:))' ,...
            z2 >= 0 , ...
            
            delta_1{p}(j) + delta_2{p}(j) + delta_3{p}(j)+ delta_4{p}(j) == 1]; 
    end   
    
   
end

Constraints = [constraints_state_input , constraints_obs ];

%% objective
objective = 0;
for p =1:P
 w_avg{p} = mean(samples{p},2);
 minimum_time =  - ( States_free{p}(N*nx+1)+ Gamma_H{p}(N*nx+1,:)*w_avg{p} ) ;
 
 objective = objective + weight{p}*( minimum_time ); 
end
 

 %% Solve the problem
    ops = sdpsettings;
    ops.solver = selected_solver;
    ops.verbose = 3;
    ops.debug = 1;
    DIAGNOSTIC = optimize(Constraints, objective, ops);
    Solvertime = DIAGNOSTIC.solvertime;
Optimal_Cost = value(objective);

%% Save Results
for p= 1:P
Optimizers.H{p} = value(H{p});
Optimizers.h{p} = value(h{p});
Optimizers.delta_1{p} = value(delta_1{p});
Optimizers.delta_2{p} = value(delta_2{p});
Optimizers.delta_3{p} = value(delta_3{p});
Optimizers.delta_4{p} = value(delta_4{p});
end

end

