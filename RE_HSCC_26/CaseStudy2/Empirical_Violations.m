function [ Truck_states, States, Inputs, part , violations ] = Empirical_Violations( H_star, h_star, params, w_realiz, tol)

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
xmax_bold = params.xmax_bold;
xmin_bold = params.xmin_bold;
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
big_M = params.big_M;
wmax_bold = params.wmax_bold;
wmin_bold = params.wmin_bold;

%%
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

%% 
violations =  0;

for i = 1:size(w_realiz,2)
    for p = 1:P
    if all(W{p}*w_realiz(:,i) <= v{p})
        States{i} = States_free{p} + Gamma_H{p}*w_realiz(:,i);
        States{i} =  reshape(States{i}, nx,[]);
        Inputs{i} =  reshape( h_star{p} + H_star{p}*w_realiz(:,i), nu,[]);
        part{i} = p;
        break
    end
    end
    
constraints_state_input = [ States{i}(1,2:end)' <= xmax_bold(nx+1:nx:end) + tol;...
                           -States{i}(1,2:end)' <= -xmin_bold(nx+1:nx:end)+ tol;...
                            States{i}(2,2:end)' <= xmax_bold(nx+2:nx:end)+ tol;...
                           -States{i}(2,2:end)' <= -xmin_bold(nx+2:nx:end)+ tol;...
%%%%%%% 3rd and 4th states have COUPLED CONSTRAINTS: ...
%       x4 - c1*x3 <= - c3
%       x4 - c1*x3 >= - c2
%       x4 + c1*x3 <= c2
%       x4 + c1*x3 >= c3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            States{i}(4,2:end)'- c1*States{i}(3,2:end)' <= -c3+ tol;...
                            States{i}(4,2:end)'- c1*States{i}(3,2:end)' >= -c2 - tol;...
                            States{i}(4,2:end)'+ c1*States{i}(3,2:end)' <= c2 + tol;...
                            States{i}(4,2:end)'+ c1*States{i}(3,2:end)' >= c3 - tol;...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       u2 - t1*u1 <= - t3
%       u2 - t1*u1 >= - t2
%       u2 + t1*u1 <= t2
%       u2 + t1*u1 >= t3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                            
                            Inputs{i}(2,:)' - t1*Inputs{i}(1,:)' <= -t3 + tol;...
                            Inputs{i}(2,:)' - t1*Inputs{i}(1,:)'  >= -t2 - tol;...
                            Inputs{i}(2,:)' + t1*Inputs{i}(1,:)'  <= t2 + tol;...
                            Inputs{i}(2,:)' + t1*Inputs{i}(1,:)'  >= t3 - tol];

constraints_obs = [ States{i}(1,2:end)' <= truck.start + w_realiz(1,i) - car.d(1)/2 + tol ,...
                    States{i}(1,2:end)' >=  truck.start + w_realiz(1,i) + truck.d(1) + car.d(1)/2 - tol,...
                    States{i}(2,2:end)' <=  middle + Matrix(2:end,:)*w_realiz(:,i) - truck.d(2)/2 - car.d(2)/2 + tol,...
                    States{i}(2,2:end)' >=  middle + Matrix(2:end,:)*w_realiz(:,i) + truck.d(2)/2 + car.d(2)/2 - tol];
aux = sum(constraints_obs,2);

if any(aux == 0) || any(constraints_state_input == 0)
    violations = violations + 1;
end 
     
    Truck_states{i} = [repmat(truck.start + w_realiz(1,i), 1, N+1); [middle +  Matrix*w_realiz(:,i)]'];
end


end

