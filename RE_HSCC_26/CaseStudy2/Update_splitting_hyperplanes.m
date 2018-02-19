function [g, l, time_splitting] = Update_splitting_hyperplanes (selected_solver, g,l, Optimizers, P, W, v, params)
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

%%
for p = 2:P

% load('Workspace')  % H_star, h_star , delta_1 , delta_2, delta_3, delta_4 + all the problem parameters 
H_star = Optimizers.H{p};
h_star = Optimizers.h{p};
delta_1_star = Optimizers.delta_1{p};
delta_2_star = Optimizers.delta_2{p};
delta_3_star = Optimizers.delta_3{p};
delta_4_star = Optimizers.delta_4{p};
%%

w = sdpvar(nw*(N+1),1,'full');


Constraints = [   W{p}*w <= v{p}  ];

Abar = eye(N * nx) - [zeros(nx, N * nx); kron(eye(N-1),A) zeros((N-1)*nx,nx)];
Bbar = kron(eye(N),B);
Gamma = [zeros(nx, N*nu); Abar \ Bbar]; % Abar*Gamma = Bbar except for first rows

States_free_init  = [];
for k = 1:N+1
    States_free_init = [ States_free_init ; (A^(k-1))*x0 ];
end
States_free = States_free_init +  Gamma*h_star;
Gamma_H = Gamma*H_star ;

x1_bold = States_free(nx+1:nx:end) + Gamma_H(nx+1:nx:end,:)* w ;
x2_bold = States_free(nx+2:nx:end) + Gamma_H(nx+2:nx:end,:)* w ;
x3_bold = States_free(nx+3:nx:end) + Gamma_H(nx+3:nx:end,:)* w ;
x4_bold = States_free(nx+4:nx:end) + Gamma_H(nx+4:nx:end,:)* w ;


% state constraints
state_constr = [  x1_bold - xmax_bold(nx+1:nx:end) ;...
                - x1_bold + xmin_bold(nx+1:nx:end);...
                  x2_bold - xmax_bold(nx+2:nx:end) ;...
                - x2_bold + xmin_bold(nx+2:nx:end);...
               % 3rd and 4th states have coupled constraints 
                
                  x4_bold - c1*x3_bold + c3 ;...
                  - x4_bold + c1*x3_bold - c2 ;...
                  x4_bold + c1*x3_bold - c2   ;...
                  - x4_bold - c1*x3_bold + c3  ];

% input constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       u2 - t1*u1 <= - t3
%       u2 - t1*u1 >= - t2
%       u2 + t1*u1 <= t2
%       u2 + t1*u1 >= t3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

u1_bold = h_star(1:nu:end) + H_star(1:nu:end,:)* w; 
u2_bold = h_star(2:nu:end) + H_star(2:nu:end,:)* w;

input_constr = [ u2_bold - t1*u1_bold + t3 ;...
                 - u2_bold + t1*u1_bold  - t2   ;...
                 u2_bold + t1*u1_bold - t2    ;...
                - u2_bold - t1*u1_bold + t3   ];
% obstacle constraints 

obs_constr = [  x1_bold - ( truck.start + w(1) - car.d(1)/2 + big_M*(1- delta_1_star) ) ;...
               - x1_bold + truck.start + w(1) + truck.d(1) + car.d(1)/2 - big_M*(1-delta_2_star);...             
                 x2_bold - ( middle + Matrix(2:end,:)*w - truck.d(2)/2 - car.d(2)/2 + big_M*(1-delta_3_star));...
                 -x2_bold + middle + Matrix(2:end,:)*w + truck.d(2)/2 + car.d(2)/2 - big_M*(1- delta_4_star) ];
  
%% Objective   
tol = 1e-4;
robust_constraints = [state_constr; input_constr; obs_constr] ;           
objective = - robust_constraints - tol;

%% Solve the problem 
ops = sdpsettings;
ops.solver = selected_solver;
ops.verbose = 0;
ops.debug = 1;

Obj = ones(length(robust_constraints),1);
W_star = zeros(nw*(N+1), length(robust_constraints));
time_splitting = 0;
for i = 1:length(robust_constraints)
DIAGNOSTIC = optimize(Constraints, objective(i), ops);
time_splitting = time_splitting + DIAGNOSTIC.solvertime;
Obj(i) = value(objective(i));
W_star(:,i) = value(w);
end
W_star = W_star(:,Obj < 0);

%% Find w_a and w_b, the most distant (using 2-norm) vectors in W_star(1:2,:)
Vectors = W_star(1:2,:);
dist = 0;
for i = 1:size(Vectors,2)/2+1
    for j = 1:size(Vectors,2)
       dist_curr = sum((Vectors(:,i)-Vectors(:,j)).^2);
       if dist_curr > dist 
           w_a = Vectors(:,i);
           w_b = Vectors(:,j);
           dist = dist_curr;
       end
    end
end
% w_a = [1.98476811042264;0];
% w_b = [0.958655172484004;0.416559715528461];
% 
new_P = 1 + 2*(P-1);
g{new_P}{1+ (p-2)*2 + 1} = [ g{P}{p}, [w_a- w_b; zeros(nw*N,1)]];
l{new_P}{1+ (p-2)*2 + 1} = [l{P}{p} ; 0.5*[w_a- w_b]'*[w_a + w_b]];
g{new_P}{1+ (p-2)*2 + 2} = [ g{P}{p}, -[w_a- w_b; zeros(nw*N,1)]];
l{new_P}{1+ (p-2)*2 + 2} = [l{P}{p} ; -0.5*[w_a- w_b]'*[w_a + w_b]];
end


end
