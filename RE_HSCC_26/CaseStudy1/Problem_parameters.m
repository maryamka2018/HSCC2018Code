%% This Script loads the needed problem parameters

% System dynamics
Ts = 0.4;
sys = c2d(ss([zeros(2) eye(2); zeros(2,4)], [zeros(2); eye(2)],...
    [eye(2) zeros(2)], zeros(2)),Ts); % Double integrator model
A = sys.A;
B = sys.B;
[nx,nu] = size(B);

%% Problem parameters
lane = 3.5;
% Car and truck dimensions
car.d = [4.5; 2.5];
car.v = 50/3.6; % initial forward velocity in m/s
truck.d = [9; 2.5];
%truck.v = -20/3.6;  %truck initial velocity
truck.theta0 = 0 ; % truck initial orientation
truck.start = [44; lane/2];
% Road dimensions
xmin = [0; -lane]; % lower bounds
xmax = [100+4*car.d(1); 0]; % upper bounds
% Vehicle dynamics restrictions
vmin = [0/3.6 ; -20/3.6]; % lower bounds
vmax = [80/3.6 ; 20/3.6] ; % upper bounds
umin = [-10; -5];
umax = [3; 5];
N = 10; % prediction horizon
middle = 0; % vertical center of the road
x0 = [ xmin(1)+car.d(1) middle-lane/2 car.v 0]'; % initial state
%% Disturbance
nw = 3; % Truck's position and orientation

%% State and input constraints (Box constraints)
xmin_bold = repmat([xmin + car.d/2; vmin], N+1,1);
xmax_bold = repmat([xmax - car.d/2; vmax], N+1,1);
umin_bold = repmat( umin, N,1);
umax_bold = repmat( umax, N,1);
%% 3rd and 4th states have COUPLED CONSTRAINTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       x4 - c1*x3 <= - c3
%       x4 - c1*x3 >= - c2
%       x4 + c1*x3 <= c2
%       x4 + c1*x3 >= c3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients for the coupled constraints for 3rd and 4th states
c1 = vmax(2)/(0.5*(vmax(1)-vmin(1)));  % they work only assuming  vmax(2) = -vmin(2)
c2 = c1*vmax(1);
c3 = c1*vmin(1);
%% Inputs have COUPLED CONSTRAINTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       u2 - t1*u1 <= - t3
%       u2 - t1*u1 >= - t2
%       u2 + t1*u1 <= t2
%       u2 + t1*u1 >= t3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients for the coupled constraints for 3rd and 4th states
t1 = umax(2)/(0.5*(umax(1)-umin(1)));  % they work only assuming  vmax(2) = -vmin(2)
t2 = t1*umax(1);
t3 = t1*umin(1);

%% 
Abar = eye(N * nx) - [zeros(nx, N * nx); kron(eye(N-1),A) zeros((N-1)*nx,nx)];
Bbar = kron(eye(N),B);
Gamma = [zeros(nx, N*nu); Abar \ Bbar]; % Abar*Gamma = Bbar except for first rows

States_free_init  = [];
for k = 1:N+1
    States_free_init = [ States_free_init ; (A^(k-1))*x0 ];
end
%%
diag = sqrt(car.d(2)^2 + car.d(1)^2)/2;
%%
