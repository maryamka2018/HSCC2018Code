%%
% System dynamics
Ts = 0.4;
cr = 0.3;
sys = c2d(ss([zeros(2) eye(2); zeros(2,4)], [zeros(2); eye(2)],...
    [eye(2) zeros(2)], zeros(2)),Ts); % Double integrator model
A = sys.A;
B = sys.B;
[nx,nu] = size(B);

%% Problem parameters
lane = 3.5;
% Car and truck dimensions
car.d = [4.5; 2];
car.v = 100/3.6; % forward velocity of the moving frame in m/s
truck.d = [9; 2.5];
truck.v = 100/3.6;
truck.start = 11;   % +- uncertainty
% Road dimensions
xmin = [0; -2*lane]; % lower bounds
xmax = [20+4*car.d(1); lane]; % upper bounds
% Vehicle dynamics restrictions
vmin = [90/3.6-car.v; -20/3.6]; % lower bounds
vmax = [130/3.6-car.v; 20/3.6] ; % upper bounds
umin = [-10; -5];
umax = [3; 5];
N = 10; % prediction horizon

%%
middle = -2*lane+3*lane/2;
x0 = [xmin(1)+car.d(1) middle+lane 110/3.6-car.v 0]'; % initial state


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
t1 = umax(2)/(0.5*(umax(1)-umin(1)));  % they work only assuming  umax(2) = -umin(2)
t2 = t1*umax(1);
t3 = t1*umin(1);

%% Uncertainties 
alpha = 1.5/3.6; % maximum lateral velocity for the obstacle truck (normalized for Ts s)
offset_max = 2;  % w.r.t. truck.start(1)
wmin = [ - offset_max ; -alpha ];
wmax = [ offset_max ; alpha];
nw = length(wmin);

wmin_bold = [ wmin ; repmat([ 0; wmin(2)], N,1)];
wmax_bold = [ wmax ; repmat([ 0;  wmax(2)], N,1)];
%%
big_M = 100;
Matrix = kron(tril(ones(N+1,N+1)) , [0, 1]);
Matrix = [zeros(1,nw*(N+1)); Matrix(1:end-1,:)];

%% 
params.lane = lane;
params.car.d  = car.d;
params.car.v = car.v;
params.x0 = x0;
params.truck.d = truck.d;
params.truck.start = truck.start;
params.xmin = xmin;
params.xmax = xmax;
params.xmax_bold = xmax_bold;
params.xmin_bold = xmin_bold;
params.vmin = vmin;
params.vmax = vmax;
params.umin = umin;
params.umax = umax;
params.N = N;
params.middle = middle;
params.nx = nx;
params.nu = nu;
params.nw = nw;
params.A = A;
params.B = B;
params.c1 = c1;
params.c2 = c2;
params.c3 = c3;
params.t1 = t1;
params.t2 = t2;
params.t3 = t3;
params.Ts = Ts;
params.wmax_bold = wmax_bold ;
params.wmin_bold = wmin_bold ;
params.Matrix = Matrix;
params.big_M = big_M;
