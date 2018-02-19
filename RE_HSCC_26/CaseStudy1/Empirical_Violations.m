function [Car_pos, violations_LTL, violations_input_state  ] = Empirical_Violations( runs, tol, H_star, h_star , params, truck_poses)

N = params.N;
nx = params.nx;
nu = params.nu;
Gamma = params.Gamma;
States_free_init = params.States_free_init;
truck.theta0 = params.truck.theta0;
lane = params.lane;
car.d = params.car.d ;
car.v = params.car.v ;
x0 = params.x0 ;
truck.d = params.truck.d ;
truck.theta0 = params.truck.theta0 ;
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
xmax_bold = params.xmax_bold;
xmin_bold = params.xmin_bold;
diag = params.diag;
Ts = params.Ts;
nw = 3;

%%    
    violations_input_state = 0;
    violations_LTL = 0;
    
    
    
    % compute empyrical violations
    w_realiz_vec = zeros(nw*(N+1), runs);
    
    for i = 1:runs
 
        w_realiz_vec(1:nw:end,i) = truck_poses{i}(1,:)';
        w_realiz_vec(2:nw:end,i) = truck_poses{i}(2,:)';
        w_realiz_vec(3:nw:end,i) = truck_poses{i}(3,:)';
        
        Inputs =  h_star + H_star*w_realiz_vec(:,i);
        States{i} = States_free_init + Gamma*Inputs;
        %             States{i} =  reshape(States{i}, nx,[]);
        %             Inputs =  reshape(Inputs, nu,[]);
        
        % Input-State violations
        %%%%%%% 3rd and 4th states have COUPLED CONSTRAINTS: ...
        %       x4 - c1*x3 <= - c3
        %       x4 - c1*x3 >= - c2
        %       x4 + c1*x3 <= c2
        %       x4 + c1*x3 >= c3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inputs have COUPLED CONSTRAINTS:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       u2 - t1*u1 <= - t3
        %       u2 - t1*u1 >= - t2
        %       u2 + t1*u1 <= t2
        %       u2 + t1*u1 >= t3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aux = [...,
            States{i}(nx+1:nx:end) <= xmax_bold(nx+1:nx:end) + tol ,...
            States{i}(nx+1:nx:end) >= xmin_bold(nx+1:nx:end) - tol,...
            States{i}(nx+2:nx:end) <= xmax_bold(nx+2:nx:end) + tol,...
            States{i}(nx+2:nx:end) >= xmin_bold(nx+2:nx:end) - tol,...
            
            States{i}(nx+4:nx:end) - c1* States{i}(nx+3:nx:end) <= - c3 + tol,...
            States{i}(nx+4:nx:end) - c1* States{i}(nx+3:nx:end) >= - c2 - tol,...
            States{i}(nx+4:nx:end) + c1* States{i}(nx+3:nx:end) <=  c2 + tol,...
            States{i}(nx+4:nx:end) + c1* States{i}(nx+3:nx:end) >=  c3 - tol,...
            
            Inputs(2:nu:end) - t1*Inputs(1:nu:end) <= -t3 + tol,...
            Inputs(2:nu:end) - t1*Inputs(1:nu:end) >= -t2 - tol,...
            Inputs(2:nu:end) + t1*Inputs(1:nu:end) <=  t2 + tol,...
            Inputs(2:nu:end) + t1*Inputs(1:nu:end) >=  t3 - tol ];
        
        if any(any((aux == 0)))
            violations_input_state = violations_input_state + 1 ;
        end
        
        % LTL violations
       truck_y1 =  w_realiz_vec(nw+1:nw:end,i);
       truck_y2 =  w_realiz_vec(nw+2:nw:end,i);
       truck_theta =  w_realiz_vec(nw+3:nw:end,i);
        
        aux = [  ...
            (States{i}(nx+2:nx:end)- truck_y2 ).*cos(truck_theta) + diag <= (sin(truck_theta).*(States{i}(nx+1:nx:end) - truck_y1  ) - truck.d(2)/2 ) + tol , ...
            (- States{i}(nx+2:nx:end) + truck_y2 ).*cos(truck_theta) + diag <= ( - sin(truck_theta).*(States{i}(nx+1:nx:end) - truck_y1  ) - truck.d(2)/2) + tol , ...
            (States{i}(nx+2:nx:end)- truck_y2 ).*sin(truck_theta) + diag <= ( - cos(truck_theta).*(States{i}(nx+1:nx:end) - truck_y1  ) - truck.d(1)/2)+ tol , ...
            (- States{i}(nx+2:nx:end) + truck_y2 ).*sin(truck_theta) + diag <= ( cos(truck_theta).*(States{i}(nx+1:nx:end) - truck_y1  ) - truck.d(1)/2)  + tol];
        aux = sum(aux,2);
        
        if any(aux == 0)
            violations_LTL = violations_LTL  + 1;
        end
        
        Car_pos{i} = reshape(States{i},nx,[]);
    end
    
    
%     violations_LTL
%     violations_input_state
    
end

