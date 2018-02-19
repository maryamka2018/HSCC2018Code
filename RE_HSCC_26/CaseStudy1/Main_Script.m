%% script header
clear all %#ok
close all

%% seed used to generate realizations
% [train, veri] = RandStream.create('mrg32k3a','NumStreams', 2, 'Seed',17);
%% 
fprintf(['-------------------------------------------------------------------------\n',...
        'This Script reproduces the data related to \nCase Study 5.1: Dealing with a turning truck.\n']);
fprintf(' Loading problem parameters...\n');
run Problem_parameters

%% decide number of samples according to Theorem 4.6
% create a matrix with the same sparsity of H, to count the number of free parameters in H  
H_sparsity = ones(nu*(N+1) , nw*(N+1));
H_sparsity(logical(kron(triu(ones(N+1)), ones(nu,nw)))) = 0;
H_sparsity = H_sparsity(nu+1:end,:);
H_sparsity(7:20,1:3)= 0;
H_sparsity(9:20,4:6)= 0;
H_sparsity(11:20,7:9)= 0;
H_sparsity(13:20,10:12)= 0;
H_sparsity(15:20,13:15)= 0;
H_sparsity(17:20,16:18)= 0;
H_sparsity(19:20,19:21)= 0;

eps = [0.3 , 0.05]; % violation parameters
% Find corresponding number of samples using a binary search on the implicit criterion 
K_tot = zeros(1,length(eps));
for idx=1:length(eps)
    beta = 1e-3;

    dim_H = sum(nonzeros(H_sparsity)); % nu*nw*N*(N+1)/2
    d  = (nu*N + dim_H );
    m = 4^(N); % The number of feasible binary configurations is << 2^(P*(4*N));
    K_tot(idx) = round( exp(1)/(exp(1)-1) * 1/eps(idx)*(log(m/beta) + d  - 1 ) ) ; % # of samples
    
    % binary search
    up = K_tot(idx);
    low = 0;
    while up-low > 1
        mid =  round( low + (up-low)/2 );
        if binocdf(d-1, mid, eps(idx)) > beta/m
            low = mid;
        else
            up = mid;
        end
    end
    K_tot(idx) = up;
end
K = sum(K_tot);

%% Set problem parameters
params.lane = lane;
params.car.d  = car.d;
params.car.v = car.v;
params.x0 = x0;
params.truck.d = truck.d;
params.truck.theta0 = truck.theta0;
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
params.Gamma = Gamma;
params.States_free_init = States_free_init;
params.c1 = c1;
params.c2 = c2;
params.c3 = c3;
params.t1 = t1;
params.t2 = t2;
params.t3 = t3;
params.diag = diag;
params.Ts = Ts;
params.K_tot = K_tot;
params.K = K;

fprintf(' Problem parameters loaded.\n');
%%
fprintf(['To quickly reproduce the results found in the paper, select option 1) for all the following questions.\n' ...
         '-------------------------------------------------------------------------\n']);

%% solver selection
ok = 0;
while ok == 0
    prompt = ['Select one of the following solvers:\n', ...
        ' 1) CPLEX\n',...
        ' 2) GUROBI (Using GUROBI the results may be slightly different)\n',...
        'Note: The decision is not relevant if option 1) is selected for all the future questions,\n',...
        'because no optimization is performed in such a case. \n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end

if selected == '1' 
    selected_solver = 'cplex';
elseif selected == '2'
    selected_solver = 'gurobi';
end

%% Solve the scenario program
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'For obtaining the optimal policy, the following three options are available:\n', ...
        ' 1) Load the optimizer of the Scenario Program used in the paper.\n',...
        ' 2) Solve the Scenario Program for the same realization of the multisample used in the paper. \n',...
        ' 3) Solve the Scenario Program for a random realization of the multisample.\n',... 
        'Note: When option 2) or 3) are selected, the Scenario Program is solved in\n' ...
        '       approximately 1 hour on a 8 Gb machine at 2.8 GHz using CPLEX.\n',...
        '\nChoose (1,2 or 3):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2' || selected == '3'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end
    
if selected == '1' 
    load('Saved_data/optimizers_15_sims') % LOAD 1x15 cells 
    H_star{1} = H_star_loaded{1};
    h_star{1} = h_star_loaded{1};
    DIAGNOSTIC{1} = DIAGNOSTIC_loaded{1};
    
elseif selected == '2' 
load('Saved_data/w_train_15_sims')
[H_star{1} , h_star{1}, DIAGNOSTIC{1}] = Solve_Scenario_Program(  params ,w_train{1} , selected_solver);

elseif selected == '3'
    % draw a new multisample
    w = cell(1,K);
    for k = 1:K
        w{k}.theta = zeros(N+1,1);
        for i = 1:N
            w{k}.theta(i) = pi/2/(N+1) - 0.03  +  min(0.06, pi/2-sum(w{k}.theta)- pi/2/(N+1) + 0.03)*rand;
        end
        w{k}.theta(N+1) = pi/2 - sum(w{k}.theta);
        w{k}.v = -22/3.6*ones(N+1,1);    % the truck has constant linear velocity at -20 Km/h
    end
 
    [H_star{1} , h_star{1}, DIAGNOSTIC{1}] = Solve_Scenario_Program( params, w , selected_solver);
end


%% Test the policy with a new realization 
fprintf(['-------------------------------------------------------------------------\n' ...
    'The generated policy will be tested against a new realization of the uncertainty, and Figure 1) will be plotted.\n']);

ok = 0;
while ok == 0
    prompt = ['The following two options are available:\n', ...
        ' 1) Evaluate the policy over the sample realization chosen in the paper.\n',...
        ' 2) Evaluate the policy over a random sample realization. \n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end

if selected == '1'
    load('Saved_data/w_verification')
    best  = Inf;
    for i=1:10000
        if sum(w_veri{i}.theta(1:5))< best 
            best = sum(w_veri{i}.theta(1:5));
            idx = i;
        end
    end
    w_new = w_veri{idx};  % corresponds to a very slow maneuver (for illustrative purposes)
elseif selected == '2'
    w_new.theta = zeros(N+1,1);
    for i = 1:N
        w_new.theta(i) = pi/2/(N+1) - 0.03  +  min(0.06, pi/2-sum(w_new.theta)- pi/2/(N+1) + 0.03)*rand;
    end
    w_new.theta(N+1) = pi/2 - sum(w_new.theta);
    w_new.v = -22/3.6*ones(N+1,1) ;    % the truck has constant linear velocity at -20 Km/h
end


[car_states, truck_states]= task_for_new_realization ( H_star{1}, h_star{1}, params, w_new);

%% Plots
h = figure();
idx = 1;
imgwidth = 3.75*0.8;
objects{1}.type = 'car';
objects{1}.d = params.car.d;
objects{2}.type = 'truck';
objects{2}.d = params.truck.d;

k = 4;
objects{1}.k = k;
objects{2}.k = k;
subplot(3,1,1);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                          'objects', objects, ...
                          'highlight', 1, ...
                          'imgwidth', imgwidth, ...
                          'h', h, ...
                          'tikz', false, 'matlab', true, ...
                          'dims', [2.5 45 -15 params.lane+1], ...
                          'xlabel', false, ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087-0.06), ...
                          'name', ['turning_traj_' num2str(k)]);
title({'Figure 1: The controlled car ({\color{blue}solid}) reacting to a new truck trajectory ({\color{red}dashed})', 'using the feedback policy.', ' ', [' time k = ',num2str(k)]});

k = 8;
objects{1}.k = k;
objects{2}.k = k;
spy.color = 'blue';
spy.magnification = 6;
spy.size = 0.195*imgwidth; % In inches?
spy.pos = ([28.33,-2.589]-[2.5,-15])*1.73/10/3.75*imgwidth;
spy.viewer = ([39,-3.5]-[2.5,-15])*1.73/10/3.75*imgwidth;
subplot(3,1,2);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                          'objects', objects, ...
                          'highlight', 1, ...
                          'collisionBoxes', true, ...
                          'imgwidth', imgwidth, ...
                          'h', h, ...   
                          'tikz', false, 'matlab', true, ...
                          'spy', spy, ...
                          'dims', [2.5 45 -15 params.lane+1], ...
                          'xlabel', false, ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087), ...
                          'name', ['turning_traj_' num2str(k)]);
title([' time k = ',num2str(k)]);
k = 10;
objects{1}.k = k;
objects{2}.k = k;
subplot(3,1,3);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                          'objects', objects, ...
                          'highlight', 1, ...
                          'imgwidth', imgwidth, ...
                          'h', h, ...
                          'tikz', false, 'matlab', true, ...
                          'dims', [2.5 45 -15 params.lane+1], ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087-0.06), ...
                          'name', ['turning_traj_' num2str(k)]);
title([' time k = ',num2str(k)]);

fprintf(['-------------------------------------------------------------------------\n' ...
    'Figure 1) shows the car performing the task for the new realization of the uncertainty.\n\n']);
savefig(h, 'Figures/Figure_1.fig')

%% TIKZ generation         
if false
    k = 4; %#ok
    objects{1}.k = k;
    objects{2}.k = k;
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                              'objects', objects, ...
                              'highlight', 1, ...
                              'imgwidth', imgwidth, ...
                              'tikz', true, 'matlab', false, ...
                              'dims', [2.5 45 -params.lane-1 params.lane+1], ...
                              'xlabel', false, ...
                              'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087-0.06), ...
                              'name', ['turning_traj_' num2str(k)]);
    k = 8;
    objects{1}.k = k;
    objects{2}.k = k;
    spy.color = 'blue';
    spy.magnification = 12;
    spy.size = 0.195*imgwidth; % In inches?
    spy.pos = ([28.7,-2.589]-[2.5,-10.5])*1.73/10/3.75*imgwidth;
    spy.viewer = ([50,-3.5]-[2.5,-10.5])*1.73/10/3.75*imgwidth;
    spy.viewerDesired = [];
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                              'objects', objects, ...
                              'highlight', 1, ...
                              'collisionBoxes', true, ...
                              'imgwidth', imgwidth, ...
                              'tikz', true, 'matlab', false, ...
                              'spy', spy, ...
                              'dims', [2.5 45 -10.5 params.lane+1], ...
                              'xlabel', false, ...
                              'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087), ...
                              'name', ['turning_traj_' num2str(k)]);
    k = 10;
    objects{1}.k = k;
    objects{2}.k = k;
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
                              'objects', objects, ...
                              'highlight', 1, ...
                              'imgwidth', imgwidth, ...
                              'tikz', true, 'matlab', false, ...
                              'dims', [2.5 45 -14 params.lane+1], ...
                              'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',0.087), ...
                              'name', ['turning_traj_' num2str(k)]);
end

%% Test the policy with 10000 new realizations 
fprintf(['-------------------------------------------------------------------------\n' ...
    'The policy is now tested against 10000 new realizations of the uncertainty, and Figure 2.a) will be plotted.\n']);
ok = 0;
while ok == 0
    prompt = ['The following two options are available:\n', ...
        ' 1) Evaluate the policy over the 10000 sample realizations used in the paper.\n',...
        ' 2) Evaluate the policy over 10000 random sample realizations. \n',...
        ' Note: For each sample realization, the truck poses are obtained through numerical integration,\n',...
        ' therefore option 2) should take around 30 seconds.\n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end
         
tol = 1e-6; % tolerance selected for violating the constraints
runs  = 10000;

if selected == '1'
    
    load('Saved_data/Truck_poses')  % 10'000 truck_poses
    [Car_pos{1}, violations_LTL{1}, violations_input_state{1}] = Empirical_Violations( runs, tol, H_star{1}, h_star{1}, params, Truck_poses) ;
    
elseif selected == '2'
    w_veri = cell(1,runs);
    Truck_poses = cell(1,runs);
    for i = 1:runs
        w_veri{i}.theta = zeros(N+1,1);
        for j = 1:N
            w_veri{i}.theta(j) = pi/2/(N+1) - 0.03  +  min(0.06, pi/2-sum(w_veri{i}.theta)- pi/2/(N+1) + 0.03)*rand;
        end
        w_veri{i}.theta(N+1) = pi/2 - sum(w_veri{i}.theta);
        w_veri{i}.v = -22/3.6*ones(N+1,1);    % the truck has constant linear velocity at -22 Km/h
    
        pose_init = [truck.start; truck.theta0];
        [tout, pose_out] = ode45( @(tout, pose_out) truck_dyn(tout, pose_out, w_veri{i}.v, w_veri{i}.theta/Ts, Ts, N) , 0:Ts:N*Ts, pose_init);
        Truck_poses{i} = pose_out';
   end
  
    [Car_pos{1}, violations_LTL{1}, violations_input_state{1}] = Empirical_Violations( runs, tol, H_star{1}, h_star{1}, params, Truck_poses) ;
   
end

%% Plot trajectories (only 300 significant ones)
objects{1}.type = 'car';
objects{1}.d = params.car.d;
objects{1}.k = 0;
objects{2}.type = 'truck';
objects{2}.d = params.truck.d;
objects{2}.k = 0;
policyIdx = 1;

imgwidth = 3.75*0.73;
realizations = 1:length(Car_pos{policyIdx});
h = generateTrajectoriesPlot(struct('car', {Car_pos{policyIdx}(realizations)}, 'truck', {Truck_poses(realizations)}), params, params.N, ...
                          'objects', objects, ...
                          'imgwidth', imgwidth, ...
                          'tikz', true, 'matlab', true, ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
                          'name', 'turning_feedback');
all_states = cell2mat(Car_pos{1}(:));
last_pos = all_states(1:nx:end,end);  % 10'000 final horizontal positions of the car
if false % For tikz generation
    Nsubs = 300; %#ok
    alltruckstates = cell2mat(Truck_poses(:));
    lasttruckpos = alltruckstates(2:3:end,end);
    range = linspace(min(lasttruckpos),max(lasttruckpos),Nsubs-2);
    realizations = 1:length(range);
    [Y, I] = sort(lasttruckpos);
    for i=1:length(range)
        realizations(i) = I(find(Y >= range(i), 1, 'first'));
    end
    realizations = [realizations, find(last_pos == min(last_pos)) , find(last_pos == max(last_pos))];
    h = generateTrajectoriesPlot(struct('car', {Car_pos{policyIdx}(realizations)}, 'truck', {Truck_poses(realizations)}), params, params.N, ...
                          'objects', objects, ...
                          'imgwidth', imgwidth, ...
                          'tikz', true, 'matlab', true, ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
                          'name', 'turning_feedback');
    realizations = round(linspace(1,length(Car_pos{policyIdx}),Nsubs));
    h = generateTrajectoriesPlot(struct('car', {Car_pos{policyIdx}(realizations)}, 'truck', {Truck_poses(realizations)}), params, params.N, ...
                          'objects', objects, ...
                          'imgwidth', imgwidth, ...
                          'tikz', true, 'matlab', true, ...
                          'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
                          'name', 'turning_feedback');
end

                      
%% compute maximum and minimum horizontal position
maximum_horizontal_position = max(last_pos);
minimum_horizontal_position = min(last_pos);

title({'Figure 2.a) Reactions to 10''000 new random truck trajectories.',...
    ['Terminal position of the car lies in the interval [', num2str(minimum_horizontal_position),' m, ', num2str(maximum_horizontal_position),' m] ']});

fprintf(['-------------------------------------------------------------------------\n' ...
    'Figure 2.a) shows the different car trajectories using the policy for new realizations of the uncertainty.\n']);
fprintf('Saving the figure ...\n')
savefig(h,'Figures/Figure_2a.fig','compact')

%% Solve 14 new instances of the scenario program 
fprintf(['The obtained empirical violations are ', num2str(violations_LTL{1}/runs) ,' and ', num2str(violations_input_state{1}/runs),...
         ' for LTL specification and input-state constraints, respectively.\n\n']);
fprintf(['-------------------------------------------------------------------------\n' ...
         'To obtain a rough estimate of the distribution of the empirical violations, ' ...
         'we can solve 14 additional instances of the Scenario Program\n',...
         'and evaluate each solution against 10''000 new realizations of the uncertainty.\n']);
ok = 0;
while ok == 0
    prompt = ['The following three options are available:\n', ...
        ' 1) Load the optimizers used in the paper.\n',...
        ' 2) Solve the Scenario Programs for the same 14 realizations of the multisample used in the paper. \n',...
        ' 3) Solve the Scenario Programs for 14 random realizations of the multisample.\n',...
        'Note: When option 2) or 3) are selected, each Scenario Program is solved in\n' ...
        '       approximately 1 hour on a 8 Gb machine at 2.8 GHz using CPLEX.\n',...
        '      Therefore option 2) and 3) take approximately 14 hours to complete.' ...
        '\nChoose (1,2 or 3):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2' || selected == '3'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end

simulations = 15;

if selected == '1' 
    fprintf('Loading optimizers...\n\n');
    load('Saved_data/optimizers_15_sims') % LOAD 1x15 cells 
    for i = 2:simulations
        H_star{i} = H_star_loaded{i}; %#ok
        h_star{i} = h_star_loaded{i}; %#ok
        DIAGNOSTIC{i} = DIAGNOSTIC_loaded{i}; %#ok
    end
elseif selected == '2' 
    load('Saved_data/w_train_15_sims')
    for i = 2:simulations
        fprintf(['\nSolving Scenario program #',num2str(i-1),'...\n'])
        [H_star{i} , h_star{i}, DIAGNOSTIC{i}] = Solve_Scenario_Program(  params ,w_train{i}, selected_solver); %#ok
    end
elseif selected == '3'
    for j=2:simulations
    % draw a new multisample
    fprintf('\nDrawing new multisample...\n');
    for k = 1:K
        w{k}.theta = zeros(N+1,1);
        for i = 1:N
            w{k}.theta(i) = pi/2/(N+1) - 0.03  +  min(0.06, pi/2-sum(w{k}.theta)- pi/2/(N+1) + 0.03)*rand;
        end
        w{k}.theta(N+1) = pi/2 - sum(w{k}.theta);
        w{k}.v = -22/3.6*ones(N+1,1) ;    % the truck has constant linear velocity at -22 Km/h
    end
    fprintf(['Solving Scenario program #',num2str(j-1),'...\n'])
    [H_star{j} , h_star{j}, DIAGNOSTIC{j}] = Solve_Scenario_Program( params, w , selected_solver); %#ok
    end
end

% Generate empirical violations
for i=2:simulations
    [Car_pos{i}, violations_LTL{i}, violations_input_state{i}] = Empirical_Violations(runs, tol, H_star{i}, h_star{i}, params, Truck_poses) ; %#ok
end

%%
average_empirical_eps_phi = mean(cell2mat(violations_LTL))/runs;
average_empirical_eps_s = mean(cell2mat(violations_input_state))/runs;

solver_times = zeros(1,simulations);
for i = 1:simulations
    solver_times(i) = DIAGNOSTIC{i}.solvertime;
end
avg_solvertime = mean(solver_times);

fprintf(['-------------------------------------------------------------------------\n' ...
    'The average solver time is: ',num2str(avg_solvertime),' seconds\n\n']);
fprintf('Evaluating each solution over 10''000 realizations of the uncertainty...\n');
fprintf(['The average empirical violation probability of the LTL spec. is: ',num2str(average_empirical_eps_phi),'\n',...
        'The average empirical violation probability of state-input constraints is: ',num2str(average_empirical_eps_s),...
        '\n\n']);
     
%% Plot empirical violations
[h,ax] = generateBoxPlot({cell2mat(violations_LTL)/runs, cell2mat(violations_input_state)/runs}, {'$\hat{\epsilon}_\varphi$', '$\hat{\epsilon}_s$'}, ...
                'bounds', [0.05, 0.3], 'boundLabels', {'$\epsilon_\varphi$', '$\epsilon_s$'}, ...
                'log', true, ...
                'tikz', false, 'matlab', true, ...
                'dims', [min([cell2mat(violations_LTL)/10000 cell2mat(violations_input_state)/10000])/2, 1], ...
                'imgwidth', 0.5, 'imgheight', 1, ...
                'name', 'turning_empiricalViolations');
title(ax, {'Figure 2.b) Empirical violation probabilities over 15 instances of the scenario program:',...
    '{\color{red}Median}, 25-th to 75-th {\color{blue}percentile}, guaranteed violation probability levels ({\color{gray}dotted})'});
ax.XAxis.TickLabelInterpreter = 'latex';

fprintf(['Figure 2.b) shows the distribution of the empirical violation probabilities over the 15 different instances.\n',...
         'As expected, they satisfy the chosen levels 0.05 and 0.3.\n',...
         '----------------------------------------------------------------------------------------------------------\n']);
savefig(h, 'Figures/Figure_2b.fig')
