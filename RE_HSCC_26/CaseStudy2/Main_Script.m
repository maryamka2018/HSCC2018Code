%% script header
clear all %#ok
close all

%% seed used to generate realizations
% [train, veri] = RandStream.create('mrg32k3a','NumStreams', 2, 'Seed',3);
%%
fprintf(['-------------------------------------------------------------------------\n',...
        'This Script reproduces the data related to \nCase Study 5.2: Safe overtake.\n']);
%% load problem parameters
fprintf('Loading problem parameters...\n');
run Problem_parameters

%% compute (through binary search) the number of samples K_w needed 
eps = 0.05;
beta = 1e-3;
% classic scenario
d  = 4*nw*(N+1);  % total number of 'free' continuous variables for boxes
m = 1;
const = exp(1) /( exp(1)-1);
K = round( const/eps*(log(m/beta) + d  - 1 ) ) ; % # of samples
% binary search
up = K;
low = 0;
while up-low > 1
    mid =  round( low + (up-low)/2 );
    if binocdf( d - 1  , mid, eps) > beta/m
        low = mid;
    else
        up = mid;
    end
end
K = up;
%%
fprintf('Problem parameters loaded.\n');
fprintf('\nTo quickly reproduce the results found in the paper, select option 1) for all the following questions.\n');

%% solver selection
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'Select one of the following solvers:\n', ...
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

%% Scenario Program for obtaining estimate of W 
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'Given samples of w, an estimate of the support set W is obtained using the Scenario Approach.\n'...
        'The following three options are available:\n', ...
        ' 1) Load the estimate W[K_w] used in the paper.\n',...
        ' 2) Solve the Scenario Program for the same realization of the multisample W^{K_w} used in the paper. \n',...
        ' 3) Solve the Scenario Program for a random realization of the multisample W^{K_w}.\n',... 
        'Note: If option 2) or 3) is selected, the Scenario Program is solved in\n' ...
        '       approximately 3 seconds on a 2.8 GHz machine with 4 cores and 8 GB of RAM, using CPLEX.\n',...
        '\nChoose (1,2 or 3):  '];
    str = input(prompt,'s');
    selected = str(1);
    if selected == '1' || selected == '2' || selected == '3'
        ok = 1;
    else
        fprintf('\nWrong character inserted: Please select a valid option.\n\n')
    end
end

%% %%  draw samples for Scenario Program, to identify W_1[K_w] and W_2[K_w]

if selected == '1' 
    load('Saved_data/estimated_W_bounds') % lower1_star , upper1_star, lower2_star, upper2_star
elseif selected == '2'
    load('Saved_data/w_train') % w
    fprintf('Building the Scenario Program...\n');
    [lower1_star , upper1_star, lower2_star, upper2_star, time_SP] = Find_boxes_with_Scenario_Approach(selected_solver, w , nw, N);
    fprintf(['Solved in ', num2str(time_SP), ' seconds.\n\n']);
elseif selected == '3'
    w = zeros(nw*(N+1),K);
    for i = 1:K
        if rand > 0.5
            w(:,i) = wmax_bold.*rand(nw*(N+1),1);
        else
            w(:,i) = wmin_bold.*rand(nw*(N+1),1);
        end
        w(1,i) = -offset_max + 2*offset_max*rand;
    end

    fprintf('Building the Scenario Program...\n');
    [lower1_star , upper1_star, lower2_star, upper2_star, time_SP] = Find_boxes_with_Scenario_Approach(selected_solver, w , nw, N);
    fprintf(['Solved in ', num2str(time_SP), ' seconds.\n\n']);
end

%% Solve the Robust Program, for different P (number of elements in partition)
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'We can now find an optimal input policy solving the Robust Program, using linear programming duality.\n'...
        'We choose u() to be piecewise affine, and \delta to be piecewise constant, defined over the same partition of P elements.\n',...
        'The following three options are available:\n', ...
        ' 1) Load the optimal policies used in the paper, for P = 1, 2, 3, 5, 9.\n',...
        ' 2) Solve five Robust Programs to obtain policies for P = 1, 2, 3, 5, 9. \n',...
        ' 3) Solve a sequence of Robust Programs to obtain policies for P = 1, 2, 3, 5, 9, ..., P_max.\n',...
        '    The number N of Robust Programs is selected by the user. Then P_max = 1 + 2^N.\n',... 
        'Note: Options 2) and 3) consist in solving a sequence of Robust Programs, while iteratively partitioning the support W[K_w].\n',...
        '      Option 2) requires approximately 5 minutes on a 2.8 GHz machine with 4 cores and 8 GB of RAM, using CPLEX.\n',...
        '      Option 3) may require substantially more time than option 2) if N > 5.\n' ...
        'After the policies are loaded or computed, the user can select one of these policies\n' ...
        'to reproduce the results of Figure 3) and 4).\n',...
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
    load('Saved_data/optimal_policies') % Optimizers, Solvertime, Time_splitting, W, v, g, l
    P_max = 9; 
    set_P = [1,2,3,5,9];
    fprintf('Optimal policies loaded.\n');
    
elseif selected == '2' || selected == '3'
    if selected == '2'
        P_max = 9; % by default
        set_P = [1,2,3,5,9];
    elseif selected == '3'
        ok = 0;
        while ok == 0
            prompt = 'Choose an integer N >= 4 : ';
            str = input(prompt,'s');
            x = str;
            if str2double(x) >= 4
                P_max = 1 + 2^str2double(x);
                ok = 1;
            else
                fprintf('\nPlease select an integer N >= 4.\n\n')
            end
        end
        set_P = [1, 2, 1+2.^(1:str2double(x))];
    end
    
% initialize empty splitting hyperplanes for P = 1 and P = 2
g{1}{1} = [];
l{1}{1} = [];
g{2}{1} = [];
g{2}{2} = [];
l{2}{1} = [];
l{2}{2} = [];
for P = set_P
    if P == 1
        wmin_bold_new{1} = lower1_star;
        wmax_bold_new{1} = upper2_star;
    else
        % Partition 1. That is, Truck goes down (trivial manouver)
        wmin_bold_new{1} = lower1_star;
        wmax_bold_new{1} = upper1_star;
        % Partition 2. That is, Truck goes up (more complex overtake)
        %% partition further partition 2 (initial position of the truck)
        for p =2:P
            wmin_bold_new{p} = lower2_star; %#ok
            wmax_bold_new{p} = upper2_star; %#ok
        end
    end
    %% create polytopic regions  W{p}*w <= v{p}
    W = cell(1,P);
    v = cell(1,P);
    for p=1:P
        W{p} = [ eye(nw*(N+1)); -eye(nw*(N+1)) ; g{P}{p}'];
        v{p} = [wmax_bold_new{p}; - wmin_bold_new{p}; l{P}{p}];
    end
    
    %% plot partitions
    if false
        figure(); %#ok
        for i=2:P
            union = PolyUnion();
            pol{1} = Polyhedron('A', W{1}, 'b', v{1});
            union.add(projection(pol{1},1:2));
            pol{i} = Polyhedron('A', W{i}, 'b', v{i});
            union.add(projection(pol{i},1:2));
            subplot((P-1)/2,(P-1)/2,i-1)
            plot(union);
        end
    end
    
    %% Compute approximate (relative) areas of polytopes
    weight = cell(1,P);
    samples = cell(1,P);
    for p=1:P
        weight{p} = 0;
        samples{p} = [];
    end
    load('Saved_data/samples_for_cost_function') % w_realiz
    for i = 1:size(w_realiz,2)
        for p=1:P
            if all(W{p}*w_realiz(:,i) <= v{p})
                weight{p} = weight{p}+1;
                samples{p} = [ samples{p}, w_realiz(:,i)];
                break
            end
        end
    end
    tot = sum(cell2mat(weight(:)));
    for p=1:P
        weight{p} = weight{p}/tot;
    end
    fprintf(['-------------------------------------------------------------------------\n' ...
        'Building Robust Program for P = ',num2str(P),' ...\n']);
    [Optimizers{P}, DIAGNOSTIC{P}, Solvertime{P}] = Solve_Robust_Program(selected_solver, P, W, v, weight, samples, params); %#ok
    fprintf(['Solved in ',num2str(Solvertime{P}),' seconds.\n']);

    if P > 1 && P < P_max
        fprintf(['-------------------------------------------------------------------------\n' ...
            'Find partition for P = ',num2str(set_P(find(set_P==P) + 1)),' ...\n']);
        [g, l, Time_splitting{P}] = Update_splitting_hyperplanes (selected_solver, g, l, Optimizers{P}, P, W, v, params); %#ok
        fprintf('Found.\n');

    end
end

end

%% Choose one of the found policies and plot Figure 3 and 4
fprintf(['-------------------------------------------------------------------------\n' ...
    'One of the found policies is now chosen and evaluated against a new realization of w, shown in Figure 3),\n',...
        'and against 10''000 new realizations, shown in Figure 4).\n']);
ok = 0;
while ok == 0
    prompt = ['\nThe following options are available: \n',...
              ' 1) Choose policy for P = 5, as it was done in the paper.\n'...
              ' 2) Choose a different policy.\n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    sel = str;
    if sel == '1' || sel == '2'
        ok = 1;
    else
        fprintf('\nPlease select a valid option.\n\n');
    end
end

if sel == '1' 
    P_selected = 5;
    
elseif sel == '2'
    ok = 0;
    while ok == 0
        prompt = ['\nPlease choose one the the found policies, for P = 1, 2, 3, 5, ... ', num2str(P_max),'.\n',...
            '\nChoose P:  '];
        str = input(prompt,'s');
        P_selected = str2double(str);
        if ismember(P_selected, set_P)
            ok = 1;
        else
            fprintf(['\nPlease select a P in,' num2str(set_P),'.\n\n']);
        end
    end
end

%% redefine polythopes to make them a proper partition of R^{n_w*(N+1)}
%% (for simplicity we make them just larger hyperrectangular, since we know that which values w can take.)
if P_selected == 1
    W_larger{1} = [ eye(nw*(N+1)); -eye(nw*(N+1))];
    v_larger{1} = [wmax_bold; - wmin_bold];
else
    % Partition 1. That is, Truck goes down (trivial manouver)
    W_larger = cell(1,P_selected);
    v_larger = cell(1,P_selected);
    W_larger{1} = [ eye(nw*(N+1)); -eye(nw*(N+1))];
    v_larger{1} = [wmax_bold(1) ;zeros(nw*(N+1)-1,1); - wmin_bold];
    % Partition 2. That is, Truck goes up (more complex overtake)
    for p=2:P_selected
        W_larger{p} = [eye(nw*(N+1)); -eye(nw*(N+1)) ; g{P_selected}{p}'];
        v_larger{p} = [wmax_bold; - wmin_bold(1) ;zeros(nw*(N+1)-1,1); l{P_selected}{p}];
    end
end

params.P = P_selected;
params.W = W_larger;
params.v = v_larger;

%% Car performing the task for a sampled w
%% draw samples for verification
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'The chosen policy will be tested against a new realization of w.\n',...
              ' Then, Figure 3) will be plotted.\n',...
              ' The following options are available: \n',...
              ' 1) Choose the realization of w used in the paper.\n'...
              ' 2) Choose a random realization of w.\n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    selected = str;
    if selected == '1' || selected == '2'
        ok = 1;
    else
        fprintf('\nPlease select a valid option.\n\n');
    end
end

if selected == '1'
    load('Saved_data/w_verification')
elseif selected == '2'
    w_veri = zeros(nw*(N+1),1);
    for i = 1
        if rand > 0.5  % truck up
            w_veri(:,i) =  wmax_bold.*rand(nw*(N+1),1);
        else  % truck down
            w_veri(:,i) =  wmin_bold.*rand(nw*(N+1),1);
        end
        w_veri(1,i) = - offset_max + 2*offset_max*rand;
    end
end

w_realiz = w_veri(:,1);
H_star = Optimizers{P_selected}.H;
h_star = Optimizers{P_selected}.h;
[ car_states, truck_states ] = task_for_new_realization( H_star, h_star, params, w_realiz);

%% Plot Figure 3
h = figure();
clear objects;
objects{1}.type = 'car';
objects{1}.d = params.car.d;
objects{2}.type = 'truck';
objects{2}.d = params.truck.d;
idx = 1;
imgwidth = 2.5;

% k = 3
k = 3;
objects{1}.k = k;
objects{2}.k = k;
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
    'case', 'overtake', ...
    'objects', objects, ...
    'highlight', 1, ...
    'imgwidth', imgwidth, ...
    'h', h, ...
    'tikz', false, 'matlab', false, ...
    'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
    'name', ['overtake_traj_' num2str(k)]);

k = 5;
objects{1}.k = k;
objects{2}.k = k;
subplot(3,1,1);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
    'case', 'overtake', ...
    'objects', objects, ...
    'highlight', 1, ...
    'imgwidth', imgwidth, ...
    'h', h, ...
    'tikz', false, 'matlab', true, ...
    'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
    'name', ['overtake_traj_' num2str(k)]);
title({'Figure 3): The controlled car ({\color{blue}solid}) reacting to a new truck trajectory ({\color{red}dashed})', 'using the chosen feedback policy', ' ', [' time k = ',num2str(k)]});

k = 7;
objects{1}.k = k;
objects{2}.k = k;
subplot(3,1,2);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
    'case', 'overtake', ...
    'objects', objects, ...
    'highlight', 1, ...
    'imgwidth', imgwidth, ...
    'h', h, ...
    'ylabel', false, ...
    'tikz', false, 'matlab', true, ...
    'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
    'name', ['overtake_traj_' num2str(k)]);
title([' time k = ',num2str(k)]);

k = 9;
objects{1}.k = k;
objects{2}.k = k;
subplot(3,1,3);
generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
    'case', 'overtake', ...
    'objects', objects, ...
    'highlight', 1, ...
    'imgwidth', imgwidth, ...
    'h', h, ...
    'ylabel', false, ...
    'tikz', false, 'matlab', true, ...
    'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
    'name', ['overtake_traj_' num2str(k)]);
title([' time k = ',num2str(k)]);

fprintf('Figure 3) shows the controlled car reacting to a new trajectory of the truck. \n');
savefig('Figures/Figure_3.fig')
%% TIKZ generation
if false
    k = 5; %#ok
    objects{1}.k = k;
    objects{2}.k = k;
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
        'case', 'overtake', ...
        'objects', objects, ...
        'highlight', 1, ...
        'imgwidth', imgwidth, ...
        'tikz', true, 'matlab', false, ...
        'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
        'name', ['overtake_traj_' num2str(k)]);

    k = 7;
    objects{1}.k = k;
    objects{2}.k = k;
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
        'case', 'overtake', ...
        'objects', objects, ...
        'highlight', 1, ...
        'imgwidth', imgwidth, ...
        'ylabel', false, ...
        'tikz', true, 'matlab', false, ...
        'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
        'name', ['overtake_traj_' num2str(k)]);

    k = 9;
    objects{1}.k = k;
    objects{2}.k = k;
    generateTrajectoriesPlot(struct('car', {{car_states}}, 'truck', {{truck_states}}), params, k, ...
        'case', 'overtake', ...
        'objects', objects, ...
        'highlight', 1, ...
        'imgwidth', imgwidth, ...
        'ylabel', false, ...
        'tikz', true, 'matlab', false, ...
        'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
        'name', ['overtake_traj_' num2str(k)]);
end

%% Empirical violations
ok = 0;
while ok == 0
    prompt = ['-------------------------------------------------------------------------\n' ...
        'The chosen policy will be tested against 10''000 new realizations of w.\n',...
              ' Then, Figure 4) will be plotted.\n',...
              ' The following options are available: \n',...
              ' 1) Choose the same 10''000 realizations of w used in the paper.\n'...
              ' 2) Choose 10''000 random realizations of w.\n',...
        '\nChoose (1 or 2):  '];
    str = input(prompt,'s');
    selected = str;
    if selected == '1' || selected == '2'
        ok = 1;
    else
        fprintf('\nPlease select a valid option.\n\n');
    end
end


runs = 10000;
if selected == '1'
    load('Saved_data/w_verification')
elseif selected == '2'
    w_veri = zeros(nw*(N+1),10000);
    for i = 1:10000
        if rand > 0.5  % truck up
            w_veri(:,i) =  wmax_bold.*rand(nw*(N+1),1);
        else  % truck down
            w_veri(:,i) =  wmin_bold.*rand(nw*(N+1),1);
        end
        w_veri(1,i) = - offset_max + 2*offset_max*rand;
    end
end

fprintf('Testing the policy ...\n');

w_realiz = w_veri(:,1:runs);
tol = 1e-6;
[Truck_states, Car_states, Inputs, used_part, violations] = Empirical_Violations( H_star, h_star, params, w_realiz, tol);
%% Plot Figure 4
clear objects;
policyIdx = P_selected;
realizations = 1:10000;
imgwidth = 3.75*0.7;
h = figure();
colors = parula(policyIdx);
tikzcolors = cell(1,policyIdx);
for p=1:policyIdx
    tikzcolors{p} = ['color' num2str(p)];
end
data = struct();
objects = cell(1,policyIdx);
for p=1:policyIdx
    objects{p}.type = ['car', num2str(p)];
    objects{p}.d = params.car.d;
    objects{p}.k = [];
    objects{p}.color = colors(p,:);
    objects{p}.tikzcolor = tikzcolors{p};
    Car_states_truncated = Car_states(realizations);
    data.(objects{p}.type) = Car_states_truncated(cell2mat(used_part(realizations))==p);
end
generateTrajectoriesPlot(data, params, params.N, ...
        'case', 'overtake', ...
        'objects', objects, ...
        'highlight', [], ...
        'dims', [0 35 -2*params.lane-1 params.lane+1], ...
        'imgwidth', imgwidth, ...
        'h', h, ...
        'tikz', false, 'matlab', true, ...
        'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
        'name', 'overtake_feedback');

title({'Figure 4 a): Car trajectories for 10''000 new truck trajectories.',...
    'The colors correspond to the different pieces of the policies.'});

fprintf('Figure 4.a) the shows the trajectories of the car for 10''000 new truck trajectories.\n');
fprintf('Saving the figure ...\n')
savefig(h, 'Figures/Figure_4a.fig', 'compact')

if false % For tikz generation
    extremeBox = [inf inf -inf -inf]; %#ok
    for p=1:policyIdx
        all_states = cell2mat(Car_states(cell2mat(used_part(1:10000))==p));
        last_pos = all_states(2,N+1:N+1:end);  % 10'000 final horizontal positions of the car
        Nsubs = 300;
        data = struct();
        objects = cell(1,1);
        objects{1}.type = ['car', num2str(p)];
        objects{1}.d = params.car.d;
        objects{1}.k = [];
        objects{1}.color = colors(p,:);
        objects{1}.tikzcolor = tikzcolors{p};
        Car_states_truncated = Car_states(cell2mat(used_part(1:10000))==p);
        allcarstates = cell2mat(Car_states_truncated);
        lastcarpos = allcarstates(2,N+1:N+1:end);
        range = linspace(min(lastcarpos),max(lastcarpos),Nsubs-2);
        realizations = 1:length(range);
        [Y, I] = sort(lastcarpos);
        for i=1:length(range)
            realizations(i) = I(find(Y >= range(i), 1, 'first'));
        end
        realizations = [realizations, find(last_pos == min(last_pos), 1, 'first') , find(last_pos == max(last_pos), 1, 'first')];
        data.(objects{1}.type) = Car_states_truncated(unique(realizations));
        objects{1}.k = 1;
        figure();
        generateTrajectoriesPlot(data, params, params.N, ...
            'case', 'overtake', ...
            'objects', objects, ...
            'highlight', [], ...
            'dims', [0 35 -2*params.lane-1 params.lane+1], ...
            'imgwidth', imgwidth, ...
            'h', h, ...
            'tikz', false, 'matlab', true, ...
            'axisOptions', sprintf('ylabel style={yshift=-%0.4gin}',imgwidth*0.035), ...
            'name', 'overtake_feedback');
        Car_states_truncated = cell2mat(Car_states_truncated(unique(realizations)));
        extremeBox = [min(extremeBox(1:2), min(Car_states_truncated(1:2,:),[],2)') max(extremeBox(3:4), max(Car_states_truncated(1:2,:),[],2)')];
    end
    extremeBox
end

%% Plot partition
figure()
clear P;
hold on 
for p = 1:policyIdx 
  v_transformed = v_larger{p} + W_larger{p}*repmat([truck.start + truck.d(1)/2;middle], N+1,1);
  Plot2DPolytope(W_larger{p}(:,1:2) , v_transformed , colors(p,:));
end
xlabel('w_{1,0}')
ylabel('w_{2,0}')
title('Figure 4 b): Partition of w_0 for the piecewise-affine policy');
fprintf('Figure 4.b) shows the partition of w_0 characterizing the piecewise-affine policy. \n');
fprintf('Saving the figure ...\n')
savefig('Figures/Figure_4b.fig')
%% Tikz generation
if false
    for p=1:policyIdx %#ok
        v_transformed = v_larger{p} + W_larger{p}*repmat([truck.start + truck.d(1)/2;middle], N+1,1);
        P{p} = projection(Polyhedron('A', W_larger{p}, 'b', v_transformed), 1:2);
        P{p}.minVRep();
    end
    imgwidth = 3.75*0.25;
    tikz = '';
    for p=1:policyIdx
        tikz = [tikz, sprintf(['\\definecolor{' tikzcolors{p} '}{rgb}{%0.3g,%0.3g,%0.3g}' '\n'], colors(p,:))];
    end
    tikz = [tikz, sprintf('\\begin{tikzpicture}[')];
    tikz = [tikz, sprintf(['scale=1]' '\n'])];
    tikz = [tikz, sprintf(['  ' '\\begin{axis}[' '\n'])];
    tikz = [tikz, sprintf(['    ' 'every outer x axis line/.append style={white!20!black},' '\n'])];
    tikz = [tikz, sprintf(['    ' 'every x tick label/.append style={font=\\color{white!20!black}},' '\n'])];
    tikz = [tikz, sprintf(['    ' 'xmin=%0.4g,' '\n'], -2.2+ 15.5)];
    tikz = [tikz, sprintf(['    ' 'xmax=%0.4g,' '\n'], 2.2 + 15.5)];
    tikz = [tikz, sprintf(['    ' 'ymin=%0.4g,' '\n'], -0.45 -1.75)];
    tikz = [tikz, sprintf(['    ' 'ymax=%0.4g,' '\n'], 0.45 -1.75)];
    tikz = [tikz, sprintf(['    ' 'xlabel={$w_{1,0}$},' '\n'])];
    tikz = [tikz, sprintf(['    ' 'ylabel={$w_{2,0}$},' '\n'])];
    tikz = [tikz, sprintf(['    ' 'axis background/.style={fill=white},' '\n'])];
    tikz = [tikz, sprintf(['    ' 'axis x line*=bottom,' '\n'])];
    tikz = [tikz, sprintf(['    ' 'axis y line=left,' '\n'])];
    tikz = [tikz, sprintf(['    ' 'ylabel style={yshift=-%0.4gin},' '\n'],imgwidth*0.035)];
    tikz = [tikz, sprintf(['    ' 'width = %0.4gin, height = %0.4gin, scale only axis]' '\n'], imgwidth, imgwidth)];
    for p=1:policyIdx
        tikz = [tikz, sprintf(['  ' '\\fill[' tikzcolors{p} '] '])];
        V = P{p}.V;
        if p==1 % HACK
            V = V([1 2 4 3],:);
        end
        for i=1:size(P{p}.V,1)

            tikz = [tikz, sprintf('(axis cs:%0.4g,%0.4g) -- ', V(i,1), V(i,2))];
        end
        tikz = [tikz, sprintf(['cycle;' '\n'])];
    end
    tikz = [tikz, sprintf(['  ' '\\end{axis}' '\n'])];
    tikz = [tikz, sprintf(['\\end{tikzpicture}' '\n'])];
    f = fopen('partition.tikz', 'w');
    fprintf(f, '%s', tikz);
    fclose(f);
end

%% plot different trajectories of the car
% num = 1000;  % less than length(States_multiple)
% plot_multiple_trajectories(num, Car_states,used_part, params)

all_States = cell2mat(Car_states(:));  % 40000 x (N+1)
last_positions = all_States(1:nx:end,end);   % 10000 x 1

min_last_pos_1 = min(last_positions(cell2mat(used_part)==1));
max_last_pos_1 = max(last_positions(cell2mat(used_part)==1));
if P_selected > 1
min_last_pos_2 = min(last_positions(cell2mat(used_part)~=1));
max_last_pos_2 = max(last_positions(cell2mat(used_part)~=1));
else
min_last_pos_2 = min_last_pos_1;
max_last_pos_2 =  max_last_pos_1;
end

fprintf(['When the truck moves right:\n',...
         ' Terminal horizontal position of the car varies between ',num2str(min_last_pos_1),' m, and ',num2str(max_last_pos_1) ' m.\n']);
fprintf(['When the truck moves left:\n',...
         ' Terminal horizontal position of the car varies between ',num2str(min_last_pos_2),' m, and ',num2str(max_last_pos_2) ' m.\n']);
fprintf(['The constraints were violated ',num2str(violations), ' times, leading to\n',...
    ' Empirical violation probability = ', num2str(violations/10000),'.\n']);

%% Iterate over all the different P
fprintf(['-------------------------------------------------------------------------\n' ...
    'To illustrate the effect of partitioning, the found policies for P = [', num2str(set_P),'] are now evaluated against\n' ...
    '10''000 realizations of w, where the truck moves left.\n' ...
    'Then, Figure 5) will be plotted.\n'])
%% draw 10'000 realizations of w, for the truck moving left
ok = 0;
while ok == 0
prompt = ['\nThe following options are available: \n',...
          ' 1) Choose the same 10''000 realizations of w used in the paper.\n'...
          ' 2) Choose 10''000 random realizations of w.\n',...
    '\nChoose (1 or 2):  '];
str = input(prompt,'s');
selected = str;
if selected == '1' || selected == '2'
    ok = 1;
else
    fprintf('\nPlease select a valid option.\n\n');
end
end

if selected == '1' 
    load('Saved_data/w_verification_left') % w_realiz
elseif selected == '2'    
w_realiz = zeros(nw*(N+1),10000);
for i = 1:10000
    if true  % truck up
        w_realiz(:,i) = wmax_bold.*rand(nw*(N+1),1);
    else  % truck down
        w_realiz(:,i) = wmin_bold.*rand(nw*(N+1),1); %#ok
    end
    w_realiz(1,i) = - offset_max + 2*offset_max*rand;
end
end

fprintf('Testing the policies ...\n');

for P = set_P
    if P == 1
        W_larger{1} = [ eye(nw*(N+1)); -eye(nw*(N+1))];
        v_larger{1} = [wmax_bold; - wmin_bold];
    else
        W_larger = cell(1,P);
        v_larger = cell(1,P);
        % Partition 1. That is, Truck goes down (trivial manouver)
        W_larger{1} = [ eye(nw*(N+1)); -eye(nw*(N+1))];
        v_larger{1} = [wmax_bold(1) ;zeros(nw*(N+1)-1,1); - wmin_bold];
        % Partition 2. That is, Truck goes up (more complex overtake)
        for p=2:P
            W_larger{p} = [ eye(nw*(N+1)); -eye(nw*(N+1)) ; g{P}{p}'];
            v_larger{p} = [wmax_bold; - wmin_bold(1) ;zeros(nw*(N+1)-1,1); l{P}{p}];
        end
    end
    
    params.P = P;
    params.W = W_larger;
    params.v = v_larger;
    
    H_star = Optimizers{P}.H;
    h_star = Optimizers{P}.h;
    
    runs = 10000;
    tol = 1e-6;
    [~, Car_states_temp, ~, ~, ~] = Empirical_Violations( H_star, h_star, params, w_realiz, tol);
    
    all_States_temp = cell2mat(Car_states_temp(:));  % 40000 x (N+1)
    last_pos{P} = all_States_temp(1:nx:end,end);   %#ok 10000 x 1
end


%%
%% Plot objective (overtaking)
clear data;
for i=1:length(set_P)
    data{i} = last_pos{set_P(i)}(last_pos{set_P(i)} <= 30); %#ok
end
[h,ax] = generateBoxPlot(data, arrayfun(@(x)(num2str(x)), set_P, 'UniformOutput', false), ...
    'percentiles', [5 95], ...
    'connected', 'mean', ...
    'median', false, ...
    'xlabel', '# of elements of partition', ...
    'ylabel', '$x_{1,N}$', ...
    'log', false, ...
    'tikz', false, 'matlab', true, ...
    'dims', [6 25], ...
    'imgwidth', 0.25*length(set_P), 'imgheight', 1, ...
    'name', 'overtake_objective');
title(ax, {'Figure 5: Distribution of terminal forward position of the car',...
    'for 10''000 new trajectories of the truck moving left.',...
    'The {\color{red}mean} and 5-th to 95-th {\color{blue}percentile}.'});
ax.YLabel.Interpreter = 'latex';

savefig(h, 'Figures/Figure_5.fig')

fprintf(['Figure 5) the shows distribution of terminal forward position of the car\n',...
    'for 10''000 new trajectories of the truck moving left.\n']);
fprintf('-------------------------------------------------------------------------------------------------------------------------------\n')



