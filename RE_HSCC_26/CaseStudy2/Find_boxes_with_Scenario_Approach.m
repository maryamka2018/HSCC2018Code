function [lower1_star , upper1_star, lower2_star, upper2_star, time] =  Find_boxes_with_Scenario_Approach(selected_solver, w, nw, N)


%%
% Clustering, or simply looking at the samples, 
% we decide to identify W as union of two disjoint boxes ...

%%
upper1 = sdpvar(nw*(N+1),1,'full');
lower1 = sdpvar(nw*(N+1),1,'full');
upper2 = sdpvar(nw*(N+1),1,'full');
lower2 = sdpvar(nw*(N+1),1,'full');

objective = sum(upper1-lower1) + sum(upper2 - lower2);

Constraints = [];
for k=1:size(w,2)
    if w(2,k) < 0
        Constraints = [Constraints, lower1 <= w(:,k) <= upper1];
    else
        Constraints = [Constraints, lower2 <= w(:,k) <= upper2];
    end
end

%% solve the problem 
ops = sdpsettings;
ops.solver = selected_solver;
ops.verbose = 0;
ops.debug = 1;
DIAGNOSTIC = optimize(Constraints, objective, ops);

upper1_star = value(upper1);
lower1_star = value(lower1);
upper2_star = value(upper2);
lower2_star = value(lower2);
time = DIAGNOSTIC.solvertime;
end
















