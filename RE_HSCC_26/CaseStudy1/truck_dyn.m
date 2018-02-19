function [ dposedt ] = truck_dyn(t, pose, v,w , Ts, N)
idx_vec = Ts*[0:1:N];
idx = find(t>= idx_vec , 1, 'last' );

dposedt = [v(idx)* cos(pose(3));...
           v(idx)* sin(pose(3));...
           w(idx)  ];
end

