function [ output_args ] = Plot2DPolytope( A, b , color )

idx = any(A,2);
A = A( idx , :);
b = b( idx);


Intersection_points = [];
for i = 1:size(A,1)
    for j=i+1:size(A,1)
        %       if A(i,1) == 0 && A(j,1) ~= 0
        %           x2 = b(i) / A(i,2);
        %           x1 = b(j)/A(j,1) - A(j,2)/A(j,1)*x2;
        %           Intersection_points = [Intersection_points ;  [x1 , x2] ];
        %
        %       elseif A(i,2) == 0 && A(j,2) ~= 0
        %           x1 = b(i) / A(i,1);
        %           x2 = b(j)/A(j,2) - A(j,1)/A(j,2)*x1;
        %           Intersection_points = [Intersection_points ;  [x1 , x2] ];
        %
        %       elseif A(i,1) ~= 0 && A(j,1) == 0
        %           x2 = b(j)/A(j,2);
        %           x1 = b(i)/A(i,1) - A(i,2)/A(i,1)*x2;
        %           Intersection_points = [Intersection_points ;  [x1 , x2] ];
        %
        %       elseif A(i,2) ~= 0 && A(j,2) == 0
        %           x1 = b(j) / A(j,1);
        %           x2 = b(i)/A(i,2) - A(i,1)/A(i,2)*x1;
        %           Intersection_points = [Intersection_points ;  [x1 , x2] ];
        %
        %       elseif A(i,1) ~= 0 && A(i,2) ~= 0 && A(j,1) ~= 0 && A(j,2) ~= 0
        %       denumerator = A(i,2)/(A(i,1)) - A(j,2)/(A(j,1)) + 1e-6;
        %       x2 = (b(i)/(A(i,1)) - b(j)/(A(j,1))) / denumerator;
        %       x1 = b(i)/(A(i,1)) - A(i,2)/(A(i,1))*x2;
        %       Intersection_points = [Intersection_points ;  [x1 , x2] ];
        %       end
        if rank(A([i j],:)) >= 2
            Intersection_points = [Intersection_points ;  (A([i j],:)\b([i j]))']; %#ok
        end
    end
end

Inner_points = [];
for i = 1:size(Intersection_points,1)
   if A*Intersection_points(i,:)' <=  b + 1e-6 
       Inner_points = [Inner_points; Intersection_points(i,:)];
   end
end
Inner_points  = unique(Inner_points, 'rows');

k = convhull(Inner_points);
V_ordered = Inner_points(k,:);

patch(V_ordered(:,1), V_ordered(:,2), color);
end

