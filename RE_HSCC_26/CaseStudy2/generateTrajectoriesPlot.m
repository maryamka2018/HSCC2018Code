function h = generateTrajectoriesPlot(varargin)
    cases = {'turning', 'overtake'};

    p = inputParser;
    p.addRequired('trajectories');
    p.addRequired('params'); % Parameter struct
    p.addRequired('N'); % Planning horizon
    p.addParameter('objects', {});
    p.addParameter('case', cases{1}, @(s)(ischar(s) && any(strcmp(s, cases))));
    p.addParameter('highlight', []); 
    p.addParameter('collisionBoxes', false, @islogical); % Enable collision boxes in tikz
    p.addParameter('background', true, @islogical);
    p.addParameter('imgwidth', 4.5); % In inches
    p.addParameter('matlab', true, @islogical); % Enable generation of MATLAB figure
    p.addParameter('h', []);
    p.addParameter('tikz', false, @islogical); % Enable generation of tikz figure
    p.addParameter('spy', []);
    p.addParameter('legend', false, @islogical); % Enable labelling
    p.addParameter('dims', []);
    p.addParameter('xlabel', true, @islogical);
    p.addParameter('ylabel', true, @islogical);
    p.addParameter('axisOptions', '', @(c)(ischar(c) || (iscell(c) && all(cellfun(@ischar, c)))));
    p.addParameter('name', 'test', @ischar);
    p.parse(varargin{:});
    options = p.Results;
    trajectories = options.trajectories;
    params = options.params;
    N = options.N;
    objects = options.objects;
    
    %% Define some colors
    colors.road = [1 1 1];
    colors.offroad = [0.9, 0.9, 0.9];
    colordefines = '';
    for j=1:length(objects)
        if ~isfield(objects{j}, 'color')
            if ~isempty(strfind(objects{j}.type, 'car'))
                colors.(objects{j}.type) = [0 0 1];
            elseif ~isempty(strfind(objects{j}.type, 'truck'))
                colors.(objects{j}.type) = [1 0 0];
            end
        else
            colors.(objects{j}.type) = objects{j}.color;
        end
        if ~isfield(objects{j}, 'tikzcolor')
            if ~isempty(strfind(objects{j}.type, 'car'))
                tikzcolors.(objects{j}.type) = 'blue';
            elseif ~isempty(strfind(objects{j}.type, 'truck'))
                tikzcolors.(objects{j}.type) = 'red';
            end
        else
            tikzcolors.(objects{j}.type) = objects{j}.tikzcolor;
            colordefines = [colordefines, sprintf(['\\definecolor{' tikzcolors.(objects{j}.type) '}{rgb}{%0.3g,%0.3g,%0.3g}' '\n'], colors.(objects{j}.type))];
        end
    end
    markers.car = 'x';
    markers.truck = '*';
    linestyle.car = '-';
    linestyle.truck = '--';
    tikzmarkers.car = 'x';
    tikzmarkers.truck = 'star';
    
    % Figure dimensions
    if isempty(options.dims)
        switch options.case
            case 'turning'
                options.dims = [0 50 -15 params.lane+1];
            case 'overtake'
                options.dims = [0 30 -2*params.lane-1 params.lane+1];
        end
    end
    imgwidth = options.imgwidth;
    
    if options.matlab && isempty(options.h)
        h = figure();
    end
    if options.tikz
        tikz = '';
        if ~isempty(colordefines)
            tikz = [tikz, colordefines];
        end
        tikz = [tikz, sprintf('\\begin{tikzpicture}[')];
        if ~isempty(options.spy)
            tikz = [tikz, sprintf(['spy using outlines={circle,' options.spy.color ',magnification=%0.4g,connect spies},'], options.spy.magnification)];
        end
        tikz = [tikz, sprintf(['scale=1,cap=round]' '\n'])];
        tikz = [tikz, sprintf(['  ' '\\tikzstyle{dashed}=[dash pattern=on 4pt off 8pt]' '\n'])];
        tikz = [tikz, sprintf(['  ' '\\begin{axis}[' '\n'])];
        tikz = [tikz, sprintf(['    ' 'every outer x axis line/.append style={white!20!black},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'every x tick label/.append style={font=\\color{white!20!black}},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'xmin=%0.4g,' '\n'], options.dims(1))];
        tikz = [tikz, sprintf(['    ' 'xmax=%0.4g,' '\n'], options.dims(2))];
        tikz = [tikz, sprintf(['    ' 'ymin=%0.4g,' '\n'], options.dims(3))];
        tikz = [tikz, sprintf(['    ' 'ymax=%0.4g,' '\n'], options.dims(4))];
        if options.xlabel
            tikz = [tikz, sprintf(['    ' 'xlabel={$x_1$},' '\n'])];
        else
            tikz = [tikz, sprintf(['    ' 'xticklabels={,,},' '\n'])];
        end
        if options.ylabel
            tikz = [tikz, sprintf(['    ' 'ylabel={$x_2$},' '\n'])];
        else
            tikz = [tikz, sprintf(['    ' 'yticklabels={,,},' '\n'])];
        end
        tikz = [tikz, sprintf(['    ' 'axis equal image,' '\n'])];
        tikz = [tikz, sprintf(['    ' 'axis background/.style={fill=white},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'axis x line*=bottom,' '\n'])];
        tikz = [tikz, sprintf(['    ' 'axis y line=left,' '\n'])];
        % Add additional axis options
        if iscell(options.axisOptions)
            for i=1:length(options.axisOptions)
                tikz = [tikz, sprintf(['    ' options.axisOptions{i} ',' '\n'])]; %#ok
            end
        elseif ~isempty(options.axisOptions)
            tikz = [tikz, sprintf(['    ' options.axisOptions ',' '\n'])];
        end
        tikz = [tikz, sprintf(['    ' 'width = %0.4gin, height = %0.4gin, scale only axis]' '\n'], imgwidth, (options.dims(4)-options.dims(3))*imgwidth/(options.dims(2)-options.dims(1)))];
    end
    
    %% Generate road
    if options.background
        switch options.case
            case 'turning'
                % Road
                road2.x = [21 34];
                % TODO: Turn magic numbers into parameters
                if options.matlab
                    % Offroad
                    set(gca,'color', colors.offroad);
                    set(gcf,'inverthardcopy','off');
                    % Road
                    patch([params.xmin(1) params.xmin(1) params.xmax(1) params.xmax(1)],[params.xmin(2), params.xmax(2)+ params.lane, params.xmax(2)+ params.lane , params.xmin(2)], colors.road, 'EdgeColor', 'none'); hold on;
                    patch(kron(road2.x,[1 1]),[params.xmin(2)-30 , params.xmin(2), params.xmin(2) , params.xmin(2)-30] , colors.road, 'EdgeColor', 'none');
                    % Borders
                    plot([0, params.xmax(1)], [params.lane, params.lane], 'k-', 'LineWidth', 2);
                    plot([0, road2.x(1), road2.x(1)], [-params.lane, -params.lane, params.xmin(2)-30], 'k-', 'LineWidth', 2);
                    plot([road2.x(2), road2.x(2), params.xmax(1)], [params.xmin(2)-30, -params.lane, -params.lane], 'k-', 'LineWidth', 2);
                    % Lane markation
                    plot([0, params.xmax(1)], [0, 0], 'k--', 'LineWidth', 2);
                    axis equal;
                    axis(options.dims);
                end
                if options.tikz
                    % Borders
                    tikz = [tikz, sprintf(['  ' '\\draw[black,line width = 2pt] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], 0, params.lane, params.xmax(1), params.lane)];
                    tikz = [tikz, sprintf(['  ' '\\draw[black,line width = 2pt] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], 0, -params.lane, road2.x(1), -params.lane, road2.x(1), params.xmin(2)-30)];
                    tikz = [tikz, sprintf(['  ' '\\draw[black,line width = 2pt] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], road2.x(2), params.xmin(2)-30, road2.x(2), -params.lane, params.xmax(1), -params.lane)];
                    % Lane markation
                    tikz = [tikz, sprintf(['  ' '\\addplot[domain=%0.4g:%0.4g,color=black,line width = 3pt, dashed] {0};' '\n'], 0, params.xmax(1))];
                end
            case 'overtake'
                if options.matlab
                    % Offroad
                    set(gca,'color', colors.offroad);
                    set(gcf,'inverthardcopy','off');
                    % Road
                    patch([0 50 50 0],[-2*params.lane -2*params.lane params.lane params.lane], colors.road, 'EdgeColor', 'none'); hold on;
                    % Borders
                    plot([0 50], [params.lane, params.lane], 'k-', 'LineWidth', 2);
                    plot([0 50], [-2*params.lane, -2*params.lane], 'k-', 'LineWidth', 2);
                    % Lane markations
                    plot([0, 50], [0, 0], 'k--', 'LineWidth', 2);
                    plot([0, 50], [-params.lane, -params.lane], 'k--', 'LineWidth', 2);
                    axis equal;
                    axis(options.dims);
                end
                if options.tikz
                    % Borders
                    tikz = [tikz, sprintf(['  ' '\\draw[black,line width = 2pt] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], 0, params.lane, 50, params.lane)];
                    tikz = [tikz, sprintf(['  ' '\\draw[black,line width = 2pt] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], 0, -2*params.lane, 50, -2*params.lane)];
                    % Lane markation
                    tikz = [tikz, sprintf(['  ' '\\addplot[domain=%0.4g:%0.4g,color=black,line width = 3pt, dashed] {0};' '\n'], 0, 50)];
                    tikz = [tikz, sprintf(['  ' '\\addplot[domain=%0.4g:%0.4g,color=black,line width = 3pt, dashed] {%0.4g};' '\n'], 0, 50, -params.lane)];
                end
            otherwise
        end
    end
    
    %% Pre-process trajectories
    if strcmp(options.case, 'overtake')
        for j=length(objects):-1:1
            if ~isempty(strfind(objects{j}.type, 'truck'))
                for i=1:length(trajectories.(objects{j}.type))
                    trajectories.(objects{j}.type){i} = trajectories.(objects{j}.type){i} + repmat([objects{j}.d(1)/2; 0],1, size(trajectories.(objects{j}.type){i},2));
                end
            end
        end
    end
    
    %% Generate trajectories
    for j=length(objects):-1:1
        for i=1:length(trajectories.(objects{j}.type))
            if isempty(options.highlight) || all(i ~= options.highlight) % Only plot non-highlighted trajectories
                if options.matlab
                    plot(trajectories.(objects{j}.type){i}(1,1:N+1), trajectories.(objects{j}.type){i}(2,1:N+1), 'Color', colors.(objects{j}.type),  'LineWidth', 1);
                end
                if options.tikz
                    opacity = 10/length(trajectories.(objects{j}.type));
                    opacity = 0.1;
                    if ~isempty(options.highlight)
                        tikz = [tikz, sprintf(['  ' '\\addplot[color=light' tikzcolors.(objects{j}.type) ',opacity=%0.4g,solid,line width=1pt,forget plot]' '\n'], opacity)];
                    else
                        tikz = [tikz, sprintf(['  ' '\\addplot[color=' tikzcolors.(objects{j}.type) ',opacity=%0.4g,solid,line width=1pt,forget plot]' '\n'], opacity)];
                    end
                    tikz = [tikz, sprintf(['    ' 'table[row sep=crcr]{' '\n'])];
                    for k=0:N
                        tikz = [tikz, sprintf(['      ' '%g' '    ' '%g' '\\\\' '\n'], trajectories.(objects{j}.type){i}(1,k+1), trajectories.(objects{j}.type){i}(2,k+1))];
                    end
                    tikz = [tikz, sprintf(['    ' '};' '\n'])];
                end
            end
        end
    end
    
    %% Generate objects
    if ~isempty(options.highlight)
        for k=0:N
            for j=length(objects):-1:1
                if any(objects{j}.k == k) % If objects should be plotted at step k
                    x0 = trajectories.(objects{j}.type){options.highlight}(1:2,k+1);
                    if ~isempty(strfind(objects{j}.type, 'car'))
                        if strcmp(options.case, 'turning')
                            theta0 = atan2(trajectories.(objects{j}.type){options.highlight}(4,k+1),trajectories.(objects{j}.type){options.highlight}(3,k+1));
                        else
                            theta0 = 0;
                        end
                    elseif ~isempty(strfind(objects{j}.type, 'truck'))
                        if strcmp(options.case, 'turning')
                            theta0 = trajectories.(objects{j}.type){options.highlight}(3,k+1);
                        else
                            theta0 = pi;
                        end
                    end
                    box = [objects{j}.d(1)/2*cos(theta0) + objects{j}.d(2)/2*sin(theta0) ; objects{j}.d(1)/2*sin(theta0) - objects{j}.d(2)/2*cos(theta0)];
                    box = [box [objects{j}.d(1)/2*cos(theta0) - objects{j}.d(2)/2*sin(theta0) ; objects{j}.d(1)/2*sin(theta0) + objects{j}.d(2)/2*cos(theta0)]]; %#ok
                    box = [box -box]; %#ok
                    if options.matlab
                        patch(x0(1)+box(1,:), x0(2)+box(2,:), colors.(objects{j}.type));
                    end
                    if options.tikz
                        if ~isempty(strfind(objects{j}.type, 'car'))
                            tikz = [tikz, sprintf(['  ' '\\node[rotate=%0.4g] at (axis cs:%0.4g,%0.4g) {\\includegraphics[width=%0.4gin,height=%0.4gin]{car.png}};' '\n'], theta0/2/pi*360, x0(1), x0(2), imgwidth/(options.dims(2)-options.dims(1))*objects{j}.d(1), imgwidth/(options.dims(2)-options.dims(1))*objects{j}.d(2))]; %#ok
                        elseif ~isempty(strfind(objects{j}.type, 'truck'))
                            tikz = [tikz, sprintf(['  ' '\\node[rotate=%0.4g] at (axis cs:%0.4g,%0.4g) {\\includegraphics[width=%0.4gin,height=%0.4gin]{truck.png}};' '\n'], theta0/2/pi*360+180, x0(1), x0(2), imgwidth/(options.dims(2)-options.dims(1))*objects{j}.d(1), imgwidth/(options.dims(2)-options.dims(1))*objects{j}.d(2))]; %#ok
                        end
                    end
                    if options.collisionBoxes && options.tikz
                        % Reset angle to zero
                        if ~isempty(strfind(objects{j}.type, 'car'))
                            theta0 = 0;
                            box = [objects{j}.d(1)/2*cos(theta0) + objects{j}.d(2)/2*sin(theta0) ; objects{j}.d(1)/2*sin(theta0) - objects{j}.d(2)/2*cos(theta0)];
                            box = [box [objects{j}.d(1)/2*cos(theta0) - objects{j}.d(2)/2*sin(theta0) ; objects{j}.d(1)/2*sin(theta0) + objects{j}.d(2)/2*cos(theta0)]]; %#ok
                            box = [box -box]; %#ok
                        end
                        % preaction/postaction ensures that the border is drawn inside the box (to not enlarge it)
                        tikz = [tikz, sprintf(['\t' '\\fill[pattern color=verylight' tikzcolors.(objects{j}.type) ',pattern=north west lines,' ...
                                                           'preaction={clip,postaction={draw=light' tikzcolors.(objects{j}.type) ',line width=2pt}}] ' ...
                                                           '(axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g) -- cycle;' '\n'], x0(1)+box(1,1), x0(2)+box(2,1), x0(1)+box(1,2), x0(2)+box(2,2), x0(1)+box(1,3), x0(2)+box(2,3), x0(1)+box(1,4), x0(2)+box(2,4))]; %#ok
                    end
                end
            end
        end
    end
    
    %% Generate highlighted trajectory
    if ~isempty(options.highlight)
        for j=length(objects):-1:1
            if options.matlab
                plot(trajectories.(objects{j}.type){options.highlight}(1,1:N+1), trajectories.(objects{j}.type){options.highlight}(2,1:N+1), 'Color', colors.(objects{j}.type), 'Marker', markers.(objects{j}.type), 'LineStyle', linestyle.(objects{j}.type), 'LineWidth', 0.5);
            end
            if options.tikz
                tikz = [tikz, sprintf(['  ' '\\addplot[color=' tikzcolors.(objects{j}.type) ',solid,mark=' tikzmarkers.(objects{j}.type) ',line width=2pt,forget plot]' '\n'])];
                tikz = [tikz, sprintf(['    ' 'table[row sep=crcr]{' '\n'])];
                for k=0:N
                    tikz = [tikz, sprintf(['      ' '%g' '    ' '%g' '\\\\' '\n'], trajectories.(objects{j}.type){options.highlight}(1,k+1), trajectories.(objects{j}.type){options.highlight}(2,k+1))];
                end
                tikz = [tikz, sprintf(['    ' '};' '\n'])];
            end
        end
    end
    
    %% Finish up
    if options.tikz
        if ~isempty(options.spy)
            tikz = [tikz, sprintf(['  ' '\\spy[' options.spy.color ',size=%0.4gin] on (%0.4g,%0.4g) in node at (%0.4g,%0.4g);' '\n'], options.spy.size, options.spy.pos(1), options.spy.pos(2), options.spy.viewer(1), options.spy.viewer(2))];
        end
        tikz = [tikz, sprintf(['  ' '\\end{axis}' '\n'])];
        tikz = [tikz, sprintf(['\\end{tikzpicture}' '\n'])];
        f = fopen([options.name '.tikz'], 'w');
        fprintf(f, '%s', tikz);
        fclose(f);
    end
    
end