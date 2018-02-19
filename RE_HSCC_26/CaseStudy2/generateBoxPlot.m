function [h, ax] = generateBoxPlot(varargin)
    p = inputParser;
    
    p.addRequired('data');
    p.addRequired('labels', @iscell);
    p.addParameter('bounds', []);
    p.addParameter('boundLabels', @iscell);
    p.addParameter('percentiles', [25 75]);
    p.addParameter('imgwidth', 4.5); % In inches
    p.addParameter('imgheight', 4.5); % In inches
    p.addParameter('log', false, @islogical);
    p.addParameter('connected', []);
    p.addParameter('median', true, @islogical);
    p.addParameter('matlab', true, @islogical); % Enable generation of MATLAB figure
    p.addParameter('tikz', false, @islogical); % Enable generation of tikz figure
    p.addParameter('dims', []);
    p.addParameter('xlabel', '', @ischar);
    p.addParameter('ylabel', '', @ischar);
    p.addParameter('axisOptions', '', @(c)(ischar(c) || (iscell(c) && all(cellfun(@ischar, c)))));
    p.addParameter('name', 'test', @ischar);
    p.parse(varargin{:});
    options = p.Results;
    data = options.data;
    labels = options.labels;
    
    tikzcolors.inner = 'blue';
    tikzcolors.median = 'red';
    tikzcolors.mean = 'lightred';
    tikzcolors.range = 'black';
    tikzcolors.bounds = 'gray';
    
    %% Sanity/integrity checks
    if length(data) ~= length(labels)
        error('Data does not match.');
    end
    for i=1:length(labels) % Ensure they are row vectors
        data{i} = data{i}(:)';
    end
    
    %% Extract statistics
    statistics = {};
    if isempty(options.dims)
        if ~isempty(options.bounds)
            options.dims = [inf, max(options.bounds)];
        else
            options.dims = [inf, -inf];
        end
    end
    for i=1:length(labels)
        statistics{i}.median = prctile(data{i},50);
        statistics{i}.mean = mean(data{i});
        statistics{i}.innerPercentile = prctile(data{i},options.percentiles);
        statistics{i}.range = [min(data{i}), max(data{i})];
        options.dims(1) = min(options.dims(1),min(data{i}));
        options.dims(2) = max(options.dims(2),max(data{i}));
    end
    
    %% Start building figure
    if options.matlab
        h = figure();
    end
    if options.tikz
        tikz = sprintf(['\\begin{tikzpicture}[scale=1,cap=round]' '\n']);
        tikz = [tikz, sprintf(['  ' '\\begin{axis}[' '\n'])];
        tikz = [tikz, sprintf(['    ' 'every outer x axis line/.append style={white!20!black},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'every x tick label/.append style={font=\\color{white!20!black}},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'xmin=%0.4g,' '\n'], 0)];
        tikz = [tikz, sprintf(['    ' 'xmax=%0.4g,' '\n'], length(labels)+1)];
        if options.log
            tikz = [tikz, sprintf(['    ' 'ymode=log,' '\n'])];
        end
        tikz = [tikz, sprintf(['    ' 'ymin=%0.4g,' '\n'], options.dims(1))];
        tikz = [tikz, sprintf(['    ' 'ymax=%0.4g,' '\n'], options.dims(2))];
        tikz = [tikz, sprintf(['    ' 'xtick={1,...,' num2str(length(labels)) '},' '\n'])];
        tikz = [tikz, sprintf(['    ' 'xticklabels={' strjoin(strrep(labels, '\', '\\'), ', ') '},' '\n'])];
        if ~isempty(options.xlabel)
            tikz = [tikz, sprintf(['    ' 'xlabel=' options.xlabel ',' '\n'])];
        end
        if ~isempty(options.ylabel)
            tikz = [tikz, sprintf(['    ' 'ylabel=' options.ylabel ',' '\n'])];
        end
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
        tikz = [tikz, sprintf(['    ' 'width = %0.4gin, height = %0.4gin, scale only axis]' '\n'], options.imgwidth, options.imgheight)];
    end
    
    %% Draw
    if options.matlab
        hb = boxplot(cell2mat(data), cell2mat(cellfun(@(c)(c(2)*ones(1,c(1))), mat2cell([cellfun(@length, data); 1:length(data)]', ones(1,length(labels)), 2), 'UniformOutput', false)'), 'labels', labels); hold on;
        if options.median
            set(hb(6,:), {'Ydata'}, cellfun(@(c)([c.median c.median]), statistics, 'UniformOutput', false)');
        else
            set(hb(6,:), {'Ydata'}, cellfun(@(c)([c.mean c.mean]), statistics, 'UniformOutput', false)');
        end
        % Adjust percentile
        set(hb(1,:), {'Ydata'}, cellfun(@(c)([c.innerPercentile(2) c.range(2)]), statistics, 'UniformOutput', false)'); % Upper wisker line
        set(hb(2,:), {'Ydata'}, cellfun(@(c)([c.range(1) c.innerPercentile(1)]), statistics, 'UniformOutput', false)'); % Lower wisker line
        set(hb(3,:), {'Ydata'}, cellfun(@(c)([c.range(2) c.range(2)]), statistics, 'UniformOutput', false)'); % Upper wisker end bar
        set(hb(4,:), {'Ydata'}, cellfun(@(c)([c.range(1) c.range(1)]), statistics, 'UniformOutput', false)'); % Lower wisker end bar
        set(hb(5,:), {'Ydata'}, cellfun(@(c)([c.innerPercentile(1) c.innerPercentile(2) c.innerPercentile(2) c.innerPercentile(1) c.innerPercentile(1)]), statistics, 'UniformOutput', false)');
        ax = findall(h, 'type', 'axes');
        hold(ax, 'on');
        % Plot mean
        if ~isempty(options.connected)
            plot(ax, 1:length(labels), cellfun(@(c)(c.(options.connected)), statistics), '-*', 'Color', [1 0.35 0.35]);
        end
        % Add bounds
        if ~isempty(options.bounds)
            for i=1:length(labels)
                plot(ax, [i-0.45 i+0.45],[options.bounds(i) options.bounds(i)], ':', 'Color', [0.5 0.5 0.5]);
            end
        end
        if options.log
            set(ax,'yscale','log');
        end
        ylim(ax, options.dims);
        if ~isempty(options.xlabel)
            xlabel(ax, options.xlabel);
        end
        if ~isempty(options.ylabel)
            ylabel(ax, options.ylabel);
        end
    end
    if options.tikz
        boxwidth = 0.4;
        for i=1:length(labels)
            if ~isempty(options.bounds)
                tikz = [tikz, sprintf(['    ' '\\draw[' tikzcolors.bounds ',thick,dotted] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g)'], i-0.45, options.bounds(i), i+0.45, options.bounds(i))];
                if ~isempty(options.boundLabels)
                    tikz = [tikz, sprintf([' node[midway,above] {' strrep(options.boundLabels{i}, '\', '\\') '}'])];
                end
                tikz = [tikz, sprintf([';' '\n'])];
            end
            tikz = [tikz, sprintf(['    ' '\\draw[' tikzcolors.range ',thin,-|] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], i, statistics{i}.innerPercentile(2), i, statistics{i}.range(2))];
            tikz = [tikz, sprintf(['    ' '\\draw[' tikzcolors.range ',thin,-|] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], i, statistics{i}.innerPercentile(1), i, statistics{i}.range(1))];
            tikz = [tikz, sprintf(['    ' '\\filldraw[' tikzcolors.inner ',fill=white,thick] (axis cs:%0.4g,%0.4g) rectangle (axis cs:%0.4g,%0.4g);' '\n'], i-boxwidth/2, statistics{i}.innerPercentile(1), i+boxwidth/2, statistics{i}.innerPercentile(2))];
            if options.median
                tikz = [tikz, sprintf(['    ' '\\draw[' tikzcolors.median ',very thick] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], i-boxwidth/2, statistics{i}.median, i+boxwidth/2, statistics{i}.median)];
            else
                tikz = [tikz, sprintf(['    ' '\\draw[' tikzcolors.mean ',very thick] (axis cs:%0.4g,%0.4g) -- (axis cs:%0.4g,%0.4g);' '\n'], i-boxwidth/2, statistics{i}.mean, i+boxwidth/2, statistics{i}.mean)];
            end
        end
        % Plot mean
        if ~isempty(options.connected)
            tikz = [tikz, sprintf(['  ' '\\addplot[color=' tikzcolors.(options.connected) ',solid,thick,mark=*,mark size=1pt,forget plot]' '\n'])];
            tikz = [tikz, sprintf(['    ' 'table[row sep=crcr]{' '\n'])];
            for i=1:length(labels)
                tikz = [tikz, sprintf(['      ' '%g' '    ' '%g' '\\\\' '\n'], i, statistics{i}.(options.connected))];
            end
            tikz = [tikz, sprintf(['    ' '};' '\n'])];
        end
    end
    
    %% Finish up
    if options.tikz
        tikz = [tikz, sprintf(['  ' '\\end{axis}' '\n'])];
        tikz = [tikz, sprintf(['\\end{tikzpicture}' '\n'])];
        f = fopen([options.name '.tikz'], 'w');
        fprintf(f, '%s', tikz);
        fclose(f);
    end
    
end