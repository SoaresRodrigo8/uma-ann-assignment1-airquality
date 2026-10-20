function L = aqlib()
% AQLIB  Function library for the UCI Air Quality feature analysis.
%   L = aqlib() returns a struct of function handles.

    % --- Import and exploration ---
    L.loadData          = @loadData;
    L.markMissing       = @markMissing;
    L.checkTimeIndex    = @checkTimeIndex;
    L.missingSummary    = @missingSummary;
    L.describeVars      = @describeVars;
    L.targetCorrelation = @targetCorrelation;

    % --- EDA plots ---
    L.plotTargetDistribution = @plotTargetDistribution;
    L.plotHistograms         = @plotHistograms;
    L.plotTargetCorrelation  = @plotTargetCorrelation;
end

%% ===================== Import and exploration =====================

function T = loadData(filePath)
% Reads the raw Air Quality CSV (';' separator, ',' decimal).
% Date and Time are merged into a single DateTime column.
% Input:  filePath - path to AirQualityUCI.csv
% Output: T        - table with DateTime followed by the numeric variables
    opts = detectImportOptions(filePath, 'Delimiter', ';', ...
        'DecimalSeparator', ',', 'VariableNamingRule', 'preserve');
    opts = setvartype(opts, {'Date', 'Time'}, 'char');
    T = readtable(filePath, opts);

    % drop trailing empty rows (no date)
    T = T(~cellfun(@isempty, T.Date), :);

    dt = datetime(strcat(T.Date, {' '}, T.Time), ...
        'InputFormat', 'dd/MM/yyyy HH.mm.ss');
    T = [table(dt, 'VariableNames', {'DateTime'}), removevars(T, {'Date', 'Time'})];
end

function T = markMissing(T, code)
% Replaces the missing-value code (e.g. -200) by NaN in every numeric column.
% Input:  T    - data table
%         code - numeric code used for missing values
% Output: T    - same table with NaN in place of the code
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    for k = 1:numel(vars)
        x = T.(vars{k});
        x(x == code) = NaN;
        T.(vars{k}) = x;
    end
end

function S = checkTimeIndex(dt)
% Checks that the timestamps form a regular hourly series.
% Input:  dt - datetime vector
% Output: S  - struct with start, stop, duplicates, sorted flag,
%              expected and missing hours
    S.start      = min(dt);
    S.stop       = max(dt);
    S.duplicates = numel(dt) - numel(unique(dt));
    S.sorted     = issorted(dt);
    fullRange    = (S.start:hours(1):S.stop)';
    S.expected   = numel(fullRange);
    S.missing    = numel(setdiff(fullRange, dt));
end

function M = missingSummary(T)
% Counts missing values (NaN) per numeric variable.
% Input:  T - data table (missing values already set to NaN)
% Output: M - table with Missing count and Percent, sorted descending
    vars  = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    X     = T{:, vars};
    count = sum(isnan(X))';
    pct   = 100 * count / size(X, 1);
    M = table(count, round(pct, 2), 'VariableNames', {'Missing', 'Percent'}, ...
        'RowNames', vars);
    M = sortrows(M, 'Missing', 'descend');
end

function D = describeVars(T)
% Descriptive statistics of every numeric variable, ignoring NaN.
% Kurtosis is reported as excess kurtosis (0 for a normal distribution).
% Input:  T - data table
% Output: D - table, one row per variable
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    n = numel(vars);
    st = zeros(n, 10);
    for k = 1:n
        x = T.(vars{k});
        x = x(~isnan(x));
        q = prctile(x, [25 50 75]);
        st(k, :) = [numel(x), mean(x), std(x), min(x), q, max(x), ...
                    skewness(x), kurtosis(x) - 3];
    end
    D = array2table(round(st, 3), 'RowNames', vars, 'VariableNames', ...
        {'N', 'Mean', 'Std', 'Min', 'Q1', 'Median', 'Q3', 'Max', 'Skewness', 'Kurtosis'});
end

function C = targetCorrelation(T, target)
% Pearson and Spearman correlation of every numeric feature with the target,
% using pairwise complete rows. Sorted by absolute Spearman correlation.
% Input:  T      - data table
%         target - name of the target variable
% Output: C      - table with Pearson, Spearman and N (pairs used)
    vars  = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    feats = setdiff(vars, target, 'stable');
    y = T.(target);
    n = numel(feats);
    r = zeros(n, 3);
    for k = 1:n
        x  = T.(feats{k});
        ok = ~isnan(x) & ~isnan(y);
        r(k, :) = [corr(x(ok), y(ok), 'Type', 'Pearson'), ...
                   corr(x(ok), y(ok), 'Type', 'Spearman'), sum(ok)];
    end
    C = array2table(r, 'RowNames', feats, 'VariableNames', {'Pearson', 'Spearman', 'N'});
    C.Pearson  = round(C.Pearson, 3);
    C.Spearman = round(C.Spearman, 3);
    [~, idx] = sort(abs(C.Spearman), 'descend');
    C = C(idx, :);
end

%% ===================== EDA plots =====================

function fig = plotTargetDistribution(y, name)
% Histogram (with mean and median) and boxplot of the target.
% Input:  y    - target vector (NaN ignored)
%         name - label of the target
% Output: fig  - figure handle
    y = y(~isnan(y));
    fig = newFigure([100 100 800 500]);

    % bin edges at x.x5 so that the 0.1-resolution values never sit on an edge
    edges = 0.05:0.2:max(y) + 0.2;
    subplot(4, 1, 1:3);
    histogram(y, edges, 'FaceColor', [0.18 0.44 0.62], 'EdgeColor', 'w');
    hold on;
    xline(mean(y),   '--r', sprintf('Mean = %.2f', mean(y)), 'LineWidth', 1.5);
    xline(median(y), ':k',  sprintf('Median = %.2f', median(y)), 'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left');
    hold off;
    ylabel('Count');
    title(['Distribution of ' name]);
    xl = xlim;

    subplot(4, 1, 4);
    boxplot(y, 'Orientation', 'horizontal', 'Symbol', '.');
    xlim(xl);
    set(gca, 'YTick', []);
    xlabel([name ' (mg/m^3)']);
end

function fig = plotHistograms(T, vars)
% Grid of histograms, one per variable.
% Input:  T    - data table
%         vars - cell array of variable names
% Output: fig  - figure handle
    n = numel(vars);
    nc = 4;
    nr = ceil(n / nc);
    fig = newFigure([50 50 1400 260 * nr]);
    for k = 1:n
        subplot(nr, nc, k);
        x = T.(vars{k});
        histogram(x(~isnan(x)), 40, 'FaceColor', [0.18 0.44 0.62], 'EdgeColor', 'w');
        title(vars{k}, 'Interpreter', 'none');
    end
end

function fig = plotTargetCorrelation(C, target)
% Horizontal bar chart of Pearson and Spearman correlation with the target.
% Input:  C      - table from targetCorrelation
%         target - name of the target (for the title)
% Output: fig    - figure handle
    fig = newFigure([100 100 800 550]);
    barh([C.Pearson, C.Spearman]);
    set(gca, 'YDir', 'reverse', 'YTick', 1:height(C), ...
        'YTickLabel', C.Properties.RowNames, 'TickLabelInterpreter', 'none');
    xline(0, 'k');
    xlim([-1 1]);
    legend({'Pearson', 'Spearman'}, 'Location', 'southeast');
    xlabel(['Correlation with ' target]);
    title(['Feature correlation with ' target], 'Interpreter', 'none');
    grid on;
end

function fig = newFigure(pos)
% Creates a white figure (light theme, also when MATLAB uses a dark theme).
% Input:  pos - figure position [left bottom width height]
% Output: fig - figure handle
    fig = figure('Color', 'w', 'Position', pos);
    if exist('theme', 'file')
        theme(fig, 'light');
    end
end
