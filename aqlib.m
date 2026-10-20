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
    L.monthlyAvailability = @monthlyAvailability;

    % --- Cleaning and preprocessing ---
    L.gapSummary        = @gapSummary;
    L.gapHours          = @gapHours;
    L.dropSparseVars    = @dropSparseVars;
    L.interpolateGaps   = @interpolateGaps;

    % --- EDA plots ---
    L.plotTargetDistribution = @plotTargetDistribution;
    L.plotHistograms         = @plotHistograms;
    L.plotTargetCorrelation  = @plotTargetCorrelation;
    L.plotMissingMap         = @plotMissingMap;
    L.plotInterpolation      = @plotInterpolation;
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

function A = monthlyAvailability(T)
% Percentage of hours with a valid value, per month and per variable.
% Shows when each instrument was working (e.g. NMHC(GT) only early on).
% Input:  T - data table with a DateTime column (missing values as NaN)
% Output: A - table with Month and one % column per numeric variable
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    [grp, month] = findgroups(dateshift(T.DateTime, 'start', 'month'));
    month.Format = 'MMM yyyy';
    A = table(month, 'VariableNames', {'Month'});
    for k = 1:numel(vars)
        A.(vars{k}) = round(splitapply(@(x) 100 * mean(~isnan(x)), T.(vars{k}), grp), 1);
    end
end

%% ===================== Cleaning and preprocessing =====================

function G = gapSummary(T, maxGap)
% Describes the gaps (runs of consecutive missing hours) of each numeric variable.
% Input:  T      - data table (missing values as NaN, one row per hour)
%         maxGap - longest gap, in hours, that is considered short
% Output: G      - table with number of gaps, median and max gap length,
%                  number of short gaps and % of missing cells in short gaps
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    n = numel(vars);
    g = zeros(n, 5);
    for k = 1:n
        r = nanRuns(T.(vars{k}));
        if isempty(r), continue; end
        short = r <= maxGap;
        g(k, :) = [numel(r), median(r), max(r), sum(short), ...
                   100 * sum(r(short)) / sum(r)];
    end
    G = array2table(g, 'RowNames', vars, 'VariableNames', ...
        {'Gaps', 'MedianLen', 'MaxLen', 'ShortGaps', 'ShortCellsPct'});
    G.ShortCellsPct = round(G.ShortCellsPct, 1);
    G = sortrows(G, 'Gaps', 'descend');
end

function H = gapHours(T, vars, maxLen)
% Counts at which hour of the day the short gaps start, to check whether
% they are random or systematic (e.g. a daily calibration of the analyser).
% Input:  T      - data table with a DateTime column (missing values as NaN)
%         vars   - cell array of variable names
%         maxLen - only gaps with at most this many hours are counted
% Output: H      - table with Hour (0-23) and one count column per variable
    H = table((0:23)', 'VariableNames', {'Hour'});
    for k = 1:numel(vars)
        [r, s] = nanRuns(T.(vars{k}));
        h = hour(T.DateTime(s(r <= maxLen)));
        H.(vars{k}) = histcounts(h, -0.5:1:23.5)';
    end
end

function [T, dropped] = dropSparseVars(T, maxPct, keep)
% Removes the numeric variables with more than maxPct % of missing values.
% Input:  T       - data table (missing values as NaN)
%         maxPct  - maximum % of missing values allowed
%         keep    - cell array of variables never removed (e.g. the target)
% Output: T       - table without the sparse variables
%         dropped - cell array with the names of the removed variables
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    pct = 100 * mean(isnan(T{:, vars}));
    dropped = vars(pct > maxPct & ~ismember(vars, keep));
    T = removevars(T, dropped);
end

function [T, filled] = interpolateGaps(T, vars, maxGap)
% Fills gaps of at most maxGap consecutive missing hours by linear
% interpolation in time. Longer gaps and gaps at the ends of the series
% are left as NaN.
% Input:  T      - data table with a DateTime column (missing values as NaN)
%         vars   - cell array of variables to fill (never the target)
%         maxGap - longest gap, in missing hours, that is filled
% Output: T      - table with the short gaps filled
%         filled - logical matrix (rows x vars), true where a value was interpolated
    filled = false(height(T), numel(vars));
    for k = 1:numel(vars)
        x = T.(vars{k});

        % mark the cells that belong to short gaps
        [r, s] = nanRuns(x);
        short = false(size(x));
        for g = find(r <= maxGap)'
            short(s(g):s(g) + r(g) - 1) = true;
        end

        xi = fillmissing(x, 'linear', 'SamplePoints', T.DateTime, 'EndValues', 'none');
        x(short) = xi(short);
        filled(:, k) = short & ~isnan(x);
        T.(vars{k}) = x;
    end
end

function [r, s] = nanRuns(x)
% Lengths and start positions of the runs of consecutive NaN in a vector.
% Input:  x - numeric vector
% Output: r - column vector with the length of each run (empty if none)
%         s - index of the first NaN of each run
    d = diff([0; isnan(x(:)); 0]);
    s = find(d == 1);
    r = find(d == -1) - s;
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

function fig = plotMissingMap(T)
% Map of missing values over time: one row per variable, dark = missing.
% Input:  T   - data table with a DateTime column (missing values as NaN)
% Output: fig - figure handle
    vars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    M = isnan(T{:, vars})';
    fig = newFigure([100 100 1100 450]);
    imagesc(M);
    colormap([1 1 1; 0.75 0.15 0.15]);
    clim([0 1]);

    % one tick at the start of each month
    ticks = find(T.DateTime == dateshift(T.DateTime, 'start', 'month'));
    set(gca, 'XTick', ticks, 'XTickLabel', cellstr(T.DateTime(ticks), 'MMM yy'), ...
        'YTick', 1:numel(vars), 'YTickLabel', vars, 'TickLabelInterpreter', 'none');
    xlabel('Time');
    title('Missing values over time (red = missing)');
end

function fig = plotInterpolation(t, Xb, Xa, F, names)
% Before/after view of the interpolation: measured values as a line,
% interpolated values as red markers, long gaps left as breaks.
% Input:  t     - datetime vector of the window to show
%         Xb    - values before interpolation (rows x variables)
%         Xa    - values after interpolation (rows x variables)
%         F     - logical matrix, true where a value was interpolated
%         names - cell array with the variable names
% Output: fig   - figure handle
    n = numel(names);
    fig = newFigure([100 100 1100 280 * n]);
    for k = 1:n
        subplot(n, 1, k);
        plot(t, Xa(:, k), ':', 'Color', [0.85 0.15 0.15], 'LineWidth', 1);
        hold on;
        plot(t, Xb(:, k), '-', 'Color', [0.18 0.44 0.62], 'LineWidth', 1.2);
        plot(t(F(:, k)), Xa(F(:, k), k), 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', 'none');
        hold off;
        ylabel(names{k}, 'Interpreter', 'none');
        grid on;
        if k == 1
            legend({'After interpolation', 'Measured', 'Interpolated value'}, ...
                'Location', 'northoutside', 'Orientation', 'horizontal');
            title('Short-gap interpolation (long gaps are left missing)');
        end
    end
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
