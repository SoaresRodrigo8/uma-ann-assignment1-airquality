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
    L.latexTable        = @latexTable;

    % --- Cleaning and preprocessing ---
    L.iqrBounds            = @iqrBounds;
    L.outlierSummary       = @outlierSummary;
    L.clipToBounds         = @clipToBounds;
    L.hourlyIndex          = @hourlyIndex;
    L.gapStats             = @gapStats;
    L.interpolateShortGaps = @interpolateShortGaps;
    L.knnImpute            = @knnImpute;
    L.zscoreFit            = @zscoreFit;
    L.zscoreApply          = @zscoreApply;
    L.compareStats         = @compareStats;

    % --- EDA plots ---
    L.plotTargetDistribution = @plotTargetDistribution;
    L.plotHistograms         = @plotHistograms;
    L.plotSkewGroups         = @plotSkewGroups;
    L.plotTargetCorrelation  = @plotTargetCorrelation;
    L.plotScalingBeforeAfter = @plotScalingBeforeAfter;
    L.plotFeaturesVsTarget   = @plotFeaturesVsTarget;
    L.plotScatterByTime      = @plotScatterByTime;
    L.plotTargetOverTime     = @plotTargetOverTime;
    L.topCorrelatedPairs     = @topCorrelatedPairs;
    L.plotPairs              = @plotPairs;

    % --- Feature selection ---
    L.featureRelevance = @featureRelevance;
end

%% ===================== Import and exploration =====================

function df = loadData(filePath, numVars)
% Reads the raw Air Quality CSV (';' separator, ',' decimal mark), keeping the
% original column names. The empty trailing columns and the blank rows at the
% end of the file are dropped, and Date + Time are parsed into Datetime.
% Input:  filePath - path to AirQualityUCI.csv
%         numVars  - string array with the numeric columns to read
% Output: df       - table with Date, Time, the numeric columns and Datetime
    opts = detectImportOptions(filePath, "Delimiter", ";", "VariableNamingRule", "preserve");
    opts.SelectedVariableNames = ["Date", "Time", numVars];   % drops the empty trailing columns
    opts = setvartype(opts, ["Date", "Time"], "string");
    opts = setvartype(opts, numVars, "double");
    opts = setvaropts(opts, numVars, "DecimalSeparator", ",");
    df = readtable(filePath, opts);

    df = df(strlength(df.Date) > 0, :);    % remove blank rows at the end of the file

    df.Datetime = datetime(df.Date + " " + df.Time, "InputFormat", "dd/MM/yyyy HH.mm.ss");
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

function S = checkTimeIndex(ts)
% Checks that the timestamps form a regular hourly series.
% Input:  ts - datetime vector
% Output: S  - struct with start, stop, duplicates, sorted flag, expected,
%              present and missing hours, and a table of the steps between rows
    h = hours(ts - min(ts));
    S.start      = min(ts);
    S.stop       = max(ts);
    S.duplicates = numel(ts) - numel(unique(ts));
    S.sorted     = issorted(ts);
    S.expected   = round(max(h)) + 1;
    S.present    = numel(unique(round(h)));
    S.missing    = S.expected - S.present;
    [step_h, ~, ic] = unique(hours(diff(ts)));
    S.steps = table(step_h, accumarray(ic, 1), 'VariableNames', {'Step_hours', 'Count'});
end

function M = missingSummary(T, vars)
% Counts missing values (NaN) per variable, sorted from most to least missing.
% Input:  T    - data table (missing values already set to NaN)
%         vars - string array of variable names
% Output: M    - table with Missing count and PctMissing (2 decimals)
    miss = sum(isnan(T{:, vars}))';
    M = table(miss, round(miss / height(T) * 100, 2), ...
        'RowNames', cellstr(vars), 'VariableNames', {'Missing', 'PctMissing'});
    M = sortrows(M, 'Missing', 'descend');
end

function D = describeVars(T, vars)
% Descriptive statistics of each variable, ignoring NaN. The skewness is the
% bias-corrected one (flag 0), the same as pandas.
% Input:  T    - data table
%         vars - string array of variable names
% Output: D    - table with mean, std, min, median, max and skew per variable
    S = T{:, vars};
    stats = [mean(S, "omitnan")', std(S, "omitnan")', min(S, [], "omitnan")', ...
             median(S, "omitnan")', max(S, [], "omitnan")', skewness(S, 0)'];
    D = array2table(stats, 'RowNames', cellstr(vars), ...
        'VariableNames', {'mean', 'std', 'min', 'median', 'max', 'skew'});
end

function rho = targetCorrelation(T, feats, target)
% Spearman correlation of each feature with the target, using for each
% feature all rows where both values are present (pairwise).
% Input:  T      - data table
%         feats  - string array of feature names
%         target - name of the target
% Output: rho    - column vector, one correlation per feature
    rho = corr(T{:, feats}, T.(target), "Type", "Spearman", "Rows", "pairwise");
end

function s = latexTable(rowNames, data, header, caption, label)
% Builds a booktabs LaTeX table (2 decimals), ready to paste into the report.
% Input:  rowNames - string array with one name per row
%         data     - numeric matrix (rows x columns)
%         header   - string (array) with the header line(s)
%         caption  - table caption
%         label    - LaTeX label
% Output: s        - string with the full table environment
    s = sprintf("\n\\begin{table}[htbp]\n\\centering\n");
    s = s + sprintf("\\caption{%s}\n\\label{%s}\n", caption, label);
    s = s + sprintf("\\begin{tabular}{l%s}\n\\toprule\n", repmat('r', 1, size(data, 2)));
    s = s + sprintf("%s\n", header);
    s = s + sprintf("\\midrule\n");
    for i = 1:size(data, 1)
        s = s + sprintf("%s & %s \\\\\n", rowNames(i), strjoin(compose("%.2f", data(i, :)), " & "));
    end
    s = s + sprintf("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
end

%% ===================== Cleaning and preprocessing =====================

function [lo, hi] = iqrBounds(X, k)
% Tukey fences of each column: Q1 - k*IQR and Q3 + k*IQR (NaN ignored).
% Input:  X  - numeric matrix (rows x variables)
%         k  - fence multiplier (1.5 is the usual choice)
% Output: lo - row vector with the lower bound of each column
%         hi - row vector with the upper bound of each column
    q1 = quantile7(X, 0.25);
    q3 = quantile7(X, 0.75);
    iqrv = q3 - q1;
    lo = q1 - k * iqrv;
    hi = q3 + k * iqrv;
end

function S = outlierSummary(T, vars, k)
% Counts the outliers of each variable with the IQR rule (Tukey fences).
% Input:  T    - data table (missing values as NaN)
%         vars - string array of variable names
%         k    - fence multiplier
% Output: S    - table with number of outliers, % of the non-missing values
%                (2 decimals), and how many are above / below the fences
    X = T{:, vars};
    [lo, hi] = iqrBounds(X, k);
    n_upper = sum(X > hi)';
    n_lower = sum(X < lo)';
    n_valid = sum(~isnan(X))';
    S = table(n_upper + n_lower, round((n_upper + n_lower) ./ n_valid * 100, 2), n_upper, n_lower, ...
        'RowNames', cellstr(vars), 'VariableNames', {'Outliers', 'Pct', 'Upper', 'Lower'});
end

function [X, nClip] = clipToBounds(X, lo, hi)
% Clips each column to [lo, hi]; missing values stay missing.
% Input:  X     - numeric matrix (rows x variables)
%         lo/hi - row vectors with the bounds of each column
% Output: X     - clipped matrix
%         nClip - column vector with the number of clipped values per column
    nClip = sum(X < lo | X > hi)';
    % "includenan" keeps missing values missing (by default min/max would replace them)
    X = min(max(X, lo, "includenan"), hi, "includenan");
end

function [hidx, inData] = hourlyIndex(dt)
% Position of each row on a complete hourly grid from the first to the last
% timestamp, so that gaps can be measured in hours even after rows were removed.
% Input:  dt     - sorted datetime vector
% Output: hidx   - grid position of each row (1 = first hour)
%         inData - logical vector over the grid, true for hours present in dt
    hidx = round(hours(dt - dt(1))) + 1;
    inData = false(hidx(end), 1);
    inData(hidx) = true;
end

function G = gapStats(x, hidx, inData, maxGap)
% Describes the gaps (runs of missing hours) of one variable on the hourly
% grid. Removed rows count as part of a gap. Only gaps that contain at least
% one row still in the data ("relevant" gaps) need filling.
% Input:  x      - values of the variable for the rows in the data
%         hidx   - grid position of each row (from hourlyIndex)
%         inData - logical grid, true for hours present in the data
%         maxGap - longest gap, in hours, that is considered short
% Output: G      - struct with nGaps (relevant), longest and median length,
%                  and the number of values to fill in short and long gaps
    v = NaN(numel(inData), 1);
    v(hidx) = x;
    [starts, ends] = nanRuns(isnan(v));
    run_len = ends - starts + 1;
    needs_fill = arrayfun(@(s, e) sum(inData(s:e)), starts, ends);
    rel = needs_fill > 0;
    short = run_len <= maxGap;
    G.nGaps   = sum(rel);
    G.longest = max(run_len(rel));
    G.median  = median(run_len(rel));
    G.inShort = sum(needs_fill(rel & short));
    G.inLong  = sum(needs_fill(rel & ~short));
end

function x = interpolateShortGaps(x, hidx, maxGap)
% Linear interpolation over time, only for gaps of at most maxGap hours on the
% hourly grid (removed rows count as part of a gap). Longer gaps and gaps at
% the ends of the series are left missing.
% Input:  x      - values of the variable for the rows in the data
%         hidx   - grid position of each row (from hourlyIndex)
%         maxGap - longest gap, in hours, that is filled
% Output: x      - values with the short gaps filled
    v = NaN(hidx(end), 1);
    v(hidx) = x;
    na = isnan(v);
    [starts, ends] = nanRuns(na);
    run_len = zeros(hidx(end), 1);
    for i = 1:numel(starts)
        run_len(starts(i):ends(i)) = ends(i) - starts(i) + 1;
    end
    v_interp = fillmissing(v, "linear", "EndValues", "none");   % only gaps with values on both sides
    fill_mask = na & run_len <= maxGap;
    v(fill_mask) = v_interp(fill_mask);
    x = v(hidx);
end

function X = knnImpute(X, k)
% Fills the missing values with the mean of the k nearest rows (KNN), on
% standardised columns so that no variable dominates the distance. Only the
% missing values are replaced. The target must not be among the columns.
% Input:  X - numeric matrix (rows x variables) with NaN
%         k - number of neighbours
% Output: X - matrix without missing values
    mu = mean(X, "omitnan");
    sd = std(X, 1, "omitnan");
    Xs = knnFill((X - mu) ./ sd, k);
    X_knn = Xs .* sd + mu;
    X(isnan(X)) = X_knn(isnan(X));          % only missing values are replaced
end

function P = zscoreFit(T, vars)
% Mean and population standard deviation of each variable (as sklearn's
% StandardScaler). Kept separate from zscoreApply so that, in a future model,
% the parameters can be computed on the training set only.
% Input:  T    - data table (no missing values)
%         vars - string array of variable names
% Output: P    - table with Mean and Std, one row per variable
    X = T{:, vars};
    P = table(mean(X)', std(X, 1)', 'RowNames', cellstr(vars), ...
        'VariableNames', {'Mean', 'Std'});
end

function T = zscoreApply(T, P)
% Standardises the variables listed in P: (x - Mean) / Std.
% Input:  T - data table
%         P - table from zscoreFit
% Output: T - table with the standardised variables
    vars = P.Properties.RowNames;
    T{:, vars} = (T{:, vars} - P.Mean') ./ P.Std';
end

function C = compareStats(Xb, Xa, names)
% Standard deviation, maximum and skewness of each variable before and after
% cleaning (NaN ignored in the data before).
% Input:  Xb    - values before (rows x variables)
%         Xa    - values after (rows x variables)
%         names - string array with the variable names
% Output: C     - table with std, max and skew, before and after
    cmp = [std(Xb, "omitnan")', std(Xa)', ...
           max(Xb, [], "omitnan")', max(Xa)', ...
           skewness(Xb, 0)', skewness(Xa, 0)'];
    C = array2table(cmp, 'RowNames', cellstr(names), 'VariableNames', ...
        {'std_before', 'std_after', 'max_before', 'max_after', 'skew_before', 'skew_after'});
end

function q = quantile7(X, p)
% Column-wise quantile with linear interpolation, ignoring NaN.
% Same method as pandas/numpy, so the IQR bounds match the Python version
% (MATLAB's built-in quantile uses a slightly different method).
    q = NaN(1, size(X, 2));
    for j = 1:size(X, 2)
        x = sort(X(~isnan(X(:, j)), j));
        n = numel(x);
        if n == 0, continue; end
        h = (n - 1) * p + 1;
        lo = floor(h);
        hi = min(lo + 1, n);
        q(j) = x(lo) + (h - lo) * (x(hi) - x(lo));
    end
end

function [starts, ends] = nanRuns(na)
% Start and end positions of each run of consecutive true values.
    d = diff([false; na(:); false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
end

function Xs = knnFill(Xs, k)
% K-nearest-neighbours imputation with NaN-aware Euclidean distance, the same
% approach as sklearn's KNNImputer: for each missing value, the k most
% similar rows that have that column observed are found (using only the
% columns observed in both rows), and the mean of their values is used.
    X = Xs;                                  % distances always use the original values
    n_feat = size(X, 2);
    [rows, cols_miss] = find(isnan(X));
    for m = 1:numel(rows)
        r = rows(m);
        c = cols_miss(m);
        donors = find(~isnan(X(:, c)));
        D = X(donors, :) - X(r, :);
        present = ~isnan(D);
        D(~present) = 0;
        n_present = sum(present, 2);
        d = sqrt(n_feat ./ max(n_present, 1) .* sum(D .^ 2, 2));
        d(n_present == 0) = Inf;
        [~, nearest] = mink(d, k);
        Xs(r, c) = mean(X(donors(nearest), c));
    end
end

%% ===================== EDA plots =====================

function fig = plotTargetDistribution(y, name, xLabel)
% Histogram with a density curve, mean and median, and a boxplot below
% sharing the same x axis.
% Input:  y      - target vector (NaN ignored)
%         name   - name of the target (for the title)
%         xLabel - label of the x axis
% Output: fig    - figure handle
    C = palette();
    y = y(~isnan(y));
    mu = mean(y);
    med = median(y);

    fig = newFigure([100 100 900 560]);
    tl = tiledlayout(fig, 5, 1, "TileSpacing", "none", "Padding", "compact");

    axH = nexttile(tl, [4 1]); hold(axH, "on");
    binw = 0.2;   % bin edges at x.x5 so each bin holds exactly two 0.1-step values
    histogram(axH, y, 0.05:binw:(max(y) + binw), ...
        "FaceColor", C.blue, "EdgeColor", "w", "FaceAlpha", 0.8);
    xi = linspace(min(y), max(y), 400);
    plot(axH, xi, ksdensity(y, xi) * numel(y) * binw, "Color", C.navy, "LineWidth", 2);
    hm = xline(axH, mu, "LineStyle", "--", "Color", C.red, "LineWidth", 1.5, ...
        "DisplayName", sprintf("Mean = %.2f", mu));
    hd = xline(axH, med, "LineStyle", ":", "Color", C.orange, "LineWidth", 2, ...
        "DisplayName", sprintf("Median = %.2f", med));
    legend(axH, [hm hd], "Box", "off");
    ylabel(axH, "Count");
    title(axH, "Distribution of " + name, "FontWeight", "bold");
    axH.XTickLabel = {};

    axB = nexttile(tl); hold(axB, "on");
    boxchart(axB, y, "Orientation", "horizontal", "BoxFaceColor", C.blue, ...
        "BoxFaceAlpha", 0.35, "BoxWidth", 0.5, "WhiskerLineColor", C.grey, ...
        "MarkerStyle", ".", "MarkerColor", C.blue);
    axB.YTick = [];
    axB.YColor = "none";
    axB.YGrid = "off";
    xlabel(axB, xLabel);
    linkaxes([axH axB], "x");
end

function fig = plotHistograms(T, vars)
% Grid of histograms (40 bins) with a density curve, one per variable.
% Input:  T    - data table (NaN ignored)
%         vars - string array of variable names (at most 12)
% Output: fig  - figure handle
    C = palette();
    fig = newFigure([50 50 1600 900]);
    tl = tiledlayout(fig, 3, 4, "TileSpacing", "compact", "Padding", "compact");
    for c = vars
        ax = nexttile(tl); hold(ax, "on");
        x = T.(c);
        x = x(~isnan(x));
        hh = histogram(ax, x, 40, "FaceColor", C.blue, "EdgeColor", "w", "FaceAlpha", 0.75);
        xi = linspace(min(x), max(x), 300);
        plot(ax, xi, ksdensity(x, xi) * numel(x) * hh.BinWidth, "Color", C.blue, "LineWidth", 1.5);
        title(ax, c, "FontWeight", "bold");
    end
end

function fig = plotSkewGroups(T, groups)
% One panel per group of variables with the density curve of each variable,
% standardised (z-score) for display only, so that only the shapes are
% compared. The skewness of each variable is shown in the legend.
% Input:  T      - data table (NaN ignored)
%         groups - n x 2 cell array: {panel title, string array of variables}
% Output: fig    - figure handle
    C = palette();
    colors = [C.blue; C.orange; C.green; C.red];
    styles = ["-", "--", "-.", ":"];
    n = size(groups, 1);

    fig = newFigure([50 50 1600 480]);
    tl = tiledlayout(fig, 1, n, "TileSpacing", "compact", "Padding", "compact");
    kde_axes = gobjects(1, n);
    for g = 1:n
        ax = nexttile(tl); hold(ax, "on");
        kde_axes(g) = ax;
        feats = groups{g, 2};
        for i = 1:numel(feats)
            x = T.(feats(i));
            x = x(~isnan(x));
            z = (x - mean(x)) / std(x);
            xi = linspace(min(z), max(z), 400);
            f = ksdensity(z, xi);
            fill(ax, [xi fliplr(xi)], [f zeros(size(f))], colors(i, :), ...
                "FaceAlpha", 0.08, "EdgeColor", "none", "HandleVisibility", "off");
            plot(ax, xi, f, "LineStyle", styles(i), "Color", colors(i, :), "LineWidth", 2.2, ...
                "DisplayName", sprintf("%s (skew = %.2f)", feats(i), skewness(x, 0)));
        end
        xline(ax, 0, "LineStyle", "--", "Color", [0.5 0.5 0.5], "HandleVisibility", "off");
        xlim(ax, [-3 6]);
        title(ax, groups{g, 1}, "FontWeight", "bold");
        xlabel(ax, "Standardised value (z-score)");
        legend(ax, "Box", "off", "Location", "northeast", "FontSize", 10);
    end
    ylabel(kde_axes(1), "Density");
    linkaxes(kde_axes, "y");
end

function fig = plotTargetCorrelation(rho, names, target)
% Horizontal bars of the correlation of each feature with the target,
% sorted by absolute value (strongest on top).
% Input:  rho    - vector of correlations (from targetCorrelation)
%         names  - string array with the feature names
%         target - name of the target (for the axis label)
% Output: fig    - figure handle
    C = palette();
    [~, idx] = sort(abs(rho), "descend");

    fig = newFigure([100 100 800 600]);
    ax = axes(fig); hold(ax, "on");
    barh(ax, rho(idx), 0.7, "FaceColor", C.blue, "EdgeColor", "none");
    xline(ax, 0, "Color", "k", "LineWidth", 0.8);
    ax.YTick = 1:numel(idx);
    ax.YTickLabel = cellstr(names(idx));
    ax.YDir = "reverse";
    ax.YGrid = "off";
    xlabel(ax, "Spearman correlation with " + target);
    title(ax, "Feature correlation with the target", "FontWeight", "bold");
end

function fig = plotScalingBeforeAfter(Xb, Xa, names)
% Boxplots of every feature before (original units) and after standardisation,
% side by side, with the +/-1 band marked on the standardised panel.
% Input:  Xb    - values before scaling (rows x variables)
%         Xa    - values after scaling (rows x variables)
%         names - string array with the variable names
% Output: fig   - figure handle
    C = palette();
    fig = newFigure([50 50 1500 650]);
    tl = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
    panels = {Xb, "Before scaling",        "Original value";
              Xa, "After standardisation", "Standardised value (z-score)"};
    box_axes = gobjects(1, 2);
    for p = 1:2
        ax = nexttile(tl); hold(ax, "on");
        box_axes(p) = ax;
        boxchart(ax, panels{p, 1}, "Orientation", "horizontal", "BoxFaceColor", C.blue, ...
            "BoxFaceAlpha", 0.35, "BoxWidth", 0.65, "LineWidth", 1.1, ...
            "WhiskerLineColor", C.grey, "MarkerStyle", ".", "MarkerColor", C.blue);
        ax.YTick = 1:numel(names);
        ax.YDir = "reverse";
        ax.YGrid = "off";
        if p == 1
            ax.YTickLabel = cellstr(names);
        else
            ax.YTickLabel = {};
        end
        title(ax, panels{p, 2}, "FontWeight", "bold");
        xlabel(ax, panels{p, 3});
    end
    linkaxes(box_axes, "y");

    ax = box_axes(2);
    yl = ylim(ax);
    band = patch(ax, [-1 1 1 -1], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
        "FaceAlpha", 0.08, "EdgeColor", "none");
    uistack(band, "bottom");
    xline(ax, 0, "LineStyle", "--", "Color", [0.33 0.33 0.33]);

    title(tl, "Effect of standardisation on the feature distributions", ...
        "FontWeight", "bold", "FontSize", 17);
end

function fig = plotFeaturesVsTarget(T, feats, target)
% Scatter of each feature against the target, with the mean of the target
% within 20 equal-sized (quantile) bins of the feature, which shows the shape
% of the relationship without assuming a straight line.
% Input:  T      - data table (original units)
%         feats  - string array with the features, in plotting order (at most 11)
%         target - name of the target
% Output: fig    - figure handle
    C = palette();
    y = T.(target);

    fig = newFigure([50 50 1600 950]);
    tl = tiledlayout(fig, 3, 4, "TileSpacing", "compact", "Padding", "compact");
    for i = 1:numel(feats)
        ax = nexttile(tl); hold(ax, "on");
        x = T.(feats(i));
        s = scatter(ax, x, y, 4, C.blue, "filled", "MarkerFaceAlpha", 0.15, "MarkerEdgeColor", "none");

        edges = unique(quantile(x, 0:0.05:1));      % unique: clipped values create repeated edges
        bin = discretize(x, edges);
        nb = numel(edges) - 1;
        xm = accumarray(bin, x, [nb 1], @mean, NaN);
        ym = accumarray(bin, y, [nb 1], @mean, NaN);
        ok = ~isnan(xm);
        t = plot(ax, xm(ok), ym(ok), "-o", "Color", C.red, "LineWidth", 2.2, ...
            "MarkerSize", 3, "MarkerFaceColor", C.red);

        title(ax, feats(i), "FontWeight", "bold", "FontSize", 12);
        if mod(i - 1, 4) == 0
            ylabel(ax, target + " (mg/m³)");
        end
        if i == 1
            ax1 = ax; hs = s; ht = t;
        end
    end
    lg = legend(ax1, [hs ht], ["Hourly observations", "Mean " + target + " per bin of the feature"], ...
        "Box", "off", "FontSize", 12);
    lg.Layout.Tile = 12;                % the empty 12th tile
    title(tl, "Features against the target", "FontWeight", "bold", "FontSize", 17);
end

function fig = plotScatterByTime(T, xVar, target)
% Scatter of one feature against the target, coloured by the number of days
% since the first measurement, to show whether the relationship changes over time.
% Input:  T      - data table with a Datetime column (original units)
%         xVar   - name of the feature
%         target - name of the target
% Output: fig    - figure handle
    fig = newFigure([100 100 800 550]);
    ax = axes(fig);
    scatter(ax, T.(xVar), T.(target), 4, days(T.Datetime - T.Datetime(1)), ...
        "filled", "MarkerFaceAlpha", 0.3);
    cb = colorbar(ax);
    cb.Label.String = "Days since the first measurement";
    xlabel(ax, xVar);
    ylabel(ax, target + " (mg/m³)");
    title(ax, xVar + " against " + target + ", coloured by date", "FontWeight", "bold");
end

function fig = plotTargetOverTime(T, target)
% Target over the whole period: hourly values on a complete hourly grid (so
% removed rows show as gaps instead of being joined by a line) and daily mean.
% Input:  T      - data table with a Datetime column
%         target - name of the target
% Output: fig    - figure handle
    C = palette();
    tt = table2timetable(T(:, ["Datetime", target]));
    hourly = retime(tt, "hourly", "fillwithmissing");
    daily = retime(tt, "daily", "mean");

    light_blue = C.blue + (1 - C.blue) * 0.65;

    fig = newFigure([100 100 1400 450]);
    ax = axes(fig); hold(ax, "on");
    plot(ax, hourly.Datetime, hourly.(target), "Color", light_blue, "LineWidth", 0.5, ...
        "DisplayName", "Hourly value");
    plot(ax, daily.Datetime, daily.(target), "Color", C.navy, "LineWidth", 1.8, ...
        "DisplayName", "Daily mean");
    ylabel(ax, target + " (mg/m³)");
    title(ax, target + " over the measurement period", "FontWeight", "bold");
    legend(ax, "Box", "off", "Location", "northwest");
    xlim(ax, [min(hourly.Datetime) max(hourly.Datetime)]);
end

function pairs = topCorrelatedPairs(T, vars, n)
% The n pairs of variables with the highest absolute Spearman correlation.
% Input:  T     - data table
%         vars  - string array of variable names
%         n     - number of pairs
% Output: pairs - n x 2 string array, one pair per row
    R = corr(T{:, vars}, "Type", "Spearman");
    R(tril(true(size(R)))) = 0;                 % keep each pair once
    [~, order] = sort(abs(R(:)), "descend");
    [ia, ib] = ind2sub(size(R), order(1:n));
    pairs = [vars(ia)', vars(ib)'];
end

function fig = plotPairs(T, pairs)
% Scatter plot of each pair of variables, one below the other, without axis
% values: a narrow band of points means one variable almost determines the other.
% Input:  T     - data table
%         pairs - n x 2 string array (from topCorrelatedPairs)
% Output: fig   - figure handle
    C = palette();
    fig = newFigure([100 50 550 1300]);
    tl = tiledlayout(fig, size(pairs, 1), 1, "TileSpacing", "compact", "Padding", "compact");
    for p = 1:size(pairs, 1)
        ax = nexttile(tl);
        scatter(ax, T.(pairs(p, 1)), T.(pairs(p, 2)), 4, C.blue, "filled", ...
            "MarkerFaceAlpha", 0.15, "MarkerEdgeColor", "none");
        xticks(ax, []);                 % no values on the axes (also removes the grid)
        yticks(ax, []);
        ax.FontSize = 16;               % axis titles scale with this
        xlabel(ax, pairs(p, 1));
        ylabel(ax, pairs(p, 2));
        title(ax, pairs(p, 1) + " vs " + pairs(p, 2), "FontWeight", "bold", "FontSize", 16);
    end
end

function C = palette()
% Colours used in every figure.
    C.blue   = [46 111 158] / 255;   % #2E6F9E
    C.navy   = [27  63  92] / 255;   % #1B3F5C
    C.red    = [192 57  43] / 255;   % #C0392B
    C.orange = [230 126 34] / 255;   % #E67E22
    C.green  = [39 174  96] / 255;   % #27AE60
    C.grey   = [0.27 0.27 0.27];
end

function fig = newFigure(pos)
% Creates a white figure with the plot defaults used in the project (light
% theme also when MATLAB uses a dark theme).
% Input:  pos - figure position [left bottom width height]
% Output: fig - figure handle
    fig = figure("Position", pos, "Color", "w", ...
        "defaultAxesFontSize", 14, "defaultAxesTickDir", "out", ...
        "defaultAxesBox", "off", "defaultAxesXGrid", "on", "defaultAxesYGrid", "on", ...
        "defaultAxesGridAlpha", 0.15);
    if exist("theme", "file")
        theme(fig, "light");
    end
end

%% ===================== Feature selection =====================

function R = featureRelevance(T, feats, target, nBins)
% Relevance of each feature for the target, measured three ways:
% Pearson (linear), Spearman (monotonic) and mutual information (any
% dependence, also non-monotonic). Sorted by absolute Spearman correlation.
% Input:  T      - data table without missing values
%         feats  - string array of feature names
%         target - name of the target
%         nBins  - number of equal-frequency bins for the mutual information
% Output: R      - table with Pearson, Spearman and MI (bits), one row per feature
    X = T{:, feats};
    y = T.(target);
    mi = zeros(numel(feats), 1);
    for j = 1:numel(feats)
        mi(j) = mutualInfo(X(:, j), y, nBins);
    end
    R = table(corr(X, y, "Type", "Pearson"), corr(X, y, "Type", "Spearman"), mi, ...
        'RowNames', cellstr(feats), 'VariableNames', {'Pearson', 'Spearman', 'MI'});
    [~, order] = sort(abs(R.Spearman), "descend");
    R = R(order, :);
end

function mi = mutualInfo(x, y, nBins)
% Mutual information (in bits) between two variables, estimated from a 2-D
% histogram with nBins equal-frequency bins per variable. Bins are built
% from the ranks, so tied values always fall in the same bin.
    n = numel(x);
    bx = min(nBins, ceil(nBins * tiedrank(x) / n));
    by = min(nBins, ceil(nBins * tiedrank(y) / n));
    pxy = accumarray([bx by], 1, [nBins nBins]) / n;   % joint probabilities
    pind = sum(pxy, 2) * sum(pxy, 1);                  % product of the marginals
    nz = pxy > 0;
    mi = sum(pxy(nz) .* log2(pxy(nz) ./ pind(nz)));
end
