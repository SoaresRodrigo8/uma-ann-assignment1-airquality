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
    L.plotRelevance    = @plotRelevance;
    L.correlationClusters = @correlationClusters;
    L.plotCorrClusters    = @plotCorrClusters;
    L.pcaAnalysis         = @pcaAnalysis;
    L.plotPCA             = @plotPCA;
    L.reliefWeights       = @reliefWeights;
    L.plotRelief          = @plotRelief;
    L.periodCorrelation   = @periodCorrelation;
    L.sameSeasonShift     = @sameSeasonShift;
    L.plotStability       = @plotStability;
    L.selectFeatures      = @selectFeatures;

    % --- Feature engineering ---
    L.inputGrid         = @inputGrid;
    L.hourlyProfile     = @hourlyProfile;
    L.plotHourlyProfile = @plotHourlyProfile;
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

function fig = plotRelevance(R, target)
% Relevance of each feature in two panels sharing the feature axis: Pearson
% and Spearman correlation with the target (left) and mutual information (right).
% Input:  R      - table from featureRelevance
%         target - name of the target
% Output: fig    - figure handle
    C = palette();
    names = R.Properties.RowNames;
    n = numel(names);

    fig = newFigure([100 100 1300 600]);
    tl = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");

    ax1 = nexttile(tl); hold(ax1, "on");
    b = barh(ax1, [R.Pearson, R.Spearman], "grouped", "EdgeColor", "none");
    b(1).FaceColor = C.orange;
    b(2).FaceColor = C.blue;
    xline(ax1, 0, "Color", "k", "LineWidth", 0.8, "HandleVisibility", "off");
    ax1.YTick = 1:n;
    ax1.YTickLabel = names;
    ax1.YDir = "reverse";
    ax1.YGrid = "off";
    xlim(ax1, [-1 1]);
    xlabel(ax1, "Correlation with " + target);
    legend(ax1, ["Pearson", "Spearman"], "Box", "off", "Location", "southwest");
    title(ax1, "Linear and monotonic relationship", "FontWeight", "bold");

    ax2 = nexttile(tl); hold(ax2, "on");
    barh(ax2, R.MI, 0.6, "FaceColor", C.green, "EdgeColor", "none");
    ax2.YTick = 1:n;
    ax2.YTickLabel = {};
    ax2.YDir = "reverse";
    ax2.YGrid = "off";
    xlabel(ax2, "Mutual information (bits)");
    title(ax2, "Any dependence (mutual information)", "FontWeight", "bold");
    linkaxes([ax1 ax2], "y");

    title(tl, "Relevance of each feature for " + target, "FontWeight", "bold", "FontSize", 17);
end

function C = correlationClusters(T, feats, minCorr)
% Hierarchical clustering of the features by redundancy, with distance
% 1 - |Spearman| (the sign does not matter: an inverse relation is also
% redundant) and complete linkage, so that inside each group every pair of
% features has |rho| >= minCorr.
% Input:  T       - data table without missing values
%         feats   - string array of feature names
%         minCorr - lowest |rho| between any two features of the same group
% Output: C       - struct with R (Spearman matrix), Z (linkage tree),
%                   order (leaf order, similar features next to each other)
%                   and groups (group number of each feature)
    R = corr(T{:, feats}, "Type", "Spearman");
    Dm = 1 - abs(R);
    Dm = (Dm + Dm') / 2;                    % exactly symmetric
    Dm(1:size(Dm, 1) + 1:end) = 0;          % exactly zero diagonal
    D = squareform(Dm);
    C.R = R;
    C.Z = linkage(D, "complete");
    C.order = optimalleaforder(C.Z, D);
    C.groups = cluster(C.Z, "Cutoff", 1 - minCorr, "Criterion", "distance");
end

function fig = plotCorrClusters(C, names, minCorr)
% Dendrogram of the feature clustering (top) aligned with the heatmap of the
% Spearman matrix in the same order (bottom). The cut that defines the
% groups of redundant features is marked on the dendrogram.
% Input:  C       - struct from correlationClusters
%         names   - string array with the feature names
%         minCorr - |rho| used to cut the tree
% Output: fig     - figure handle
    P = palette();
    n = numel(names);
    names = names(C.order);
    R = C.R(C.order, C.order);

    fig = newFigure([100 50 1000 1050]);
    tl = tiledlayout(fig, 4, 1, "TileSpacing", "compact", "Padding", "compact");

    % Dendrogram, with the height shown as |rho| (1 - distance)
    axD = nexttile(tl);
    h = dendrogram(C.Z, 0, "Reorder", C.order);
    set(h, "Color", P.navy, "LineWidth", 1.5);
    hold(axD, "on");
    yline(axD, 1 - minCorr, "--", sprintf("|\\rho| = %.2f", minCorr), "Color", P.red, ...
        "LineWidth", 1.2, "LabelHorizontalAlignment", "left", "FontSize", 11);
    xlim(axD, [0.5 n + 0.5]);
    axD.XTick = [];
    axD.XGrid = "off";
    yt = 0:0.2:1;
    axD.YTick = yt;
    axD.YTickLabel = compose("%.1f", 1 - yt);
    ylabel(axD, "|\rho| at which groups join");
    title(axD, "Redundancy between features (Spearman, complete linkage)", "FontWeight", "bold");

    % Heatmap of the correlation matrix, values written in each cell
    axH = nexttile(tl, [3 1]);
    imagesc(axH, R);
    t = linspace(0, 1, 128)';
    cmap = [P.blue + (1 - P.blue) .* t; 1 - (1 - P.red) .* t];   % blue - white - red
    colormap(axH, cmap);
    clim(axH, [-1 1]);
    cb = colorbar(axH);
    cb.Layout.Tile = "east";
    cb.Label.String = "Spearman \rho";
    for i = 1:n
        for j = 1:n
            if abs(R(i, j)) > 0.6, tc = "w"; else, tc = "k"; end
            text(axH, j, i, sprintf("%.2f", R(i, j)), "HorizontalAlignment", "center", ...
                "FontSize", 10, "Color", tc);
        end
    end
    axH.XTick = 1:n;
    axH.YTick = 1:n;
    axH.XTickLabel = cellstr(names);
    axH.YTickLabel = cellstr(names);
    axH.XTickLabelRotation = 45;
    axH.XGrid = "off";
    axH.YGrid = "off";
    axH.TickLength = [0 0];
end

function P = pcaAnalysis(T, feats)
% Principal component analysis, step by step: (1) standardise each feature
% (mean 0, variance 1), (2) covariance matrix, (3) eigendecomposition,
% (4) eigenvalues sorted from largest to smallest and explained variance.
% Input:  T     - data table without missing values
%         feats - string array of feature names
% Output: P     - struct with eigval, explained and cumulative (% of the
%                 total variance), coeff (loadings: one column per
%                 component, one row per feature) and score (data projected
%                 on the components)
    X = T{:, feats};
    Z = (X - mean(X)) ./ std(X);              % (1) standardise
    S = cov(Z);                               % (2) covariance matrix
    [V, D] = eig(S);                          % (3) eigenvectors and eigenvalues
    [P.eigval, order] = sort(diag(D), "descend");
    V = V(:, order);
    % the sign of an eigenvector is arbitrary: make the largest loading positive
    [~, imax] = max(abs(V));
    V = V .* sign(V(sub2ind(size(V), imax, 1:size(V, 2))));
    P.coeff = V;
    P.explained = 100 * P.eigval / sum(P.eigval);   % (4) % of the total variance
    P.cumulative = cumsum(P.explained);
    P.score = Z * V;
end

function fig = plotPCA(P, names, threshold)
% Two panels: scree plot (variance explained by each component, with the
% cumulative curve and the threshold used to choose the number of
% components) and loadings of each feature on the first two components.
% Input:  P         - struct from pcaAnalysis
%         names     - string array with the feature names
%         threshold - cumulative explained variance (%) to mark
% Output: fig       - figure handle
    C = palette();
    n = numel(P.eigval);

    fig = newFigure([100 100 1400 600]);
    tl = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");

    % Scree plot
    ax1 = nexttile(tl); hold(ax1, "on");
    bar(ax1, 1:n, P.explained, 0.6, "FaceColor", C.blue, "EdgeColor", "none", ...
        "DisplayName", "Explained by the component");
    plot(ax1, 1:n, P.cumulative, "-o", "Color", C.red, "LineWidth", 2, ...
        "MarkerFaceColor", C.red, "MarkerSize", 5, "DisplayName", "Cumulative");
    yline(ax1, threshold, "--", sprintf("%d%%", threshold), "Color", C.grey, ...
        "LineWidth", 1.2, "HandleVisibility", "off", "LabelHorizontalAlignment", "right");
    xlim(ax1, [0.4 n + 0.6]);
    ylim(ax1, [0 105]);
    ax1.XTick = 1:n;
    ax1.XGrid = "off";
    xlabel(ax1, "Principal component");
    ylabel(ax1, "Variance explained (%)");
    legend(ax1, "Box", "off", "Location", "east");
    title(ax1, "Scree plot", "FontWeight", "bold");

    % Loadings on PC1 and PC2: one arrow per feature
    ax2 = nexttile(tl); hold(ax2, "on");
    xy = P.coeff(:, 1:2);
    for j = 1:n
        plot(ax2, [0 xy(j, 1)], [0 xy(j, 2)], "-", "Color", C.blue, "LineWidth", 1.5);
        plot(ax2, xy(j, 1), xy(j, 2), "o", "Color", C.blue, ...
            "MarkerFaceColor", C.blue, "MarkerSize", 5);
        % two almost equal arrows: label the higher one above, the lower one below
        d = vecnorm(xy - xy(j, :), 2, 2);
        d(j) = Inf;
        [dmin, k] = min(d);
        va = "middle";
        if dmin < 0.06
            if xy(j, 2) >= xy(k, 2), va = "bottom"; else, va = "top"; end
        end
        text(ax2, xy(j, 1), xy(j, 2), "  " + names(j), "FontSize", 11, ...
            "HorizontalAlignment", "left", "VerticalAlignment", va);
    end
    xline(ax2, 0, "Color", [0.6 0.6 0.6]);
    yline(ax2, 0, "Color", [0.6 0.6 0.6]);
    axis(ax2, "equal");
    xlim(ax2, [-0.5 0.75]);
    ylim(ax2, [-0.4 0.7]);
    xlabel(ax2, sprintf("PC1 (%.1f%% of the variance)", P.explained(1)));
    ylabel(ax2, sprintf("PC2 (%.1f%% of the variance)", P.explained(2)));
    title(ax2, "Loadings on the first two components", "FontWeight", "bold");

    title(tl, "Principal component analysis of the features", "FontWeight", "bold", "FontSize", 17);
end

function W = reliefWeights(T, feats, target, k, seed)
% RReliefF weights (the regression version of Relief): for each hour, its k
% nearest hours in the space of ALL the features are found, and a feature
% gains weight when its differences between neighbours go together with
% differences in the target. Because it compares neighbours in all the
% features at once, it can detect features that only help in combination
% with others. A column of random noise is added as a reference: a useful
% feature must get a clearly higher weight than noise.
% Input:  T      - data table without missing values
%         feats  - string array of feature names
%         target - name of the target
%         k      - number of nearest neighbours
%         seed   - seed of the random noise column (reproducible results)
% Output: W      - table with the Weight of each feature and of the noise
%                  column ("Random noise"), sorted from highest to lowest
    s = RandStream("mt19937ar", "Seed", seed);     % local generator, global rng untouched
    noise = randn(s, height(T), 1);
    [~, w] = relieff([T{:, feats}, noise], T.(target), k);
    W = table(w', 'RowNames', [cellstr(feats), {'Random noise'}], 'VariableNames', {'Weight'});
    W = sortrows(W, "Weight", "descend");
end

function fig = plotRelief(W, k)
% Horizontal bars of the RReliefF weights, with the random-noise reference
% in grey and marked by a dashed line.
% Input:  W   - table from reliefWeights
%         k   - number of neighbours used (for the title)
% Output: fig - figure handle
    C = palette();
    names = W.Properties.RowNames;
    n = height(W);
    is_noise = strcmp(names, 'Random noise');

    fig = newFigure([100 100 900 600]);
    ax = axes(fig); hold(ax, "on");
    b = barh(ax, W.Weight, 0.65, "FaceColor", "flat", "EdgeColor", "none");
    b.CData = repmat(C.blue, n, 1);
    b.CData(is_noise, :) = [0.65 0.65 0.65];
    xline(ax, W.Weight(is_noise), "--", "Noise level", "Color", C.grey, "LineWidth", 1.2, ...
        "LabelVerticalAlignment", "bottom", "FontSize", 11);
    xline(ax, 0, "Color", "k", "LineWidth", 0.8);
    xlim(ax, [1.8 * min(0, min(W.Weight)), 1.08 * max(W.Weight)]);   % room for the noise label
    ax.YTick = 1:n;
    ax.YTickLabel = names;
    ax.YDir = "reverse";
    ax.YGrid = "off";
    xlabel(ax, "RReliefF weight");
    title(ax, sprintf("RReliefF weights (k = %d neighbours)", k), "FontWeight", "bold");
end

function P = periodCorrelation(T, feats, target, splitDate)
% Spearman correlation of each feature with the target before and after a
% date, to see which relationships change over time.
% Input:  T         - data table with a Datetime column
%         feats     - string array of feature names
%         target    - name of the target
%         splitDate - first datetime of the second period
% Output: P         - table with Before, After and Change (After - Before)
    late = T.Datetime >= splitDate;
    X = T{:, feats};
    y = T.(target);
    before = corr(X(~late, :), y(~late), "Type", "Spearman");
    after  = corr(X(late, :),  y(late),  "Type", "Spearman");
    P = table(before, after, after - before, 'RowNames', cellstr(feats), ...
        'VariableNames', {'Before', 'After', 'Change'});
end

function S = sameSeasonShift(T, feats, target, months, years, targetRange)
% Median of each feature in the same months of two different years, using
% only the hours with the target inside the same range. Same season and same
% pollution level: a feature that still changes has changed its relationship
% with the target (drift of the sensor, recalibration or a real change).
% Input:  T           - data table with a Datetime column
%         feats       - string array of feature names
%         target      - name of the target
%         months      - months to compare (e.g. [3 4] for March-April)
%         years       - the two years, [first second]
%         targetRange - [min max] of the target in the hours used
% Output: S           - table with the median in each year and the % change;
%                       S.Properties.UserData has the number of hours used
    y = T.(target);
    mid = y >= targetRange(1) & y <= targetRange(2) & ismember(T.Datetime.Month, months);
    in1 = mid & T.Datetime.Year == years(1);
    in2 = mid & T.Datetime.Year == years(2);
    m1 = median(T{in1, feats})';
    m2 = median(T{in2, feats})';
    S = table(m1, m2, 100 * (m2 - m1) ./ m1, 'RowNames', cellstr(feats), ...
        'VariableNames', {sprintf('Median%d', years(1)), sprintf('Median%d', years(2)), 'PctChange'});
    S.Properties.UserData = [sum(in1), sum(in2)];
end

function fig = plotStability(S, maxShift)
% Horizontal bars of the % change of each feature between the two years of
% sameSeasonShift, with the band of changes considered stable shaded.
% Input:  S        - table from sameSeasonShift
%         maxShift - largest |% change| considered stable
% Output: fig      - figure handle
    C = palette();
    [~, order] = sort(abs(S.PctChange), "descend");
    S = S(order, :);
    names = S.Properties.RowNames;
    n = height(S);
    unstable = abs(S.PctChange) > maxShift;

    fig = newFigure([100 100 950 600]);
    ax = axes(fig); hold(ax, "on");
    patch(ax, [-maxShift maxShift maxShift -maxShift], [0.4 0.4 n + 0.6 n + 0.6], C.green, ...
        "FaceAlpha", 0.12, "EdgeColor", "none");
    b = barh(ax, S.PctChange, 0.65, "FaceColor", "flat", "EdgeColor", "none");
    b.CData = repmat(C.blue, n, 1);
    b.CData(unstable, :) = repmat(C.red, sum(unstable), 1);
    xline(ax, 0, "Color", "k", "LineWidth", 0.8);
    for j = 1:n
        v = round(S.PctChange(j));
        if v == 0, lbl = "0%"; else, lbl = sprintf("%+d%%", v); end   % no "-0%"
        if S.PctChange(j) >= 0, ha = "left"; else, ha = "right"; end
        text(ax, S.PctChange(j), j, "  " + lbl + "  ", "HorizontalAlignment", ha, "FontSize", 11);
    end
    ax.YTick = 1:n;
    ax.YTickLabel = names;
    ax.YDir = "reverse";
    ax.YGrid = "off";
    xlim(ax, [min(-60, 1.3 * min(S.PctChange)), 1.2 * max(S.PctChange)]);
    ylim(ax, [0.4 n + 0.6]);
    xlabel(ax, "Change of the median, same season and same CO level (%)");
    title(ax, sprintf("Stability of the features over one year (stable band: \\pm%d%%)", maxShift), ...
        "FontWeight", "bold");
end

function D = selectFeatures(D, minCorr, maxShift)
% Applies the selection rules to a table of criteria, in this order:
%   1. unstable: |% change| in the same-season test above maxShift;
%   2. no added information: RReliefF weight not above the random-noise
%      weight for every number of neighbours tested;
%   3. redundant: |rho| >= minCorr with another feature that is more stable
%      (smaller |% change|), which is the one kept.
% A feature is kept only if no rule applies; every rule that applies is
% listed in Reason.
% Input:  D        - table, one row per feature, with variables Shift (%),
%                    AboveNoise (logical), MaxCorr (|rho|) and MostSimilar
%                    (name of the feature with that |rho|)
%         minCorr  - |rho| from which two features are redundant
%         maxShift - largest |% change| considered stable
% Output: D        - same table with Keep (logical) and Reason (string)
    names = string(D.Properties.RowNames);
    n = height(D);
    reason = strings(n, 1);
    for j = 1:n
        r = strings(0);
        if abs(D.Shift(j)) > maxShift
            r(end + 1) = sprintf("unstable (%+.0f%%)", D.Shift(j));
        end
        if ~D.AboveNoise(j)
            r(end + 1) = "no information above noise (RReliefF)";
        end
        k = find(names == D.MostSimilar(j));
        if D.MaxCorr(j) >= minCorr && abs(D.Shift(k)) < abs(D.Shift(j))
            r(end + 1) = sprintf("redundant with %s (|rho| = %.3f, more stable)", names(k), D.MaxCorr(j));
        end
        reason(j) = strjoin(r, "; ");
    end
    D.Keep = reason == "";
    D.Reason = reason;
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

%% ===================== Feature engineering =====================

function [Xg, gridTime] = inputGrid(T, vars, lo, hi)
% Places the inputs of every measured hour on a complete hourly grid,
% clipped to the given bounds, so that windows and hour-to-hour changes use
% real hours. Hours without a row, or with a missing value, stay NaN.
% Input:  T        - data table with a sorted Datetime column
%         vars     - string array of input names
%         lo/hi    - row vectors with the clipping bounds of each input
% Output: Xg       - matrix (hours x inputs) on the grid
%         gridTime - datetime of each hour of the grid
    [hidx, ~] = hourlyIndex(T.Datetime);
    gridTime = T.Datetime(1) + hours(0:hidx(end) - 1)';
    Xg = NaN(hidx(end), numel(vars));
    Xg(hidx, :) = clipToBounds(T{:, vars}, lo, hi);
end

function P = hourlyProfile(T, target)
% Median and interquartile range of the target for each hour of the day,
% separately for weekdays and weekends.
% Input:  T      - data table with a Datetime column
%         target - name of the target
% Output: P      - table, one row per hour (0-23), with the median, Q1 and Q3
%                  of the target on weekdays and on weekends
    h = T.Datetime.Hour;
    we = isweekend(T.Datetime);
    y = T.(target);
    S = NaN(24, 6);
    for k = 0:23
        yd = y(h == k & ~we);
        ye = y(h == k & we);
        S(k + 1, :) = [quantile7(yd, 0.5), quantile7(yd, 0.25), quantile7(yd, 0.75), ...
                       quantile7(ye, 0.5), quantile7(ye, 0.25), quantile7(ye, 0.75)];
    end
    P = array2table([(0:23)', S], 'VariableNames', {'Hour', 'WeekdayMedian', 'WeekdayQ1', ...
        'WeekdayQ3', 'WeekendMedian', 'WeekendQ1', 'WeekendQ3'});
end

function fig = plotHourlyProfile(P, target)
% Daily profile of the target: median per hour as a line and the
% interquartile range as a shaded band, for weekdays and weekends.
% Input:  P      - table from hourlyProfile
%         target - name of the target
% Output: fig    - figure handle
    C = palette();
    fig = newFigure([100 100 1100 650]);
    ax = axes(fig); hold(ax, "on");
    groups = ["Weekday", "Weekend"];
    colors = [C.blue; C.orange];
    for g = 1:2
        q1 = P.(groups(g) + "Q1");
        q3 = P.(groups(g) + "Q3");
        fill(ax, [P.Hour; flipud(P.Hour)], [q1; flipud(q3)], colors(g, :), ...
            "FaceAlpha", 0.15, "EdgeColor", "none", "HandleVisibility", "off");
        plot(ax, P.Hour, P.(groups(g) + "Median"), "-o", "Color", colors(g, :), ...
            "LineWidth", 2, "MarkerFaceColor", colors(g, :), "MarkerSize", 6, ...
            "DisplayName", groups(g) + " (median)");
    end
    xlim(ax, [0 23]);
    xticks(ax, 0:2:22);
    ylim(ax, [0 ceil(max([P.WeekdayQ3; P.WeekendQ3]))]);
    xlabel(ax, "Hour of day");
    ylabel(ax, target + " (mg/m³)");
    title(ax, "Daily profile of " + target + ": weekdays vs weekends", "FontWeight", "bold");
    subtitle(ax, "Lines: median; shaded bands: interquartile range");
    legend(ax, "Box", "off", "Location", "northwest");
end
