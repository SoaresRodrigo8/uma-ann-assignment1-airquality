%% Air Quality (UCI): cleaning, preprocessing and plots
% MATLAB version of main.py
% Requires MATLAB R2021a or newer and the Statistics and Machine Learning
% Toolbox (ksdensity, skewness, corr).

clear; clc; close all;

% -------------------------------------------------------------------------
% Paths, relative to the folder this script is in (the repository root)
% -------------------------------------------------------------------------
root = fileparts(mfilename("fullpath"));
if root == "" || startsWith(root, tempdir)          % e.g. when running a single section
    try
        root = fileparts(matlab.desktop.editor.getActiveFilename);
    catch
        root = pwd;
    end
end

raw_path = fullfile(root, "data", "Raw", "AirQualityUCI.csv");
out_dir  = fullfile(root, "data", "Processed");
fig_dir  = fullfile(root, "Figures");
if ~isfolder(out_dir), mkdir(out_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end

% -------------------------------------------------------------------------
% Colours and plot defaults
% -------------------------------------------------------------------------
blue   = [46 111 158] / 255;   % #2E6F9E
navy   = [27  63  92] / 255;   % #1B3F5C
red    = [192 57  43] / 255;   % #C0392B
orange = [230 126 34] / 255;   % #E67E22
green  = [39 174  96] / 255;   % #27AE60
grey   = [0.27 0.27 0.27];

set(groot, "defaultAxesFontSize", 14, "defaultAxesTickDir", "out", ...
    "defaultAxesBox", "off", "defaultAxesXGrid", "on", "defaultAxesYGrid", "on", ...
    "defaultAxesGridAlpha", 0.15, "defaultFigureColor", "w");

% -------------------------------------------------------------------------
% Column lists
% -------------------------------------------------------------------------
target_col = "CO(GT)";
all_num = ["CO(GT)", "PT08.S1(CO)", "NMHC(GT)", "C6H6(GT)", "PT08.S2(NMHC)", ...
           "NOx(GT)", "PT08.S3(NOx)", "NO2(GT)", "PT08.S4(NO2)", "PT08.S5(O3)", ...
           "T", "RH", "AH"];
feature_cols = ["PT08.S1(CO)", "C6H6(GT)", "PT08.S2(NMHC)", "NOx(GT)", "PT08.S3(NOx)", ...
                "NO2(GT)", "PT08.S4(NO2)", "PT08.S5(O3)", "T", "RH", "AH"];
cols = [target_col, feature_cols];
sensor_cols = ["PT08.S1(CO)", "PT08.S2(NMHC)", "C6H6(GT)", "PT08.S3(NOx)", ...
               "PT08.S4(NO2)", "PT08.S5(O3)", "T", "RH", "AH"];
targets = ["NOx(GT)", "NO2(GT)"];   % columns that still need imputation
max_gap = 3;                        % gaps up to 3 hours are interpolated

%% 1. Load the data
% Semicolon separator, comma as decimal mark, original column names kept.
opts = detectImportOptions(raw_path, "Delimiter", ";", "VariableNamingRule", "preserve");
opts.SelectedVariableNames = ["Date", "Time", all_num];   % drops the empty trailing columns
opts = setvartype(opts, ["Date", "Time"], "string");
opts = setvartype(opts, all_num, "double");
opts = setvaropts(opts, all_num, "DecimalSeparator", ",");
df = readtable(raw_path, opts);

df = df(strlength(df.Date) > 0, :);    % remove blank rows at the end of the file

X = df{:, all_num};
X(X == -200) = NaN;                    % -200 is the missing-value code
df{:, all_num} = X;
df_before = df;

df.Datetime = datetime(df.Date + " " + df.Time, "InputFormat", "dd/MM/yyyy HH.mm.ss");

%% 2. Timestamp checks
ts = df.Datetime;
h = hours(ts - min(ts));
fprintf("Duplicates: %d\n", numel(ts) - numel(unique(ts)));
fprintf("Sorted: %d\n", issorted(ts));
fprintf("Range: %s to %s\n", string(min(ts)), string(max(ts)));
fprintf("Expected hours: %d, present: %d, missing: %d\n", ...
    round(max(h)) + 1, numel(unique(round(h))), round(max(h)) + 1 - numel(unique(round(h))));
[step_h, ~, ic] = unique(hours(diff(ts)));
disp(table(step_h, accumarray(ic, 1), 'VariableNames', {'Step_hours', 'Count'}))

%% 3. Missing values per column
miss = sum(isnan(df{:, all_num}))';
missing_tbl = table(miss, round(miss / height(df) * 100, 2), ...
    'RowNames', cellstr(all_num), 'VariableNames', {'Missing', 'PctMissing'});
missing_tbl = sortrows(missing_tbl, 'Missing', 'descend');
fprintf("\nTotal missing values: %d\n", sum(miss));
disp(missing_tbl)

%% 4. Distribution of the target
target = df.(target_col);
target = target(~isnan(target));
mu = mean(target);
med = median(target);

fig = figure("Position", [100 100 900 560]);
tl = tiledlayout(fig, 5, 1, "TileSpacing", "none", "Padding", "compact");

axH = nexttile(tl, [4 1]); hold(axH, "on");
binw = 0.2;   % bin edges at x.x5 so each bin holds exactly two 0.1-step values
histogram(axH, target, 0.05:binw:(max(target) + binw), ...
    "FaceColor", blue, "EdgeColor", "w", "FaceAlpha", 0.8);
xi = linspace(min(target), max(target), 400);
plot(axH, xi, ksdensity(target, xi) * numel(target) * binw, "Color", navy, "LineWidth", 2);
hm = xline(axH, mu, "LineStyle", "--", "Color", red, "LineWidth", 1.5, ...
    "DisplayName", sprintf("Mean = %.2f", mu));
hd = xline(axH, med, "LineStyle", ":", "Color", orange, "LineWidth", 2, ...
    "DisplayName", sprintf("Median = %.2f", med));
legend(axH, [hm hd], "Box", "off");
ylabel(axH, "Count");
title(axH, "Distribution of CO(GT)", "FontWeight", "bold");
axH.XTickLabel = {};

axB = nexttile(tl); hold(axB, "on");
boxchart(axB, target, "Orientation", "horizontal", "BoxFaceColor", blue, ...
    "BoxFaceAlpha", 0.35, "BoxWidth", 0.5, "WhiskerLineColor", grey, ...
    "MarkerStyle", ".", "MarkerColor", blue);
axB.YTick = [];
axB.YColor = "none";
axB.YGrid = "off";
xlabel(axB, "CO concentration (mg/m³)");
linkaxes([axH axB], "x");

exportgraphics(fig, fullfile(fig_dir, "co_distribution.pdf"), "ContentType", "vector");

%% 5. Descriptive statistics (Table II)
S = df{:, cols};
stats = [mean(S, "omitnan")', std(S, "omitnan")', min(S, [], "omitnan")', ...
         median(S, "omitnan")', max(S, [], "omitnan")', skewness(S, 0)'];
disp(array2table(round(stats, 2), 'RowNames', cellstr(cols), ...
    'VariableNames', {'mean', 'std', 'min', 'median', 'max', 'skew'}))

printLatex(cols, stats, " & mean & std & min & median & max & skew \\", ...
    "Descriptive statistics of the target and features.", "tab:descriptive_stats");

%% 6. Histograms of the target and every feature
fig = figure("Position", [50 50 1600 900]);
tl = tiledlayout(fig, 3, 4, "TileSpacing", "compact", "Padding", "compact");
for c = cols
    ax = nexttile(tl); hold(ax, "on");
    x = df.(c);
    x = x(~isnan(x));
    hh = histogram(ax, x, 40, "FaceColor", blue, "EdgeColor", "w", "FaceAlpha", 0.75);
    xi = linspace(min(x), max(x), 300);
    plot(ax, xi, ksdensity(x, xi) * numel(x) * hh.BinWidth, "Color", blue, "LineWidth", 1.5);
    title(ax, c, "FontWeight", "bold");
end
exportgraphics(fig, fullfile(fig_dir, "feature_distributions.pdf"), "ContentType", "vector");

%% 7. Features grouped by distribution shape (standardised KDEs)
groups = {"Strongly right-skewed",   ["NOx(GT)", "C6H6(GT)", "PT08.S3(NOx)"];
          "Moderately right-skewed", ["PT08.S1(CO)", "PT08.S2(NMHC)", "PT08.S5(O3)", "NO2(GT)"];
          "Approximately symmetric", ["PT08.S4(NO2)", "T", "RH", "AH"]};
palette = [blue; orange; green; red];
styles = ["-", "--", "-.", ":"];

fig = figure("Position", [50 50 1600 480]);
tl = tiledlayout(fig, 1, 3, "TileSpacing", "compact", "Padding", "compact");
kde_axes = gobjects(1, 3);
for g = 1:3
    ax = nexttile(tl); hold(ax, "on");
    kde_axes(g) = ax;
    feats = groups{g, 2};
    for i = 1:numel(feats)
        x = df.(feats(i));
        x = x(~isnan(x));
        z = (x - mean(x)) / std(x);
        xi = linspace(min(z), max(z), 400);
        f = ksdensity(z, xi);
        fill(ax, [xi fliplr(xi)], [f zeros(size(f))], palette(i, :), ...
            "FaceAlpha", 0.08, "EdgeColor", "none", "HandleVisibility", "off");
        plot(ax, xi, f, "LineStyle", styles(i), "Color", palette(i, :), "LineWidth", 2.2, ...
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
exportgraphics(fig, fullfile(fig_dir, "feature_distributions_grouped.pdf"), "ContentType", "vector");


%% 8. Spearman correlation with the target
rho = corr(df{:, feature_cols}, df.(target_col), "Type", "Spearman", "Rows", "pairwise");
[~, idx] = sort(abs(rho), "descend");

fig = figure("Position", [100 100 800 600]);
ax = axes(fig); hold(ax, "on");
barh(ax, rho(idx), 0.7, "FaceColor", blue, "EdgeColor", "none");
xline(ax, 0, "Color", "k", "LineWidth", 0.8);
ax.YTick = 1:numel(idx);
ax.YTickLabel = cellstr(feature_cols(idx));
ax.YDir = "reverse";
ax.YGrid = "off";
xlabel(ax, "Spearman correlation with CO(GT)");
title(ax, "Feature correlation with the target", "FontWeight", "bold");
exportgraphics(fig, fullfile(fig_dir, "corr_with_target.pdf"), "ContentType", "vector");

%% 9. Outliers (IQR rule) on the data before cleaning
X = df{:, all_num};
q1 = quantile7(X, 0.25);
q3 = quantile7(X, 0.75);
iqrv = q3 - q1;
n_upper = sum(X > q3 + 1.5 * iqrv)';
n_lower = sum(X < q1 - 1.5 * iqrv)';
n_valid = sum(~isnan(X))';
disp(table(n_upper + n_lower, round((n_upper + n_lower) ./ n_valid * 100, 2), n_upper, n_lower, ...
    'RowNames', cellstr(all_num), 'VariableNames', {'Outliers', 'Pct', 'Upper', 'Lower'}))

%% 10. Remove NMHC(GT), sensor-offline rows and rows with a missing target
df = removevars(df, "NMHC(GT)");

offline = all(isnan(df{:, sensor_cols}), 2);
fprintf("Sensor-offline rows dropped: %d\n", sum(offline));
df = df(~offline, :);

no_target = isnan(df.(target_col));
fprintf("Missing-target rows dropped: %d\n", sum(no_target));
df = df(~no_target, :);
fprintf("Rows remaining: %d\n", height(df));

df = sortrows(df, "Datetime");

%% 11. Clipping to the IQR bounds (features only, not the target)
X = df{:, feature_cols};
q1 = quantile7(X, 0.25);
q3 = quantile7(X, 0.75);
iqrv = q3 - q1;
lo = q1 - 1.5 * iqrv;
hi = q3 + 1.5 * iqrv;

n_clip = sum(X < lo | X > hi)';
disp(table(n_clip, round(n_clip ./ sum(~isnan(X))' * 100, 2), ...
    'RowNames', cellstr(feature_cols), 'VariableNames', {'Clipped', 'Pct'}))
fprintf("%.2f%% of all feature values clipped\n", sum(n_clip) / sum(~isnan(X), "all") * 100);

% "includenan" keeps missing values missing (by default min/max would replace them)
df{:, feature_cols} = min(max(X, lo, "includenan"), hi, "includenan");

%% 12. Gap analysis on a complete hourly grid
hidx = round(hours(df.Datetime - df.Datetime(1))) + 1;   % position of each row on the grid
n_grid = hidx(end);
in_df = false(n_grid, 1);
in_df(hidx) = true;

fprintf("\nMissing before imputation:\n");
disp(array2table(sum(isnan(df{:, targets})), 'VariableNames', cellstr(targets)))

for c = targets
    v = NaN(n_grid, 1);
    v(hidx) = df.(c);
    [starts, ends] = nanRuns(isnan(v));
    run_len = ends - starts + 1;
    needs_fill = arrayfun(@(s, e) sum(in_df(s:e)), starts, ends);
    rel = needs_fill > 0;
    short = run_len <= max_gap;
    fprintf("%s: %d relevant gaps, longest = %d h, median = %g h | " + ...
            "values in short gaps: %d, in long gaps: %d\n", c, sum(rel), ...
            max(run_len(rel)), median(run_len(rel)), ...
            sum(needs_fill(rel & short)), sum(needs_fill(rel & ~short)));
end

%% 13. Imputation
% Step 1: linear interpolation over time for gaps of up to max_gap hours
for c = targets
    v = NaN(n_grid, 1);
    v(hidx) = df.(c);
    na = isnan(v);
    [starts, ends] = nanRuns(na);
    run_len = zeros(n_grid, 1);
    for i = 1:numel(starts)
        run_len(starts(i):ends(i)) = ends(i) - starts(i) + 1;
    end
    v_interp = fillmissing(v, "linear", "EndValues", "none");   % only gaps with values on both sides
    fill_mask = na & run_len <= max_gap;
    v(fill_mask) = v_interp(fill_mask);
    df.(c) = v(hidx);
end
fprintf("\nLeft for KNN:\n");
disp(array2table(sum(isnan(df{:, targets})), 'VariableNames', cellstr(targets)))

% Step 2: KNN (k = 5) for the remaining long-gap values, target not used
X = df{:, feature_cols};
mu_knn = mean(X, "omitnan");
sd_knn = std(X, 1, "omitnan");
Xs = knnImpute((X - mu_knn) ./ sd_knn, 5);
X_knn = Xs .* sd_knn + mu_knn;
X(isnan(X)) = X_knn(isnan(X));          % only missing values are replaced
df{:, feature_cols} = X;
fprintf("Missing after imputation: %d\n", sum(isnan(X), "all"));

%% 14. Standardisation
df_before_scaling = df;

X = df{:, feature_cols};
scale_mu = mean(X);
scale_sd = std(X, 1);                   % population std, as in sklearn's StandardScaler
df{:, feature_cols} = (X - scale_mu) ./ scale_sd;

Z = df{:, feature_cols};
disp(table(round(mean(Z), 3)', round(std(Z), 3)', ...
    'RowNames', cellstr(feature_cols), 'VariableNames', {'mean', 'std'}))

%% 15. Before/after standardisation plot
fig = figure("Position", [50 50 1500 650]);
tl = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
panels = {df_before_scaling{:, feature_cols}, "Before scaling",        "Original value";
          df{:, feature_cols},                "After standardisation", "Standardised value (z-score)"};
box_axes = gobjects(1, 2);
for p = 1:2
    ax = nexttile(tl); hold(ax, "on");
    box_axes(p) = ax;
    boxchart(ax, panels{p, 1}, "Orientation", "horizontal", "BoxFaceColor", blue, ...
        "BoxFaceAlpha", 0.35, "BoxWidth", 0.65, "LineWidth", 1.1, ...
        "WhiskerLineColor", grey, "MarkerStyle", ".", "MarkerColor", blue);
    ax.YTick = 1:numel(feature_cols);
    ax.YDir = "reverse";
    ax.YGrid = "off";
    if p == 1
        ax.YTickLabel = cellstr(feature_cols);
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
exportgraphics(fig, fullfile(fig_dir, "scaling_before_after.pdf"), "ContentType", "vector");

%% 16. Before/after cleaning comparison table
before = df_before{:, cols};
after = df_before_scaling{:, cols};
cmp = [std(before, "omitnan")', std(after)', ...
       max(before, [], "omitnan")', max(after)', ...
       skewness(before, 0)', skewness(after, 0)'];
disp(array2table(round(cmp, 2), 'RowNames', cellstr(cols), 'VariableNames', ...
    {'std_before', 'std_after', 'max_before', 'max_after', 'skew_before', 'skew_after'}))

header = [" & \multicolumn{2}{c}{std} & \multicolumn{2}{c}{max} & \multicolumn{2}{c}{skew} \\"
          "\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}"
          " & before & after & before & after & before & after \\"];
printLatex(cols, cmp, header, ...
    "Standard deviation, maximum and skewness before and after cleaning (row removal, clipping and imputation).", ...
    "tab:stats_comparison");

%% 17. Save the processed files
save_cols = ["Datetime", target_col, feature_cols];

clean_out = df_before_scaling(:, save_cols);    % cleaned, original units
scaled_out = df(:, save_cols);                   % features standardised, target unchanged
clean_out.Datetime.Format = "yyyy-MM-dd HH:mm:ss";
scaled_out.Datetime.Format = "yyyy-MM-dd HH:mm:ss";

writetable(clean_out, fullfile(out_dir, "airquality_clean.csv"));
writetable(scaled_out, fullfile(out_dir, "airquality_scaled.csv"));

% Mean and std of each feature, needed to scale new data the same way
% or to convert scaled values back to original units
scaling_params = table(feature_cols', scale_mu', scale_sd', ...
    'VariableNames', {'feature', 'mean', 'std'});
writetable(scaling_params, fullfile(out_dir, "scaling_params.csv"));

fprintf("Saved %d rows to %s\n", height(clean_out), out_dir);

%% 18. EDA: each feature against the target
% Uses the cleaned but unscaled data, so the axes are in original units.
% The red line is the mean CO(GT) within 20 equal-sized bins of the feature,
% which shows the shape of the relationship without assuming a straight line.
eda = df_before_scaling;
y = eda.(target_col);
feats = feature_cols(idx);          % same order as the correlation bar chart (section 8)

fig = figure("Position", [50 50 1600 950]);
tl = tiledlayout(fig, 3, 4, "TileSpacing", "compact", "Padding", "compact");
for i = 1:numel(feats)
    ax = nexttile(tl); hold(ax, "on");
    x = eda.(feats(i));
    s = scatter(ax, x, y, 4, blue, "filled", "MarkerFaceAlpha", 0.15, "MarkerEdgeColor", "none");

    edges = unique(quantile(x, 0:0.05:1));      % unique: clipped values create repeated edges
    bin = discretize(x, edges);
    nb = numel(edges) - 1;
    xm = accumarray(bin, x, [nb 1], @mean, NaN);
    ym = accumarray(bin, y, [nb 1], @mean, NaN);
    ok = ~isnan(xm);
    t = plot(ax, xm(ok), ym(ok), "-o", "Color", red, "LineWidth", 2.2, ...
        "MarkerSize", 3, "MarkerFaceColor", red);

    title(ax, feats(i), "FontWeight", "bold", "FontSize", 12);
    if mod(i - 1, 4) == 0
        ylabel(ax, "CO(GT) (mg/m³)");
    end
    if i == 1
        ax1 = ax; hs = s; ht = t;
    end
end
lg = legend(ax1, [hs ht], ["Hourly observations", "Mean CO(GT) per bin of the feature"], ...
    "Box", "off", "FontSize", 12);
lg.Layout.Tile = 12;                % the empty 12th tile
title(tl, "Features against the target", "FontWeight", "bold", "FontSize", 17);

% Rasterised: about 90,000 points would make a vector PDF very slow to open
exportgraphics(fig, fullfile(fig_dir, "features_vs_target.pdf"), ...
    "ContentType", "image", "Resolution", 300);

fig = figure("Position", [100 100 800 550]);
ax = axes(fig);
scatter(ax, eda.("NOx(GT)"), y, 4, days(eda.Datetime - eda.Datetime(1)), ...
    "filled", "MarkerFaceAlpha", 0.3);
cb = colorbar(ax);
cb.Label.String = "Days since the first measurement";
xlabel(ax, "NOx(GT)");
ylabel(ax, "CO(GT) (mg/m³)");
title(ax, "NOx(GT) against CO(GT), coloured by date", "FontWeight", "bold");

exportgraphics(fig, fullfile(fig_dir, "NOx_bands.pdf"), ...
    "ContentType", "image", "Resolution", 300);


%% 19. EDA: CO(GT) over the measurement period
% Hourly values on a complete hourly grid, so removed rows show as gaps
% instead of being joined by a line, plus the daily mean to show the trend.
tt = table2timetable(eda(:, ["Datetime", target_col]));
hourly = retime(tt, "hourly", "fillwithmissing");
daily = retime(tt, "daily", "mean");

light_blue = blue + (1 - blue) * 0.65;

fig = figure("Position", [100 100 1400 450]);
ax = axes(fig); hold(ax, "on");
plot(ax, hourly.Datetime, hourly.(target_col), "Color", light_blue, "LineWidth", 0.5, ...
    "DisplayName", "Hourly value");
plot(ax, daily.Datetime, daily.(target_col), "Color", navy, "LineWidth", 1.8, ...
    "DisplayName", "Daily mean");
ylabel(ax, "CO(GT) (mg/m³)");
title(ax, "CO(GT) over the measurement period", "FontWeight", "bold");
legend(ax, "Box", "off", "Location", "northwest");
xlim(ax, [min(hourly.Datetime) max(hourly.Datetime)]);
exportgraphics(fig, fullfile(fig_dir, "co_over_time.pdf"), ...
    "ContentType", "image", "Resolution", 300);

%% 20. EDA: feature redundancy
% Scatter plots of the three most similar feature pairs. A narrow band of
% points means one feature almost determines the other. The correlation
% matrix is only used to pick the pairs; its values are not shown.
R = corr(eda{:, feature_cols}, "Type", "Spearman");
R(tril(true(size(R)))) = 0;                 % keep each pair once
[~, order] = sort(abs(R(:)), "descend");
[ia, ib] = ind2sub(size(R), order(1:3));
pairs = [feature_cols(ia)', feature_cols(ib)'];
% To choose the pairs yourself instead:
% pairs = ["C6H6(GT)", "PT08.S2(NMHC)"; "PT08.S1(CO)", "PT08.S5(O3)"; "NOx(GT)", "NO2(GT)"];

fig = figure("Position", [100 50 550 1300]);
tl = tiledlayout(fig, size(pairs, 1), 1, "TileSpacing", "compact", "Padding", "compact");
for p = 1:size(pairs, 1)
    ax = nexttile(tl);
    scatter(ax, eda.(pairs(p, 1)), eda.(pairs(p, 2)), 4, blue, "filled", ...
        "MarkerFaceAlpha", 0.15, "MarkerEdgeColor", "none");
    xticks(ax, []);                 % no values on the axes (also removes the grid)
    yticks(ax, []);
    ax.FontSize = 16;               % axis titles scale with this
    xlabel(ax, pairs(p, 1));
    ylabel(ax, pairs(p, 2));
    title(ax, pairs(p, 1) + " vs " + pairs(p, 2), "FontWeight", "bold", "FontSize", 16);
end
exportgraphics(fig, fullfile(fig_dir, "redundant_pairs.pdf"), ...
    "ContentType", "image", "Resolution", 300);

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

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

function Xs = knnImpute(Xs, k)
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

function printLatex(row_names, data, header, caption, label)
% Prints a booktabs LaTeX table to the Command Window.
fprintf("\n\\begin{table}[htbp]\n\\centering\n");
fprintf("\\caption{%s}\n\\label{%s}\n", caption, label);
fprintf("\\begin{tabular}{l%s}\n\\toprule\n", repmat('r', 1, size(data, 2)));
fprintf("%s\n", header);
fprintf("\\midrule\n");
for i = 1:size(data, 1)
    fprintf("%s & %s \\\\\n", row_names(i), strjoin(compose("%.2f", data(i, :)), " & "));
end
fprintf("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
end