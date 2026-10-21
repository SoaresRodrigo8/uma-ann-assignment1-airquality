%% Air Quality (UCI): cleaning, preprocessing and plots
% MATLAB version of main.py. All the functions are in aqlib.m.
% Requires MATLAB R2023a or newer and the Statistics and Machine Learning
% Toolbox (ksdensity, skewness, corr).

clear; clc; close all;
L = aqlib();

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
% -200 is the missing-value code; df_before keeps the raw data for section 16.
df = L.loadData(raw_path, all_num);
df = L.markMissing(df, -200);
df_before = df;

%% 2. Timestamp checks
S = L.checkTimeIndex(df.Datetime);
fprintf("Duplicates: %d\n", S.duplicates);
fprintf("Sorted: %d\n", S.sorted);
fprintf("Range: %s to %s\n", string(S.start), string(S.stop));
fprintf("Expected hours: %d, present: %d, missing: %d\n", S.expected, S.present, S.missing);
disp(S.steps)

%% 3. Missing values per column
missing_tbl = L.missingSummary(df, all_num);
fprintf("\nTotal missing values: %d\n", sum(missing_tbl.Missing));
disp(missing_tbl)

%% 4. Distribution of the target
fig = L.plotTargetDistribution(df.(target_col), target_col, "CO concentration (mg/m³)");
exportgraphics(fig, fullfile(fig_dir, "co_distribution.pdf"), "ContentType", "vector");

%% 5. Descriptive statistics (Table II)
D = L.describeVars(df, cols);
disp(round(D, 2))

fprintf("%s", L.latexTable(cols, D{:, :}, " & mean & std & min & median & max & skew \\", ...
    "Descriptive statistics of the target and features.", "tab:descriptive_stats"));

%% 6. Histograms of the target and every feature
fig = L.plotHistograms(df, cols);
exportgraphics(fig, fullfile(fig_dir, "feature_distributions.pdf"), "ContentType", "vector");

%% 7. Features grouped by distribution shape (standardised KDEs)
groups = {"Strongly right-skewed",   ["NOx(GT)", "C6H6(GT)", "PT08.S3(NOx)"];
          "Moderately right-skewed", ["PT08.S1(CO)", "PT08.S2(NMHC)", "PT08.S5(O3)", "NO2(GT)"];
          "Approximately symmetric", ["PT08.S4(NO2)", "T", "RH", "AH"]};
fig = L.plotSkewGroups(df, groups);
exportgraphics(fig, fullfile(fig_dir, "feature_distributions_grouped.pdf"), "ContentType", "vector");

%% 8. Spearman correlation with the target
rho = L.targetCorrelation(df, feature_cols, target_col);
[~, idx] = sort(abs(rho), "descend");      % order reused in section 18

fig = L.plotTargetCorrelation(rho, feature_cols, target_col);
exportgraphics(fig, fullfile(fig_dir, "corr_with_target.pdf"), "ContentType", "vector");

%% 9. Outliers (IQR rule) on the data before cleaning
disp(L.outlierSummary(df, all_num, 1.5))

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
[lo, hi] = L.iqrBounds(X, 1.5);
[X_clip, n_clip] = L.clipToBounds(X, lo, hi);

disp(table(n_clip, round(n_clip ./ sum(~isnan(X))' * 100, 2), ...
    'RowNames', cellstr(feature_cols), 'VariableNames', {'Clipped', 'Pct'}))
fprintf("%.2f%% of all feature values clipped\n", sum(n_clip) / sum(~isnan(X), "all") * 100);

df{:, feature_cols} = X_clip;

%% 12. Gap analysis on a complete hourly grid
% Removed rows count as part of a gap, so gap lengths are real hours.
[hidx, in_df] = L.hourlyIndex(df.Datetime);

fprintf("\nMissing before imputation:\n");
disp(array2table(sum(isnan(df{:, targets})), 'VariableNames', cellstr(targets)))

for c = targets
    G = L.gapStats(df.(c), hidx, in_df, max_gap);
    fprintf("%s: %d relevant gaps, longest = %d h, median = %g h | " + ...
            "values in short gaps: %d, in long gaps: %d\n", c, G.nGaps, ...
            G.longest, G.median, G.inShort, G.inLong);
end

%% 13. Imputation
% Step 1: linear interpolation over time for gaps of up to max_gap hours
for c = targets
    df.(c) = L.interpolateShortGaps(df.(c), hidx, max_gap);
end
fprintf("\nLeft for KNN:\n");
disp(array2table(sum(isnan(df{:, targets})), 'VariableNames', cellstr(targets)))

% Step 2: KNN (k = 5) for the remaining long-gap values, target not used
X = L.knnImpute(df{:, feature_cols}, 5);
df{:, feature_cols} = X;
fprintf("Missing after imputation: %d\n", sum(isnan(X), "all"));

%% 14. Standardisation
df_before_scaling = df;

P = L.zscoreFit(df, feature_cols);        % population std, as in sklearn's StandardScaler
df = L.zscoreApply(df, P);

Z = df{:, feature_cols};
disp(table(round(mean(Z), 3)', round(std(Z), 3)', ...
    'RowNames', cellstr(feature_cols), 'VariableNames', {'mean', 'std'}))

%% 15. Before/after standardisation plot
fig = L.plotScalingBeforeAfter(df_before_scaling{:, feature_cols}, df{:, feature_cols}, feature_cols);
exportgraphics(fig, fullfile(fig_dir, "scaling_before_after.pdf"), "ContentType", "vector");

%% 16. Before/after cleaning comparison table
cmp = L.compareStats(df_before{:, cols}, df_before_scaling{:, cols}, cols);
disp(round(cmp, 2))

header = [" & \multicolumn{2}{c}{std} & \multicolumn{2}{c}{max} & \multicolumn{2}{c}{skew} \\"
          "\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}"
          " & before & after & before & after & before & after \\"];
fprintf("%s", L.latexTable(cols, cmp{:, :}, header, ...
    "Standard deviation, maximum and skewness before and after cleaning (row removal, clipping and imputation).", ...
    "tab:stats_comparison"));

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
scaling_params = table(feature_cols', P.Mean, P.Std, ...
    'VariableNames', {'feature', 'mean', 'std'});
writetable(scaling_params, fullfile(out_dir, "scaling_params.csv"));

fprintf("Saved %d rows to %s\n", height(clean_out), out_dir);

%% 18. EDA: each feature against the target
% Uses the cleaned but unscaled data, so the axes are in original units.
% The red line is the mean CO(GT) within 20 equal-sized bins of the feature,
% which shows the shape of the relationship without assuming a straight line.
eda = df_before_scaling;
feats = feature_cols(idx);          % same order as the correlation bar chart (section 8)

fig = L.plotFeaturesVsTarget(eda, feats, target_col);
% Rasterised: about 90,000 points would make a vector PDF very slow to open
exportgraphics(fig, fullfile(fig_dir, "features_vs_target.pdf"), ...
    "ContentType", "image", "Resolution", 300);

fig = L.plotScatterByTime(eda, "NOx(GT)", target_col);
exportgraphics(fig, fullfile(fig_dir, "NOx_bands.pdf"), ...
    "ContentType", "image", "Resolution", 300);

%% 19. EDA: CO(GT) over the measurement period
% Hourly values on a complete hourly grid, so removed rows show as gaps
% instead of being joined by a line, plus the daily mean to show the trend.
fig = L.plotTargetOverTime(eda, target_col);
exportgraphics(fig, fullfile(fig_dir, "co_over_time.pdf"), ...
    "ContentType", "image", "Resolution", 300);

%% 20. EDA: feature redundancy
% Scatter plots of the three most similar feature pairs. A narrow band of
% points means one feature almost determines the other. The correlation
% matrix is only used to pick the pairs; its values are not shown.
pairs = L.topCorrelatedPairs(eda, feature_cols, 3);
% To choose the pairs yourself instead:
% pairs = ["C6H6(GT)", "PT08.S2(NMHC)"; "PT08.S1(CO)", "PT08.S5(O3)"; "NOx(GT)", "NO2(GT)"];

fig = L.plotPairs(eda, pairs);
exportgraphics(fig, fullfile(fig_dir, "redundant_pairs.pdf"), ...
    "ContentType", "image", "Resolution", 300);
