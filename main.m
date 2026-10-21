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

%% 21. Feature selection: relevance of each feature for the target
% Cleaned data in original units (eda, section 18). Three measures, because
% each one sees a different kind of relationship: Pearson only linear,
% Spearman any monotonic one, mutual information any dependence at all.
% MI uses 10 equal-frequency bins (deciles): about 73 rows per cell of the
% 10 x 10 joint histogram, enough for a stable estimate.
n_bins = 10;
rel = L.featureRelevance(eda, feature_cols, target_col, n_bins);
disp(round(rel, 3))

fig = L.plotRelevance(rel, target_col);
exportgraphics(fig, fullfile(fig_dir, "feature_relevance.pdf"), "ContentType", "vector");

%% 22. Feature selection: redundancy between features (correlation clusters)
% Features are grouped by hierarchical clustering with distance 1 - |rho|
% (Spearman). Complete linkage: a group only forms if EVERY pair inside it
% has |rho| >= min_corr, so each group is a set of mutually redundant features.
min_corr = 0.9;
cc = L.correlationClusters(eda, feature_cols, min_corr);

fig = L.plotCorrClusters(cc, feature_cols, min_corr);
exportgraphics(fig, fullfile(fig_dir, "correlation_clusters.pdf"), "ContentType", "vector");

ord = feature_cols(cc.order);               % similar features next to each other
disp(array2table(round(cc.R(cc.order, cc.order), 2), ...
    'RowNames', cellstr(ord), 'VariableNames', cellstr(ord)))

% Sensitivity of the groups to the threshold: which groups form as the
% required |rho| is lowered (min_corr = 0.9 is the one used in the decision)
for t = [0.95 0.90 0.85 0.80]
    ct = L.correlationClusters(eda, feature_cols, t);
    fprintf("\nGroups with |rho| >= %.2f between every pair:\n", t);
    for g = unique(ct.groups)'
        in_g = ct.groups == g;
        if sum(in_g) > 1
            Rg = abs(ct.R(in_g, in_g));
            fprintf("  %s (lowest |rho| inside: %.3f)\n", ...
                strjoin(feature_cols(in_g), ", "), min(Rg(:)));
        end
    end
end

%% 23. Feature selection: principal component analysis (PCA)
% Steps of the lecture: standardise, covariance matrix, eigendecomposition,
% components sorted by eigenvalue. Used as a diagnostic of redundancy: if a
% few components explain most of the variance, the 11 features carry much
% less than 11 independent pieces of information. The features themselves
% are kept (a subset of the original features is selected, not replaced).
pc = L.pcaAnalysis(eda, feature_cols);
n_pc = numel(pc.eigval);
disp(table((1:n_pc)', round(pc.eigval, 3), round(pc.explained, 1), round(pc.cumulative, 1), ...
    'VariableNames', {'PC', 'Eigenvalue', 'ExplainedPct', 'CumulativePct'}))

% Number of components needed for the usual explained-variance thresholds
for t = [80 90 95]
    fprintf("Components for %d%% of the variance: %d\n", t, find(pc.cumulative >= t, 1));
end

fig = L.plotPCA(pc, feature_cols, 90);
exportgraphics(fig, fullfile(fig_dir, "pca.pdf"), "ContentType", "vector");

% Loadings: which features each component is made of. The first three
% carry most of the variance; the last one (eigenvalue close to 0) is a
% combination of features that is almost constant, i.e. a redundancy.
disp(array2table(round(pc.coeff(:, [1:3 n_pc]), 2), 'RowNames', cellstr(feature_cols), ...
    'VariableNames', {'PC1', 'PC2', 'PC3', sprintf('PC%d', n_pc)}))

%% 24. Feature selection: RReliefF (relevance in the context of all features)
% Each hour is compared with its k most similar hours in ALL the features.
% A feature gains weight when its small differences between similar hours
% go together with differences in CO. Unlike section 21, which only
% look at pairs of variables, this can reveal features that only help in
% combination with others (e.g. temperature correcting the sensors).
% k = 10 is the usual choice for RReliefF; a column of random noise is added
% as a reference, since the weights have no absolute scale.
k_relief = 10;
relief = L.reliefWeights(eda, feature_cols, target_col, k_relief, 1);
disp(round(relief, 4))

fig = L.plotRelief(relief, k_relief);
exportgraphics(fig, fullfile(fig_dir, "relieff.pdf"), "ContentType", "vector");

% Sensitivity to the number of neighbours and to the noise column (seed)
ks = [5 10 20 50];
sens = table('Size', [numel(feature_cols) + 1, numel(ks)], 'VariableTypes', repmat("double", 1, numel(ks)), ...
    'VariableNames', compose("k%d", ks), 'RowNames', relief.Properties.RowNames);
for i = 1:numel(ks)
    w = L.reliefWeights(eda, feature_cols, target_col, ks(i), 1);
    sens{:, i} = w{sens.Properties.RowNames, "Weight"};
end
disp(round(sens, 4))
for seed = 2:5
    w = L.reliefWeights(eda, feature_cols, target_col, k_relief, seed);
    fprintf("Seed %d: noise weight %.4f\n", seed, w{"Random noise", "Weight"});
end

% Relief does not penalise redundancy: two features carrying the same
% information share the weight. Removing C6H6(GT) shows how much of the
% weight of PT08.S2(NMHC) it was taking.
no_c6h6 = L.reliefWeights(eda, feature_cols(feature_cols ~= "C6H6(GT)"), target_col, k_relief, 1);
fprintf("PT08.S2(NMHC) weight: %.4f with C6H6(GT), %.4f without it\n", ...
    relief{"PT08.S2(NMHC)", "Weight"}, no_c6h6{"PT08.S2(NMHC)", "Weight"});

%% 25. Feature selection: stability of the relationships over time
% (a) Spearman correlation with CO before and after September 2004, when the
%     level of NOx(GT) changes (NOx_bands.pdf, section 18).
split_date = datetime(2004, 9, 1);
by_period = L.periodCorrelation(eda, feature_cols, target_col, split_date);
disp(round(by_period, 3))
fprintf("Rows: before %d | after %d\n", sum(eda.Datetime < split_date), sum(eda.Datetime >= split_date));

% (b) The two periods are also different seasons, so a change between them
%     may just be seasonal. Control for season and pollution level: compare
%     March-April 2004 with March-April 2005, only in the hours with CO
%     between 1.5 and 2.5 mg/m3 (around the median, 1.8). T, RH and AH should
%     barely change (same season): they show how well the comparison is
%     controlled. A feature that still changes has changed its relationship
%     with CO.
max_shift = 20;                     % largest |% change| considered stable
shift = L.sameSeasonShift(eda, feature_cols, target_col, [3 4], [2004 2005], [1.5 2.5]);
disp(round(shift, 1))
fprintf("Hours used: %d in 2004, %d in 2005\n", shift.Properties.UserData);

% Which features are unstable for other thresholds
for t = [10 20 30 50]
    fprintf("|change| > %d%%: %s\n", t, strjoin(feature_cols(abs(shift.PctChange') > t), ", "));
end

fig = L.plotStability(shift, max_shift);
exportgraphics(fig, fullfile(fig_dir, "stability.pdf"), "ContentType", "vector");

%% 26. Feature selection: decision and final subset
% One row per feature with every criterion. The rules (selectFeatures) were
% fixed before looking at the result:
%   1. unstable over time (|change| > max_shift)                 -> remove
%   2. RReliefF weight not above the noise for every k tested    -> remove
%   3. |rho| >= min_corr with a more stable feature (redundant)  -> remove
% Low relevance on its own is NOT a reason to remove a feature: section 24
% showed that T, RH and AH only help together with the sensors.
Rabs = abs(cc.R);
Rabs(1:numel(feature_cols) + 1:end) = 0;            % ignore each feature with itself
[max_corr, j_sim] = max(Rabs, [], 2);
above_noise = all(sens{cellstr(feature_cols), :} > sens{'Random noise', :}, 2);

crit = table(rel{cellstr(feature_cols), "Spearman"}, rel{cellstr(feature_cols), "MI"}, ...
    max_corr, feature_cols(j_sim)', relief{cellstr(feature_cols), "Weight"}, above_noise, ...
    shift.PctChange, 'RowNames', cellstr(feature_cols), 'VariableNames', ...
    {'Spearman', 'MI', 'MaxCorr', 'MostSimilar', 'Relief', 'AboveNoise', 'Shift'});
crit = L.selectFeatures(crit, min_corr, max_shift);

show = crit;                                        % rounded copy for display
show.Spearman = round(show.Spearman, 3);
show.MI = round(show.MI, 3);
show.MaxCorr = round(show.MaxCorr, 3);
show.Relief = round(1000 * show.Relief, 2);         % x 1e-3
show.Shift = round(show.Shift, 1);
disp(show)

selected = feature_cols(crit.Keep);
fprintf("Selected features (%d of %d): %s\n", numel(selected), numel(feature_cols), strjoin(selected, ", "));
fprintf("Removed: %s\n", strjoin(feature_cols(~crit.Keep), ", "));

% Save the selected subset (same rows as section 17)
sel_cols = ["Datetime", target_col, selected];
writetable(clean_out(:, sel_cols), fullfile(out_dir, "airquality_selected.csv"));
writetable(scaled_out(:, sel_cols), fullfile(out_dir, "airquality_selected_scaled.csv"));
