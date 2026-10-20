% MAIN  Feature analysis and engineering of the UCI Air Quality dataset.
%   Target: CO(GT). Each section is one step of the assignment (Run Section).

L = aqlib();

rawFile    = fullfile('data', 'Raw', 'AirQualityUCI.csv');
figDir     = 'Figures';
target     = 'CO(GT)';
missCode   = -200;
maxGap     = 3;     % longest gap (hours) filled by interpolation
maxMissPct = 50;    % variables with more missing values than this are removed
skewThr    = 1;     % inputs with |skewness| above this are log-transformed
procDir    = fullfile('data', 'Processed');
shiftDate  = datetime(2004, 9, 1);   % NOx(GT) level shift found in 2.5

if ~exist(figDir, 'dir'), mkdir(figDir); end

%% 1.1 Load the dataset
% Missing values are encoded as -200 in the raw file; they become NaN here
% so that they are not counted as real measurements in the statistics.
Traw = L.loadData(rawFile);
T    = L.markMissing(Traw, missCode);

numVars = T.Properties.VariableNames(2:end);
fprintf('Observations: %d\n', height(T));
fprintf('Variables:    %d (DateTime + %d numeric)\n', width(T), numel(numVars));
fprintf('Missing values (total): %d of %d cells (%.2f%%)\n', ...
    sum(isnan(T{:, numVars}), 'all'), numel(T{:, numVars}), ...
    100 * mean(isnan(T{:, numVars}), 'all'));

%% 1.2 Time index
% The data are hourly; check for gaps, duplicates and ordering, since this
% matters later for time-based imputation and lag features.
S = L.checkTimeIndex(T.DateTime);
fprintf('Range: %s to %s\n', S.start, S.stop);
fprintf('Expected hours: %d | missing hours: %d | duplicates: %d | sorted: %d\n', ...
    S.expected, S.missing, S.duplicates, S.sorted);

%% 1.3 Missing values per variable
M = L.missingSummary(T);
disp(M);

% When was each instrument working? (% of valid hours per month)
A = L.monthlyAvailability(T);
disp(A);

fig = L.plotMissingMap(T);
exportgraphics(fig, fullfile(figDir, 'missing_map.pdf'));

%% 1.4 Distribution of the target
fig = L.plotTargetDistribution(T.(target), target);
exportgraphics(fig, fullfile(figDir, 'target_distribution.pdf'));

%% 1.5 Descriptive statistics (summary table)
D = L.describeVars(T);
disp(D);

fig = L.plotHistograms(T, numVars);
exportgraphics(fig, fullfile(figDir, 'feature_histograms.pdf'));

%% 1.6 Correlation with the target
% Spearman is reported alongside Pearson because several variables are
% skewed and the relation with CO may be monotonic but not linear.
C = L.targetCorrelation(T, target);
disp(C);

fig = L.plotTargetCorrelation(C, target);
exportgraphics(fig, fullfile(figDir, 'target_correlation.pdf'));

%% 2.1 Missing-value pattern
% Before choosing an imputation method: do the gaps coincide across
% variables, and how long are they? Short gaps can be interpolated in time;
% long gaps (days) cannot be filled credibly without erasing the daily cycle.
sensorVars = {'PT08.S1(CO)', 'PT08.S2(NMHC)', 'PT08.S3(NOx)', 'PT08.S4(NO2)', ...
              'PT08.S5(O3)', 'C6H6(GT)', 'T', 'RH', 'AH'};
sensorMiss = isnan(T{:, sensorVars});
fprintf('Sensor variables missing on the same rows: %d\n', ...
    isequal(all(sensorMiss, 2), any(sensorMiss, 2)));
fprintf('Rows without %s that also miss NOx(GT): %d of %d\n', target, ...
    sum(isnan(T.(target)) & isnan(T.('NOx(GT)'))), sum(isnan(T.(target))));

G = L.gapSummary(T, maxGap);
disp(G);

% Sensitivity of the maxGap choice: % of missing cells that would be
% interpolated for different limits (inputs only matter; CO is the target).
limits = [1 3 6 12 24];
S = table();
for k = 1:numel(limits)
    Gk = L.gapSummary(T, limits(k));
    S.(sprintf('maxGap_%dh', limits(k))) = Gk.ShortCellsPct;
end
S.Properties.RowNames = Gk.Properties.RowNames;
disp(S);

% Where do the 1-hour gaps of the analysers fall during the day?
% A single dominant hour points to a systematic cause (daily calibration).
H = L.gapHours(T, {target, 'NOx(GT)', 'NO2(GT)'}, 1);
disp(H(any(H{:, 2:end} > 0, 2), :));

%% 2.2 Remove variables with too many missing values
% A variable that is mostly missing cannot be imputed without inventing most
% of it. NMHC(GT) is the only one above the limit (90%); the next worst input
% has ~18%, so any limit between ~20% and 90% gives the same result.
[Tc, dropped] = L.dropSparseVars(T, maxMissPct, {target});
fprintf('Removed (> %d%% missing): %s\n', maxMissPct, strjoin(dropped, ', '));
fprintf('Variables before: %d | after: %d\n', width(T) - 1, width(Tc) - 1);

inputPct = M.Percent(~ismember(M.Properties.RowNames, [dropped, {target}]));
fprintf('Highest %% missing among the remaining inputs: %.2f%%\n', max(inputPct));

%% 2.3 Interpolate short gaps in the input features
% Linear interpolation in time, only for gaps of at most maxGap hours
% (mostly the daily 1-hour calibration gaps of NOx/NO2). Long gaps stay NaN.
% The target is never interpolated: missing CO values are not invented.
inputVars = setdiff(Tc.Properties.VariableNames(2:end), {target}, 'stable');
[Ti, filled] = L.interpolateGaps(Tc, inputVars, maxGap);

% Before/after: missing values per variable
Mb = L.missingSummary(Tc);
Ma = L.missingSummary(Ti);
Ma = Ma(Mb.Properties.RowNames, :);
B = table(Mb.Missing, Ma.Missing, Mb.Missing - Ma.Missing, ...
    'RowNames', Mb.Properties.RowNames, ...
    'VariableNames', {'MissingBefore', 'MissingAfter', 'Interpolated'});
disp(B);
fprintf('Cells interpolated: %d | target unchanged: %d\n', ...
    sum(filled, 'all'), isequaln(Tc.(target), Ti.(target)));

% Before/after: the statistics of the inputs should barely change, since
% only a few isolated hours were filled
Db = L.describeVars(Tc(:, inputVars));
Da = L.describeVars(Ti(:, inputVars));
disp(table(Db.Mean, Da.Mean, Db.Std, Da.Std, 'RowNames', inputVars, ...
    'VariableNames', {'MeanBefore', 'MeanAfter', 'StdBefore', 'StdAfter'}));

% Why does the NOx/NO2 mean drop? The interpolated hours are mostly at
% 03:00 (low traffic), so before interpolation night hours were missing
% more often than day hours and the mean was biased upwards.
for v = {'NOx(GT)', 'NO2(GT)'}
    k = find(strcmp(inputVars, v{1}));
    fprintf('%s: mean of interpolated values %.1f vs mean of measured values %.1f\n', ...
        v{1}, mean(Ti.(v{1})(filled(:, k))), mean(Tc.(v{1}), 'omitnan'));
end

% Before/after plot for one week (Mon 8 - Sun 14 Nov 2004), chosen because it
% contains both cases: daily 1-hour gaps that are filled and a long gap that
% is deliberately left missing.
win  = Tc.DateTime >= datetime(2004, 11, 8) & Tc.DateTime < datetime(2004, 11, 15);
show = {'NOx(GT)', 'NO2(GT)'};
[~, cols] = ismember(show, inputVars);   % column of each shown variable in 'filled'
fig = L.plotInterpolation(Tc.DateTime(win), Tc{win, show}, Ti{win, show}, ...
    filled(win, cols), show);
exportgraphics(fig, fullfile(figDir, 'interpolation_example.pdf'));

% Whole-period view: valid values per hour of the day, before and after.
% Shows what the interpolation corrects (the systematic 03:00 gap) using all
% the data, not only the example week above.
Vb = L.validByHour(Tc, show);
Va = L.validByHour(Ti, show);
disp(table(Vb.Hour, Vb.('NOx(GT)'), Va.('NOx(GT)'), Vb.('NO2(GT)'), Va.('NO2(GT)'), ...
    'VariableNames', {'Hour', 'NOx_before', 'NOx_after', 'NO2_before', 'NO2_after'}));

fig = L.plotHourlyCoverage(Vb, Va, show);
exportgraphics(fig, fullfile(figDir, 'interpolation_hourly_coverage.pdf'));

%% 2.4 Keep only complete rows (target and all inputs present)
% Rows without the target cannot be used to relate inputs to CO, and the
% target is never invented. Rows with a long gap in an input are also
% excluded instead of filling days of data. Rows are marked, not deleted:
% Ti stays hourly (needed for lag features); analyses use Tk = Ti(keep, :).
keep = L.completeRows(Ti, [{target}, inputVars]);
Tk   = Ti(keep, :);

noTarget = isnan(Ti.(target));
noInput  = any(isnan(Ti{:, inputVars}), 2);
fprintf('Rows before: %d | kept: %d (%.1f%%) | excluded: %d\n', ...
    height(Ti), sum(keep), 100 * mean(keep), sum(~keep));
fprintf('  only target missing: %d | only inputs missing: %d | both: %d\n', ...
    sum(noTarget & ~noInput), sum(~noTarget & noInput), sum(noTarget & noInput));

% Is the kept subset still representative of the target?
% Compare CO over all hours where it was measured vs the kept rows only.
Dall  = L.describeVars(Ti(~noTarget, {target}));
Dkeep = L.describeVars(Tk(:, {target}));
Dall.Properties.RowNames  = {'CO, all measured hours'};
Dkeep.Properties.RowNames = {'CO, kept rows'};
disp([Dall; Dkeep]);

% Rows kept per month: does the exclusion remove whole periods?
[g, month] = findgroups(dateshift(Ti.DateTime, 'start', 'month'));
month.Format = 'MMM yyyy';
disp(table(month, splitapply(@sum, keep, g), round(100 * splitapply(@mean, keep, g), 1), ...
    'VariableNames', {'Month', 'RowsKept', 'PctKept'}));

%% 2.5 Outliers
% Identification with Tukey's rule (1.5 IQR) on the complete rows. An
% outlier is only removed if it is implausible; pollution peaks are real
% events and are exactly the hours a model of CO must be able to explain.
allVars = [{target}, inputVars];
[So, O] = L.outlierSummary(Tk, allVars, 1.5);
disp(So);

% 1) Physical plausibility: no negative concentrations, RH within 0-100 %
fprintf('Min of all variables except T: %.3f | RH range: %.1f to %.1f %%\n', ...
    min(Tk{:, setdiff(allVars, {'T'})}, [], 'all'), min(Tk.RH), max(Tk.RH));

% 2) Do the outliers of one pollutant coincide with outliers of others?
% (meteorological variables are left out: they are not pollutants)
polVars = setdiff(allVars, {'T', 'RH', 'AH'}, 'stable');
[~, iPol] = ismember(polVars, allVars);
Co = L.outlierCooccurrence(O(:, iPol), polVars);

% Median CO during the outlier hours of each variable (overall median: 1.8).
% Unlike the co-occurrence, this does not depend on the sign of the relation:
% PT08.S3(NOx) responds inversely, so its high values mean LOW pollution.
Co.MedianCO = arrayfun(@(j) median(Tk.(target)(O(:, j))), iPol)';
disp(Co);
fprintf('Median CO over all kept rows: %.2f\n', median(Tk.(target)));

% 3) When do the outliers happen? (number per month, CO and NOx)
[g, month] = findgroups(dateshift(Tk.DateTime, 'start', 'month'));
month.Format = 'MMM yyyy';
iNOx = find(strcmp(allVars, 'NOx(GT)'));
disp(table(month, splitapply(@sum, keep(keep), g), splitapply(@sum, O(:, 1), g), ...
    splitapply(@sum, O(:, iNOx), g), ...
    'VariableNames', {'Month', 'Rows', 'CO_outliers', 'NOx_outliers'}));

% NOx has no outliers from Mar to Aug 2004 but many in Mar 2005 (same
% month, other year): check whether the NOx level itself changed over time,
% using monthly medians of NOx, CO and their ratio.
medNOx = splitapply(@median, Tk.('NOx(GT)'), g);
medCO  = splitapply(@median, Tk.(target), g);
disp(table(month, medCO, medNOx, round(medNOx ./ medCO), ...
    'VariableNames', {'Month', 'MedianCO', 'MedianNOx', 'NOx_per_CO'}));

fig = L.plotBoxplots(Tk, allVars);
exportgraphics(fig, fullfile(figDir, 'outliers_boxplots.pdf'));

%% 2.6 Transformation and scaling of the input features
% (a) Log transform of the strongly skewed inputs (|skewness| > skewThr):
%     compresses the long right tail without removing any value (the
%     outliers are kept, see 2.5). The target stays in mg/m^3: there is no
%     model yet, and transforming the target is a modelling decision.
Ds = L.describeVars(Tk(:, inputVars));
logVars = inputVars(abs(Ds.Skewness') > skewThr);
[Tl, logNames] = L.logTransform(Tk, logVars);

Dl = L.describeVars(Tl(:, logNames));
disp(table(Ds{logVars, 'Skewness'}, Dl.Skewness, Ds{logVars, 'Kurtosis'}, Dl.Kurtosis, ...
    'RowNames', logVars, ...
    'VariableNames', {'SkewBefore', 'SkewAfter', 'KurtBefore', 'KurtAfter'}));

fig = L.plotTransformHist(Tk{:, logVars}, Tl{:, logNames}, logVars, logNames);
exportgraphics(fig, fullfile(figDir, 'log_transform.pdf'));

% (b) Z-score of all inputs: the inputs have very different scales (AH ~1,
%     PT08.S4 ~1500), which would let large-valued variables dominate
%     distance- and gradient-based methods. Parameters are computed once
%     (zscoreFit) and applied separately (zscoreApply), so that a future model
%     can fit them on the training set only.
modelVars = setdiff(Tl.Properties.VariableNames(2:end), {target}, 'stable');
P  = L.zscoreFit(Tl, modelVars);
Tz = L.zscoreApply(Tl, P);

Dz = L.describeVars(Tz(:, modelVars));
disp(table(round(P.Mean, 3), round(P.Std, 3), Dz.Mean, Dz.Std, Dz.Min, Dz.Max, ...
    'RowNames', modelVars, ...
    'VariableNames', {'MeanBefore', 'StdBefore', 'MeanAfter', 'StdAfter', 'MinAfter', 'MaxAfter'}));

%% 2.7 Save the processed data
% Tk: complete rows, original units (for EDA and interpretation)
% Tz: complete rows, log-transformed and standardised inputs, CO in mg/m^3
% P:  scaling parameters, needed to apply the same scaling to new data
if ~exist(procDir, 'dir'), mkdir(procDir); end
writetable(Tk, fullfile(procDir, 'airquality_clean.csv'));
writetable(Tz, fullfile(procDir, 'airquality_scaled.csv'));
writetable(P,  fullfile(procDir, 'scaling_params.csv'), 'WriteRowNames', true);
fprintf('Saved %d rows to %s\n', height(Tk), procDir);

%% 3.1 Each input against the target
% Original units (Tk), so that the relations can be read physically. The
% points are split in the two periods found in 2.5 (before/after the NOx
% level shift) to see which relations with CO stay stable over time. The
% line is the median CO in 15 quantile bins of the input: it shows the shape
% of the relation (linear, curved, flat) without assuming any model.
late = Tk.DateTime >= shiftDate;
periodNames = {'Mar-Aug 2004', 'Sep 2004-Apr 2005'};
fig = L.plotFeatureVsTarget(Tk, inputVars, target, late, periodNames);
exportgraphics(fig, fullfile(figDir, 'features_vs_target.pdf'), ...
    'ContentType', 'image', 'Resolution', 200);

% The same relations as numbers: Spearman correlation with CO per period
Ce = L.targetCorrelation(Tk(~late, [{'DateTime', target}, inputVars]), target);
Cl = L.targetCorrelation(Tk(late,  [{'DateTime', target}, inputVars]), target);
Cl = Cl(Ce.Properties.RowNames, :);
disp(table(Ce.Spearman, Cl.Spearman, Cl.Spearman - Ce.Spearman, Ce.Pearson, Cl.Pearson, ...
    'RowNames', Ce.Properties.RowNames, ...
    'VariableNames', {'SpearmanEarly', 'SpearmanLate', 'Change', 'PearsonEarly', 'PearsonLate'}));
fprintf('Rows: early %d | late %d\n', sum(~late), sum(late));

% The two periods are also different seasons, so a change between them may
% be seasonal. Control for both season and CO level: compare spring 2004 with
% spring 2005 (Mar-Apr) at the same CO level (1.5-2.5 mg/m^3). An input that
% still changes points to a change in the instrument (drift/recalibration).
mo  = Tk.DateTime.Month;   % (property, because 'month' is a variable above)
yr  = Tk.DateTime.Year;
mid = Tk.(target) >= 1.5 & Tk.(target) <= 2.5;
s04 = ismember(mo, [3 4]) & yr == 2004 & mid;
s05 = ismember(mo, [3 4]) & yr == 2005 & mid;
a = varfun(@median, Tk(s04, inputVars));
b = varfun(@median, Tk(s05, inputVars));
disp(table(a{1, :}', b{1, :}', round(100 * (b{1, :} - a{1, :}) ./ a{1, :})', ...
    'RowNames', inputVars, 'VariableNames', {'Spring2004', 'Spring2005', 'PctChange'}));
fprintf('Rows with CO in [1.5, 2.5]: spring 2004 %d | spring 2005 %d\n', sum(s04), sum(s05));
