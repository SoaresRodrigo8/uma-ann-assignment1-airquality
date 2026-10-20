% MAIN  Feature analysis and engineering of the UCI Air Quality dataset.
%   Target: CO(GT). Each section is one step of the assignment (Run Section).

L = aqlib();

rawFile    = fullfile('data', 'Raw', 'AirQualityUCI.csv');
figDir     = 'Figures';
target     = 'CO(GT)';
missCode   = -200;
maxGap     = 3;     % longest gap (hours) filled by interpolation
maxMissPct = 50;    % variables with more missing values than this are removed

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
