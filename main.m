% MAIN  Feature analysis and engineering of the UCI Air Quality dataset.
%   Target: CO(GT). Each section is one step of the assignment (Run Section).

L = aqlib();

rawFile   = fullfile('data', 'Raw', 'AirQualityUCI.csv');
figDir    = 'Figures';
target    = 'CO(GT)';
missCode  = -200;
maxGap    = 3;      % longest gap (hours) filled by interpolation

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
