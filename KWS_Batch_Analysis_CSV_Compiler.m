function out = kws_results_aggregate(baseFolder)
% KWS Results Aggregator (robust, single file)
% ------------------------------------------------------------
% - Accepts either a parent results folder (recommended) OR a det_out_* folder.
% - Merges all ablation_summary_*.csv in the parent folder.
% - Merges all det_out_*/det_*__*.csv per-file detection tables.
% - Produces tidy tables + simple pivots and a heatmap.
% - Saves a MAT and XLSX (multi-sheet) plus a few CSVs.
%
% Usage:
%   kws_results_aggregate(pwd)
%   kws_results_aggregate('C:\path\to\results')                 % parent folder
%   kws_results_aggregate('C:\path\to\results\det_out_YYYYMMDD')% det_out folder
%
% Output struct fields:
%   SummaryAll, DetAll, Counts, TotLen, Peak, Tidy, outMat, outXlsx
%
% Notes:
% - "Summary" files live in the parent folder (ablation_summary_*.csv).
% - "Detection" files live inside det_out_* subfolders.

if nargin < 1 || isempty(baseFolder), baseFolder = pwd; end

% If a det_out_* folder was passed, promote to parent where the summary csv lives
[pf, nm] = fileparts(baseFolder);
if startsWith(string(nm), "det_out_")
    baseFolder = pf;
end

ts = datestr(now,'yyyymmdd_HHMMSS');

%% 1) Load all summary CSVs (parent folder)
sumFiles = dir(fullfile(baseFolder,'ablation_summary_*.csv'));
assert(~isempty(sumFiles), 'No ablation_summary_*.csv found in %s', baseFolder);

SummaryAll = table();
for k = 1:numel(sumFiles)
    p = fullfile(sumFiles(k).folder, sumFiles(k).name);
    T = readtable(p);

    % skip empty tables safely
    if height(T) == 0
        continue
    end

    % Tolerate older schemas (ensure essential columns exist)
    if ~ismember('exp', T.Properties.VariableNames),  T.exp  = strings(height(T),1); end
    if ~ismember('file',T.Properties.VariableNames),  T.file = strings(height(T),1); end
    if ~ismember('num_events',T.Properties.VariableNames)
        T.num_events = zeros(height(T),1);
    end

    % Add source filename column (match height)
    T.src_summary = repmat(string(sumFiles(k).name), height(T), 1);

    SummaryAll = [SummaryAll; T]; %#ok<AGROW>
end

% Normalize file names to just the basename.ext for grouping
if height(SummaryAll) > 0
    if ~isstring(SummaryAll.file), SummaryAll.file = string(SummaryAll.file); end
    [~,baseNames,exts] = cellfun(@fileparts, cellstr(SummaryAll.file), 'UniformOutput', false);
    SummaryAll.file_base = string(baseNames) + string(exts);
else
    warning('No rows loaded from summary CSVs in %s', baseFolder);
end

%% 2) Load all per-file detections (from det_out_* subfolders)
detDirs = dir(fullfile(baseFolder,'det_out_*'));
DetAll = table();

for d = 1:numel(detDirs)
    dd = fullfile(detDirs(d).folder, detDirs(d).name);
    detFiles = dir(fullfile(dd,'det_*__*.csv'));   % pattern: det_<wav>__<exp>.csv
    for j = 1:numel(detFiles)
        p = fullfile(detFiles(j).folder, detFiles(j).name);
        D = readtable(p);

        % Parse file/exp from filename
        nm = erase(detFiles(j).name, '.csv');   % det_<file>__<exp>
        nm = extractAfter(nm, 'det_');
        parts = split(nm, '__');
        if numel(parts) ~= 2
            fileBase = string(detFiles(j).name);
            expName  = "unknown";
        else
            fileBase = parts(1);
            expName  = parts(2);
        end

        if height(D) == 0
            tmp = table(fileBase, expName, [], [], [], ...
                'VariableNames', {'file_base','exp','start_s','end_s','dur_s'});
        else
            if ~ismember('start_s', D.Properties.VariableNames) || ~ismember('end_s', D.Properties.VariableNames)
                warning('Detection file missing start_s/end_s: %s', detFiles(j).name);
                continue
            end
            dur = D.end_s - D.start_s;
            tmp = table( repmat(fileBase,height(D),1), repmat(expName,height(D),1), ...
                         D.start_s, D.end_s, dur, ...
                        'VariableNames', {'file_base','exp','start_s','end_s','dur_s'});
        end
        tmp.src_detdir = repmat(string(detDirs(d).name), height(tmp), 1);
        DetAll = [DetAll; tmp]; %#ok<AGROW>
    end
end

%% 3) Build pivots / tidy views  (CLEANED)
% counts (files × experiments)
GCounts = groupsummary(SummaryAll, {'file_base','exp'}, 'sum', 'num_events');

% Drop the GroupCount column so it can't leak into labels later
if ismember('GroupCount', GCounts.Properties.VariableNames)
    GCounts.GroupCount = [];
end

Counts  = unstack(GCounts, 'sum_num_events', 'exp', 'AggregationFunction', @sum);

% Ensure table exists and fill gaps
if ~isempty(Counts)
    Counts{:,2:end} = fillmissing(Counts{:,2:end}, 'constant', 0);
end

% total event duration (files × experiments)
if ismember('total_event_s', SummaryAll.Properties.VariableNames)
    GTot   = groupsummary(SummaryAll, {'file_base','exp'}, 'sum', 'total_event_s');
    if ismember('GroupCount', GTot.Properties.VariableNames); GTot.GroupCount = []; end
    TotLen = unstack(GTot, 'sum_total_event_s', 'exp', 'AggregationFunction', @sum);
    if ~isempty(TotLen)
        TotLen{:,2:end} = fillmissing(TotLen{:,2:end}, 'constant', 0);
    end
else
    TotLen = table(SummaryAll.file_base(1:0));
end

% mean posterior peak per file×exp (if column exists)
if ismember('peak_max', SummaryAll.Properties.VariableNames)
    GPeak = groupsummary(SummaryAll, {'file_base','exp'}, 'mean', 'peak_max');
    if ismember('GroupCount', GPeak.Properties.VariableNames); GPeak.GroupCount = []; end
    Peak  = unstack(GPeak, 'mean_peak_max', 'exp', 'AggregationFunction', @mean);
else
    Peak = table(SummaryAll.file_base(1:0));
end

% A tidy unique per (file,exp) row (first occurrence)
if height(SummaryAll) > 0
    [~,idx] = unique(SummaryAll(:,{'file_base','exp'}), 'rows', 'stable');
    keepCols = intersect({'file_base','exp','num_events','total_event_s','avg_event_s', ...
                          'pct_frames_above_on','thr_on','thr_off','runtime_s','src_summary'}, ...
                          SummaryAll.Properties.VariableNames, 'stable');
    Tidy = SummaryAll(idx, keepCols);
else
    Tidy = SummaryAll;
end

%% 4) Visual quick-look: heatmap of counts (clean labels, no GroupCount)
if ~isempty(Counts)
    % Raw labels from table var names (2:end are experiments)
    exps_raw = string(Counts.Properties.VariableNames(2:end));

    % Just in case: drop any stray "GroupCount" column if it slipped in
    keepMask = exps_raw ~= "GroupCount";
    exps_raw = exps_raw(keepMask);

    Z = Counts{:, [false keepMask]};     % numeric matrix (files × kept exps)

    % Make nice display labels (no underscores/commas)
    %   "A_baseline,cmvn,3,vad5,92" → "A baseline · cmvn · 3 · vad5 · 92"
    exps_disp = regexprep(exps_raw, {'(^[A-Za-z])_', '_', ','}, {'$1 ', ' ', ' · '});

    files_disp = string(Counts.file_base);
    files_disp = regexprep(files_disp, '_', ' ');  % readability

    figure('Color','w','Name','KWS Ablations: num_events heatmap');
    h = heatmap(files_disp, exps_disp, Z');  % experiments on Y, files on X
    xlabel('File'); ylabel('Experiment');
    title('# Detections per file × experiment');

    cmax = max(1, max(Z(:)));
    h.ColorLimits = [0 cmax];

    % Optional: reduce clutter on long labels
    h.FontSize = 10;
end

%% 5) Save compiled results
outMat  = fullfile(baseFolder, ['compiled_results_' ts '.mat']);
outXlsx = fullfile(baseFolder, ['compiled_results_' ts '.xlsx']);
save(outMat, 'SummaryAll','DetAll','Counts','TotLen','Peak','Tidy');

% Also save handy CSVs
csvSummaryAll = fullfile(baseFolder, ['ALL_ablation_summary_' ts '.csv']);
csvDetAll     = fullfile(baseFolder, ['ALL_detections_' ts '.csv']);
csvCounts     = fullfile(baseFolder, ['pivot_counts_' ts '.csv']);
csvTotLen     = fullfile(baseFolder, ['pivot_total_len_' ts '.csv']);
csvTidy       = fullfile(baseFolder, ['tidy_unique_' ts '.csv']);

writetable(SummaryAll, csvSummaryAll);
if ~isempty(DetAll), writetable(DetAll, csvDetAll); end
if ~isempty(Counts),  writetable(Counts,  csvCounts); end
if ~isempty(TotLen),  writetable(TotLen,  csvTotLen); end
writetable(Tidy, csvTidy);

% Excel workbook with multiple sheets (best-effort)
try
    writetable(SummaryAll, outXlsx, 'Sheet','summary_all');
    if ~isempty(DetAll), writetable(DetAll, outXlsx, 'Sheet','detections_all'); end
    if ~isempty(Counts),  writetable(Counts,  outXlsx, 'Sheet','pivot_counts'); end
    if ~isempty(TotLen),  writetable(TotLen,  outXlsx, 'Sheet','pivot_total_len'); end
    if ~isempty(Peak),    writetable(Peak,    outXlsx, 'Sheet','pivot_peak'); end
    writetable(Tidy, outXlsx, 'Sheet','tidy_unique');
catch ME
    warning('Could not write Excel workbook: %s', ME.message);
end

fprintf('Saved:\n  %s\n  %s\n', outMat, outXlsx);
fprintf('  %s\n', csvSummaryAll);
if ~isempty(DetAll), fprintf('  %s\n', csvDetAll); end
if ~isempty(Counts),  fprintf('  %s\n', csvCounts); end
if ~isempty(TotLen),  fprintf('  %s\n', csvTotLen); end
fprintf('  %s\n', csvTidy);

% Return struct (optional)
out = struct('SummaryAll',SummaryAll,'DetAll',DetAll, ...
             'Counts',Counts,'TotLen',TotLen,'Peak',Peak,'Tidy',Tidy, ...
             'outMat',outMat,'outXlsx',outXlsx, ...
             'csvSummaryAll',csvSummaryAll,'csvDetAll',csvDetAll, ...
             'csvCounts',csvCounts,'csvTotLen',csvTotLen,'csvTidy',csvTidy);
end
