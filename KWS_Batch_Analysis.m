function kws_ablation_batch(wavFolder)
% KWS ablation runner (no retraining).
% - Loads KWSNet from KWSBaseline1.mat
% - Runs many front-end / post-proc configs across a folder of WAVs
% - Exports a summary CSV + per-file detection CSVs
%
% Usage:
%   kws_ablation_batch(pwd)
%   kws_ablation_batch('my_wavs')

%% -------------------- GLOBAL CONFIG --------------------
NET.mat_path   = fullfile(pwd,'KWSBaseline1.mat');
NET.var_name   = 'KWSNet';          % Variable name inside MAT
CLASSES        = ["0","1"];         % Binary model: background/keyword
TARGET_IDX     = 2;                 % class "1" is the keyword

% Front-end
FE.fs          = 16e3;
FE.num_mfcc    = 13;                % MFCC base count
FE.win_ms      = 25;
FE.hop_ms      = 10;

% Output folders/files
ts         = datestr(now,'yyyymmdd_HHMMSS');
out_dir    = fullfile(wavFolder, ['det_out_' ts]);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
summary_csv = fullfile(wavFolder, ['ablation_summary_' ts '.csv']);

%% -------------------- LOAD NET -------------------------
S = load(NET.mat_path);
assert(isfield(S,NET.var_name),'Could not find %s in %s', NET.var_name, NET.mat_path);
net = S.(NET.var_name);
fprintf('[NET] Loaded "%s" (%s)\n', NET.var_name, class(net));

%% -------------------- DEFINE EXPERIMENTS ---------------
% Each struct toggles a combination of filters & post-proc.
% Feel free to add/remove/duplicate and tweak parameters.

EXPS = [
struct('name',"A_baseline_cmvn_s3_vad-65_p92", ...
    'FIL',fil(false,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"B_no_filters_cmvn_off_vad-65_p92", ...
    'FIL',fil(false,false,[ ],false,20,0.15,false,-20,0.97), ...
    'CMVN',false,'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"C_preemph_bandpass_p92", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"D_spectral_subtract_p92", ...
    'FIL',fil(true,true,[150 4000],true,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"E_smooth5_p92", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',5, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.14, 'MIN_GAP',0.20 )

struct('name',"F_lower_threshold_p88", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.88, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"G_no_VAD_p92", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.0, 'SMOOTH',3, 'VAD_BYPASS',true,  'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"H_temp_0p7_sharper_p92", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',0.7, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )

struct('name',"I_temp_1p5_softer_p92", ...
    'FIL',fil(true,true,[150 4000],false,20,0.15,true,-20,0.97), ...
    'CMVN',true, 'TEMP',1.5, 'SMOOTH',3, 'VAD_BYPASS',false, 'VAD_DB',-65, ...
    'HYST_P',0.92, 'MIN_DUR',0.12, 'MIN_GAP',0.18 )
];

%% -------------------- GATHER FILES -----------------------
if nargin<1 || isempty(wavFolder), wavFolder = pwd; end
files = dir(fullfile(wavFolder,'*.wav'));
assert(~isempty(files),'No .wav files found in %s', wavFolder);

%% -------------------- RUN -------------------------------
all = table();
for iExp = 1:numel(EXPS)
    E = EXPS(iExp);
    fprintf('\n=== EXP: %s ===\n', E.name);
    for k = 1:numel(files)
        wavPath = fullfile(files(k).folder, files(k).name);
        [row, detTbl] = run_one(net, wavPath, NET, FE, CLASSES, TARGET_IDX, E, out_dir);
        row.exp = string(E.name);
        all = [all; row]; %#ok<AGROW>

        % save per-file detection list
        det_name = sprintf('det_%s__%s.csv', erase(files(k).name,'.wav'), E.name);
        writetable(detTbl, fullfile(out_dir, det_name));
    end
end

% Write summary
writetable(all, summary_csv);
fprintf('\nSaved SUMMARY: %s\n', summary_csv);
fprintf('Saved DETECTIONS per file in: %s\n', out_dir);

end % function kws_ablation_batch

%% =========================================================
%% ===================== helpers ===========================
%% =========================================================

function F = fil(preemph, bandpass, bp_hz, specsub, ss_nf, ss_floor, norm_rms, tgt_db, preemph_a)
% Pack filtering options into a struct
F.enable_preemph  = preemph;  F.preemph_coef = preemph_a;
F.enable_bandpass = bandpass; F.bp_hz        = bp_hz;
F.enable_specsub  = specsub;  F.ss_noise_frames = ss_nf; F.ss_floor = ss_floor;
F.normalize_rms   = norm_rms; F.target_rms   = tgt_db;
end

function [row, detTbl] = run_one(net, wavPath, NET, FE, CLASSES, TARGET_IDX, E, out_dir)
[x,fs] = audioread(wavPath);
if size(x,2)>1, x = mean(x,2); end
if fs ~= FE.fs, x = resample(x, FE.fs, fs); fs = FE.fs; end
dur = numel(x)/fs;

% ----- filtering -----
x = apply_preproc(x, fs, E.FIL);
x = max(-1,min(1,x));

% ----- features (MFCC + Δ + ΔΔ) -----
win = hamming(round(FE.win_ms*1e-3*fs),'periodic');
hop = round(FE.hop_ms*1e-3*fs);
M0  = mfcc(x,fs,"NumCoeffs",FE.num_mfcc,"LogEnergy","Ignore", ...
           "Window",win,"OverlapLength",numel(win)-hop).';
D   = diff([M0(:,1) M0],1,2);
DD  = diff([D(:,1) D],1,2);
feats = [M0; D; DD];    % 39 x T
if E.CMVN
    feats = feats - mean(feats,2);
    s = std(feats,0,2); s(s<1e-6) = 1; feats = feats ./ s;
end
T = size(feats,2); tfrm = (0:T-1)*FE.hop_ms*1e-3;

% ----- inference -----
tic;
Y = predict(net,{feats});
runtime_s = toc;
if iscell(Y), Y = Y{1}; end

% normalize to C×T
ncls = numel(CLASSES);
if size(Y,2) == ncls && size(Y,1) ~= ncls, Y = Y.'; end
if size(Y,1) == ncls && size(Y,2) == 1, Y = repmat(Y,1,T); end
if size(Y,2) == ncls && size(Y,1) == 1, Y = repmat(Y.',1,T); end
if size(Y,1) == T && size(Y,2) ~= ncls, Y = Y.'; end

% temperature scaling (on logits/scores)
Y = Y ./ max(E.TEMP, eps);

% softmax
Y = Y - max(Y,[],1);
P = exp(Y) ./ max(sum(exp(Y),1), eps);

p0 = P(1,:); p1 = P(2,:);            % diagnostics
p   = P(TARGET_IDX,:);               % keyword = class "1" (row 2)

% ----- post-proc -----
% smoothing
if E.SMOOTH>1, kf = ones(1,E.SMOOTH)/E.SMOOTH; p = conv(p,kf,'same'); end

% VAD
if E.VAD_BYPASS
    vmask = true(1,T);
else
    en = movmean(x.^2, round(FE.win_ms*1e-3*fs), 'Endpoints','shrink');
    en = en(1:hop:end);
    db = 10*log10(en + 1e-12);
    vmask = (db > E.VAD_DB).';
end
p_eff = p; p_eff(~vmask) = 0;

% adaptive thresholds from percentile
on  = max(0.12, quantile(p_eff, E.HYST_P));
off = max(0.08, 0.6*on);

% events
ev = hyst_events(p_eff, tfrm, on, off, E.MIN_GAP, E.MIN_DUR);

% metrics (no GT): simple descriptive
num_ev = size(ev,1);
tot_len = sum(ev(:,2)-ev(:,1));
avg_len = (num_ev>0) * mean(ev(:,2)-ev(:,1));
pct_active = 100*sum(p_eff > on)/numel(p_eff);
pk = max(p_eff); mn = mean(p_eff);

% per-file detection table
if num_ev>0
    detTbl = table(ev(:,1), ev(:,2), 'VariableNames',{'start_s','end_s'});
else
    detTbl = table([],[],'VariableNames',{'start_s','end_s'});
end

% one-row summary for this file × experiment
row = table( ...
    string(wavPath), string(E.name), dur, runtime_s, ...
    num_ev, tot_len, avg_len, pct_active, on, off, ...
    pk, mn, E.SMOOTH, E.VAD_BYPASS, E.VAD_DB, E.CMVN, ...
    E.FIL.enable_preemph, E.FIL.enable_bandpass, E.FIL.enable_specsub, E.TEMP, ...
    'VariableNames', {'file','exp','dur_s','runtime_s','num_events','total_event_s','avg_event_s', ...
                      'pct_frames_above_on','thr_on','thr_off','peak_max','mean_prob', ...
                      'smooth_win','vad_bypass','vad_db','cmvn', ...
                      'preemph','bandpass','specsub','temp'});
end

function ev = hyst_events(p, t, onThr, offThr, minGap, minDur)
    % Robust to small length mismatches between p and t.
    L = min(numel(p), numel(t));
    p = p(1:L); t = t(1:L);

    st = false; t0 = NaN; ev = [];
    for i = 1:L
        if ~st && p(i) >= onThr
            st = true; t0 = t(i);
        end
        if st && p(i) < offThr
            t1 = t(i);
            if (t1 - t0) >= minDur
                ev = [ev; t0 t1]; %#ok<AGROW>
            end
            st = false; t0 = NaN;
        end
    end
    if st
        t1 = t(end);
        if (t1 - t0) >= minDur
            ev = [ev; t0 t1];
        end
    end
    if isempty(ev), return; end

    % Merge close events
    m = ev(1,:);
    for k = 2:size(ev,1)
        if ev(k,1) - m(end,2) < minGap
            m(end,2) = ev(k,2);
        else
            m = [m; ev(k,:)]; %#ok<AGROW>
        end
    end
    ev = m;
end


function y = apply_preproc(x, fs, F)
y = x(:);

% pre-emphasis
if F.enable_preemph
    y = filter([1 -F.preemph_coef], 1, y);
end

% band-pass FIR (zero-phase)
if F.enable_bandpass
    bp = F.bp_hz;
    if isempty(bp), bp=[150 4000]; end
    bp(1) = max(20, bp(1));
    bp(2) = min(0.45*fs, bp(2));
    ord = max(64, round(fs/100));
    d = designfilt('bandpassfir','FilterOrder',ord, ...
                   'CutoffFrequency1',bp(1), ...
                   'CutoffFrequency2',bp(2), ...
                   'SampleRate',fs);
    y = filtfilt(d.Coefficients,1,y);
end

% spectral subtraction (simple)
if F.enable_specsub
    N=1024; H=256; w=hann(N,'periodic');
    X = stft(y, "Window",w, "OverlapLength",N-H, "Centered",false);
    nf = min(F.ss_noise_frames, size(X,2));
    noise = median(abs(X(:,1:nf)),2);
    Ymag = max(abs(X) - noise, F.ss_floor*noise);
    y = istft(Ymag .* exp(1j*angle(X)), "Window",w, "OverlapLength",N-H, "Centered",false);
    if numel(y)<numel(x), y(end+1:numel(x))=0; end
    y = y(1:numel(x));
end

% loudness normalization to target RMS (dBFS)
if F.normalize_rms
    rms_now = sqrt(mean(y.^2)+eps);
    tgt_lin = 10^(F.target_rms/20);
    if rms_now>0, y = y * (tgt_lin/rms_now); end
    m = max(abs(y)); hr = 10^(-1/20); if m>hr, y = y*(hr/m); end
end
end
