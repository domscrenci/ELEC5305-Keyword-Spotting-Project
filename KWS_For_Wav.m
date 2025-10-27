% Keyword Spotting: simple, robust, and tuned for your binary model.
% Files:  KWSBaseline1.mat  (var: KWSNet)
%         keywordTestSignal.wav
% Target: class "1" (keyword present)

clear; clc; close all;

%% ---------------- CONFIG ----------------
CFG.wav_path    = fullfile(pwd,'yes__01.wav');
CFG.net_path    = fullfile(pwd,'KWSBaseline1.mat');
CFG.net_var     = 'KWSNet';           % variable name in MAT

CFG.class_list  = ["0","1"];          % binary: 0 = background, 1 = keyword
CFG.target_keyword = "1";             % force keyword = class "1"

% Front-end (matches 39-dim input: 13 MFCC + Δ + ΔΔ)
CFG.fs_target   = 16e3;
CFG.num_mfcc    = 13;
CFG.win_ms      = 25;
CFG.hop_ms      = 10;
CFG.add_deltas  = true;

% Feature normalization
CFG.use_cmvn    = true;               % cepstral mean/var norm (utterance)

% Post-processing
PP.smooth_win   = 1;                  % low smoothing to keep peaks
PP.vad_db       = -65;                % looser VAD (only used if BYPASS_VAD=false)
PP.min_gap_s    = 0.15;               % merge close events
PP.min_dur_s    = 0.10;               % prune very short

% Controls
DBG.plot_both_classes = true;         % overlay class0 & class1 tracks
CTRL.BYPASS_VAD       = false;         % <<< start with VAD OFF to confirm detections
CTRL.FORCE_CLASS1     = true;         % <<< force class "1" as keyword

% Output
CFG.save_csv    = true;

%% ---------------- LOAD ----------------
S = load(CFG.net_path);
assert(isfield(S,CFG.net_var),"Network var '%s' not found in MAT.",CFG.net_var);
net = S.(CFG.net_var);
fprintf('[NET] Loaded "%s" (%s)\n',CFG.net_var,class(net));

[x,fs] = audioread(CFG.wav_path);
if size(x,2)>1, x = mean(x,2); end
if fs~=CFG.fs_target, x = resample(x,CFG.fs_target,fs); fs = CFG.fs_target; end
x = max(-1,min(1,x));
tSig = (0:numel(x)-1)/fs;

%% ---------------- FEATURE EXTRACTION (39-dim) ----------------
win = hamming(round(CFG.win_ms*1e-3*fs),'periodic');
hop = round(CFG.hop_ms*1e-3*fs);
M0  = mfcc(x,fs,"NumCoeffs",CFG.num_mfcc,"LogEnergy","Ignore", ...
           "Window",win,"OverlapLength",numel(win)-hop).';
if CFG.add_deltas
    D  = diff([M0(:,1) M0],1,2);
    DD = diff([D(:,1) D],1,2);
    feats = [M0; D; DD];   % 39 x T
else
    feats = M0;
end
if CFG.use_cmvn
    feats = feats - mean(feats,2);
    s = std(feats,0,2); s(s<1e-6) = 1;
    feats = feats ./ s;
end
T = size(feats,2);
tfrm = (0:T-1)*CFG.hop_ms*1e-3;

%% ---------------- INFERENCE (robust C×T) ----------------
tic;
Y = predict(net,{feats});
runtime_s = toc;
if iscell(Y), Y = Y{1}; end

% Normalize to C×T
ncls = numel(CFG.class_list);
if size(Y,2) == ncls && size(Y,1) ~= ncls
    Y = Y.';                         % T×C -> C×T
elseif size(Y,1) == ncls && size(Y,2) == 1
    Y = repmat(Y,1,T);               % clip-level -> broadcast
elseif size(Y,2) == ncls && size(Y,1) == 1
    Y = repmat(Y.',1,T);             % clip-level row -> broadcast
elseif size(Y,1) == T && size(Y,2) ~= ncls
    Y = Y.';                         % assume T×? -> transpose
end

% Softmax along classes
Y = Y - max(Y,[],1);
P = exp(Y) ./ max(sum(exp(Y),1), eps);

fprintf('[INFO] Posterior size = %d×%d (C×T)\n', size(P,1), size(P,2));
fprintf('[TIME] %.2f s audio processed in %.2f s (x%.2f real-time)\n', ...
        tSig(end), runtime_s, runtime_s/tSig(end));

% Diagnostics
p0 = P(1,:); p1 = P(2,:);
fprintf('[DIAG] class "0": mean=%.3f max=%.3f | class "1": mean=%.3f max=%.3f\n', ...
        mean(p0), max(p0), mean(p1), max(p1));

% Choose target class
if CTRL.FORCE_CLASS1
    idxTarget = 2;    % "1"
else
    % auto-pick class with higher variability
    [~, idxTarget] = max([std(p0) std(p1)]);
end
p_raw = P(idxTarget,:);

%% ---------------- POSTPROCESS ----------------
% Smoothing
if PP.smooth_win>1
    k = ones(1,PP.smooth_win)/PP.smooth_win;
    p_raw = conv(p_raw,k,'same');
end

% VAD
if CTRL.BYPASS_VAD
    vmask = true(1, size(P,2));
else
    vmask = vad_mask(x,fs,PP.vad_db,CFG);
    vmask = vmask(1:min(end,size(P,2)));
end
p_eff = p_raw; p_eff(~vmask) = 0;

% Adaptive hysteresis from clip stats (scale-aware)
on_th  = max(0.15, quantile(p_eff, 0.90));   % top % frames
off_th = max(0.10, 0.6*on_th);
fprintf('[THR] on=%.3f  off=%.3f  (adaptive)\n', on_th, off_th);

events = hyst_events(p_eff,tfrm,on_th,off_th,PP.min_gap_s,PP.min_dur_s);

% Peaks report (top 5)
[pk,loc] = findpeaks(p_eff, 'MinPeakDistance', max(1,round(0.25/(CFG.hop_ms*1e-3))));
[pk,ord] = sort(pk,'descend'); loc = loc(ord);
tpk = (loc-1)*CFG.hop_ms*1e-3;
fprintf('[PEAKS] '); fprintf('%.2fs(%.3f) ', tpk(1:min(5,end)), pk(1:min(5,end))); fprintf('\n');

%% ---------------- METRICS + CSV ----------------
dur = numel(x)/fs;
num_ev = size(events,1);
pct_frames = 100*sum(p_eff>on_th)/numel(p_eff);
fprintf('[DETECTIONS] %d events | %.1f%% frames > on_th\n', num_ev, pct_frames);

tbl = table(string(CFG.wav_path), num_ev, pct_frames, runtime_s, ...
    CTRL.BYPASS_VAD, on_th, off_th, ...
    'VariableNames',{'file','num_events','pct_frames_active','runtime_s','vad_bypassed','thr_on','thr_off'});
if CFG.save_csv
    fname = sprintf('results_kw1_%s.csv', datestr(now,'yyyymmdd_HHMMSS'));
    writetable(tbl,fname);
    fprintf('Saved: %s\n', fname);
end

%% ---------------- PLOTS ----------------
figure('Color','w','Name','Timeline');
subplot(3,1,1); plot(tSig,x); ylabel('Audio'); title('Waveform');

subplot(3,1,2); 
plot(tfrm,p_raw,'-','LineWidth',1.2); hold on;
plot(tfrm,p_eff,':','LineWidth',1.1);
yline(on_th,'--'); yline(off_th,':');
if DBG.plot_both_classes
    plot(tfrm,p0,'--'); plot(tfrm,p1,'-.');
    legend('raw','eff','on','off','class0','class1','Location','best');
else
    legend('raw','eff','on','off','Location','best');
end
ylabel('Posterior'); xlim([0 tSig(end)]);

subplot(3,1,3); hold on;
for i=1:num_ev, plot([events(i,1) events(i,2)],[0.5 0.5],'LineWidth',6); end
ylim([0 1]); xlim([0 tSig(end)]); ylabel('Events'); xlabel('Time (s)');
title('Detected keyword "1"');

%% ---------------- Local helpers ----------------
function vmask = vad_mask(x,fs,thr,cfg)
    win = round(cfg.win_ms*1e-3*fs);
    hop = round(cfg.hop_ms*1e-3*fs);
    en  = movmean(x.^2, win, 'Endpoints','shrink');
    en  = en(1:hop:end);
    db  = 10*log10(en + 1e-12);
    vmask = (db > thr).';
end

function ev = hyst_events(p,t,onThr,offThr,minGap,minDur)
    st=false; t0=NaN; ev=[];
    for i=1:numel(p)
        if ~st && p(i)>=onThr, st=true; t0=t(i); end
        if st && p(i)<offThr
            t1=t(i); if (t1-t0)>=minDur, ev=[ev; t0 t1]; end %#ok<AGROW>
            st=false; t0=NaN;
        end
    end
    if st, t1=t(end); if (t1-t0)>=minDur, ev=[ev; t0 t1]; end, end
    if isempty(ev), return; end
    m=ev(1,:);
    for k=2:size(ev,1)
        if ev(k,1)-m(end,2)<minGap, m(end,2)=ev(k,2);
        else, m=[m; ev(k,:)]; end %#ok<AGROW>
    end
    ev=m;
end