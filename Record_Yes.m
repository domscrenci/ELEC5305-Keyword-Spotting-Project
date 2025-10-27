% record_yes_wavs.m
% Record several "YES" utterances and save them as .wav files.

fs       = 16e3;           % sample rate (Hz)
nBits    = 16;             % bit depth
nChan    = 1;              % mono
numClips = 5;              % how many samples to record
dur_s    = 5;            % seconds per recording
pause_s  = 5;            % pause between takes

recObj = audiorecorder(fs, nBits, nChan);

disp('--------------------------------------------------');
disp('Recording YES clips...');
disp('You will be prompted before each take.');
disp('Press Ctrl+C to cancel early.');
disp('--------------------------------------------------');

for k = 1:numClips
    fprintf('\nClip %d/%d: Say "YES" after the beep...\n', k, numClips);
    beep;                                  % make a short beep
    pause(0.3);
    recordblocking(recObj, dur_s);         % record for dur_s seconds
    y = getaudiodata(recObj);
    y = y / max(abs(y)+eps);               % normalize amplitude
    fname = sprintf('yes__%02d.wav', k);
    audiowrite(fname, y, fs);
    fprintf('Saved: %s  (%.1f s)\n', fname, dur_s);
    pause(pause_s);
end

disp('All recordings complete!');