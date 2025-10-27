clear; clc;

% 1) List workspace nets
disp('--- Networks currently in workspace ---');
whos SeriesNetwork DAGNetwork dlnetwork

% 2) List variables inside your MAT file
disp('--- Variables inside KWSBaseline.mat ---');
whos -file KWSBaseline1.mat

% 3) Try to load it
S = load('KWSBaseline1.mat');
disp('--- Variable classes in KWSBaseline.mat ---');
for f = fieldnames(S)'
    fprintf('%s : %s\n', f{1}, class(S.(f{1})));
end
