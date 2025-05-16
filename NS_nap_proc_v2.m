%% Awareness Cessation Project
% Mar 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria


%%%%%% STEP 2: CESSATION SESSION (NAP + WAKE RS + NS) DATA PROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

mymachine = input('Select local (1) or server (0):');
if mymachine == 1
    addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0\'); eeglab;
    addpath('C:\Users\Leila\Local\MATLAB\toolboxes\fieldtrip-20240110\'); ft_defaults;
    addpath(genpath(codepath))
    datapath = 'C:\Users\Leila\Local\NSproject\data\';
    derivpath = 'C:\Users\Leila\Local\NSproject\data\derivatives\';
    if isempty(gcp('nocreate'))
    parpool;  % Initializes default parallel pool with all available cores
end
elseif mymachine == 0
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/eeglab2025.0.0/'); eeglab;
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/fieldtrip-20240110/'); ft_defaults;
    addpath(genpath(codepath))
    datapath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data/';
    derivpath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data/derivatives/';
end

mydate = datetime('now');
myoutpath = [derivpath, codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'), '\'];
mkdir(myoutpath);

% Define frequency ranges
freq_bands = [0.5 2; 2 4; 4 8; 8 12; 12 16; 16 30; 30 45; 45 80];
slope_ranges = [0.5 30; 30 45; 0.5 45; 1 20; 20 40];

%% Load data

clc; disp('Loading data...');

total_time = tic;  % Start timing

% Load nap session data
myeeg_prep = pop_loadset(fullfile(derivpath, 'NS_nap_preproc_v9', 'prep_data', 'DA_fullNAP_PREP.set'));
myeeg_prep = eeg_checkset(myeeg_prep);
eeglab redraw;

%% Define sleep stages segments within nap session 
% (N1 vs Wake)


% Step 1: Get all events
all_events = {EEG.event.type};

% Step 2: Prepare new event list marking durations
eventList = EEG.event;
N = length(eventList);
%%%% N1

n1_bouts = find_stage_segments(myeeg_prep.preproc.stages, 1);
EEG_n1_bouts = extract_eeg_segments(myeeg_prep, n1_bouts);

pow_n1_bouts = NaN(length(EEG_n1_bouts), myeeg_prep.nbchan, length(freq_bands));
osc_n1_bouts = NaN(length(EEG_n1_bouts), myeeg_prep.nbchan, length(freq_bands));
slp_n1_bouts = NaN(length(EEG_n1_bouts), myeeg_prep.nbchan, length(slope_ranges));
itc_n1_bouts = NaN(length(EEG_n1_bouts), myeeg_prep.nbchan, length(slope_ranges));

pow_n1_epochs = cell(length(EEG_n1_bouts), 1); % Initialize cell arrays to store raw segment data before averaging
osc_n1_epochs = cell(length(EEG_n1_bouts), 1);
slp_n1_epochs = cell(length(EEG_n1_bouts), 1); % Initialize cell arrays to store raw segment data before averaging
itc_n1_epochs = cell(length(EEG_n1_bouts), 1);

for n1b = 1:length(EEG_n1_bouts)

    fprintf('Loading bout %d / %d\n', n1b, length(EEG_n1_bouts));
    
    [orig_bout, osc_bout, slp_bout, itc_bout] = analyze_eeg_segments(EEG_n1_bouts{n1b}, freq_bands, slope_ranges);
    
    %Store epoch-level data before averaging
    pow_n1_epochs{n1b} = orig_bout;  
    osc_n1_epochs{n1b} = osc_bout;  
    slp_n1_epochs{n1b} = slp_bout;
    itc_n1_epochs{n1b} = itc_bout;  

    pow_n1_bouts(n1b,:,:) = squeeze(mean(orig_bout, 1)); % average across 5s epochs
    osc_n1_bouts(n1b,:,:) = squeeze(mean(osc_bout, 1)); % average across 5s epochs
    slp_n1_bouts(n1b,:,:) = squeeze(mean(slp_bout, 1));
    itc_n1_bouts(n1b,:,:) = squeeze(mean(itc_bout, 1));
end

% Save averaged results AND epoch-level data
save(fullfile(myoutpath, 'nap_freqslp_analysis_n1_output.mat'), ...
    'pow_n1_bouts', 'osc_n1_bouts', 'slp_n1_bouts', 'itc_n1_bouts', ...
    'pow_n1_epochs', 'osc_n1_epochs', 'slp_n1_epochs', 'itc_n1_epochs')


elapsed_time = toc(total_time);  % End timing
fprintf('Total running time : %.2f minutes (%.2f hours)\n', elapsed_time/60, elapsed_time/3600);

disp('Analysis results and epoch-level data saved successfully.');

% Total running time: 26.81 minutes (0.45 hours)


%%%% Wake

wake_bouts = find_stage_segments(myeeg_prep.preproc.stages, 0);
EEG_wake_bouts = extract_eeg_segments(myeeg_prep, wake_bouts);

pow_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(freq_bands));
osc_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(freq_bands));
slp_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(slope_ranges));
itc_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(slope_ranges));

pow_wake_epochs = cell(length(EEG_wake_bouts), 1);
osc_wake_epochs = cell(length(EEG_wake_bouts), 1); % Initialize cell arrays to store raw segment data before averaging
slp_wake_epochs = cell(length(EEG_wake_bouts), 1);
itc_wake_epochs = cell(length(EEG_wake_bouts), 1);


for wkb = 1:length(EEG_wake_bouts)
    fprintf('Loading bout %d / %d\n', wkb, length(EEG_wake_bouts));
    
    [orig_bout, osc_bout, slp_bout, itc_bout] = analyze_eeg_segments(EEG_wake_bouts{wkb}, freq_bands, slope_ranges);
    
    % Store epoch-level data before averaging
    pow_wake_epochs{wkb} = orig_bout;
    osc_wake_epochs{wkb} = osc_bout;  
    slp_wake_epochs{wkb} = slp_bout;
    itc_wake_epochs{wkb} = itc_bout;

    pow_wake_bouts(wkb,:,:) = squeeze(mean(orig_bout, 1)); % average across 5s epochs
    osc_wake_bouts(wkb,:,:) = squeeze(mean(osc_bout, 1)); 
    slp_wake_bouts(wkb,:,:) = squeeze(mean(slp_bout, 1));
    itc_wake_bouts(wkb,:,:) = squeeze(mean(itc_bout, 1));
end

% Save averaged results AND epoch-level data
save(fullfile(myoutpath, 'nap_freqslp_analysis_wk_output.mat'), ...
    'osc_wake_bouts', 'slp_wake_bouts', 'osc_wake_epochs', 'slp_wake_epochs', ...
    '-v7.3'); 

%%Last run: Total running time: 4957.42 minutes (82.62 hours) !!! 
%%(Pb: no pow/orig nor itc saved; nested parfor..)


elapsed_time = toc(total_time);  % End timing
fprintf('Total running time: %.2f minutes (%.2f hours)\n', elapsed_time/60, elapsed_time/3600);

disp('Analysis results and epoch-level data saved successfully.');