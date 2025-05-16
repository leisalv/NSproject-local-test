%% Awareness Cessation Project
% Apr 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria

%%%%%% CESSATION SESSION (NAP + WAKE RS + NS) DATA PREPROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

mymachine = input('Select local (1) or server (0):');
if mymachine == 1
    addpath(genpath(codepath));
    addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0'); eeglab;
    addpath('C:\Users\Leila\Local\MATLAB\toolboxes\fieldtrip-20240110'); ft_defaults;
    datapath = 'C:\Users\Leila\Local\NSproject\data';
    derivpath = 'C:\Users\Leila\Local\NSproject\data\derivatives';
elseif mymachine == 0
    addpath(genpath(codepath))
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/eeglab2025.0.0'); eeglab;
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/fieldtrip-20240110'); ft_defaults;
    datapath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data';
    derivpath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data/derivatives';
end

mydate = datetime('now');
myoutpath = fullfile(derivpath, strcat(codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'))); mkdir(myoutpath);
myoutdata = fullfile(myoutpath, 'prep_data'); mkdir(myoutdata);
myoutreps = fullfile(myoutpath, 'prep_reps'); mkdir(myoutreps);

% if isempty(gcp('nocreate'))
%     parpool;  % Initializes default parallel pool with all available cores
% end

%% Parameter definition
params.evaluationChannels = 1:64; % Exclude non EEG channels
params.recreference       = [];   % Recording reference (excluded from bad-channel detection)   %!!!??? CHANGE (currently none)
params.newreference       = [];   % New reference for data re-referencing
params.newsr              = 256;  % Resampling rate (Hz)
params.heogchannel        = 65;   % Vertical EOG channel (excluded from bad-channel detection)
params.veogchannel        = 66;   % Vertical EOG channel (excluded from bad-channel detection)
params.emgchannel         = 67;   % Vertical EOG channel (excluded from bad-channel detection)
params.highPassFilter     = 0.5;  % High-pass filter for ICA (0 = no filter)
params.lowPassFilter      = 45;   % Low-pass filter for ICA (0 = no filter)
params.resampleICA        = 0;  % Resample the data for ICA saves computation time (0 = no resampling)
params.applyPCAtoICA      = 0;   % Dimensionality reduction prior to ICA (0 = no PCA)
params.brainThrICA        = 0.25; % Minimum IClabel threshold for defining an IC as brain
params.brainNrtICA        = 0.80; % Minimum ratio between brain and first noise IC
params.badchanmaxICs      = 1;   % Number of ICs to consider when looking for bad channels (proportion; between 0=none and 1=all)
params.badchanthrICA      = 0.80; % Threshold for detecting bad channels based on ICA (0 = no detection)
params.useICLoverride     = 0;    % Override ICL classification to keep/reject specific components (0 = no)

%% Bad channel detection parameters
params.badchandet.robustDeviationThreshold    = []; % Z score cutoff for robust channel deviation (default = 5)
params.badchandet.highFrequencyNoiseThreshold = []; % Z score cutoff for SNR (signal above 50 Hz; default = 5)
params.badchandet.correlationWindowSeconds    = []; % Correlation window size in seconds (default = 1 sec)
params.badchandet.correlationThreshold        = []; % Correlation below which window is bad (default = 0.4)
params.badchandet.badTimeThreshold            = []; % Cutoff fraction of bad corr windows (default = 0.01)
params.badchandet.ransacSampleSize            = []; % Samples for computing ransac (default = 50)
params.badchandet.ransacChannelFraction       = []; % Fraction of channels for robust reconstruction (default = 0.25)
params.badchandet.ransacCorrelationThreshold  = []; % Cutoff correlation for abnormal wrt neighbors (default = 0.75)
params.badchandet.ransacUnbrokenTime          = []; % Cutoff fraction of time channel can have poor ransac predictability (default = 0.4)
params.badchandet.ransacWindowSeconds         = []; % Correlation window for ransac (default = 5 sec)

%% Load data
[rawdatafile, rawdatapath] = uigetfile({'*.set', 'Raw EEG Dataset Files (*.set)'; '*.*', 'All Files (*.*)'}, 'Select Raw EEG Dataset File');
if isequal(rawdatafile, 0)
    error('No file selected. Exiting...');
else
    myfilepath = fullfile(rawdatapath, rawdatafile);
    disp(['Selected file: ', myfilepath]);
end

% Load nap data
clc; disp('Loading data...');
myeeg_raw = pop_loadset(myfilepath);
myeeg_raw = eeg_checkset(myeeg_raw);
eeglab redraw;

% Load sleep scoring data
[scdatafile, scdatapath] = uigetfile({'*.txt', 'Sleep Scoring Files (*.txt)'; '*.*', 'All Files (*.*)'}, 'Select Sleep Scoring File');
if isequal(scdatafile, 0)
    error('No file selected. Exiting...');
else
    myscfilepath = fullfile(scdatapath, scdatafile);
    disp(['Selected file: ', myscfilepath]);
end     

% Store sleep stages into EEG structure as events
myeeg_rawsc = pop_importevent(myeeg_raw, 'event', myscfilepath, 'fields', {'latency' 'type'}, 'timeunit', 1, 'skipline', 1);
for i = 1:length(myeeg_rawsc.event) % Add 30s duration (in points)
    myeeg_rawsc.event(i).duration = 30 * myeeg_rawsc.srate;
end
myeeg_rawsc = eeg_checkset(myeeg_rawsc);
myeeg_rawsc = pop_saveset(myeeg_rawsc, 'filename', 'myeeg_rawsc.set', 'filepath', myoutdata);
% eeglab redraw

%% Preprocessing 
runtime = datetime;

% Select scalp EEG data
myeeg_in = pop_select(myeeg_rawsc, 'channel', params.evaluationChannels);  

% Apply high-pass filter

% Downsampling
disp('>> Downsample');
myeeg_rs = pop_resample(myeeg_in, params.newsr);                            
myeeg_rs = eeg_checkset(myeeg_rs); % Check data
myeeg_rs = pop_saveset(myeeg_rs, 'filename', 'myeeg_rs.set', 'filepath', myoutdata);
eeglab redraw

% Remove Line noise
disp('>> Remove Line Noise');
myeeg_zl = myzapcleaner(myeeg_rs, myoutreps);
myeeg_zl = eeg_checkset(myeeg_zl); % Check data
myeeg_zl = pop_saveset(myeeg_zl, 'filename', 'myeeg_cl.set', 'filepath', myoutdata);

% Detect Bad Channels 
disp('>> Bad Channel Detection');
badchans = mybadchanfinder(myeeg_zl, params);   
indbadelec = badchans.all_marked;
goodchans = setxor(params.evaluationChannels, indbadelec);
mybadchansrep(badchans, params, myeeg_zl, myoutreps);

% Detect & Reject artifacted data
disp('>> Artifact Rejection');
myeeg_cl = pop_clean_rawdata(myeeg_zl); % Select high-pass (needed for ASR) + ASR without reconstruction + 

%% Run Independent Component Analysis
disp('>> Independent Component Analysis');
myeeg_ica = myicarun(myeeg_cl, goodchans, params);
[myeeg_ica, reject, keep, classcomps, badchans.bc_ica] = myicreject(myeeg_ica, params, myoutreps);

if ~isempty(badchans.bc_ica)
    warning('Bad channels found during ICA-based artifact reduction!');
    warning('I will exclude these channels and try again...');
    clear goodchans indbadelec myeeg_ica reject keep classcomps;
    indbadelec = sort(unique([badchans.all_marked, badchans.bc_ica]));
    goodchans = setxor(params.evaluationChannels, indbadelec);
    mybadchansrep(badchans, params, myeeg_cl, myoutreps);
    myeeg_ica = myicarun(myeeg_cl, goodchans, params);
    [myeeg_ica, reject, keep, classcomps, ~] = myicreject(myeeg_ica, params, myoutreps);
end

if params.useICLoverride == 1 % If needed, use to correct based on visual detection creating txt files
    [reject, keep] = ICLoverride(reject, keep, myoutreps, myoutreps);
end

report_iccleaning(myeeg_ica, reject, keep, classcomps, goodchans, myoutreps);

myeeg_cl.icawinv     = myeeg_ica.icawinv;
myeeg_cl.icasphere   = myeeg_ica.icasphere;
myeeg_cl.icaweights  = myeeg_ica.icaweights;
myeeg_cl.icachansind = myeeg_ica.icachansind;
myeeg_cl = eeg_checkset(myeeg_cl);  % Check data structure
clear myeeg_ica; % This structure was needed only for IC-label application

myeeg_cl.reject.gcompreject(reject) = true;
myeeg = pop_subcomp(myeeg_cl, reject, 0);
myeeg = eeg_checkset(myeeg); % Check data

% Rereference to avg
%myeeg_reref = pop_reref(myeeg, params.newreference);                   %%Rereference ???!!!!!

% Interpolate bad channels
disp('>> Bad Channel Interpolation');
myeeg = eeg_interp(myeeg, indbadelec, 'spherical'); % Use spherical spline for interpolation
myeeg = eeg_checkset(myeeg); % Check data

% Save output
disp('Saving Preprocessed Data');
myeeg.preproc.myparams = params;
myeeg.preproc.badchans = badchans;
myeeg.preproc.icreject = reject;
myeeg.preproc.icccomps = classcomps;
myeeg.preproc.sruntime = runtime;
pop_saveset(myeeg, 'filename', 'DA_fullNAP_PREP.set', 'filepath', myoutdata, 'version','7.3');
close all;

disp('Program Completed!');