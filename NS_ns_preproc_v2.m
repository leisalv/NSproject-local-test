%% Awareness Cessation Project
% Mar 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria


%%%%%% CESSATION SESSION (NAP + WAKE RS + NS) DATA PREPROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0\'); eeglab;
addpath('C:\Users\Leila\Local\MATLAB\toolboxes\zapline-plus-main\');
addpath(genpath(codepath))

inpath = 'C:\Users\Leila\Local\NSproject\data\';
outpath = 'C:\Users\Leila\Local\NSproject\data\derivatives\';

mydate = datetime('now');
myoutpath = [outpath, codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'), '\'];
mkdir(myoutpath);
myoutdata = strcat(myoutpath, 'prep_data', '\');
myoutreps = strcat(myoutpath, 'prep_reps', '\');
myoutsets = strcat(myoutpath, 'prep_sets', '\');
mkdir(myoutdata); mkdir(myoutreps); mkdir(myoutsets);

%% Load data

clc; disp('Loading data...');
% Load NS data
myeeg_raw = pop_loadset([inpath, 'raw\med\DA_meditation_RAW.set']);
recsr = myeeg_raw.srate;
myeeg_raw = eeg_checkset(myeeg_raw);
eeglab redraw;


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
params.resampleICA        = 250;  % Resample the data for ICA saves computation time (0 = no resampling)
params.applyPCAtoICA      = 30;   % Dimensionality reduction prior to ICA (0 = no PCA)
params.brainThrICA        = 0.25; % Minimum IClabel threshold for defining an IC as brain
params.brainNrtICA        = 0.80; % Minimum ratio between brain and first noise IC
params.badchanmaxICs      = 30;   % Number of ICs to consider when looking for bad channels (between 0=none and 1=all)
params.badchanthrICA      = 0.80; % Threshold for detecting bad channels based on ICA (0 = no detection)
params.useICLoverride     = 1;    % Override ICL classification to keep/reject specific components (0 = no)


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

%% Preprocessing 
runtime = datetime;

% % Save RNG seed for replicability
% seedpath = ([myoutpath 'RNGseed.mat']);
% if exist([seedpath],'file')>0 && replace_seed==0
%     temp = load([seedpath]);
%     warning('The program is using a previously saved RNG seed!');
%     rng(temp.seed);
% else
%     rng('shuffle');
%     seed = rng;
%     cdate = datetime;
%     save([seedpath],'seed','cdate'); %removed filename
% end; clear temp seed cdate seedpath;


%myeeg = pop_reref(myeeg_raw, params.newreference);                   %%Rereference !!!!!
myeeg = pop_select(myeeg_raw, 'channel', params.evaluationChannels);  % Include only EEG channels
myeeg = pop_resample(myeeg, params.newsr);                            % Downsampling
myeeg = eeg_checkset(myeeg); % Check data

% Remove Line noise
disp('>> Remove Line Noise');
myeeg = zapclean(myeeg, params, myoutreps);
myeeg = eeg_checkset(myeeg); % Check data

% Detect Bad Channels
disp('>> Bad Channel Detection');
badchans = findbadchans(myeeg, params);                                      % Kept Dami's version
indbadelec = sort(unique([badchans.all_marked{:}]));
goodchannels = setxor(params.evaluationChannels, indbadelec);
report_badchannels(badchans, params, myeeg.chanlocs, myoutreps);
myeeg = eeg_checkset(myeeg); % Check data

%% Run Independent Component Analysis
disp('>> Independent Component Analysis');

% epochs = find(myeegstg==0|myeegstg==1|myeegstg==2|myeegstg==3|myeegstg==5);  % modify to split into NREM1 / wake??
% myeeg_preica=myeeg;
% myeeg_preica.data=myeeg_preica.data(:,epochs);
% myeeg_preica.times=myeeg_preica.times(epochs);
% myeeg_preica.pnts=length(myeeg_preica.times);
% myeeg_ica = icarunnight(myeeg_preica, goodchannels, params);

myeeg_ica = icarunnight(myeeg, goodchannels, params);
[myeeg_ica, reject, keep, classcomps, badchans.bc_ica] = icreject(myeeg_ica, params);
% [reject, keep, classcomps] = icreject(myeeg_ica, params);

if ~isempty(badchans.bc_ica)
    warning('Bad channels found during ICA-based artifact reduction!');
    warning('I will exclude these channels and try again...');
    clear myeeg_ica indbadelec goodchannels reject keep classcomps;
    indbadelec = sort(unique([badchans.all_marked{:}, badchans.bc_ica]));
    goodchannels = setxor(params.evaluationChannels, indbadelec);
    report_badchannels(badchans, params, myeeg.chanlocs, myoutreps);
    myeeg_ica = icarunnight(myeeg, goodchannels, params);
    [myeeg_ica, reject, keep, classcomps, ~] = icreject(myeeg_ica, params);
end

if params.useICLoverride == 1
    [reject, keep] = ICLoverride(reject, keep, myoutreps, myoutreps);
end

report_iccleaning(myeeg_ica, reject, keep, classcomps, goodchannels, myoutreps);
myeeg.icawinv     = myeeg_ica.icawinv;
myeeg.icasphere   = myeeg_ica.icasphere;
myeeg.icaweights  = myeeg_ica.icaweights;
myeeg.icachansind = myeeg_ica.icachansind;
myeeg = eeg_checkset(myeeg);  % Check data structure
clear myeeg_ica; % This structure was needed only for IC-label application

% If needed, add manually some IC after visual inspection
% manualadd = [11]; % MODIFY DEPENDING ON VISUAL INSPECTION
% reject = setxor(reject, manualadd);
% clear manualadd;

myeeg.reject.gcompreject(reject) = true;
myeeg = pop_subcomp(myeeg, reject, 0);
myeeg = eeg_checkset(myeeg); % Check data

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