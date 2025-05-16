%% Awareness Cessation Project
% Mar 2025, Leila Salvesen
% adapted from EEG03_DataProcessing [G. Bernardi]

%% Initialization

clear variables; close all; clc;

addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0\');
addpath('C:\Users\Leila\Local\MATLAB\toolboxes\zapline-plus-main\');
datapath = 'C:\Users\Leila\Local\NSproject\data\';

%% Parameter Definition
params.evaluationChannels = 1:64; % Exclude non EEG channels
params.recreference       = [];   % Recording reference (excluded from bad-channel detection)   %!!! CHANGE (currently none)
params.heogchannel       = 65;   % Vertical EOG channel (excluded from bad-channel detection)
params.veogchannel       = 66;   % Vertical EOG channel (excluded from bad-channel detection)
params.emgchannel       = 67;   % Vertical EOG channel (excluded from bad-channel detection)
params.epochDur         = 30;     % Epoch duration in seconds
params.highPassFilter     = 0.5;   % High-pass filter for ICA (0 = no filter)
params.lowPassFilter      = 45;    % Low-pass filter for ICA (0 = no filter)
params.resampleICA        = 0;     % Resample the data for ICA saves computation time (0 = no resampling)
params.applyPCAtoICA      = 0;   % Dimensionality reduction prior to ICA (0 = no PCA)
params.brainThrICA        = 0.25;  % Minimum IClabel threshold for defining an IC as brain
params.brainNrtICA        = 0.80;  % Minimum ratio between brain and first noise IC
params.badchanmaxICs      = 30;    % Number of ICs to consider when looking for bad channels (between 0=none and 1=all)
params.badchanthrICA      = 0.80;  % Threshold for detecting bad channels based on ICA (0 = no detection)
params.useICLoverride     = 1;     % Override ICL classification to keep/reject specific components (0 = no)

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

%% Data preprocessing 

runtime = datetime;

%Define paths
outproc = 'C:\Users\Leila\Local\NSproject\data\reports_preproc'
outsets = 'C:\Users\Leila\Local\NSproject\data\proc_settings'
eeglab; close all;
disp(['>> Starting pre-processing']);

% Load and Epoch Data
clc; disp('Loading data...');
EEG_nap = pop_loadset('C:\Users\Leila\Local\NSproject\data\raw\nap\DA_fullNAP_RAW.set');
EEG_nap = pop_select(EEG_nap, 'channel', 1:64);
EEG_nap = eeg_checkset(EEG_nap);

% Remove Line noise
disp('>> Remove Line Noise');
EEG_nap = zapclean(EEG_nap, params, outproc);
EEG_nap = eeg_checkset(EEG_nap); % Check data

% Detect Bad Channels
disp('>> Bad Channel Detection');
badchans = findbadchans(EEG_nap, params);
indelec = sort(unique([badchans.all_marked{:}]));
goodchannels = setxor(params.evaluationChannels,indelec);
report_badchannels(badchans, params, EEG_nap.chanlocs, outproc);
% 
% % Save EOG Channels (L1G-R1G, VEOGR-R1Z)
% disp('>> Storing EOG Channels');
% EOG = eeg_interp(EEG_nap,indelec,'spherical'); % Use spherical spline for interpolation
% h1 = find(matches({EOG.chanlocs.labels},'L1G'));
% h2 = find(matches({EOG.chanlocs.labels},'R1G'));
% v1 = find(matches({EOG.chanlocs.labels},'VEOGR'));
% v2 = find(matches({EOG.chanlocs.labels},'R1Z'));
% eog_data(1,:,:) = EOG.data(h2,:,:)-EOG.data(h1,:,:); % Right - Left
% eog_data(2,:,:) = EOG.data(v2,:,:)-EOG.data(v1,:,:); % Superior - Inferior
% clear EOG v1 v2 h1 h2;

% Run Independent Component Analysis
disp('>> Independent Component Analysis');
EEGica = icarun(EEG_nap, goodchannels, params);
[EEGica, reject, keep, classcomps, badchans.bc_ica] = icreject(EEGica, params);
% if ~isempty(badchans.bc_ica) 
%     warning('Bad channels found during ICA-based artifact reduction!');
%     warning('I will exclude these channels and try again...');
%     clear EEGica indelec goodchannels reject keep classcomps;
%     indelec = sort(unique([badchans.all_marked{:},badchans.bc_ica]));
%     goodchannels = setxor(params.evaluationChannels, indelec);
%     report_badchannels(badchans, params, EEG_nap.chanlocs, outproc);
%     EEGica = icarun(EEG_nap, goodchannels, params);
%     [EEGica, reject, keep, classcomps, ~] = icreject(EEGica, params);
% end
if params.useICLoverride == 1
    [reject, keep] = ICLoverride(reject, keep, outsets, outproc);
end
report_iccleaning(EEGica, reject, keep, classcomps, goodchannels, outproc);
EEG_nap.icawinv     = EEGica.icawinv;
EEG_nap.icasphere   = EEGica.icasphere;
EEG_nap.icaweights  = EEGica.icaweights;
EEG_nap.icachansind = EEGica.icachansind;
EEG_nap = eeg_checkset(EEG_nap);  % Check data structure
clear EEGica; % This structure was needed only for IC-label application
EEG_nap.reject.gcompreject(reject) = true;
EEG_nap = pop_subcomp(EEG_nap, reject ,0);
EEG_nap = eeg_checkset(EEG_nap); % Check data

% Interpolate bad channels
disp('>> Bad Channel Interpolation');
EEG_nap = eeg_interp(EEG_nap,indelec,'spherical'); % Use spherical spline for interpolation
EEG_nap = eeg_checkset(EEG_nap); % Check data

% Re-add EOG channels
% EEG_nap.data(EEG_nap.nbchan+1:EEG_nap.nbchan+2,:,:) = eog_data;
% EEG_nap.chanlocs(EEG_nap.nbchan+1).labels = 'HEOG';
% EEG_nap.chanlocs(EEG_nap.nbchan+2).labels = 'VEOG';
% EEG_nap = pop_chanedit(EEG_nap, 'settype',{EEG_nap.nbchan+1,'EOG'});
% EEG_nap = pop_chanedit(EEG_nap, 'settype',{EEG_nap.nbchan+2,'EOG'});
% EEG_nap = eeg_checkset(EEG_nap); % Check data

% Save output
disp('Saving Preprocessed Data');
EEG_nap.preproc.myparams = params;
EEG_nap.preproc.badchans = badchans;
EEG_nap.preproc.icreject = reject;
EEG_nap.preproc.icccomps = classcomps;
EEG_nap.preproc.sruntime = runtime;
pop_saveset(EEG_nap, 'filename','DA_fullNAP_eegpreproc', 'filepath', 'C:\Users\Leila\Local\NSproject\data\preproc\nap\', 'version','7.3');

close all;
clearvars
disp('Program Completed!');