%% Awareness Cessation Project
% Apr 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria

% v11: considers recording ref in avg + tests M1/M2

%%%%%% CESSATION SESSION (NAP + WAKE RS + NS) DATA PREPROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

mymachine = input('Select local (1) or server (0):');
if mymachine == 1
    addpath(genpath(codepath));
    addpath('C:\Users\Leila\VULocal\MATLAB\toolboxes\eeglab2025.0.0'); eeglab;
    addpath('C:\Users\Leila\VULocal\MATLAB\toolboxes\fieldtrip-20240110'); ft_defaults;
    datapath = 'C:\Users\Leila\VULocal\NSproject\data';
    derivpath = 'C:\Users\Leila\VULocal\NSproject\data\derivatives';
elseif mymachine == 0
    addpath(genpath(codepath))
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/eeglab2025.0.0'); eeglab;
    addpath('/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/MATLAB/toolboxes/fieldtrip-20240110'); ft_defaults;
    datapath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data';
    derivpath = '/specialmount/pscst001-krb5/FGB-ETP-CogPsy-CessationStudy/LEILA/NSproject/data/derivatives';
    if isempty(gcp('nocreate'))
        parpool;  % Initializes default parallel pool with all available cores
    end
end

mydate = datetime('now');
myoutpath = fullfile(derivpath, strcat(codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'))); mkdir(myoutpath);
myoutdatapath = fullfile(myoutpath, 'prep_data'); mkdir(myoutdatapath);
myoutrepspath = fullfile(myoutpath, 'prep_reps'); mkdir(myoutrepspath);

%% Parameter definition
params.recref             = [];   % Recording reference (excluded from bad-channel detection)   %!!!??? CHANGE (currently none)
params.newref             = 'Oz';   % New reference for data re-referencing
params.newsr              = 256;  % Resampling rate (Hz)
params.eegchans           = 1:64;  
params.heogchannel        = 65;   % Vertical EOG channel (excluded from bad-channel detection)
params.veogchannel        = 66;   % Vertical EOG channel (excluded from bad-channel detection)
params.emgchannel         = 67;   % Vertical EOG channel (excluded from bad-channel detection)
params.outchans           = 65:67; %13, 19, 
params.inchans            = setdiff(params.eegchans, params.outchans);
params.evaluationChannels = 1:length(params.inchans); % Exclude non EEG channels
params.highPassFilter     = 0;    % High-pass filter pre-ASR with mirroring (0 = no filter)
params.mirrorWindow       = 5;    % Length of mirrored fake data to compute hp filtering 
params.highPassFilterICA  = 0.5;  % High-pass filter for ICA (0 = no filter)
params.lowPassFilterICA   = 45;  % Low-pass filter for ICA (0 = no filter)
params.resampleICA        = 0;    % Resample the data for ICA saves computation time (0 = no resampling)
params.applyPCAtoICA      = 0;    % Dimensionality reduction prior to ICA (0 = no PCA)
params.brainThrICA        = 0.25; % Minimum IClabel threshold for defining an IC as brain
params.brainNrtICA        = 0.80; % Minimum ratio between brain and first noise IC
params.badchanmaxICs      = 1;    % Number of ICs to consider when looking for bad channels (proportion; between 0=none and 1=all)
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

% Select meditation data
clc; disp('Loading data...');
myeeg_raw = pop_loadset(myfilepath);
sr = myeeg_raw.srate;

% Define NS section within meditation session (10min rs - 90min NS - 10min rs)
med_start = 10*60*sr; %10mins * 60secs * Hz
med_end = 100*60*sr; %end = (10+90mins * 60secs * Hz)
myeeg_med = pop_select(myeeg_raw, 'point',[med_start med_end-1]);
myeeg_med = eeg_checkset(myeeg_med);

% Set 30s events to homogeneise processing with nap data
maxep = floor(myeeg_med.pnts / (30 * sr)); % Number of 30s epochs
myevents = struct('type', {}, 'latency', {}, 'duration', {}); % Create new event structure

for i = 1:maxep
    myevents(i).type = '30sEp';
    myevents(i).latency = (i - 1) * 30 * sr + 1; % latency in points
    myevents(i).duration = 30 * sr; % duration in points
end

myeeg_med.event = myevents;
myeeg_med = eeg_checkset(myeeg_med);

% Set included EEG channels
myeeg_in = pop_select(myeeg_med, 'channel', params.eegchans);  
myeeg_in.ref = params.recref; % set recording reference to avg
myeeg_in = eeg_checkset(myeeg_in);
myeeg_in = pop_saveset(myeeg_in, 'filename', 'myeeg_in.set', 'filepath', myoutdatapath);

%% Testing set to check noisy chans (M1/M2)
% 
% % Find indices of electrodes of interest
% idxM1 = find(strcmp({myeeg_in.chanlocs.labels}, 'M1'));
% idxM2 = find(strcmp({myeeg_in.chanlocs.labels}, 'M2'));
% idxCz = find(strcmp({myeeg_in.chanlocs.labels}, 'Cz'));
% 
% 
% %%% Create the new signal: M1 - M2
% myeeg_M1M2 = myeeg_in;
% newChanData = myeeg_M1M2.data(idxM1,:) - myeeg_M1M2.data(idxM2,:);
% myeeg_M1M2.data(end+1,:) = newChanData; % Add the new channel to EEG.data
% myeeg_M1M2.nbchan = size(myeeg_M1M2.data, 1); % Update number of channels
% 
% % Copy chanloc structure from an existing channel and edit label
% newChanLoc = myeeg_M1M2.chanlocs(idxM1);  % Copy M1 metadata as template
% newChanLoc.labels = 'M1-M2';
% newChanLoc.urchan = myeeg_M1M2.nbchan + 1;
% myeeg_M1M2.chanlocs(end+1) = newChanLoc; % Add to chanlocs
% 
% myeeg_M1M2 = eeg_checkset(myeeg_M1M2);
% myeeg_M1M2 = pop_saveset(myeeg_M1M2, 'filename', 'myeeg_M1-M2.set', 'filepath', myoutdatapath);
% 
% 
% %%% Create the new signal: M2-Cz / M1-Cz
% myeeg_MCz = myeeg_in;
% newChanData1 = myeeg_MCz.data(idxM1,:) - myeeg_MCz.data(idxCz,:);
% newChanData2 = myeeg_MCz.data(idxM2,:) - myeeg_MCz.data(idxCz,:);
% 
% myeeg_MCz.data(end+1,:) = newChanData1; % Add the new channels to EEG.data
% myeeg_MCz.data(end+1,:) = newChanData2;
% myeeg_MCz.nbchan = size(myeeg_MCz.data, 1); % Update number of channels
% 
% % Copy chanloc structures from an existing channel and edit labels
% newChanLoc1 = myeeg_MCz.chanlocs(idxM1);  % Copy metadata template
% newChanLoc1.labels = 'M1-Cz';
% newChanLoc1.urchan = myeeg_MCz.nbchan + 1;
% myeeg_MCz.chanlocs(end+1) = newChanLoc1; % Add to chanlocs
% 
% newChanLoc2 = myeeg_MCz.chanlocs(idxM1);  % Copy metadata template 
% newChanLoc2.labels = 'M2-Cz';
% newChanLoc2.urchan = myeeg_MCz.nbchan + 1;
% myeeg_MCz.chanlocs(end+1) = newChanLoc2; % Add to chanlocs
% 
% myeeg_MCz = eeg_checkset(myeeg_MCz);
% myeeg_MCz = pop_saveset(myeeg_MCz, 'filename', 'myeeg_M-Cz.set', 'filepath', myoutdatapath);
% 


%% Preprocessing 
runtime = datetime;

% Initial reref to Oz
idxRef = find(strcmp({myeeg_in.chanlocs.labels}, params.newref));
myeeg_reref = pop_reref(myeeg_in, idxRef, 'keepref', 'on');                
myeeg_reref = eeg_checkset(myeeg_reref);
myeeg_reref = pop_saveset(myeeg_reref, 'filename', 'myeeg_reref.set', 'filepath', myoutdatapath);

% High-pass filtering - set to 0 when doing highpass filtering during ASR
% WARNING: 'raw' data is 0.1 HZ highpass filtered - avoid further filtering
if params.highPassFilter>0  
    disp('>> High-pass filtering');
    [myeeg_hp, myeeg_mir] = myhpfilt(myeeg_in, params.highPassFilter, params.mirrorWindow); 
    myeeg_hp = pop_saveset(myeeg_hp, 'filename', 'myeeg_hp.set', 'filepath', myoutdatapath);
end

% Downsample
disp('>> Downsampling');
if exist('myeeg_hp','var')
    myeeg_rs = pop_resample(myeeg_hp, params.newsr);
else
    myeeg_rs = pop_resample(myeeg_reref, params.newsr);
end
myeeg_rs = eeg_checkset(myeeg_rs); % Check data

% Remove Line noise
disp('>> Removing Line Noise');
myeeg_zl = myzapcleaner(myeeg_rs, myoutrepspath);
myeeg_zl = eeg_checkset(myeeg_zl); % Check data
myeeg_zl = pop_saveset(myeeg_zl, 'filename', 'myeeg_preASR.set', 'filepath', myoutdatapath);

% Detect Bad Channels 
disp('>> Bad Channel Detection');
badchans = mybadchanfinder(myeeg_zl, params);   
indbadelec = badchans.all_marked;
mybadchansrep(badchans, params, myeeg_zl, myoutrepspath);
goodchans = params.inchans;
goodchans(indbadelec) = [];
[~, indgoodchans] = ismember(goodchans, [myeeg_zl.chanlocs.urchan]);

% Reject Artifacted data segments (ASR)
disp('>> Artifact Rejection');
myeeg_asr = pop_clean_rawdata(myeeg_zl, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
myeeg_asr = eeg_checkset(myeeg_asr); % Check data
myeeg_asr = pop_saveset(myeeg_asr, 'filename', 'myeeg_postASR.set', 'filepath', myoutdatapath);

%% Run Independent Component Analysis
disp('>> Independent Component Analysis');
myeeg_ica = myicarun(myeeg_asr, indgoodchans, params);
[myeeg_ica, reject, keep, classcomps, badchans.bc_ica] = myicreject(myeeg_ica, params, myoutrepspath);

if ~isempty(badchans.bc_ica)
    warning('Bad channels found during ICA-based artifact reduction!');
    warning('I will exclude these channels and try again...');
    clear goodchans indbadelec indgoodchans myeeg_ica reject keep classcomps;
    indbadelec = sort(unique([badchans.all_marked, badchans.bc_ica]));
    mybadchansrep(badchans, params, myeeg_asr, myoutrepspath);
    goodchans = params.inchans;
    goodchans(indbadelec) = [];   
    [~, indgoodchans] = ismember(goodchans, [myeeg_asr.chanlocs.urchan]);
    myeeg_ica = myicarun(myeeg_asr, indgoodchans, params);
    [myeeg_ica, reject, keep, classcomps, ~] = myicreject(myeeg_ica, params, myoutrepspath);
end

if params.useICLoverride == 1 % If needed, use to correct based on visual detection creating txt files
    [reject, keep] = ICLoverride(reject, keep, myoutrepspath, myoutrepspath);
end

myiccleanrep(myeeg_ica, reject, keep, classcomps, indgoodchans, myoutrepspath);

myeeg_asr.icawinv     = myeeg_ica.icawinv;
myeeg_asr.icasphere   = myeeg_ica.icasphere;
myeeg_asr.icaweights  = myeeg_ica.icaweights;
myeeg_asr.icachansind = myeeg_ica.icachansind;
myeeg_asr = eeg_checkset(myeeg_asr);  % Check data structure
clear myeeg_ica; % This structure was needed only for IC-label application

myeeg_asr.reject.gcompreject(reject) = true;
myeeg_postica = pop_subcomp(myeeg_asr, reject, 0);
myeeg_postica = eeg_checkset(myeeg_postica); % Check data
myeeg_postica = pop_saveset(myeeg_postica, 'filename', 'myeeg_postICA.set', 'filepath', myoutdatapath);

% Interpolate bad channels
disp('>> Bad Channel Interpolation');
myeeg_itp = eeg_interp(myeeg_postica, indbadelec, 'spherical'); % Use spherical spline for interpolation
myeeg_itp = eeg_checkset(myeeg_itp); % Check data

% Exclude mastoids & rereference to avg
idxM1 = find(strcmp({myeeg_itp.chanlocs.labels}, 'M1'));
idxM2 = find(strcmp({myeeg_itp.chanlocs.labels}, 'M2'));
myeeg_noM1M2 = pop_select(myeeg_itp, 'nochannel', [idxM1 idxM2]);
myeeg_out = pop_reref(myeeg_noM1M2, []);                   %%Rereference to avg after ICA

% Save output
disp('Saving Preprocessed Data');
myeeg_out.preproc.myparams = params;
myeeg_out.preproc.badchans = badchans;
myeeg_out.preproc.icreject = reject;
myeeg_out.preproc.icccomps = classcomps;
myeeg_out.preproc.sruntime = runtime;

% Epoch & save data
pop_saveset(myeeg_out, 'filename', 'DA_med_PREP.set', 'filepath', myoutdatapath, 'version','7.3');

myeeg_epMed = pop_epoch(myeeg_out, {'30sEp'}, [0 30], 'epochinfo', 'yes');
pop_saveset(myeeg_epMed, 'filename', 'DA_med_PREP.set', 'filepath', myoutdatapath, 'version','7.3');

close all;
disp('Program Completed!');