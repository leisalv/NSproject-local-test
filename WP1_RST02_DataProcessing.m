%% EEG data preprocessing for TweakDreams Project
% Jul 2024, Giulio Bernardi [giulio.bernardi@imtlucca.it]
%

%% Initialization

restoredefaultpath; clear variables; close all; clc;

addpath('/opt/MATLAB/TOOLBOX/eeglab2021.1');
addpath('/home/DATA/SCRIPTS/tweakdreams/');
addpath('/home/DATA/SCRIPTS/zapline-plus/');
datapath = '/home/DATA/TweakDreams/BIDS/';

replace_seed = 0; % Set to 1 to replace RNG, or 0 to use existing one (if present)

%% Parameter Definition
params.evaluationChannels = 1:257; % Exclude non EEG channels
params.recreference       = 257;   % Recording reference (excluded from bad-channel detection)   
params.veogrchannel       = 248;   % Vertical EOG channel (excluded from bad-channel detection)
params.epochStart         = 116;     % Time before alarm in seconds
params.epochEnd           = 0;     % Time after alarm in seconds
params.highPassFilter     = 0.5;   % High-pass filter for ICA (0 = no filter)
params.lowPassFilter      = 45;    % Low-pass filter for ICA (0 = no filter)
params.resampleICA        = 0;     % Resample the data for ICA saves computation time (0 = no resampling)
params.applyPCAtoICA      = 120;   % Dimensionality reduction prior to ICA (0 = no PCA)
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

%% File selection

[filenames, filepath] = uigetfile([datapath,'*trest_eeg.set'],'Select eeglab set file','MultiSelect','on');
if ~iscell(filenames); temp{1} = filenames; clear filenames; filenames = temp; clear temp; end

%% Data preprocessing loop

for nf = 1:length(filenames)

    runtime = datetime;
    
    % Define paths
    filename = filenames{nf}; eeglab; close all;
    disp(['>> Processing ',filename]);
    params.name = filename(1:end-4); % Data name
    dataset = filepath(strfind(filepath,'dataset'):strfind(filepath,'dataset')+11);
    subject = filename(strfind(filename,'sub')+4:strfind(filename,'sub')+8);
    night   = filename(strfind(filename,'ses')+4:strfind(filename,'ses')+6);
    mkdir([datapath,dataset,'/derivatives/'],['sub-',subject]);
    outpath = [datapath,dataset,'/derivatives/sub-',subject,'/'];
    outdata = [filename(1:end-4),'_clean.set'];
    mkdir(outpath,'reports_preproc_trest');
    mkdir([outpath,'reports_preproc_trest'],['sub-',subject,'_ses-',night,'_trest']);
    outproc = [outpath,'reports_preproc_trest/sub-',subject,'_ses-',night,'_trest/'];
    mkdir(filepath,'proc_settings_trest');
    mkdir([filepath,'proc_settings_trest'],['sub-',subject,'_ses-',night,'_trest']);
    outsets = [filepath,'proc_settings_trest/sub-',subject,'_ses-',night,'_trest/'];
    
    % Save RNG seed for replicability
    seedname = ['sub-',subject,'_ses-',night,'_trest_rng_seed.mat'];
    if exist([outsets,seedname],'file')>0 && replace_seed==0
        temp = load([outsets,seedname]);
        warning('The program is using a previously saved RNG seed!');
        rng(temp.seed);
    else
        rng('shuffle');
        seed = rng;
        cdate = datetime;
        save([outsets,seedname],'seed','subject','night','cdate');
    end; clear temp seed cdate seedname;

    % Load and Epoch Data
    clc; disp('Loading data...');
    EEG = pop_loadset(filename, filepath);
    EEG.epoch = NaN(1,EEG.trials); % Add epochs
    for e = 1:EEG.trials
        EEG.epoch(e) = size(EEG.data,2).*(e-1)+1;
        EEG.event(e).type = 's30';
        EEG.event(e).latency = EEG.epoch(e)+(60000-(EEG.srate*2));
        EEG.event(e).epoch = e;
    end; clear e;
    EEG = eeg_checkset(EEG);
    
    % Remove Line noise
    disp('>> Remove Line Noise');
    EEG = zapclean(EEG, params, outproc);
    EEG = eeg_checkset(EEG); % Check data

    % Detect Bad Channels
    disp('>> Bad Channel Detection');
    badchans = findbadchans(EEG, params);
    indelec = sort(unique([badchans.all_marked{:}]));
    goodchannels = setxor(params.evaluationChannels,indelec);
    report_badchannels(badchans, params, EEG.chanlocs, outproc);
    EEG.tp_epoch = EEG.epoch; EEG.epoch = [];
    EEG = eeg_checkset(EEG); % Check data
    
    % Save EOG Channels (L1G-R1G, VEOGR-R1Z)
    disp('>> Storing EOG Channels');
    EOG = eeg_interp(EEG,indelec,'spherical'); % Use spherical spline for interpolation
    h1 = find(matches({EOG.chanlocs.labels},'L1G'));
    h2 = find(matches({EOG.chanlocs.labels},'R1G'));
    v1 = find(matches({EOG.chanlocs.labels},'VEOGR'));
    v2 = find(matches({EOG.chanlocs.labels},'R1Z'));
    eog_data(1,:,:) = EOG.data(h2,:,:)-EOG.data(h1,:,:); % Right - Left
    eog_data(2,:,:) = EOG.data(v2,:,:)-EOG.data(v1,:,:); % Superior - Inferior
    clear EOG v1 v2 h1 h2;
    
    % Run Independent Component Analysis
    disp('>> Independent Component Analysis');
    EEGica = icarun(EEG, goodchannels, params);
    [EEGica, reject, keep, classcomps, badchans.bc_ica] = icreject(EEGica, params);
    if ~isempty(badchans.bc_ica) 
        warning('Bad channels found during ICA-based artifact reduction!');
        warning('I will exclude these channels and try again...');
        clear EEGica indelec goodchannels reject keep classcomps;
        indelec = sort(unique([badchans.all_marked{:},badchans.bc_ica]));
        goodchannels = setxor(params.evaluationChannels, indelec);
        report_badchannels(badchans, params, EEG.chanlocs, outproc);
        EEGica = icarun(EEG, goodchannels, params);
        [EEGica, reject, keep, classcomps, ~] = icreject(EEGica, params);
    end
    if params.useICLoverride == 1
        [reject, keep] = ICLoverride(reject, keep, outsets, outproc);
    end
    report_iccleaningrest_prova(EEGica, reject, keep, classcomps, goodchannels, outproc);
    EEG.icawinv     = EEGica.icawinv;
    EEG.icasphere   = EEGica.icasphere;
    EEG.icaweights  = EEGica.icaweights;
    EEG.icachansind = EEGica.icachansind;
    EEG = eeg_checkset(EEG);  % Check data structure
    clear EEGica; % This structure was needed only for IC-label application
    EEG.reject.gcompreject(reject) = true;
    EEG = pop_subcomp(EEG, reject ,0);
    EEG = eeg_checkset(EEG); % Check data

    % Interpolate bad channels
    disp('>> Bad Channel Interpolation');
    EEG = eeg_interp(EEG,indelec,'spherical'); % Use spherical spline for interpolation
    EEG = eeg_checkset(EEG); % Check data
    
    % Re-add EOG channels
    EEG.data(EEG.nbchan+1:EEG.nbchan+2,:,:) = eog_data;
    EEG.chanlocs(EEG.nbchan+1).labels = 'HEOG';
    EEG.chanlocs(EEG.nbchan+2).labels = 'VEOG';
    EEG = pop_chanedit(EEG, 'settype',{EEG.nbchan+1,'EOG'});
    EEG = pop_chanedit(EEG, 'settype',{EEG.nbchan+2,'EOG'});
    EEG = eeg_checkset(EEG); % Check data

    % Save output
    disp('Saving Preprocessed Data');
    EEG.preproc.myparams = params;
    EEG.preproc.badchans = badchans;
    EEG.preproc.icreject = reject;
    EEG.preproc.icccomps = classcomps;
    EEG.preproc.sruntime = runtime;
    pop_saveset(EEG, 'filename', outdata, 'filepath', outpath, 'version','7.3');
    close all;
    
    clearvars -except nf datapath params filenames filepath replace_seed; close all;
    
end; clear nf;

disp('Program Completed!');