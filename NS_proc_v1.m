%% Awareness Cessation Project
% Mar 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria


%%%%%% STEP 2: CESSATION SESSION (NAP + WAKE RS + NS) DATA PROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0\'); eeglab;
addpath(genpath(codepath))

datapath = 'C:\Users\Leila\Local\NSproject\data\';
derivpath = 'C:\Users\Leila\Local\NSproject\data\derivatives\';
mydate = datetime('now');
myoutpath = [derivpath, codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'), '\'];
mkdir(myoutpath);

%% Load data

clc; disp('Loading data...');
% Load meditation session data
myeegmed_prep = pop_loadset([derivpath, 'NS_med_preproc_v3\prep_data\DA_meditation_PREP.set']);
myeegmed_prep = eeg_checkset(myeegmed_prep);
eeglab redraw;

% Load nap session data
myeegnap_prep = pop_loadset([derivpath, 'NS_nap_preproc_v3\prep_data\DA_fullNAP_PREP.set']);
myeegnap_prep = eeg_checkset(myeegnap_prep);
eeglab redraw;


%% Define sleep stages tp within nap session 
% (N1 vs Wake)


% Identify segments for both stages without duration constraint
n1_segments_all = find_all_segments(myeegnap_prep.preproc.stages, 1);
wake_segments_all = find_all_segments(myeegnap_prep.preproc.stages, 0);

nrem1_segs = find_continuous_stage_segments(stage_vec, 1, 30, myeegnap_sr);  % NREM 1 - 30 sec minimun
wake_segs  = find_continuous_stage_segments(stage_vec, 0, 30, myeegnap_sr);  % WAKE - 30 sec min
% 
% myeeg_ns = pop_select(myeegns_prep, 'point', medtp);
% myeeg_ns = eeg_checkset(myeeg_ns);
% myeeg_rs = pop_select(myeegns_prep, 'nopoint', medtp);
% myeeg_rs = eeg_checkset(myeeg_rs);

check = [nrem1_segs(:,2) - nrem1_segs(:,1)]
check2 = round(check/256)/30
