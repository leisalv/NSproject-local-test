%% Awareness Cessation Project
% Mar 2025, Leila Salvesen
% adapted from EEG03_DataProcessing / Script_Prep_Austria


%%%%%% STEP 2: CESSATION SESSION (NAP + WAKE RS + NS) DATA PROCESSING SCRIPT

%% Initialization

clear variables; close all; clc;
restoredefaultpath;
[codepath,codefile,codeext] = fileparts(matlab.desktop.editor.getActiveFilename);

addpath('C:\Users\Leila\Local\MATLAB\toolboxes\eeglab2025.0.0\'); eeglab;
addpath('C:\Users\Leila\Local\MATLAB\toolboxes\fieldtrip-20240110\'); ft_defaults;
addpath(genpath(codepath))

datapath = 'C:\Users\Leila\Local\NSproject\data\';
derivpath = 'C:\Users\Leila\Local\NSproject\data\derivatives\';
mydate = datetime('now');
myoutpath = [derivpath, codefile, '_', datestr(mydate, 'yyyy-mm-dd_HH-MM-SS'), '\'];
mkdir(myoutpath);

% Define frequency ranges
freq_bands = [0.5 2; 2 4; 4 8; 8 12; 12 16; 16 25; 25 45; 45 80];
slope_ranges = [0.5 30; 30 45; 1 20; 20 40; 1 45];

%% Load data

clc; disp('Loading data...');

% Load nap session data
myeeg_prep = pop_loadset([derivpath, 'NS_nap_preproc_v3\prep_data\DA_fullNAP_PREP.set']);
myeeg_prep = eeg_checkset(myeeg_prep);
eeglab redraw;

%% Define sleep stages segments within nap session 
% (N1 vs Wake)

% Identify segments for both stages
n1_bouts = find_stage_segments(myeeg_prep.preproc.stages, 1); %totalSumN1 = sum(n1_bouts(:,2) - n1_bouts(:,1));
wake_bouts = find_stage_segments(myeeg_prep.preproc.stages, 0); %totalSumWk = sum(wake_bouts(:,2) - wake_bouts(:,1));

% Extract N1 and Wake segments
EEG_n1_bouts = extract_eeg_segments(myeeg_prep, n1_bouts);
pow_n1_bouts = NaN(length(EEG_n1_bouts), length(freq_bands));
slp_n1_bouts = NaN(length(EEG_n1_bouts), length(slope_ranges));

EEG_wake_bouts = extract_eeg_segments(myeeg_prep, wake_bouts);
pow_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(freq_bands));
slp_wake_bouts = NaN(length(EEG_wake_bouts), myeeg_prep.nbchan, length(slope_ranges));

for n1b = 1:length(EEG_n1_bouts)
    fprintf('Loading bout %d / %d\n', n1b, length(EEG_n1_bouts));
    [orig_bout, osc_bout, slp_bout, itc_bout] = analyze_eeg_segment(EEG_n1_bouts{n1b}, freq_bands, slope_ranges);
    pow_n1_bouts(n1b,:,:) = squeeze(mean(osc_bout, 1)); % average across 5s epochs
    slp_n1_bouts(n1b,:,:) = squeeze(mean(slp_bout, 1));
end

for wkb = 1:length(EEG_wake_bouts)
    fprintf('Loading bout %d / %d\n', wkb, length(EEG_wake_bouts));
    [orig_bout, osc_bout, slp_bout, itc_bout] = analyze_eeg_segment(EEG_wake_bouts{wkb}, freq_bands, slope_ranges);
    pow_wake_bouts(wkb,:,:) = squeeze(mean(osc_bout, 1)); % average across 5s epochs
    slp_wake_bouts(wkb,:,:) = squeeze(mean(slp_bout, 1));
end

%SAVE STRUCTURES!!!!!!!


% Define output path
save(fullfile(myoutpath, 'nap_analyses_results.mat'), ...
     'pow_n1_bouts', 'slp_n1_bouts', 'pow_wake_bouts', 'slp_wake_bouts','-v7.3');


disp('Analysis results saved successfully.');


%%% TRYING OUT TOPOPLOTS


%% Plot Power in Frequency Bands for N1 Bouts

%% Remove last row if NaNs
if all(all(isnan(pow_n1_bouts(end,:,:))))
    pow_n1_bouts = pow_n1_bouts(1:end-1,:,:);
end


%% Compute mean across bouts for N1 (channels x frequency bands)
mean_pow_n1 = squeeze(mean(pow_n1_bouts, 1)); % 64 chans x 8 bands

%% EEG bands (excluding high gamma)
freq_labels = {'0.5-2 Hz (Delta)', '2-4 Hz (Delta)', '4-8 Hz (Theta)', ...
               '8-12 Hz (Alpha)', '12-16 Hz (Sigma)', '16-25 Hz (Beta)', ...
               '25-45 Hz (Gamma)'};

mean_pow_n1 = mean_pow_n1(:,1:7); % remove high gamma (8th band)

%% Configuration for topoplots
cfg = [];
cfg.layout    = 'easycapM1.mat'; % confirm correct layout
cfg.marker    = 'on';
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
cfg.gridscale = 200;
cfg.colormap  = jet;

% Common color scale across all plots
cfg.zlim = [min(mean_pow_n1(:)), max(mean_pow_n1(:))];

%% Generate individual hidden plots first
nBands = size(mean_pow_n1,2);
plot_handles = gobjects(1, nBands);

for b = 1:nBands
    data = [];
    data.label = {myeeg_prep.chanlocs.labels};
    data.dimord = 'chan_time';
    data.time = 0;
    data.avg = mean_pow_n1(:, b);
    data.avg(isnan(data.avg)) = 0;

    plot_handles(b) = figure('visible','off'); % hidden figures
    ft_topoplotER(cfg, data);
    title(freq_labels{b}, 'FontSize',12,'FontWeight','bold');
end

%% Combine into one figure clearly
figure('Name','EEG Power (N1 Sleep)','Position',[50 50 1600 700]);

for b = 1:nBands
    subplot(2,4,b); % 2 rows, 4 columns, last subplot empty

    % Copy hidden plot axes
    ax_old = findobj(plot_handles(b),'type','axes');
    ax_new = copyobj(ax_old,gcf);
    set(ax_new,'Position',get(gca,'Position'));
    delete(gca); % remove empty axes
    close(plot_handles(b)); % close hidden plot

    title(freq_labels{b},'FontSize',12,'FontWeight','bold');
end

% Remove the empty subplot (8th plot position)
subplot(2,4,8); axis off;

sgtitle('EEG Power Topoplots during N1 Sleep','FontSize',16,'FontWeight','bold');

disp('✅ EEG topoplots successfully generated (Delta to Gamma).');


%% Plot Slope for N1 Bouts
% Similarly, check for slopes if needed

if all(all(isnan(slp_n1_bouts(end,:,:))))
    slp_n1_bouts = slp_n1_bouts(1:end-1, :, :);
    disp('Last row with NaNs removed from slope data.');
end

for s = 1:size(slope_ranges,1)
    figure;
    cfg.parameter = 'avg';  
    cfg.zlim = 'maxabs';  

    slope_label = sprintf('%.1f-%.1f Hz', slope_ranges(s,1), slope_ranges(s,2));
    cfg.title = ['N1 Slope - ', slope_label];

    data = [];
    data.label = {myeeg_prep.chanlocs.labels};
    data.dimord = 'chan_time';
    data.time = 0;
    data.avg = squeeze(mean(slp_n1_bouts(:,:,s), 1))';

    ft_topoplotER(cfg, data);
    title(['Slope in ', slope_label, ' during N1']);
end


disp('Topographical plots generated successfully.');