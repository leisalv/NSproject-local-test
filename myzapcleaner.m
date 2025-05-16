function EEG = myzapcleaner(EEG, outpath)
%
% Clean data using zapline-plus and save report / adapted for non-epoched data
% Klug and Kloosterman, Human Brain Mapping 2022

srate = EEG.srate;    
signal = EEG.data';
[cleand, ~, summary] = clean_data_with_zapline_plus(signal, srate, 'noisefreqs',50, 'winSizeCompleteSpectrum', floor(length(signal)/srate/8));

if ~isempty(summary.foundNoise)
    EEG.data = cleand'; % replace data in EEGLAB structure
end
    
set(gcf,'units','normalized','outerposition',[0.01 0.01 .99 .99]);
figname = fullfile(outpath, 'Rep_LineNoise.png'); % Use fullfile for path concatenation
disp(['Figure name: ', figname]); % Debugging statement
try
    exportfig(gcf, figname, 'Format', 'png', 'Color', 'cmyk', 'Resolution', 300, 'Renderer', 'opengl');
catch ME
    disp('Error using exportfig:');
    disp(ME.message);
    disp('Using saveas as an alternative...');
    saveas(gcf, figname, 'png');
end
close all;