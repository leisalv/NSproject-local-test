function [EEG, reject, keep, classcomps, bc_ica] = myicreject(EEG, params, outpath)
%% [EEG, reject, keep, classcomps] = icreject(EEG, params)
%
% Remove bad components
%

gchannels = EEG.icachansind;
thr_brain = params.brainThrICA;
thr_ratio = params.brainNrtICA;

EEG        = iclabel(EEG); % Classify components using iclabel
classcomps = EEG.etc.ic_classification.ICLabel.classifications;

temp_rej1 = find(sum(classcomps(:,1) < classcomps(:,2:end),2)>1); % Components where brain is first or second (excluding 'other')
temp_rej2 = find(sum(classcomps(:,1) < classcomps(:,2:end),2)==1 & (classcomps(:,1)./max(classcomps(:,2:end),[],2)) < thr_ratio); % Components where brain is second but ratio with first is high (excluding 'other')
temp_rej3 = find(classcomps(:,1) < thr_brain); % Brain-class treshold

reject    = sort(unique([temp_rej1;temp_rej2;temp_rej3]));
keep      = setxor(1:size(classcomps,1),reject);

EEG.reject.gcompreject(reject) = true;
    
% Find bad channels based on ICA

bc_ica = [];
if params.badchanthrICA > 0
    
    idx = find(classcomps(1:(params.badchanmaxICs * size(classcomps,1)),6) > params.badchanthrICA); % Look for bad channels
    components = EEG.icawinv;

    if ~isempty(idx)
        figure('Visible', 'off'); % create invisible figure (use 'on' for debugging)
        nICs = length(idx);
        for n = 1:nICs
            subplot(1, nICs, n);
            ic = abs(components(:, idx(n)))';
            topoplot(ic, EEG.chanlocs(EEG.icachansind), 'electrodes', 'numbers');
            title(['IC ', num2str(idx(n))]);
            thr = max(ic) * 0.99;
            bc_ica = [bc_ica, gchannels(ic >= thr)];  %#ok<AGROW>
        end
        bc_ica = sort(unique(bc_ica));

        % Save figure
        figname = fullfile(outpath, 'Rep_BadChansICA_Topo.png');
        saveas(gcf, figname);
        close(gcf); % close after saving
    end
    
end

