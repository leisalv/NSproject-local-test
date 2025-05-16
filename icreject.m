function [EEG, reject, keep, classcomps, bc_ica] = icreject(EEG, params)
%% [EEG, reject, keep, classcomps] = icreject(EEG, params)
%
% Remove bad components
%

gchannels = EEG.icachansind;
thr_brain = params.brainThrICA;
thr_ratio = params.brainNrtICA;

EEG         = iclabel(EEG); % Classify components using iclabel
classcomps = EEG.etc.ic_classification.ICLabel.classifications;

temp_rej1 = find(sum(classcomps(:,1) < classcomps(:,2:end),2)>1); % Components where brain is first or second
temp_rej2 = find(sum(classcomps(:,1) < classcomps(:,2:end),2)==1 & (classcomps(:,1)./max(classcomps(:,2:end),[],2)) < thr_ratio); 
temp_rej3 = find(classcomps(:,1) < thr_brain); % Brain-class tresh

reject    = sort(unique([temp_rej1;temp_rej2;temp_rej3]));
keep      = setxor(1:size(classcomps,1),reject);

EEG.reject.gcompreject(reject) = true;

    
% Find bad channels based on ICA

bc_ica = [];
if params.badchanthrICA > 0
    
    idx = find(classcomps(1:params.badchanmaxICs,6)>params.badchanthrICA); % Look for bad channels
    components = EEG.icawinv;

    if ~isempty (idx)
        for n = 1:length(idx)
            ic = abs(components(:,idx(n)))';
%             figure; topoplot(ic,EEG.chanlocs(EEG.icachansind),'electrodes','numbers');
            thr = max(ic).*0.99;
            bc_ica = [bc_ica,gchannels(ic>=thr)];  %#ok<AGROW>
        end; clear n;
        bc_ica = sort(unique(bc_ica));
    end
    
end

