function myeegstg = scorecsv2mat(scorefile, eegfile)
    
    % SCORECSV2MAT Converts sleep scoring CSV directly to sample-wise hypnogram.
    %
    % Inputs:
    %   scorefile      - table from CSV containing sleep scoring epochs
    %   eegfile        - EEG file
    %
    % Output:
    %   myeegstg       - hypnogram (numeric stage labels) in timepoints
    
    % Load scoring information
    scoring = table2cell(scorefile);
    
    % Extract epoch onset and offset (samples)
    ep_str = ([scoring{:,2}] .* eegfile.srate) + 1;  
    ep_end = ([scoring{:,3}] .* eegfile.srate);
    
    % Define stage mapping
    stage_map = containers.Map(...
        {'Wake','NREM1','NREM2','NREM3','REM'}, ...
        [0, 1, 2, 3, 5]);
    
    % Initialize hypnogram vector
    myeegstg = zeros(1, eegfile.pnts);
    
    % Assign numeric stages directly
    for ep = 1:length(ep_str)
        label = scoring{ep,4};
        if stage_map.isKey(label)
            myeegstg(ep_str(ep):min(ep_end(ep), eegfile.pnts)) = stage_map(label);
        else
            warning('Unknown stage label at epoch %d: "%s"', ep, label);
        end
    end

end
