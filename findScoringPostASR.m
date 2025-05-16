%% Load Pre-ASR Dataset (with full labels)
myeeg_preasr = pop_loadset();

% Step 1: Get original total length
original_total_length = myeeg_preasr.pnts;

% Step 2: Initialize full-length label array
original_labels = repmat({'undefined'}, 1, original_total_length);
myevents = myeeg_preasr.event;

% Step 3: Fill in original labels from events
for i = 1:length(myevents)
    start_idx = round(myevents(i).latency);
    if i < length(myevents)
        end_idx = round(myevents(i+1).latency) - 1;
    else
        end_idx = original_total_length;
    end
    start_idx = max(start_idx, 1);
    end_idx = min(end_idx, original_total_length);
    original_labels(start_idx:end_idx) = {myevents(i).type};
end

%% Load Cleaned EEG (after ASR/artifact rejection)
myeeg_prep = pop_loadset();

% Step 4: Initialize cleaned_labels from original_labels
cleaned_labels = original_labels;

% Step 5: Get boundary events
cleaned_events = myeeg_prep.event;
is_boundary = strcmp({cleaned_events.type}, 'boundary');
boundary_events = cleaned_events(is_boundary);
[~, sort_idx] = sort([boundary_events.latency]);
boundary_events = boundary_events(sort_idx);

% Step 6: Remove boundary data from labels
total_removed = 0;
for i = 1:length(boundary_events)
    b_start = round(boundary_events(i).latency - total_removed);
    b_end = b_start + boundary_events(i).duration - 1;
    b_start = max(b_start, 1);
    b_end = min(b_end, length(cleaned_labels));
    cleaned_labels(b_start:b_end) = [];
    total_removed = total_removed + boundary_events(i).duration;
end

% Step 7: Sanity check
if length(cleaned_labels) ~= myeeg_prep.pnts
    warning('⚠️ Length mismatch: cleaned_labels (%d) vs myeeg_prep.pnts (%d)', ...
        length(cleaned_labels), myeeg_prep.pnts);
else
    disp('✅ Cleaned labels match EEG length.');
end

%% SEGMENT INTO LONGEST UNINTERRUPTED BOUTS OF CONSISTENT STAGE

% Step 8: Reconstruct original boundary latencies (after cutting)
boundary_events = myeeg_prep.event(strcmp({myeeg_prep.event.type}, 'boundary'));
boundary_latencies = [];
removed_so_far = 0;

for i = 1:length(boundary_events)
    latency = round(boundary_events(i).latency - removed_so_far);
    duration = boundary_events(i).duration;
    boundary_latencies = [boundary_latencies, latency];
    removed_so_far = removed_so_far + duration;
end

% Step 9: Get continuous segments
segment_starts = [1, boundary_latencies + 1];
segment_ends = [boundary_latencies - 1, length(cleaned_labels)];

% Step 10: Initialize segmented event structure
segmented_events = [];
N = length(cleaned_labels);

% Step 11: Loop over segments and find consistent runs
for seg = 1:length(segment_starts)
    s = max(1, segment_starts(seg));
    e = min(N, segment_ends(seg));

    if e < s
        continue;
    end

    segment_labels = cleaned_labels(s:e);
    current_label = segment_labels{1};
    run_start = s;

    for idx = s+1:e
        if idx > N
            break;
        end

        if ~strcmp(cleaned_labels{idx}, current_label)
            % Save event
            event = struct();
            event.type = current_label;
            event.latency = run_start;
            event.duration = idx - run_start;
            segmented_events = [segmented_events, event];

            % Start new run
            current_label = cleaned_labels{idx};
            run_start = idx;
        end
    end

    % Add final run
    if run_start <= e
        event = struct();
        event.type = current_label;
        event.latency = run_start;
        event.duration = e - run_start + 1;
        segmented_events = [segmented_events, event];
    end
end

%% OPTIONAL: Remove segments shorter than X seconds
fs = myeeg_prep.srate;  % sampling rate
min_duration = 3 * fs;  % 3 seconds
segmented_events2 = segmented_events([segmented_events.duration] >= min_duration);

%% DONE 🎉

% You can now:
% - Add these to EEG.event if desired
% - Save to .mat for later
% - Use to extract epochs or for analysis

disp(['✅ Found ', num2str(length(segmented_events)), ' clean uniform segments.']);
%%%%%%%



mycleanevents = myeeg_prep.event;
for i=1:length(mycleanevents)
    mycleanevents(i).latency = ceil(mycleanevents(i).latency);
    mycleanevents(i).type = cleaned_labels(mycleanevents(i).latency);
end



%%%%%%
%% V2
% Initialize total_removed
total_removed = 0;

% Initialize an empty array to hold the new events
mycleanevents = struct('latency', {}, 'type', {}, 'duration', {});

% Step 4: Loop through each boundary sequentially
for i = 1:length(boundary_events)
    % Get latency and duration
    latency = boundary_events(i).latency;
    duration = boundary_events(i).duration;

    % Adjust latency for previously removed samples
    adjusted_start = ceil(latency - total_removed);
    adjusted_end = adjusted_start + duration - 1;

    % Ensure indices stay within bounds (if desired)
    % adjusted_start = max(1, adjusted_start);
    % adjusted_end = min(length(cleaned_labels), adjusted_end);

    % Remove the specified segment
    cleaned_labels(adjusted_start:adjusted_end) = [];

    % Update how many total points have been removed
    total_removed = total_removed + (adjusted_end - adjusted_start + 1);

    % Fill mycleanevents structure with the adjusted information
    mycleanevents(i).latency = adjusted_start;
    mycleanevents(i).type = cleaned_labels(adjusted_start); % 'NA' or other label
    mycleanevents(i).duration = adjusted_end - adjusted_start + 1;
end

% Display the new structure to check the output
disp(mycleanevents);


%%%% last version
% 1) Load and sort all events
all_events    = myeeg_prep.event;
[~, orderAll] = sort([all_events.latency]);
all_events    = all_events(orderAll);  % chronological

% 2) Identify boundary vs real
is_boundary  = strcmp({all_events.type}, 'boundary');
boundary_idx = find(is_boundary);

% 3) Precompute cumulative removal up to each boundary
b_lats = [all_events(boundary_idx).latency];
b_durs   = [all_events(boundary_idx).duration];
cum_durs = cumsum(b_durs);

% 4) Pre-allocate the output struct array
mycleanevents = struct('type', {}, 'latency', {}, 'duration', {});

removed_so_far = 0;

% 5) Walk through every event
for k = 1:numel(all_events)
    ev = all_events(k);
    
    if strcmp(ev.type, 'boundary')
        % --- a CUT marker ---
        % Need a next event to know how far to recover
        if k == numel(all_events)
            break; % no “next event” after the final boundary
        end
        
        nextEv = all_events(k+1);
        
        % Stage just *before* the cut, from original_labels
        idx_before = max(1, ceil(ev.latency)-1);
        recovered_stage = original_labels{idx_before};
        
        % Build the recovered event
        rec.type     = recovered_stage;
        rec.latency  = ceil(ev.latency);                  
        rec.duration = nextEv.latency - ev.latency; 
        
        % Append
        mycleanevents(end+1) = rec;  %#ok<SAGROW>
        
        % Update how many points have been removed so far
        removed_so_far = removed_so_far + ev.duration;
        
    else
        % --- a real Wake/NREM event ---
        shiftedL = ev.latency - removed_so_far;
        
        ne.type     = ev.type;
        ne.latency  = shiftedL;
        ne.duration = ev.duration;
        
        mycleanevents(end+1) = ne;
    end
end

% 6) Finally, sort by the *new* latencies
[~, ord2]      = sort([mycleanevents.latency]);
mycleanevents = mycleanevents(ord2);

% Display how many events you have now
disp(['Reconstructed events: ', num2str(numel(mycleanevents))]);




%% Really the last: 
% Step 1: Sort events chronologically
all_events = myeeg_prep.event;
[~, sort_idx] = sort([all_events.latency]);
all_events = all_events(sort_idx);

% Step 2: Initialize variables
reconstructed = struct('type', {}, 'latency', {}, 'duration', {});
removed_so_far = 0;
prev_type = '';  % to track stage transitions
reco_idx = 0;    % index for reconstructed

% Step 3: Loop through all events
for k = 1:length(all_events)
    ev = all_events(k);

    if strcmp(ev.type, 'boundary')
        % 3a. Recover the stage before the boundary
        if k < length(all_events)
            next_ev = all_events(k+1);

            % Use the original_labels just before the cut
            prev_idx = max(1, floor(ev.latency) - 1);
            stage = original_labels{prev_idx};

            % Define the recovered event duration
            duration = next_ev.latency - ev.latency;

            % Merge if same as previous
            if strcmp(stage, prev_type)
                % Extend previous event
                reconstructed(reco_idx).duration = ...
                    reconstructed(reco_idx).duration + duration;
            else
                % New recovered event
                reco_idx = reco_idx + 1;
                reconstructed(reco_idx).type = stage;
                reconstructed(reco_idx).latency = ev.latency - removed_so_far;
                reconstructed(reco_idx).duration = duration;
                prev_type = stage;
            end
        end

        % Update total removed
        removed_so_far = removed_so_far + ev.duration;

    else
        % 3b. Process a real event
        curr_type = ev.type;
        curr_lat = ev.latency - removed_so_far;
        curr_dur = ev.duration;

        % Merge with previous if same type
        if strcmp(curr_type, prev_type)
            reconstructed(reco_idx).duration = ...
                reconstructed(reco_idx).duration + curr_dur;
        else
            reco_idx = reco_idx + 1;
            reconstructed(reco_idx).type = curr_type;
            reconstructed(reco_idx).latency = curr_lat;
            reconstructed(reco_idx).duration = curr_dur;
            prev_type = curr_type;
        end
    end
end

% Optional: sort again by latency (if needed)
[~, ord] = sort([reconstructed.latency]);
mycleanevents = reconstructed(ord);

% Display result
disp(['Total reconstructed events: ', num2str(length(mycleanevents))]);
