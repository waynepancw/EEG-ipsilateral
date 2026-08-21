%% 1. INITIALIZATION & RAW .CNT LOADING
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % Open EEGLAB background workspace

filepath = 'C:\Users\cwpan\Desktop\EEG-ipsilateral position\';
raw_filename = 'SERP08_position_45-52_單_proc_convert.cdt.cnt'; % Raw Neuroscan recording in .cnt file format

fprintf('Importing raw Neuroscan data: %s...\n', raw_filename);
% Load the raw data format automatically while handling large files safely
EEG = pop_loadcnt(fullfile(filepath, raw_filename), 'dataformat', 'auto', 'memmapfile', '');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

%% 2. CREATE EVENTLIST & GENERATE BINS TEXT ON THE FLY
fprintf('Creating ERPLAB EventList...\n');
EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% --- INTELLIGENT AUTO-BIN DETECTOR ---
% Extract triggers and convert to pure numbers for safe sorting
ev_types = {EEG.event.type};
ev_types(strcmpi(ev_types, 'boundary')) = []; % Ignore boundaries
if ischar(ev_types{1}) || isstring(ev_types{1})
    ev_nums = cellfun(@str2double, ev_types);
else
    ev_nums = cell2mat(ev_types);
end

% Compute appearance frequencies
[unique_trigs, ~, idx] = unique(ev_nums);
counts = histcounts(idx, 1:length(unique_trigs)+1);

% Sort by frequency to find the top 4 real experimental triggers
[~, freq_sort_idx] = sort(counts, 'descend');
top_4_trigs = unique_trigs(freq_sort_idx(1:min(4, length(unique_trigs))));

% Sort mathematically (smallest to largest)
top_4_trigs = sort(top_4_trigs, 'ascend');
num_bins = length(top_4_trigs);

fprintf('\n=== AUTO-DETECTED & SORTED EXPERIMENTAL TRIGGERS ===\n');
detected_triggers = cell(1, num_bins);
for b = 1:num_bins
    detected_triggers{b} = num2str(top_4_trigs(b)); 
    trig_count = counts(unique_trigs == top_4_trigs(b));
    fprintf('Bin %d -> Trigger %s (Found %d in raw data)\n', b, detected_triggers{b}, trig_count);
end
fprintf('====================================================\n\n');

% DYNAMIC BINS FILE CREATION
bin_file_path = fullfile(filepath, 'autobins.txt');
fid = fopen(bin_file_path, 'w');
for b = 1:num_bins
    fprintf(fid, 'bin %d\nTrigger_%s\n.{%s}\n\n', b, detected_triggers{b}, detected_triggers{b});
end
fclose(fid);

% Assign the dynamically discovered triggers to your bins
EEG = pop_binlister(EEG, 'BDF', bin_file_path, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% Clean up text file
delete(bin_file_path);

%% 3. EXTRACT BIN-BASED EPOCHS WITH BASELINE CORRECTION
fprintf('Extracting epochs (-500 to 1000 ms) with custom baseline (-200 to 0 ms)...\n');
EEG = pop_epochbin(EEG, [-500.0  1000.0], [-200 0]);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 1);

%% 4. AUTOMATED ARTIFACT DETECTION (THRESHOLD \pm100 \muV)
fprintf('Scanning for eye blinks and high voltage motor artifacts...\n');
EEG = pop_artextval(EEG, 'Channel', 1:34, 'Flag', 1, 'LowPass', -1, 'Threshold', [-100 100], 'Twindow', [-500 1000]);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% 5. COMPUTE AVERAGED ERP & PRINT REJECTION SUMMARY
fprintf('Averaging clean epochs into ERP structure...\n');
ERP = pop_averager(EEG, 'Criterion', 'good', 'DQ_flag', 0); 

% === NEW: PRINT DETAILED BIN COUNTS ===
fprintf('\n==============================================\n');
fprintf('             FINAL PIPELINE STATUS            \n');
fprintf('==============================================\n');
fprintf('ERP Set Generated: %s\n\n', ERP.erpname);
fprintf('Trials Surviving Artifact Rejection:\n');

for b = 1:num_bins
    accepted = ERP.ntrials.accepted(b);
    rejected = ERP.ntrials.rejected(b);
    total = accepted + rejected;
    
    if total > 0
        survival_rate = (accepted / total) * 100;
    else
        survival_rate = 0;
    end
    
    fprintf('  Bin %d (Trig %s): %2d / %2d valid trials kept (%5.1f%%)\n', ...
        b, detected_triggers{b}, accepted, total, survival_rate);
end
fprintf('==============================================\n\n');

%% 6. GENERATE THE STACKED SUBPLOT GRAPH (Fz, FCz, Cz) - LIGHT MODE
fprintf('Generating final 3-channel publication figure...\n');
figure('Color', 'w', 'Position', [100, 100, 600, 800]); 
time_ms = ERP.times;
target_chans = {'Fz', 'FCz', 'Cz'};

legend_labels = cell(1, num_bins);
for b = 1:num_bins
    legend_labels{b} = sprintf('Bin %d (Trig %s)', b, detected_triggers{b});
end

for i = 1:length(target_chans)
    ch_idx = find(strcmpi({ERP.chanlocs.labels}, target_chans{i}));
    
    if ~isempty(ch_idx)
        subplot(3, 1, i); 
        hold on;
        
        plot(time_ms, squeeze(ERP.bindata(ch_idx, :, 1:num_bins)), 'LineWidth', 1.8);
        
        xlim([-200 800]);
        grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); 
        line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1); 
        
        ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
        title(['Electrode: ' target_chans{i}], 'Color', 'k', 'FontWeight', 'bold');
        
        if i == 3
            xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
            legend(legend_labels, 'Location', 'best', 'TextColor', 'k', 'Color', 'w');
        end
        hold off;
    end
end
fprintf('Done! Single-subject analysis finished.\n');