clear all, close all, clc

%% INITIALIZATION & SETTINGS
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; 

filepath = 'C:\Users\cwpan\Desktop\EEG-ipsilateral position\data';

% Define the subjects into their respective groups
control_IDs = {'08','09','11','13','23','25','29','32','33','37','39','40'};
pdcd_IDs    = {'10','16','17','18','19','20','21','26','27','28','30','36'};

groups = {'Control', 'pDCD'};
all_IDs = {control_IDs, pdcd_IDs};

% Store the final ERP structures
ERP_results = struct('Control', {cell(1, length(control_IDs))}, ...
                     'pDCD', {cell(1, length(pdcd_IDs))});

% Keep track of if a subject was plotted or skipped, if applicable
status_log = struct('Control', {cell(1, length(control_IDs))}, ...
                    'pDCD', {cell(1, length(pdcd_IDs))});


%% AUTOMATED BATCH PROCESSING LOOP WITH AUTO-BINS
for g = 1:2 % two groups (control and pDCD)
    current_group_IDs = all_IDs{g};
    current_group_name = groups{g};
    
    for s = 1:length(current_group_IDs)
        current_id = current_group_IDs{s};
        filename = sprintf('SERP%s_position_45-52_單_proc_convert.cdt.cnt', current_id);
        full_file_path = fullfile(filepath, filename);
        
        fprintf('\n==========================================\n');
        fprintf('Processing %s Subject: SERP%s\n', current_group_name, current_id);
        fprintf('==========================================\n');
        
        % Check 1: Does the file exist?
        if ~exist(full_file_path, 'file')
            warning('FILE NOT FOUND.');
            status_log.(current_group_name){s} = 'Skipped: File Not Found';
            continue; 
        end
        
        % Check 2: Is the file locked/corrupted?
        test_fid = fopen(full_file_path, 'r');
        if test_fid == -1
            warning('FILE LOCKED OR CORRUPTED.');
            status_log.(current_group_name){s} = 'Skipped: File Locked or Corrupted';
            continue;
        end
        fclose(test_fid); 
        
        try
            % A. Load Raw Data
            EEG = pop_loadcnt(full_file_path, 'dataformat', 'auto', 'memmapfile', '');
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
            
            % --- AUTOMATED-BIN DETECTOR (Sorted by Numerical Value) ---
            % Extract triggers and convert to numbers
            ev_types = {EEG.event.type};
            ev_types(strcmpi(ev_types, 'boundary')) = []; % Ignore boundaries
            if ischar(ev_types{1}) || isstring(ev_types{1})
                ev_nums = cellfun(@str2double, ev_types);
            else
                ev_nums = cell2mat(ev_types);
            end
            
            % 1. Find the top 4 most frequent triggers
            [unique_trigs, ~, idx] = unique(ev_nums);
            counts = histcounts(idx, 1:length(unique_trigs)+1);
            [~, freq_sort_idx] = sort(counts, 'descend');
            
            % Isolate the 4 most frequent triggers
            top_4_frequent = unique_trigs(freq_sort_idx(1:min(4, length(unique_trigs))));
            
            % 2. Sort the triggers numerically (smallest to largest)
            top_4_trigs = sort(top_4_frequent, 'ascend');
            
            % Now Bin 1 is the lowest trigger number, Bin 4 is the highest trigger number!
            
            % Create a custom BDF file just for this subject
            bin_file_path = fullfile(filepath, sprintf('temp_bins_%s.txt', current_id));
            fid = fopen(bin_file_path, 'w');
            for b = 1:length(top_4_trigs)
                fprintf(fid, 'bin %d\nTrigger_%d\n.{%d}\n\n', b, top_4_trigs(b), top_4_trigs(b));
            end
            fclose(fid);
            % -------------------------------------
            
            % B. Create EventList
            EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % C. Assign Bins using their custom BDF file
            EEG = pop_binlister(EEG, 'BDF', bin_file_path, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG');
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % Clean up the temporary text file
            delete(bin_file_path);
            
            % D. Extract Epochs & Baseline Correct
            EEG = pop_epochbin(EEG, [-500.0  1000.0], [-200 0]);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % E. Artifact Rejection
            EEG = pop_artextval(EEG, 'Channel', 1:34, 'Flag', 1, 'LowPass', -1, 'Threshold', [-100 100], 'Twindow', [-500 1000]);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % F. Compute Averaged ERP
            ERP = pop_averager(EEG, 'Criterion', 'good', 'DQ_flag', 0);
            ERP_results.(current_group_name){s} = ERP;
            
            % Log their survival status
            accepted_trials = ERP.ntrials.accepted(1);
            if accepted_trials == 0
                status_log.(current_group_name){s} = 'Skipped: 100% Trials Rejected (Extreme Artifacts)';
            else
                status_log.(current_group_name){s} = sprintf('Plotted (Survived with %d valid trials for Bin 1)', accepted_trials);
            end
            
        catch ME
            status_log.(current_group_name){s} = ['Skipped: Pipeline Error -> ' ME.message];
        end
    end
end

%% === PRINT STATUS REPORT TO COMMAND WINDOW ===
fprintf('\n\n=======================================================\n');
fprintf('                FINAL SUBJECT STATUS REPORT              \n');
fprintf('=======================================================\n');
for g = 1:2
    fprintf('--- %s GROUP ---\n', upper(groups{g}));
    for s = 1:length(all_IDs{g})
        fprintf('SERP%s : %s\n', all_IDs{g}{s}, status_log.(groups{g}){s});
    end
    fprintf('\n');
end
fprintf('=======================================================\n\n');

%% Ensure we have at least one valid ERP to extract time and channel maps
valid_erp = [];
for g = 1:2
    for s = 1:length(ERP_results.(groups{g}))
        if ~isempty(ERP_results.(groups{g}){s}) && ERP_results.(groups{g}){s}.ntrials.accepted(1) > 0
            valid_erp = ERP_results.(groups{g}){s};
            break;
        end
    end
    if ~isempty(valid_erp), break; end
end

if isempty(valid_erp)
    error('NO VALID DATA WAS PROCESSED. Check the status report above to see why all subjects failed.');
end

time_ms = valid_erp.times;
target_chans = {'Fz', 'FCz', 'Cz'};

%% PLOT 1: ALL INDIVIDUAL SUBJECTS
bin_to_plot = 3;


figure('Color', 'w', 'Position', [100, 100, 600, 800], 'Name', 'All Individual Subjects'); 

for i = 1:length(target_chans)
    subplot(3, 1, i); hold on;
    ch_idx = find(strcmpi({valid_erp.chanlocs.labels}, target_chans{i}));
    
    if ~isempty(ch_idx)
        h1 = []; h2 = [];
        % Plot valid Control subjects (Light Blue)
        for s = 1:length(control_IDs)
            if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(bin_to_plot) > 0
                h1 = plot(time_ms, squeeze(ERP_results.Control{s}.bindata(ch_idx, :, bin_to_plot)), ...
                    'Color', [0.6 0.8 1.0], 'LineWidth', 1);
            end
        end
        % Plot valid pDCD subjects (Light Red)
        for s = 1:length(pdcd_IDs)
            if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(bin_to_plot) > 0
                h2 = plot(time_ms, squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, bin_to_plot)), ...
                    'Color', [1.0 0.6 0.6], 'LineWidth', 1);
            end
        end
        
        xlim([-200 800]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); 
        line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1); 
        
        ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
        title(sprintf('Electrode: %s (Individual Subjects)', target_chans{i}), 'Color', 'k', 'FontWeight', 'bold');
        
        if i == 3 && ~isempty(h1) && ~isempty(h2)
            xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
            legend([h1, h2], {'Control Individuals', 'pDCD Individuals'}, 'Location', 'southwest');
        end
    end
    hold off;
end

%% PLOT 2: AVERAGE OF BOTH GROUPS (DCD and control)
fprintf('Generating Plot 2 (Group Averages)...\n');
figure('Color', 'w', 'Position', [750, 100, 600, 800], 'Name', 'Grand Average ERPs');

for i = 1:length(target_chans)
    subplot(3, 1, i); hold on;
    ch_idx = find(strcmpi({valid_erp.chanlocs.labels}, target_chans{i}));
    
    if ~isempty(ch_idx)
        valid_control_data = [];
        for s = 1:length(control_IDs)
            if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(bin_to_plot) > 0
                valid_control_data(end+1, :) = squeeze(ERP_results.Control{s}.bindata(ch_idx, :, bin_to_plot));
            end
        end
        mean_control = mean(valid_control_data, 1);
        
        valid_pdcd_data = [];
        for s = 1:length(pdcd_IDs)
            if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(bin_to_plot) > 0
                valid_pdcd_data(end+1, :) = squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, bin_to_plot));
            end
        end
        mean_pdcd = mean(valid_pdcd_data, 1);
        
        plot(time_ms, mean_control, 'b', 'LineWidth', 2.5);
        plot(time_ms, mean_pdcd, 'r', 'LineWidth', 2.5);
        
        xlim([-200 800]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); 
        line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1); 
        
        ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
        title(sprintf('Electrode: %s (Group Grand Average)', target_chans{i}), 'Color', 'k', 'FontWeight', 'bold');
        
        if i == 3
            xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
            legend({'Control Group Mean', 'pDCD Group Mean'}, 'Location', 'southwest');
        end
    end
    hold off;
end
fprintf('Done! Batch processing and plotting finished.\n');

%% FIGURES: 2x2 SUBPLOTS FOR EACH CHANNEL (Fz, FCz, Cz) and (ALL 4 TRIGGERS/BINS)
fprintf('\nGenerating Figures (2x2 Bin comparisons by Channel)...\n');

target_chans = {'Fz', 'FCz', 'Cz'};
% Give the bins descriptive names for the titles
bin_names = {'Bin 1 (Passive down)', 'Bin 2 (Passive up)', 'Bin 3 (Active down)', 'Bin 4 (Passive up)'}; 

for c = 1:length(target_chans)
    current_chan = target_chans{c};
    ch_idx = find(strcmpi({valid_erp.chanlocs.labels}, current_chan));

    if isempty(ch_idx)
        continue; % Skip if the channel name is somehow missing
    end

    % Create a new figure for this specific electrode
    fig_name = sprintf('Grand Average - %s (All 4 Conditions)', current_chan);
    figure('Color', 'w', 'Position', [200+(c*50), 150+(c*50), 900, 700], 'Name', fig_name);

    % Loop through all 4 Bins to create the 2x2 grid
    for b = 1:4 
        subplot(2, 2, b);
        hold on;

        % --- Calculate Mean for Valid Control Subjects ---
        valid_control_data = [];
        for s = 1:length(control_IDs)
            % Ensure the subject exists AND has valid trials for this specific bin
            if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(b) > 0
                valid_control_data(end+1, :) = squeeze(ERP_results.Control{s}.bindata(ch_idx, :, b));
            end
        end
        if ~isempty(valid_control_data)
            mean_control = mean(valid_control_data, 1);
            plot(time_ms, mean_control, 'b', 'LineWidth', 2.5);
        end

        % --- Calculate Mean for Valid pDCD Subjects ---
        valid_pdcd_data = [];
        for s = 1:length(pdcd_IDs)
            if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(b) > 0
                valid_pdcd_data(end+1, :) = squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, b));
            end
        end
        if ~isempty(valid_pdcd_data)
            mean_pdcd = mean(valid_pdcd_data, 1);
            plot(time_ms, mean_pdcd, 'r', 'LineWidth', 2.5);
        end

        % Formatting the individual subplot
        xlim([-200 800]); grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5); % Time Zero
        line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1.5); % Zero Voltage

        title(sprintf('%s', bin_names{b}), 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 12);

        % Only add X-labels to the bottom row to keep it visually clean
        if b == 3 || b == 4
            xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
        end
        % Only add Y-labels to the left column
        if b == 1 || b == 3
            ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
        end

        % Put the legend only in the first subplot so it doesn't block data in the others
        if b == 1
            legend({'Control Mean', 'pDCD Mean'}, 'Location', 'southwest', 'FontSize', 10);
        end

        hold off;
    end

    % Add a large Master Title at the top of the entire Figure
    sgtitle(sprintf('Electrode: %s - Group Comparison Across All 4 Triggers', current_chan), ...
        'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
end

fprintf('Done plotting all figures!\n');


%% 7. COMPREHENSIVE AUTOMATED STATISTICAL ANALYSIS (BINS 1 & 3)
fprintf('\n=======================================================\n');
fprintf('       EXTRACTING PEAKS AND RUNNING STATS FOR BINS       \n');
fprintf('=======================================================\n');

target_chans = {'Fz', 'FCz', 'Cz'};
t_window_N1 = [100, 200];
t_window_P3 = [300, 600];
bins_to_analyze = [1, 2, 3, 4]; 

for b = 1:length(bins_to_analyze)
    bin_idx = bins_to_analyze(b);
    fprintf('\n\n#######################################################\n');
    fprintf('                 STARTING ANALYSIS: BIN %d               \n', bin_idx);
    fprintf('#######################################################\n');
    
    % Prepare empty arrays
    SubjectList = {}; GroupList = {};
    N1_amp = []; N1_lat = [];
    P3_amp = []; P3_lat = [];
    
    % --- STEP 1: EXTRACT THE PEAKS ---
    for g = 1:2 
        current_group = groups{g};
        for s = 1:length(ERP_results.(current_group))
            erp = ERP_results.(current_group){s};
            
            if ~isempty(erp) && erp.ntrials.accepted(bin_idx) > 0
                time_ms = erp.times;
                idx_N1 = find(time_ms >= t_window_N1(1) & time_ms <= t_window_N1(2));
                idx_P3 = find(time_ms >= t_window_P3(1) & time_ms <= t_window_P3(2));
                
                temp_N1_amp = zeros(1, 3); temp_N1_lat = zeros(1, 3);
                temp_P3_amp = zeros(1, 3); temp_P3_lat = zeros(1, 3);
                
                for c = 1:length(target_chans)
                    ch_idx = find(strcmpi({erp.chanlocs.labels}, target_chans{c}));
                    wave = squeeze(erp.bindata(ch_idx, :, bin_idx));
                    
                    [min_val, min_loc] = min(wave(idx_N1));
                    temp_N1_amp(c) = min_val; temp_N1_lat(c) = time_ms(idx_N1(min_loc)); 
                    
                    [max_val, max_loc] = max(wave(idx_P3));
                    temp_P3_amp(c) = max_val; temp_P3_lat(c) = time_ms(idx_P3(max_loc));
                end
                
                SubjectList{end+1} = sprintf('SERP%s', all_IDs{g}{s});
                GroupList{end+1} = current_group;
                N1_amp(end+1, :) = temp_N1_amp; N1_lat(end+1, :) = temp_N1_lat;
                P3_amp(end+1, :) = temp_P3_amp; P3_lat(end+1, :) = temp_P3_lat;
            end
        end
    end
    
    % --- STEP 2: LOOP THROUGH ALL 4 METRICS FOR STATS ---
    GroupList_cat = categorical(GroupList');
    SubjectList = SubjectList';
    within_design = table(categorical({'Fz'; 'FCz'; 'Cz'}), 'VariableNames', {'Electrode'});
    
    metrics = {'N1 Amplitude', 'P3 Amplitude', 'N1 Latency', 'P3 Latency'};
    data_pack = {N1_amp, P3_amp, N1_lat, P3_lat};
    
    for m = 1:4
        curr_name = metrics{m};
        curr_data = data_pack{m};
        
        fprintf('\n=======================================================\n');
        fprintf('  METRIC: %s (BIN %d)\n', upper(curr_name), bin_idx);
        fprintf('=======================================================\n');
        
        % 1. TWO-WAY MIXED ANOVA (Group x Electrode)
        fprintf('--- 1. TWO-WAY MIXED ANOVA (Spatial Interaction) ---\n');
        tbl = table(SubjectList, GroupList_cat, curr_data(:,1), curr_data(:,2), curr_data(:,3), ...
            'VariableNames', {'Subject', 'Group', 'Fz', 'FCz', 'Cz'});
        rm = fitrm(tbl, 'Fz-Cz ~ Group', 'WithinDesign', within_design);
        disp(ranova(rm));
        
        % 2. ONE-WAY ANOVA (Averaged ROI)
        fprintf('--- 2. ONE-WAY ANOVA (Averaged Fronto-Central ROI) ---\n');
        avg_data = mean(curr_data, 2);
        [p_roi, ~] = anova1(avg_data, GroupList_cat, 'off');
        fprintf('Main Effect of Group (Averaged ROI): p = %.4f\n', p_roi);
        
        if p_roi < 0.05
            mean_ctrl = mean(avg_data(GroupList_cat == 'Control'));
            mean_pdcd = mean(avg_data(GroupList_cat == 'pDCD'));
            fprintf('  >> [*] SIGNIFICANT DIFFERENCE FOUND!\n');
            fprintf('  >> Control Mean: %.2f | pDCD Mean: %.2f\n', mean_ctrl, mean_pdcd);
            fprintf('  >> Absolute Difference: %.2f\n', abs(mean_ctrl - mean_pdcd));
        else
            fprintf('  >> No significant group difference in ROI.\n');
        end
        
        % 3. POST-HOC TESTS (Individual Electrodes)
        fprintf('\n--- 3. POST-HOC TESTS (Group Differences per Electrode) ---\n');
        for c = 1:3
            [p_elec, ~] = anova1(curr_data(:, c), GroupList_cat, 'off');
            fprintf('Electrode %s: p = %.4f\n', target_chans{c}, p_elec);
            
            if p_elec < 0.05
                mean_ctrl_e = mean(curr_data(GroupList_cat == 'Control', c));
                mean_pdcd_e = mean(curr_data(GroupList_cat == 'pDCD', c));
                fprintf('  >> [*] SIGNIFICANT at %s! Control: %.2f | pDCD: %.2f | Diff: %.2f\n', ...
                    target_chans{c}, mean_ctrl_e, mean_pdcd_e, abs(mean_ctrl_e - mean_pdcd_e));
            end
        end
        fprintf('-------------------------------------------------------\n');
    end
end
fprintf('\n=======================================================\n');
fprintf('               ALL STATISTICAL TESTS COMPLETE            \n');
fprintf('=======================================================\n');


%% BIN 3 (Active Down) — N1 & P3 amplitude/latency/onset, Group comparison
target_chans = {'Fz','FCz','Cz'};
bin_idx = 3;

t_window_N1 = [0, 130];
t_window_P3 = [50, 350];

SubjectList = {}; GroupList = {};
N1_amp = []; N1_lat = []; P3_amp = []; P3_lat = []; P3_onset = [];
edge_flags = {};

for g = 1:2
    current_group = groups{g};
    for s = 1:length(ERP_results.(current_group))
        erp = ERP_results.(current_group){s};
        if isempty(erp) || erp.ntrials.accepted(bin_idx) == 0
            continue
        end
        time_ms = erp.times;
        idx_N1 = find(time_ms >= t_window_N1(1) & time_ms <= t_window_N1(2));
        idx_P3 = find(time_ms >= t_window_P3(1) & time_ms <= t_window_P3(2));

        temp_N1_amp = nan(1,3); temp_N1_lat = nan(1,3);
        temp_P3_amp = nan(1,3); temp_P3_lat = nan(1,3);
        temp_P3_onset = nan(1,3);

        for c = 1:length(target_chans)
            ch_idx = find(strcmpi({erp.chanlocs.labels}, target_chans{c}));
            if isempty(ch_idx)
                warning('%s missing channel %s — skipped', current_group, target_chans{c});
                continue
            end
            wave = squeeze(erp.bindata(ch_idx, :, bin_idx));

            [min_val, min_loc] = min(wave(idx_N1));
            temp_N1_amp(c) = min_val; temp_N1_lat(c) = time_ms(idx_N1(min_loc));
            if min_loc == 1 || min_loc == length(idx_N1)
                edge_flags{end+1} = sprintf('%s %s N1 hit window edge', current_group, target_chans{c});
            end

            [max_val, max_loc] = max(wave(idx_P3));
            peak_idx_full = idx_P3(max_loc);
            temp_P3_amp(c) = max_val; temp_P3_lat(c) = time_ms(peak_idx_full);
            if max_loc == 1 || max_loc == length(idx_P3)
                edge_flags{end+1} = sprintf('%s %s P3 hit window edge', current_group, target_chans{c});
            end

            % NEW: 50% rise-time onset latency — time the wave first
            % reaches half the peak height, searching backward from the peak
            onset_lat = frac_peak_latency(time_ms, wave, peak_idx_full, idx_P3(1), max_val, 0.5);
            temp_P3_onset(c) = onset_lat;
            if isnan(onset_lat)
                edge_flags{end+1} = sprintf('%s %s P3 onset not found (already above 50%% at window start?)', current_group, target_chans{c});
            end
        end

        SubjectList{end+1} = sprintf('SERP%s', all_IDs{g}{s});
        GroupList{end+1} = current_group;
        N1_amp(end+1,:) = temp_N1_amp; N1_lat(end+1,:) = temp_N1_lat;
        P3_amp(end+1,:) = temp_P3_amp; P3_lat(end+1,:) = temp_P3_lat;
        P3_onset(end+1,:) = temp_P3_onset;
    end
end

if ~isempty(edge_flags)
    fprintf('WARNING: %d issue(s) flagged — double check these:\n', numel(edge_flags));
    fprintf('  %s\n', edge_flags{:});
end

%% Mixed RM-ANOVA (Group x Electrode) per metric
SubjectList = SubjectList(:);
GroupList_cat = categorical(GroupList(:));
within_design = table(categorical({'Fz';'FCz';'Cz'}), 'VariableNames', {'Electrode'});

metrics   = {'N1 Amplitude','N1 Latency','P3 Amplitude','P3 Latency','P3 Onset Latency (50%)'};
data_pack = {N1_amp, N1_lat, P3_amp, P3_lat, P3_onset};
n_metrics = numel(metrics);

p_group_manova = nan(1,n_metrics);
p_group_ttest  = nan(1,n_metrics);
p_int          = nan(1,n_metrics);
rm_models      = cell(1,n_metrics);

for m = 1:n_metrics
    good_rows = ~any(isnan(data_pack{m}), 2);
    if ~all(good_rows)
        fprintf('%s: dropping %d subject(s) with missing values before fitting.\n', ...
            metrics{m}, sum(~good_rows));
    end

    tbl = table(SubjectList(good_rows), GroupList_cat(good_rows), ...
        data_pack{m}(good_rows,1), data_pack{m}(good_rows,2), data_pack{m}(good_rows,3), ...
        'VariableNames', {'Subject','Group','Fz','FCz','Cz'});
    rm = fitrm(tbl, 'Fz-Cz ~ Group', 'WithinDesign', within_design);
    rm_models{m} = rm;

    fprintf('\n=== %s (Bin 3) ===\n', metrics{m});
    disp(mauchly(rm));

    rtbl = ranova(rm);
    disp(rtbl);
    p_int(m) = rtbl.pValue(strcmp(rtbl.Properties.RowNames, 'Group:Electrode'));

    mtbl = manova(rm);
    between_str   = strtrim(string(mtbl.Between));
    statistic_str = strtrim(string(mtbl.Statistic));
    group_rows = between_str == "Group" & statistic_str == "Wilks";
    if sum(group_rows) == 1
        p_group_manova(m) = mtbl.pValue(group_rows);
    else
        fprintf('%s: expected 1 match for Group/Wilks, found %d.\n', metrics{m}, sum(group_rows));
        disp(unique(between_str));
    end

    avg_data = mean(data_pack{m}(good_rows,:), 2);
    grp_here = GroupList_cat(good_rows);
    [~, p_t] = ttest2(avg_data(grp_here=='Control'), avg_data(grp_here=='pDCD'));
    p_group_ttest(m) = p_t;
end

%% Holm-Bonferroni correction across the 5 Group main-effect tests
[p_sorted, order] = sort(p_group_manova);
n = numel(p_sorted);
p_holm_sorted = min(cummax(p_sorted .* (n - (1:n) + 1)), 1);
p_holm = nan(1,n);
p_holm(order) = p_holm_sorted;

fprintf('\n--- Group main effect (MANOVA, Wilks), Holm-corrected across all 5 metrics ---\n');
for m = 1:n_metrics
    fprintf('%-25s MANOVA p = %.4f | Holm p = %.4f | (t-test cross-check p = %.4f)\n', ...
        metrics{m}, p_group_manova(m), p_holm(m), p_group_ttest(m));
end

%% Only chase electrode-specific effects where the interaction is real
fprintf('\n--- Group:Electrode interaction check ---\n');
for m = 1:n_metrics
    fprintf('%-25s interaction p = %.4f\n', metrics{m}, p_int(m));
    if p_int(m) < 0.05
        fprintf('  -> significant interaction, running simple effects:\n');
        disp(multcompare(rm_models{m}, 'Group', 'By', 'Electrode'));
    end
end

%% Local function — must stay at the end of the script file
function lat = frac_peak_latency(time_vec, wave_vec, peak_idx, window_start_idx, peak_val, fraction)
% Finds the time at which wave_vec first rises to fraction*peak_val,
% searching backward from the peak toward window_start_idx.
% Returns NaN if no such crossing exists in that range.
    target = fraction * peak_val;
    lat = NaN;
    for k = peak_idx:-1:(window_start_idx + 1)
        if wave_vec(k) >= target && wave_vec(k-1) < target
            t0 = time_vec(k-1); t1 = time_vec(k);
            v0 = wave_vec(k-1); v1 = wave_vec(k);
            lat = t0 + (target - v0) / (v1 - v0) * (t1 - t0);
            return
        end
    end
end