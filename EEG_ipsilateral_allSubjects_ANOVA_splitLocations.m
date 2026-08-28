clear all, close all, clc

%% INITIALIZATION & SETTINGS
fprintf('Initializing paths...\n');

% 1. Hardcode the exact data path (No more pop-ups)
filepath = '/home/cwpan/Documents/EEG-ipsilateral/data';

% 2. Add the main EEGLAB folder
eeglab_path = '/home/cwpan/Documents/MATLAB/eeglab2026.0.0';
addpath(eeglab_path);

% 3. Launch EEGLAB silently
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; 

% 4. BRUTE FORCE: Physically scan the hard drive for the ERPLAB folder
plugins_dir = fullfile(eeglab_path, 'plugins');
all_plugins = dir(plugins_dir);
erplab_found = false;

for i = 1:length(all_plugins)
    % Look for any folder containing "erplab" (ignores capitalization/version numbers)
    if all_plugins(i).isdir && contains(lower(all_plugins(i).name), 'erplab')
        erplab_root = fullfile(plugins_dir, all_plugins(i).name);
        addpath(genpath(erplab_root)); % Grab it and all subfolders (like pop_functions)
        fprintf('Force-added ERPLAB and all subfolders from: %s\n', erplab_root);
        erplab_found = true;
        break;
    end
end

if ~erplab_found
    error('Physically cannot find an ERPLAB folder inside %s.', plugins_dir);
end

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

%% PLOT: Cz, Bin 1, All Individual Subjects (Control = Blue theme, pDCD = Red theme)
fprintf('Generating Plot: Cz, Bin 1, All Individual Subjects...\n');

bin_to_plot = 1;
chan_to_plot = 'Cz';
ch_idx = find(strcmpi({valid_erp.chanlocs.labels}, chan_to_plot));

if isempty(ch_idx)
    error('Channel %s not found in channel locations.', chan_to_plot);
end

% N1 peak-picking window: the peak is taken as the most negative sample
% within this latency range (standard peak-picking, e.g. as used by
% ERPLAB's measurement tools). Adjust to match your paradigm's expected
% N1 latency.
n1_window_ms = [0 300];

% Build color gradients so individual subjects are distinguishable
% within their group's theme (Control = blues, pDCD = reds)
n_control = length(control_IDs);
n_pdcd    = length(pdcd_IDs);

control_light = [0.65 0.80 1.00];  % light blue
control_dark  = [0.00 0.10 0.55];  % dark blue
control_cmap = [linspace(control_light(1), control_dark(1), n_control)', ...
                linspace(control_light(2), control_dark(2), n_control)', ...
                linspace(control_light(3), control_dark(3), n_control)'];

pdcd_light = [1.00 0.65 0.65];  % light red
pdcd_dark  = [0.55 0.00 0.00];  % dark red
pdcd_cmap = [linspace(pdcd_light(1), pdcd_dark(1), n_pdcd)', ...
             linspace(pdcd_light(2), pdcd_dark(2), n_pdcd)', ...
             linspace(pdcd_light(3), pdcd_dark(3), n_pdcd)'];

figure('Color', 'w', 'Position', [100, 100, 800, 500], 'Name', 'Cz - Bin 1 - All Individual Subjects');
hold on;

% Plot valid Control subjects (blue theme, one shade per subject)
% Also stash each subject's trace, N1 peak, and color so the grand-average
% plot (below) can reuse them without recomputing
h1 = [];
valid_control_data = [];
n1_control_amp = [];
n1_control_lat = [];
valid_control_colors = [];
for s = 1:n_control
    if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.Control{s}.bindata(ch_idx, :, bin_to_plot))';
        h1(end+1) = plot(time_ms, this_trace, 'Color', control_cmap(s, :), 'LineWidth', 1); %#ok<SAGROW>
        valid_control_data(end+1, :) = this_trace; %#ok<SAGROW>

        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        n1_control_amp(end+1) = n1_amp; %#ok<SAGROW>
        n1_control_lat(end+1) = n1_lat; %#ok<SAGROW>
        valid_control_colors(end+1, :) = control_cmap(s, :); %#ok<SAGROW>
        plot(n1_lat, n1_amp, 'o', 'MarkerEdgeColor', control_cmap(s, :), ...
            'MarkerFaceColor', control_cmap(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
end

% Plot valid pDCD subjects (red theme, one shade per subject)
h2 = [];
valid_pdcd_data = [];
n1_pdcd_amp = [];
n1_pdcd_lat = [];
valid_pdcd_colors = [];
for s = 1:n_pdcd
    if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, bin_to_plot))';
        h2(end+1) = plot(time_ms, this_trace, 'Color', pdcd_cmap(s, :), 'LineWidth', 1); %#ok<SAGROW>
        valid_pdcd_data(end+1, :) = this_trace; %#ok<SAGROW>

        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        n1_pdcd_amp(end+1) = n1_amp; %#ok<SAGROW>
        n1_pdcd_lat(end+1) = n1_lat; %#ok<SAGROW>
        valid_pdcd_colors(end+1, :) = pdcd_cmap(s, :); %#ok<SAGROW>
        plot(n1_lat, n1_amp, 'o', 'MarkerEdgeColor', pdcd_cmap(s, :), ...
            'MarkerFaceColor', pdcd_cmap(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
end

xlim([-200 500]); grid on;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);

xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
title('Electrode: Cz - Bin 1 (Individual Subjects, o = N1 peak)', 'Color', 'k', 'FontWeight', 'bold');

if ~isempty(h1) && ~isempty(h2)
    legend([h1(1), h2(1)], {'Control Individuals', 'pDCD Individuals'}, 'Location', 'southwest');
end

hold off;
fprintf('Done! Plotting finished.\n');

%% PLOT: Cz, Bin 1, Grand Average +/- 1 SD and Range (Control = Blue theme, pDCD = Red theme)
fprintf('Generating Plot: Cz, Bin 1, Grand Average with SD and Range...\n');

figure('Color', 'w', 'Position', [950, 100, 800, 500], 'Name', 'Cz - Bin 1 - Grand Average, SD, Range');
hold on;

leg_handles = [];
leg_labels = {};
x_patch = [time_ms, fliplr(time_ms)];

if ~isempty(valid_control_data)
    mean_control = mean(valid_control_data, 1);
    sd_control   = std(valid_control_data, 0, 1);
    min_control  = min(valid_control_data, [], 1);
    max_control  = max(valid_control_data, [], 1);

    % Range (min-max across subjects) - lightest shading
    p_range_c = fill(x_patch, [max_control, fliplr(min_control)], control_light, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % +/- 1 SD band - medium shading
    p_sd_c = fill(x_patch, [mean_control + sd_control, fliplr(mean_control - sd_control)], control_light, ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');

    % Grand average - solid line
    l_mean_c = plot(time_ms, mean_control, 'Color', control_dark, 'LineWidth', 2.5);

    % Individual N1 peaks (o), each in its subject's shade
    for s = 1:length(n1_control_lat)
        plot(n1_control_lat(s), n1_control_amp(s), 'o', 'MarkerEdgeColor', valid_control_colors(s, :), ...
            'MarkerFaceColor', valid_control_colors(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    % N1 of the grand-average trace itself (x)
    [n1_amp_avg_c, n1_lat_avg_c] = find_n1_peak(mean_control, time_ms, n1_window_ms);
    l_n1_avg_c = plot(n1_lat_avg_c, n1_amp_avg_c, 'x', 'Color', control_dark, 'MarkerSize', 10, 'LineWidth', 2.5);

    leg_handles = [leg_handles, l_mean_c, p_sd_c, p_range_c, l_n1_avg_c];
    leg_labels  = [leg_labels, {'Control Mean', 'Control \pm1 SD', 'Control Range', 'Control Average N1'}];
end

if ~isempty(valid_pdcd_data)
    mean_pdcd = mean(valid_pdcd_data, 1);
    sd_pdcd   = std(valid_pdcd_data, 0, 1);
    min_pdcd  = min(valid_pdcd_data, [], 1);
    max_pdcd  = max(valid_pdcd_data, [], 1);

    % Range (min-max across subjects) - lightest shading
    p_range_p = fill(x_patch, [max_pdcd, fliplr(min_pdcd)], pdcd_light, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % +/- 1 SD band - medium shading
    p_sd_p = fill(x_patch, [mean_pdcd + sd_pdcd, fliplr(mean_pdcd - sd_pdcd)], pdcd_light, ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');

    % Grand average - solid line
    l_mean_p = plot(time_ms, mean_pdcd, 'Color', pdcd_dark, 'LineWidth', 2.5);

    % Individual N1 peaks (o), each in its subject's shade
    for s = 1:length(n1_pdcd_lat)
        plot(n1_pdcd_lat(s), n1_pdcd_amp(s), 'o', 'MarkerEdgeColor', valid_pdcd_colors(s, :), ...
            'MarkerFaceColor', valid_pdcd_colors(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    % N1 of the grand-average trace itself (x)
    [n1_amp_avg_p, n1_lat_avg_p] = find_n1_peak(mean_pdcd, time_ms, n1_window_ms);
    l_n1_avg_p = plot(n1_lat_avg_p, n1_amp_avg_p, 'x', 'Color', pdcd_dark, 'MarkerSize', 10, 'LineWidth', 2.5);

    leg_handles = [leg_handles, l_mean_p, p_sd_p, p_range_p, l_n1_avg_p];
    leg_labels  = [leg_labels, {'pDCD Mean', 'pDCD \pm1 SD', 'pDCD Range', 'pDCD Average N1'}];
end

xlim([-200 500]); grid on;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);

xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
title('Electrode: Cz - Bin 1 (Grand Average \pm1 SD, Range; o = individual N1, x = average N1)', ...
    'Color', 'k', 'FontWeight', 'bold', 'FontSize', 11);

if ~isempty(leg_handles)
    legend(leg_handles, leg_labels, 'Location', 'southwest');
end

hold off;
fprintf('Done! Plotting finished.\n');

%% LOCAL FUNCTIONS
function [peak_amp, peak_lat] = find_n1_peak(trace, t, window)
% FIND_N1_PEAK  Peak-pick the N1 as the most negative sample within a
% latency window (standard peak-picking approach; window is [start end] in ms).
%   trace  - 1 x N amplitude vector
%   t      - 1 x N time vector (ms), same length as trace
%   window - [t_start t_end] latency range to search within (ms)
    idx = find(t >= window(1) & t <= window(2));
    if isempty(idx)
        peak_amp = NaN;
        peak_lat = NaN;
        return;
    end
    [peak_amp, rel_idx] = min(trace(idx));
    peak_lat = t(idx(rel_idx));
end

%% PLOT: Cz, Bin 1, All 12 Subjects per Group (3x4 grid, individual + N1 'o')
fprintf('Generating Plot: Cz, Bin 1, 3x4 Individual Subject Grids...\n');

% --- Control group: 3x4 grid, one subplot per subject ---
figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Name', 'Control - Cz Bin 1 - All Subjects');
for s = 1:n_control
    subplot(3, 4, s); hold on;
    if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.Control{s}.bindata(ch_idx, :, bin_to_plot))';
        plot(time_ms, this_trace, 'Color', control_cmap(s, :), 'LineWidth', 1.2);

        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        plot(n1_lat, n1_amp, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', control_cmap(s, :), 'MarkerSize', 6);
    else
        text(0.5, 0.5, 'No valid data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    end

    xlim([-200 500]); grid on;
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 8);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 0.75);
    title(sprintf('SERP%s', control_IDs{s}), 'FontSize', 9, 'Color', 'k');
    hold off;
end
sgtitle('Control Group - Electrode: Cz - Bin 1 (Individual Subjects, o = N1)', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Color', 'k');

% --- pDCD group: 3x4 grid, one subplot per subject ---
figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Name', 'pDCD - Cz Bin 1 - All Subjects');
for s = 1:n_pdcd
    subplot(3, 4, s); hold on;
    if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, bin_to_plot))';
        plot(time_ms, this_trace, 'Color', pdcd_cmap(s, :), 'LineWidth', 1.2);

        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        plot(n1_lat, n1_amp, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', pdcd_cmap(s, :), 'MarkerSize', 6);
    else
        text(0.5, 0.5, 'No valid data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    end

    xlim([-200 500]); grid on;
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 8);
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
    line(xlim, [0 0], 'Color', 'k', 'LineWidth', 0.75);
    title(sprintf('SERP%s', pdcd_IDs{s}), 'FontSize', 9, 'Color', 'k');
    hold off;
end
sgtitle('pDCD Group - Electrode: Cz - Bin 1 (Individual Subjects, o = N1)', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Color', 'k');

fprintf('Done! 3x4 subject grids plotted.\n');


%% PLOT: Cz, Bin 1, Grand Average After Exclusions (Manual + N1-Outlier) with N1 Peaks
fprintf('Generating Plot: Cz, Bin 1, Grand Average after exclusions...\n');

% --- Manual exclusion: drop SERP08 and SERP11 from the Control group ---
manual_exclude_control = {'08'};

% --- Step 1: gather each remaining subject's trace + N1 (keep the ID alongside) ---
control_ID_list = {}; control_trace_list = []; control_n1amp_list = []; control_n1lat_list = []; control_color_list = [];
for s = 1:n_control
    if ismember(control_IDs{s}, manual_exclude_control)
        continue; % manual drop
    end
    if ~isempty(ERP_results.Control{s}) && ERP_results.Control{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.Control{s}.bindata(ch_idx, :, bin_to_plot))';
        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        control_ID_list{end+1} = control_IDs{s}; %#ok<SAGROW>
        control_trace_list(end+1, :) = this_trace; %#ok<SAGROW>
        control_n1amp_list(end+1) = n1_amp; %#ok<SAGROW>
        control_n1lat_list(end+1) = n1_lat; %#ok<SAGROW>
        control_color_list(end+1, :) = control_cmap(s, :); %#ok<SAGROW>
    end
end

% --- Manual exclusion: drop SERP26 and SERP27 from the pDCD group ---
manual_exclude_pdcd = {'17', '18', '21', '26', '27'};

pdcd_ID_list = {}; pdcd_trace_list = []; pdcd_n1amp_list = []; pdcd_n1lat_list = []; pdcd_color_list = [];
for s = 1:n_pdcd
    if ismember(pdcd_IDs{s}, manual_exclude_pdcd)
        continue; % manual drop
    end
    if ~isempty(ERP_results.pDCD{s}) && ERP_results.pDCD{s}.ntrials.accepted(bin_to_plot) > 0
        this_trace = squeeze(ERP_results.pDCD{s}.bindata(ch_idx, :, bin_to_plot))';
        [n1_amp, n1_lat] = find_n1_peak(this_trace, time_ms, n1_window_ms);
        pdcd_ID_list{end+1} = pdcd_IDs{s}; %#ok<SAGROW>
        pdcd_trace_list(end+1, :) = this_trace; %#ok<SAGROW>
        pdcd_n1amp_list(end+1) = n1_amp; %#ok<SAGROW>
        pdcd_n1lat_list(end+1) = n1_lat; %#ok<SAGROW>
        pdcd_color_list(end+1, :) = pdcd_cmap(s, :); %#ok<SAGROW>
    end
end

% --- Step 2: exclude N1-amplitude outliers (> 2 SD from that group's own mean) ---
control_mean_n1 = mean(control_n1amp_list);
control_sd_n1   = std(control_n1amp_list);
control_keep    = abs(control_n1amp_list - control_mean_n1) <= 2 * control_sd_n1;

pdcd_mean_n1 = mean(pdcd_n1amp_list);
pdcd_sd_n1   = std(pdcd_n1amp_list);
pdcd_keep    = abs(pdcd_n1amp_list - pdcd_mean_n1) <= 2 * pdcd_sd_n1;

fprintf('Control: manually excluded SERP%s\n', strjoin(manual_exclude_control, ', SERP'));
if any(~control_keep)
    fprintf('Control: excluded N1-outlier(s) SERP%s\n', strjoin(control_ID_list(~control_keep), ', SERP'));
else
    fprintf('Control: no N1-amplitude outliers beyond 2 SD.\n');
end
if any(~pdcd_keep)
    fprintf('pDCD: excluded N1-outlier(s) SERP%s\n', strjoin(pdcd_ID_list(~pdcd_keep), ', SERP'));
else
    fprintf('pDCD: no N1-amplitude outliers beyond 2 SD.\n');
end

control_trace_kept = control_trace_list(control_keep, :);
control_n1amp_kept = control_n1amp_list(control_keep);
control_n1lat_kept = control_n1lat_list(control_keep);
control_color_kept = control_color_list(control_keep, :);

pdcd_trace_kept = pdcd_trace_list(pdcd_keep, :);
pdcd_n1amp_kept = pdcd_n1amp_list(pdcd_keep);
pdcd_n1lat_kept = pdcd_n1lat_list(pdcd_keep);
pdcd_color_kept = pdcd_color_list(pdcd_keep, :);

% --- Step 3: plot grand average +/- 1SD, range, and N1 markers on the cleaned subject sets ---
figure('Color', 'w', 'Position', [950, 100, 800, 500], 'Name', 'Cz - Bin 1 - Grand Average (Excl. Outliers)');
hold on;

leg_handles = [];
leg_labels = {};
x_patch = [time_ms, fliplr(time_ms)];

if ~isempty(control_trace_kept)
    mean_control = mean(control_trace_kept, 1);
    sd_control   = std(control_trace_kept, 0, 1);
    min_control  = min(control_trace_kept, [], 1);
    max_control  = max(control_trace_kept, [], 1);

    p_range_c = fill(x_patch, [max_control, fliplr(min_control)], control_light, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    p_sd_c = fill(x_patch, [mean_control + sd_control, fliplr(mean_control - sd_control)], control_light, ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
    l_mean_c = plot(time_ms, mean_control, 'Color', control_dark, 'LineWidth', 2.5);

    for s = 1:length(control_n1lat_kept)
        plot(control_n1lat_kept(s), control_n1amp_kept(s), 'o', 'MarkerEdgeColor', control_color_kept(s, :), ...
            'MarkerFaceColor', control_color_kept(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    [n1_amp_avg_c, n1_lat_avg_c] = find_n1_peak(mean_control, time_ms, n1_window_ms);
    l_n1_avg_c = plot(n1_lat_avg_c, n1_amp_avg_c, 'x', 'Color', control_dark, 'MarkerSize', 10, 'LineWidth', 2.5);

    leg_handles = [leg_handles, l_mean_c, p_sd_c, p_range_c, l_n1_avg_c];
    leg_labels  = [leg_labels, {'Control Mean', 'Control \pm1 SD', 'Control Range', 'Control Average N1'}];
end

if ~isempty(pdcd_trace_kept)
    mean_pdcd = mean(pdcd_trace_kept, 1);
    sd_pdcd   = std(pdcd_trace_kept, 0, 1);
    min_pdcd  = min(pdcd_trace_kept, [], 1);
    max_pdcd  = max(pdcd_trace_kept, [], 1);

    p_range_p = fill(x_patch, [max_pdcd, fliplr(min_pdcd)], pdcd_light, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    p_sd_p = fill(x_patch, [mean_pdcd + sd_pdcd, fliplr(mean_pdcd - sd_pdcd)], pdcd_light, ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
    l_mean_p = plot(time_ms, mean_pdcd, 'Color', pdcd_dark, 'LineWidth', 2.5);

    for s = 1:length(pdcd_n1lat_kept)
        plot(pdcd_n1lat_kept(s), pdcd_n1amp_kept(s), 'o', 'MarkerEdgeColor', pdcd_color_kept(s, :), ...
            'MarkerFaceColor', pdcd_color_kept(s, :), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    [n1_amp_avg_p, n1_lat_avg_p] = find_n1_peak(mean_pdcd, time_ms, n1_window_ms);
    l_n1_avg_p = plot(n1_lat_avg_p, n1_amp_avg_p, 'x', 'Color', pdcd_dark, 'MarkerSize', 10, 'LineWidth', 2.5);

    leg_handles = [leg_handles, l_mean_p, p_sd_p, p_range_p, l_n1_avg_p];
    leg_labels  = [leg_labels, {'pDCD Mean', 'pDCD \pm1 SD', 'pDCD Range', 'pDCD Average N1'}];
end

xlim([-200 500]); grid on;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1);

xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
title('Electrode: Cz - Bin 1 - Grand Average (Excl. Manual + N1-Outliers; o = individual N1, x = average N1)', ...
    'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10);

if ~isempty(leg_handles)
    legend(leg_handles, leg_labels, 'Location', 'southwest');
end

hold off;
fprintf('Done! Outlier-excluded grand average plotted.\n');

%% STATS: Group Comparison of N1 Amplitude and Latency (After Exclusions)
fprintf('\n=======================================================\n');
fprintf('   N1 AMPLITUDE & LATENCY: CONTROL vs pDCD (Cz, Bin 1)   \n');
fprintf('=======================================================\n');

% --- N1 AMPLITUDE ---
fprintf('\n--- N1 Amplitude (uV) ---\n');
fprintf('Control: n=%d, mean=%.2f, SD=%.2f\n', length(control_n1amp_kept), mean(control_n1amp_kept), std(control_n1amp_kept));
fprintf('pDCD:    n=%d, mean=%.2f, SD=%.2f\n', length(pdcd_n1amp_kept), mean(pdcd_n1amp_kept), std(pdcd_n1amp_kept));

% Normality check (Lilliefors) - guides whether the t-test or rank-sum is more trustworthy
[~, p_norm_control_amp] = lillietest(control_n1amp_kept);
[~, p_norm_pdcd_amp]    = lillietest(pdcd_n1amp_kept);
fprintf('Normality (Lilliefors) p-values: Control=%.3f, pDCD=%.3f (p<.05 suggests non-normal)\n', ...
    p_norm_control_amp, p_norm_pdcd_amp);

% Welch's t-test (unequal variance - safer default with small/uneven N)
[~, p_amp, ci_amp, stats_amp] = ttest2(control_n1amp_kept, pdcd_n1amp_kept, 'Vartype', 'unequal');
pooled_sd_amp = sqrt((std(control_n1amp_kept)^2 + std(pdcd_n1amp_kept)^2) / 2);
cohend_amp = (mean(control_n1amp_kept) - mean(pdcd_n1amp_kept)) / pooled_sd_amp;
fprintf('Welch t-test: t(%.1f) = %.2f, p = %.4f, 95%% CI of diff = [%.2f, %.2f], Cohen''s d = %.2f\n', ...
    stats_amp.df, stats_amp.tstat, p_amp, ci_amp(1), ci_amp(2), cohend_amp);

% Nonparametric alternative (use if normality above looks violated)
[p_amp_rs] = ranksum(control_n1amp_kept, pdcd_n1amp_kept);
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_amp_rs);

% --- N1 LATENCY ---
fprintf('\n--- N1 Latency (ms) ---\n');
fprintf('Control: n=%d, mean=%.1f, SD=%.1f\n', length(control_n1lat_kept), mean(control_n1lat_kept), std(control_n1lat_kept));
fprintf('pDCD:    n=%d, mean=%.1f, SD=%.1f\n', length(pdcd_n1lat_kept), mean(pdcd_n1lat_kept), std(pdcd_n1lat_kept));

[~, p_norm_control_lat] = lillietest(control_n1lat_kept);
[~, p_norm_pdcd_lat]    = lillietest(pdcd_n1lat_kept);
fprintf('Normality (Lilliefors) p-values: Control=%.3f, pDCD=%.3f (p<.05 suggests non-normal)\n', ...
    p_norm_control_lat, p_norm_pdcd_lat);

[~, p_lat, ci_lat, stats_lat] = ttest2(control_n1lat_kept, pdcd_n1lat_kept, 'Vartype', 'unequal');
pooled_sd_lat = sqrt((std(control_n1lat_kept)^2 + std(pdcd_n1lat_kept)^2) / 2);
cohend_lat = (mean(control_n1lat_kept) - mean(pdcd_n1lat_kept)) / pooled_sd_lat;
fprintf('Welch t-test: t(%.1f) = %.2f, p = %.4f, 95%% CI of diff = [%.2f, %.2f], Cohen''s d = %.2f\n', ...
    stats_lat.df, stats_lat.tstat, p_lat, ci_lat(1), ci_lat(2), cohend_lat);

[p_lat_rs] = ranksum(control_n1lat_kept, pdcd_n1lat_kept);
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_lat_rs);

fprintf('\n=======================================================\n\n');


%% STATS (Reference-Point Method): N1 Amplitude Re-Baselined to a Pre-Movement Window
fprintf('\n=======================================================\n');
fprintf(' N1 STATS USING A PRE-MOVEMENT REFERENCE BASELINE (Toledo et al. style)\n');
fprintf('=======================================================\n');

% NOTE: Toledo et al. (2016, Clin Neurophysiol) baseline-correct N1 amplitude
% against a reference interval 1.5-0.5 s BEFORE movement onset. Your current
% epoch only spans -500 to 1000 ms, so that exact window isn't available
% without re-running pop_epochbin with a wider pre-stimulus range. As the
% closest available proxy, this block re-references amplitude to the
% earliest window your epoch actually contains: -500 to -200 ms (the slice
% immediately before ERPLAB's own -200 to 0 ms baseline).
%
% Latency is NOT touched here: Toledo et al. measure N1 latency simply as
% time-from-onset within the 0-300 ms search window (no reference
% correction needed), which is already what find_n1_peak returns. Also note
% re-referencing amplitude by a constant offset cannot change *where* the
% minimum falls in time, so the latency values below will match your
% earlier _kept values exactly - shown here for completeness.

ref_window_ms = [-500 -200];   % proxy pre-movement reference window
ref_idx = find(time_ms >= ref_window_ms(1) & time_ms <= ref_window_ms(2));

if isempty(ref_idx)
    error('Reference window %d to %d ms falls outside the current epoch.', ref_window_ms(1), ref_window_ms(2));
end

% --- Control group: re-reference each kept subject's trace, re-pick N1 ---
control_n1amp_ref = zeros(1, size(control_trace_kept, 1));
control_n1lat_ref = zeros(1, size(control_trace_kept, 1));
for s = 1:size(control_trace_kept, 1)
    ref_offset = mean(control_trace_kept(s, ref_idx));
    rereferenced_trace = control_trace_kept(s, :) - ref_offset;
    [amp, lat] = find_n1_peak(rereferenced_trace, time_ms, n1_window_ms);
    control_n1amp_ref(s) = amp;
    control_n1lat_ref(s) = lat;
end

% --- pDCD group: same procedure ---
pdcd_n1amp_ref = zeros(1, size(pdcd_trace_kept, 1));
pdcd_n1lat_ref = zeros(1, size(pdcd_trace_kept, 1));
for s = 1:size(pdcd_trace_kept, 1)
    ref_offset = mean(pdcd_trace_kept(s, ref_idx));
    rereferenced_trace = pdcd_trace_kept(s, :) - ref_offset;
    [amp, lat] = find_n1_peak(rereferenced_trace, time_ms, n1_window_ms);
    pdcd_n1amp_ref(s) = amp;
    pdcd_n1lat_ref(s) = lat;
end

% --- N1 AMPLITUDE (reference-corrected) ---
fprintf('\n--- N1 Amplitude, re-referenced to %d to %d ms baseline (uV) ---\n', ref_window_ms(1), ref_window_ms(2));
fprintf('Control: n=%d, mean=%.2f, SD=%.2f\n', length(control_n1amp_ref), mean(control_n1amp_ref), std(control_n1amp_ref));
fprintf('pDCD:    n=%d, mean=%.2f, SD=%.2f\n', length(pdcd_n1amp_ref), mean(pdcd_n1amp_ref), std(pdcd_n1amp_ref));

[~, p_amp_ref, ci_amp_ref, stats_amp_ref] = ttest2(control_n1amp_ref, pdcd_n1amp_ref, 'Vartype', 'unequal');
pooled_sd_amp_ref = sqrt((std(control_n1amp_ref)^2 + std(pdcd_n1amp_ref)^2) / 2);
cohend_amp_ref = (mean(control_n1amp_ref) - mean(pdcd_n1amp_ref)) / pooled_sd_amp_ref;
fprintf('Welch t-test: t(%.1f) = %.2f, p = %.4f, 95%% CI of diff = [%.2f, %.2f], Cohen''s d = %.2f\n', ...
    stats_amp_ref.df, stats_amp_ref.tstat, p_amp_ref, ci_amp_ref(1), ci_amp_ref(2), cohend_amp_ref);

p_amp_ref_rs = ranksum(control_n1amp_ref, pdcd_n1amp_ref);
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_amp_ref_rs);

% --- N1 LATENCY (unchanged - already trigger-referenced, matching the paper) ---
fprintf('\n--- N1 Latency, time from movement onset (ms) - unchanged ---\n');
fprintf('Control: n=%d, mean=%.1f, SD=%.1f\n', length(control_n1lat_ref), mean(control_n1lat_ref), std(control_n1lat_ref));
fprintf('pDCD:    n=%d, mean=%.1f, SD=%.1f\n', length(pdcd_n1lat_ref), mean(pdcd_n1lat_ref), std(pdcd_n1lat_ref));

[~, p_lat_ref, ci_lat_ref, stats_lat_ref] = ttest2(control_n1lat_ref, pdcd_n1lat_ref, 'Vartype', 'unequal');
pooled_sd_lat_ref = sqrt((std(control_n1lat_ref)^2 + std(pdcd_n1lat_ref)^2) / 2);
cohend_lat_ref = (mean(control_n1lat_ref) - mean(pdcd_n1lat_ref)) / pooled_sd_lat_ref;
fprintf('Welch t-test: t(%.1f) = %.2f, p = %.4f, 95%% CI of diff = [%.2f, %.2f], Cohen''s d = %.2f\n', ...
    stats_lat_ref.df, stats_lat_ref.tstat, p_lat_ref, ci_lat_ref(1), ci_lat_ref(2), cohend_lat_ref);

p_lat_ref_rs = ranksum(control_n1lat_ref, pdcd_n1lat_ref);
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_lat_ref_rs);

fprintf('\n=======================================================\n\n');