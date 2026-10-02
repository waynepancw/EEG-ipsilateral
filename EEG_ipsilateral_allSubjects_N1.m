clear all, close all, clc

%% =====================================================================
%  N1 analysis: Control vs pDCD, passive left-foot stimulation
%  Revised after PI review:
%    [PI #1] N1 search window narrowed to 0-200 ms (Tseng et al., 2024)
%    [PI #2] SERP27 re-included (narrower window resolves the misdetection)
%    [PI #3] Accepted trials, rejection rate and baseline noise reported
%            for every subject; SERP08 now excluded manually
%    [PI #4] Amplitude and latency at Fz, FCz, Cz, Pz, C3 and C4
%  Output: N1 (latency, amplitude) marked on every plot and exported to an
%  Excel file with one row per subject (no statistical tests), plus
%  per-subject scalp maps at fixed latencies.
%  =====================================================================

%% INITIALIZATION & SETTINGS
fprintf('Initializing paths...\n');

% 1. Hardcode the exact data path
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
    if all_plugins(i).isdir && contains(lower(all_plugins(i).name), 'erplab')
        erplab_root = fullfile(plugins_dir, all_plugins(i).name);
        addpath(genpath(erplab_root));
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

groups  = {'Control', 'pDCD'};
all_IDs = {control_IDs, pdcd_IDs};
n_control = length(control_IDs);
n_pdcd    = length(pdcd_IDs);

%% ANALYSIS PARAMETERS  (every decision from the PI review lives here)

bin_to_plot = 1;

% [PI #4] Electrodes to analyse (also the column order in the Excel file).
% Labels must match the cap labels (case-insensitive). A label that is not
% found gives a warning and empty (NaN) columns.
chans_to_analyze = {'Fz', 'FCz', 'Cz', 'Pz', 'C3', 'C4'};

% Electrode used for decisions made ONCE per subject (baseline-noise QC and
% the N1-amplitude outlier rule), so that all electrodes are analysed on
% exactly the same participants.
ref_chan = 'Cz';

% [PI #1] N1 search window: was [0 300] (Toledo et al., 2016).
% Now 0-200 ms (Tseng et al., 2024). Use [100 200] for the stricter option.
n1_window_ms = [0 200];

% [PI #2] SERP27: the PI expects the narrower window alone to fix the
% misdetection, so the N1 is picked with the SAME method as before: the most
% negative sample inside n1_window_ms.
% (Optional, not requested: set e.g. 10 to require a local peak within
% +/- 10 ms, like ERPLAB's "local peak" option.)
peak_neighborhood_ms = 0;

baseline_window_ms = [-200 0];    % same as the pop_epochbin baseline
pm_window_ms       = [-200 500];  % window for the plus-minus noise estimate

% [PI #3] Accepted trial count, artifact rejection rate and baseline noise
% are REPORTED for every subject. Nothing is excluded automatically on
% these numbers; the exclusion decision rests with the PI.
qc_spotlight_IDs = {'08'};   % subjects to print a detailed summary for

% Manual exclusions (visual inspection), applied to every electrode:
%  - Control: SERP08 excluded after reviewing his data-quality numbers (PI #3)
%  - pDCD:    SERP26 excluded, as in the original script
%  - SERP27 is no longer excluded (PI #2: narrower window resolves it)
% To re-include a subject, remove his ID, e.g. manual_exclude_control = {};
manual_exclude_control = {'08'};
manual_exclude_pdcd    = {'26'};

% N1-amplitude outlier rule (evaluated at ref_chan, within each group)
apply_outlier_rule = true;
outlier_sd = 3;

plot_individual_grids = true;
% false = each subject's panel is auto-scaled (small responses look noisy);
% true  = all panels of an electrode share one y-axis, so amplitudes and
%         noise can be compared by eye across subjects
grid_common_ylim = false;
output_dir = fullfile(filepath, 'N1_results');

%% STORAGE
ERP_results      = struct('Control', {cell(1, n_control)}, 'pDCD', {cell(1, n_pdcd)});
status_log       = struct('Control', {cell(1, n_control)}, 'pDCD', {cell(1, n_pdcd)});
pm_noise_results = struct('Control', {cell(1, n_control)}, 'pDCD', {cell(1, n_pdcd)});

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
            ev_types = {EEG.event.type};
            ev_types(strcmpi(ev_types, 'boundary')) = [];
            if ischar(ev_types{1}) || isstring(ev_types{1})
                ev_nums = cellfun(@str2double, ev_types);
            else
                ev_nums = cell2mat(ev_types);
            end

            [unique_trigs, ~, idx] = unique(ev_nums);
            counts = histcounts(idx, 1:length(unique_trigs)+1);
            [~, freq_sort_idx] = sort(counts, 'descend');
            top_4_frequent = unique_trigs(freq_sort_idx(1:min(4, length(unique_trigs))));
            top_4_trigs = sort(top_4_frequent, 'ascend');

            bin_file_path = fullfile(filepath, sprintf('temp_bins_%s.txt', current_id));
            fid = fopen(bin_file_path, 'w');
            for b = 1:length(top_4_trigs)
                fprintf(fid, 'bin %d\nTrigger_%d\n.{%d}\n\n', b, top_4_trigs(b), top_4_trigs(b));
            end
            fclose(fid);

            % B. Create EventList
            EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

            % C. Assign Bins using their custom BDF file
            EEG = pop_binlister(EEG, 'BDF', bin_file_path, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG');
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
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

            % G. [PI #3] Trial-level noise estimate for QC (plus-minus
            %    average: alternate trials are sign-flipped so the ERP
            %    cancels and only noise remains)
            [pm_vec, pm_n] = compute_plusminus_noise(EEG, bin_to_plot, pm_window_ms);
            pm_noise_results.(current_group_name){s} = pm_vec;
            if ~isnan(pm_n) && abs(pm_n - ERP.ntrials.accepted(bin_to_plot)) > 1
                warning('SERP%s: plus-minus estimate used %d epochs but ERPLAB accepted %d. Treat +/- RMS with caution.', ...
                    current_id, pm_n, ERP.ntrials.accepted(bin_to_plot));
            end

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

time_ms = double(valid_erp.times(:)');
nT = numel(time_ms);
nbhd_samples = round(peak_neighborhood_ms * valid_erp.srate / 1000);   % 0 = most negative sample

chan_found = true(1, numel(chans_to_analyze));
for c = 1:numel(chans_to_analyze)
    if ~any(strcmpi({valid_erp.chanlocs.labels}, chans_to_analyze{c}))
        chan_found(c) = false;
        warning('Channel %s not found (its columns will be empty). Available labels: %s', ...
            chans_to_analyze{c}, strjoin({valid_erp.chanlocs.labels}, ', '));
    end
end
if ~ismember(lower(ref_chan), lower(chans_to_analyze))
    error('ref_chan (%s) must be one of chans_to_analyze.', ref_chan);
end

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% SUBJECT BOOKKEEPING (one row per subject, both groups)
all_subj_ids   = [control_IDs, pdcd_IDs]';
all_subj_group = [repmat({'Control'}, n_control, 1); repmat({'pDCD'}, n_pdcd, 1)];
all_subj_index = [(1:n_control)'; (1:n_pdcd)'];
is_control = strcmp(all_subj_group, 'Control');
is_pdcd    = strcmp(all_subj_group, 'pDCD');
nS = numel(all_subj_ids);
nC = numel(chans_to_analyze);
rc = find(strcmpi(chans_to_analyze, ref_chan), 1);

% Colour gradients (Control = blues, pDCD = reds)
control_light = [0.65 0.80 1.00];  control_dark = [0.00 0.10 0.55];
pdcd_light    = [1.00 0.65 0.65];  pdcd_dark    = [0.55 0.00 0.00];
control_cmap = [linspace(control_light(1), control_dark(1), n_control)', ...
                linspace(control_light(2), control_dark(2), n_control)', ...
                linspace(control_light(3), control_dark(3), n_control)'];
pdcd_cmap    = [linspace(pdcd_light(1), pdcd_dark(1), n_pdcd)', ...
                linspace(pdcd_light(2), pdcd_dark(2), n_pdcd)', ...
                linspace(pdcd_light(3), pdcd_dark(3), n_pdcd)'];
subj_color = [control_cmap; pdcd_cmap];

%% EXTRACT TRACES, QC METRICS AND N1 MEASURES (all subjects x all electrodes)
has_data   = false(nS, 1);
n_acc      = nan(nS, 1);
n_rej      = nan(nS, 1);
traces     = nan(nS, nC, nT);
bl_rms     = nan(nS, nC);   % RMS of the averaged ERP in the baseline (uV)
pm_rms     = nan(nS, nC);   % RMS of the plus-minus average (uV)
snr        = nan(nS, nC);   % RMS in N1 window / baseline RMS
p2p_epoch  = nan(nS, nC);   % peak-to-peak range from -200 to 500 ms (uV)
n1_amp     = nan(nS, nC);
n1_lat     = nan(nS, nC);
n1_local   = false(nS, nC);
n1_edge    = false(nS, nC);

bl_idx   = time_ms >= baseline_window_ms(1) & time_ms <= baseline_window_ms(2);
win_idx  = time_ms >= n1_window_ms(1)       & time_ms <= n1_window_ms(2);
show_idx = time_ms >= -200                  & time_ms <= 500;

for i = 1:nS
    ERP = ERP_results.(all_subj_group{i}){all_subj_index(i)};
    if isempty(ERP) || ERP.ntrials.accepted(bin_to_plot) == 0
        continue;
    end
    if numel(ERP.times) ~= nT
        warning('SERP%s has a different epoch length; skipped.', all_subj_ids{i});
        continue;
    end
    has_data(i) = true;
    n_acc(i) = ERP.ntrials.accepted(bin_to_plot);
    if isfield(ERP.ntrials, 'rejected')
        n_rej(i) = ERP.ntrials.rejected(bin_to_plot);
    end
    pm_vec = pm_noise_results.(all_subj_group{i}){all_subj_index(i)};

    for c = 1:nC
        ch = find(strcmpi({ERP.chanlocs.labels}, chans_to_analyze{c}), 1);
        if isempty(ch), continue; end
        tr = double(ERP.bindata(ch, :, bin_to_plot));
        tr = tr(:)';
        traces(i, c, :) = tr;

        bl_rms(i, c)    = sqrt(mean(tr(bl_idx).^2));
        snr(i, c)       = sqrt(mean(tr(win_idx).^2)) / bl_rms(i, c);
        p2p_epoch(i, c) = max(tr(show_idx)) - min(tr(show_idx));
        if numel(pm_vec) >= ch, pm_rms(i, c) = pm_vec(ch); end

        [n1_amp(i, c), n1_lat(i, c), n1_local(i, c), n1_edge(i, c)] = ...
            find_n1_peak(tr, time_ms, n1_window_ms, nbhd_samples);
    end
end
rej_rate = n_rej ./ (n_acc + n_rej);

%% [PI #3] DATA-QUALITY REPORT (reported only; no automatic exclusion)
fprintf('\n=======================================================================================\n');
fprintf(' DATA QUALITY REPORT (Bin %d, noise metrics at %s)\n', bin_to_plot, ref_chan);
fprintf('=======================================================================================\n');
fprintf('%-8s %-8s %5s %5s %6s | %8s | %8s | %6s | %8s\n', ...
    'Subject', 'Group', 'Acc', 'Rej', 'Rej%', 'BL RMS', '+/- RMS', 'SNR', 'P2P');
for i = 1:nS
    if ~has_data(i)
        fprintf('SERP%-4s %-8s   (no valid data)\n', all_subj_ids{i}, all_subj_group{i});
        continue;
    end
    fprintf('SERP%-4s %-8s %5d %5d %5.1f%% | %8.2f | %8.2f | %6.2f | %8.2f\n', ...
        all_subj_ids{i}, all_subj_group{i}, n_acc(i), n_rej(i), 100*rej_rate(i), ...
        bl_rms(i, rc), pm_rms(i, rc), snr(i, rc), p2p_epoch(i, rc));
end
fprintf(['Acc/Rej = accepted/rejected epochs; BL RMS = RMS of the average from %d to %d ms;\n' ...
         '+/- RMS = plus-minus average noise; SNR = RMS in N1 window / BL RMS; P2P = peak-to-peak -200..500 ms (uV)\n'], ...
         baseline_window_ms);

% Detailed summary for specific subjects (e.g. SERP08): value, sample
% median and rank among the subjects with data (1 = lowest)
spot_names = {'Accepted trials', 'Rejection rate (%)', ['Baseline RMS @' ref_chan], ...
              ['Plus-minus RMS @' ref_chan], ['SNR @' ref_chan], ['Peak-to-peak @' ref_chan]};
spot_vals  = {n_acc, 100*rej_rate, bl_rms(:, rc), pm_rms(:, rc), snr(:, rc), p2p_epoch(:, rc)};
for k = 1:numel(qc_spotlight_IDs)
    i = find(strcmp(all_subj_ids, qc_spotlight_IDs{k}), 1);
    if isempty(i) || ~has_data(i), continue; end
    fprintf('\n--- Data quality: SERP%s (%s) vs. all %d subjects with data ---\n', ...
        all_subj_ids{i}, all_subj_group{i}, sum(has_data));
    fprintf('%-22s %8s %8s %12s\n', 'Metric', 'SERP', 'Median', 'Rank (1=low)');
    for q = 1:numel(spot_names)
        v = spot_vals{q}(has_data);
        v = v(~isnan(v));
        fprintf('%-22s %8.2f %8.2f %6d / %-3d\n', spot_names{q}, spot_vals{q}(i), median(v), ...
            sum(v < spot_vals{q}(i)) + 1, numel(v));
    end
end

%% EXCLUSIONS (applied identically to every electrode)
excl_reason = repmat({''}, nS, 1);
excl_short  = repmat({''}, nS, 1);   % short code shown in the figure titles
for i = 1:nS
    if ~has_data(i)
        excl_reason{i} = 'No valid data';
        excl_short{i}  = 'no data';
        continue;
    end
    is_manual = (is_control(i) && ismember(all_subj_ids{i}, manual_exclude_control)) || ...
                (is_pdcd(i)    && ismember(all_subj_ids{i}, manual_exclude_pdcd));
    if is_manual
        excl_reason{i} = 'Manual (visual inspection)';
        excl_short{i}  = 'manual';
        continue;
    end
end

if apply_outlier_rule
    for g = 1:2
        sel = find(strcmp(all_subj_group, groups{g}) & cellfun(@isempty, excl_reason));
        a = n1_amp(sel, rc);
        z = (a - mean(a, 'omitnan')) / std(a, 'omitnan');
        out = sel(abs(z) > outlier_sd);
        for k = 1:numel(out)
            excl_reason{out(k)} = sprintf('N1 amplitude outlier at %s (|z| > %g)', ref_chan, outlier_sd);
            excl_short{out(k)}  = 'outlier';
        end
    end
end
included = cellfun(@isempty, excl_reason);

fprintf('\n=======================================================\n');
fprintf(' EXCLUSIONS (same subject set used at every electrode)\n');
fprintf('=======================================================\n');
for g = 1:2
    sel = strcmp(all_subj_group, groups{g});
    fprintf('%s: %d included of %d\n', groups{g}, sum(included & sel), sum(sel));
    idx = find(sel & ~included);
    for k = 1:numel(idx)
        fprintf('   SERP%s excluded -> %s\n', all_subj_ids{idx(k)}, excl_reason{idx(k)});
    end
end

if peak_neighborhood_ms > 0
    fprintf('\nIncluded subjects with NO local N1 peak in %d-%d ms (absolute minimum used; please inspect):\n', n1_window_ms);
    any_flag = false;
    for i = find(included)'
        for c = 1:nC
            if ~n1_local(i, c)
                fprintf('   SERP%s (%s) at %-4s: %.2f uV at %.0f ms\n', all_subj_ids{i}, all_subj_group{i}, ...
                    chans_to_analyze{c}, n1_amp(i, c), n1_lat(i, c));
                any_flag = true;
            end
        end
    end
    if ~any_flag, fprintf('   none\n'); end
end

%% EXPORT N1 PEAKS TO EXCEL
% Sheet "N1_peaks": one row per subject; for each electrode two columns,
% N1 amplitude (uV) and N1 latency (ms). All subjects are listed; the
% second sheet says who is included/excluded and why, plus the data-quality
% numbers.
peak_tbl = table(strcat('SERP', all_subj_ids), all_subj_group, 'VariableNames', {'Subject', 'Group'});
for c = 1:nC
    cn = matlab.lang.makeValidName(chans_to_analyze{c});
    peak_tbl.([cn '_Amplitude_uV']) = n1_amp(:, c);
    peak_tbl.([cn '_Latency_ms'])   = n1_lat(:, c);
end

info_tbl = table(strcat('SERP', all_subj_ids), all_subj_group, included, excl_reason, ...
    n_acc, n_rej, 100*rej_rate, bl_rms(:, rc), pm_rms(:, rc), snr(:, rc), ...
    'VariableNames', {'Subject', 'Group', 'Included', 'ExclusionReason', 'AcceptedTrials', ...
    'RejectedTrials', 'RejectPct', ['BaselineRMS_uV_' ref_chan], ['PlusMinusRMS_uV_' ref_chan], ['SNR_' ref_chan]});

xlsx_file = fullfile(output_dir, sprintf('N1_peaks_bin%d_%d-%dms.xlsx', bin_to_plot, n1_window_ms));
try
    if exist(xlsx_file, 'file'), delete(xlsx_file); end   % no leftovers from an older run
    writetable(peak_tbl, xlsx_file, 'Sheet', 'N1_peaks');
    writetable(info_tbl, xlsx_file, 'Sheet', 'Inclusion_and_QC');
    fprintf('\nN1 peaks written to %s\n', xlsx_file);
catch ME
    warning('Could not write %s (is it open in Excel?): %s', xlsx_file, ME.message);
end

fprintf('\n=== N1 peaks (Bin %d, window %d-%d ms) ===\n', bin_to_plot, n1_window_ms);
disp(peak_tbl);

%% PLOT: 3x4 individual-subject grids for every electrode and group
if plot_individual_grids
    for c = 1:nC
        if ~chan_found(c), continue; end
        X_all = reshape(traces(has_data, c, show_idx), sum(has_data), []);
        common_yl = [floor(min(X_all(:))/5)*5, ceil(max(X_all(:))/5)*5];
        for g = 1:2
            ids = all_IDs{g};
            light_figure('Position', [100, 100, 1200, 800], ...
                'Name', sprintf('%s - %s Bin %d - All Subjects', groups{g}, chans_to_analyze{c}, bin_to_plot));
            for s = 1:numel(ids)
                i = find(strcmp(all_subj_group, groups{g}) & all_subj_index == s, 1);
                subplot(3, 4, s); hold on;
                if has_data(i)
                    tr = reshape(traces(i, c, :), 1, []);
                    plot(time_ms, tr, 'Color', subj_color(i, :), 'LineWidth', 1.2);
                    if n1_local(i, c) || peak_neighborhood_ms == 0, face = subj_color(i, :); else, face = 'w'; end
                    plot(n1_lat(i, c), n1_amp(i, c), 'o', 'MarkerEdgeColor', 'k', ...
                        'MarkerFaceColor', face, 'MarkerSize', 6);
                else
                    text(0.5, 0.5, 'No valid data', 'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
                end
                xlim([-200 500]); grid on;
                set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 8);
                if grid_common_ylim, ylim(common_yl); end
                yl = ylim;
                has_peak = has_data(i) && ~isnan(n1_amp(i, c));
                if has_peak   % make room below the valley for its label
                    yr = diff(yl);
                    yl(1) = min(yl(1), n1_amp(i, c) - 0.25 * yr);
                    yr = diff(yl);
                end
                hp = patch([n1_window_ms(1) n1_window_ms(2) n1_window_ms(2) n1_window_ms(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
                uistack(hp, 'bottom');
                line([0 0], yl, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
                line(xlim, [0 0], 'Color', 'k', 'LineWidth', 0.75);
                ylim(yl);
                if has_peak
                    text(n1_lat(i, c), n1_amp(i, c) - 0.04 * yr, ...
                        sprintf('(%.0f ms, %.1f \\muV)', n1_lat(i, c), n1_amp(i, c)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                        'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
                end
                ttl = sprintf('SERP%s', ids{s}); tcol = 'k';
                if ~included(i), ttl = sprintf('%s (%s)', ttl, excl_short{i}); tcol = [0.6 0.6 0.6]; end
                title(ttl, 'FontSize', 9, 'Color', tcol);
                hold off;
            end
            sgtitle(sprintf('%s - %s - Bin %d (shaded = N1 window %d-%d ms, o = N1 (latency, amplitude))', ...
                groups{g}, chans_to_analyze{c}, bin_to_plot, n1_window_ms), 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k');
        end
    end
end

%% PLOT: Grand averages after exclusions, all electrodes
light_figure('Position', [50, 50, 1500, 850], 'Name', 'Grand averages after exclusions - all electrodes');
light_c = {control_light, pdcd_light};
dark_c  = {control_dark,  pdcd_dark};
n_rows = ceil(nC / 3);
for c = 1:nC
    subplot(n_rows, 3, c); hold on;
    title(chans_to_analyze{c}, 'FontWeight', 'bold', 'Color', 'k');
    if ~chan_found(c)
        text(0.5, 0.5, 'Channel not found', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        axis off; continue;
    end
    hl = []; lab = {}; ga_txt = {}; ga_col = {};
    for g = 1:2
        sel = included & strcmp(all_subj_group, groups{g});
        X = reshape(traces(sel, c, :), [], nT);
        X = X(all(~isnan(X), 2), :);
        if isempty(X), continue; end
        mu = mean(X, 1); sd = std(X, 0, 1);
        fill([time_ms, fliplr(time_ms)], [mu + sd, fliplr(mu - sd)], light_c{g}, ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hl(end+1) = plot(time_ms, mu, 'Color', dark_c{g}, 'LineWidth', 2); %#ok<SAGROW>
        lab{end+1} = sprintf('%s mean \\pm1 SD (n=%d)', groups{g}, size(X, 1)); %#ok<SAGROW>
        [a, l] = find_n1_peak(mu, time_ms, n1_window_ms, nbhd_samples);
        plot(l, a, 'x', 'Color', dark_c{g}, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        ga_txt{end+1} = sprintf('%s N1: (%.0f ms, %.1f \\muV)', groups{g}, l, a); %#ok<SAGROW>
        ga_col{end+1} = dark_c{g}; %#ok<SAGROW>
    end
    xlim([-200 500]); grid on;
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 10);
    yl = ylim;
    hp = patch([n1_window_ms(1) n1_window_ms(2) n1_window_ms(2) n1_window_ms(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    uistack(hp, 'bottom');
    line([0 0], yl, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
    line(xlim, [0 0], 'Color', 'k', 'HandleVisibility', 'off');
    ylim(yl);
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    for k = 1:numel(ga_txt)   % N1 (latency, amplitude) in the lower-right corner
        text(0.98, 0.04 + 0.08 * (numel(ga_txt) - k), ga_txt{k}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', ga_col{k});
    end
    if c == 1 && ~isempty(hl), legend(hl, lab, 'Location', 'southwest', 'FontSize', 8); end
    hold off;
end
sgtitle(sprintf('Bin %d grand averages after exclusions (x = N1 of grand average, shaded = %d-%d ms)', ...
    bin_to_plot, n1_window_ms), 'FontWeight', 'bold', 'Color', 'k');

fprintf('\nDone.\n');


%% PLOT: Scalp topographies at fixed latencies (rows = subjects, cols = times)
topo_times_ms = [0 100 200];
topo_scale    = 'subject';   % 'subject' = one colour scale per subject; 'all' = one scale for the whole figure
topo_clim     = [];          % e.g. [-15 15] to force a fixed scale (overrides topo_scale)

% topoplot needs electrode coordinates; Neuroscan .cnt files often carry
% labels only, so look the positions up from the standard cap if needed.
topo_chanlocs = valid_erp.chanlocs;
if ~isfield(topo_chanlocs, 'theta') || all(cellfun(@isempty, {topo_chanlocs.theta}))
    elp = which('standard-10-5-cap385.elp');
    if isempty(elp)
        d = dir(fullfile(eeglab_path, 'plugins', '**', 'standard-10-5-cap385.elp'));
        if ~isempty(d), elp = fullfile(d(1).folder, d(1).name); end
    end
    if isempty(elp)
        warning('No electrode coordinates and no standard cap file found; skipping scalp maps.');
        topo_chanlocs = [];
    else
        topo_chanlocs = pop_chanedit(topo_chanlocs, 'lookup', elp);
    end
end

if ~isempty(topo_chanlocs)
    good_ch = find(arrayfun(@(x) ~isempty(x.theta) && ~isempty(x.radius), topo_chanlocs));
    [~, topo_idx] = arrayfun(@(tt) min(abs(time_ms - tt)), topo_times_ms);
    nCt = numel(topo_times_ms);
    blue_red = [[linspace(0, 1, 32)'; ones(32, 1)], ...
        [linspace(0, 1, 32)'; linspace(1, 0, 32)'], ...
        [ones(32, 1); linspace(1, 0, 32)']];

    for g = 1:2
        ids = all_IDs{g};
        nR  = numel(ids);
        fig = light_figure('Position', [60, 40, 190*nCt + 160, 150*nR + 60], ...
            'Name', sprintf('%s - scalp maps - Bin %d', groups{g}, bin_to_plot));
        colormap(fig, blue_red);

        % values for every subject of this group
        V = nan(numel(good_ch), nCt, nR);
        for s = 1:nR
            i = find(strcmp(all_subj_group, groups{g}) & all_subj_index == s, 1);
            if ~has_data(i), continue; end
            ERP = ERP_results.(groups{g}){s};
            V(:, :, s) = double(ERP.bindata(good_ch, topo_idx, bin_to_plot));
        end
        if strcmpi(topo_scale, 'all'), cl_all = max(abs(V(:))); end

        for s = 1:nR
            i = find(strcmp(all_subj_group, groups{g}) & all_subj_index == s, 1);
            if ~isempty(topo_clim)
                cl = topo_clim;
            elseif strcmpi(topo_scale, 'all')
                cl = [-cl_all cl_all];
            else
                m = max(abs(reshape(V(:, :, s), [], 1)));
                if isempty(m) || isnan(m) || m == 0, m = 1; end
                cl = [-m m];
            end

            for k = 1:nCt
                subplot(nR, nCt, (s - 1)*nCt + k);
                if has_data(i)
                    topoplot(V(:, k, s), topo_chanlocs(good_ch), 'maplimits', cl, ...
                        'electrodes', 'on', 'shading', 'interp', 'conv', 'on');
                else
                    axis off;
                    text(0.5, 0.5, 'no data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                end
                if s == 1
                    title(sprintf('%d ms', round(time_ms(topo_idx(k)))), 'FontSize', 10, 'Color', 'k');
                end
                if k == 1   % row label, left of the head
                    lab = sprintf('SERP%s', ids{s});
                    lcol = 'k';
                    if ~included(i), lab = sprintf('%s (%s)', lab, excl_short{i}); lcol = [0.6 0.6 0.6]; end
                    text(-0.75, 0, sprintf('%s\n\\pm%.0f \\muV', lab, cl(2)), 'Clipping', 'off', ...
                        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', lcol);
                end
            end
        end
        sgtitle(sprintf('%s - scalp maps, Bin %d (blue = negative; colour scale per %s)', ...
            groups{g}, bin_to_plot, topo_scale), 'FontWeight', 'bold', 'Color', 'k');
    end
end



%% =====================================================================
%  LOCAL FUNCTIONS (must stay at the end of the script)
%  =====================================================================
function [peak_amp, peak_lat, is_local, at_edge] = find_n1_peak(trace, t, window, nbhd)
% FIND_N1_PEAK  N1 = most negative sample within window (ms).
%   nbhd = 0 (current setting, peak_neighborhood_ms = 0): absolute minimum
%   in the window.
%   nbhd > 0 (optional): the most negative LOCAL peak is used instead. A
%   sample k is a local peak if trace(k) is lower than all samples in
%   k-nbhd..k-1 and not higher than all samples in k+1..k+nbhd (neighbours
%   may lie outside the window). If no local peak exists, the absolute
%   minimum is returned and is_local = false.
    if nargin < 4, nbhd = 0; end
    trace = trace(:)'; t = t(:)';
    peak_amp = NaN; peak_lat = NaN; is_local = false; at_edge = false;
    idx = find(t >= window(1) & t <= window(2));
    if isempty(idx) || all(isnan(trace(idx))), return; end

    cand = [];
    if nbhd > 0
        N = numel(trace);
        for k = idx
            left  = trace(max(1, k - nbhd):k - 1);
            right = trace(k + 1:min(N, k + nbhd));
            if isempty(left) || isempty(right), continue; end
            if trace(k) < min(left) && trace(k) <= min(right)
                cand(end + 1) = k; %#ok<AGROW>
            end
        end
    end

    if ~isempty(cand)
        [peak_amp, j] = min(trace(cand));
        kpk = cand(j);
        is_local = true;
    else
        [peak_amp, j] = min(trace(idx));
        kpk = idx(j);
    end
    peak_lat = t(kpk);
    at_edge = (kpk == idx(1)) || (kpk == idx(end));
end

function [noise_rms, n_used] = compute_plusminus_noise(EEG, bin_num, win_ms)
% COMPUTE_PLUSMINUS_NOISE  Residual noise of the average (uV, per channel).
%   Uses the accepted epochs of bin_num, flips the sign of every second
%   epoch and averages ("plus-minus average"), which cancels the ERP while
%   keeping noise of the same magnitude as in the real average. Returns the
%   RMS over win_ms. Returns NaN if the epoch/flag fields are unavailable.
    noise_rms = nan(EEG.nbchan, 1);
    n_used = NaN;
    try
        nep = EEG.trials;
        rej = false(1, nep);
        if isfield(EEG, 'reject') && isfield(EEG.reject, 'rejmanual') && numel(EEG.reject.rejmanual) == nep
            rej = logical(EEG.reject.rejmanual);
        end
        inbin = false(1, nep);
        for e = 1:nep
            lat  = EEG.epoch(e).eventlatency;
            bini = EEG.epoch(e).eventbini;
            if ~iscell(lat),  lat  = {lat};  end
            if ~iscell(bini), bini = {bini}; end
            for k = 1:numel(lat)
                if abs(double(lat{k})) < 1 && any(double(bini{k}) == bin_num)  % time-locking event
                    inbin(e) = true;
                    break;
                end
            end
        end
        good = find(inbin & ~rej);
        n_used = numel(good);
        n_even = n_used - mod(n_used, 2);
        if n_even < 2, return; end
        good = good(1:n_even);
        signs = ones(1, 1, n_even);
        signs(1, 1, 2:2:end) = -1;
        pm = mean(double(EEG.data(:, :, good)) .* signs, 3);
        tidx = EEG.times >= win_ms(1) & EEG.times <= win_ms(2);
        noise_rms = sqrt(mean(pm(:, tidx).^2, 2));
    catch
        noise_rms = nan(EEG.nbchan, 1);
        n_used = NaN;
    end
end

function fig = light_figure(varargin)
% LIGHT_FIGURE  figure() with a white background regardless of MATLAB's
%   dark theme (R2025a+), so black lines, labels and titles stay visible.
    fig = figure(varargin{:});
    try
        theme(fig, 'light');
    catch
        set(fig, 'Color', 'w');
    end
    set(fig, 'Color', 'w');
end