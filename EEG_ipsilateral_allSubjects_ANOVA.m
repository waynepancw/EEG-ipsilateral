clear all, close all, clc

%% =====================================================================
%  N1 analysis: Control vs pDCD, passive left-foot stimulation
%  Revised after PI review:
%    [PI #1] N1 search window narrowed to 0-200 ms (Tseng et al., 2024)
%    [PI #2] SERP27 re-included (narrower window resolves the misdetection)
%    [PI #3] Accepted trials, rejection rate and baseline noise reported
%            for every subject (SERP08 decision left to the PI)
%    [PI #4] Amplitude and latency at Cz, FCz, Fz, C3 and C4
%    [PI note] Cohen's d with 95% CI reported next to every t test
%  =====================================================================

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

% [PI #4] Electrodes to analyse. Left-foot stimulation -> C4 is
% contralateral; Cz/FCz sit over the mesial foot representation.
chans_to_analyze = {'Cz', 'FCz', 'Fz'};

% Electrode used for decisions made ONCE per subject (baseline-noise QC and
% the N1-amplitude outlier rule), so that all electrodes are analysed on
% exactly the same participants. Required for the Group x Electrode ANOVA.
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
ref_window_ms      = [-500 -200]; % Toledo-style proxy reference (secondary amplitude measure)
pm_window_ms       = [-200 500];  % window for the plus-minus noise estimate

% [PI #3] Accepted trial count, artifact rejection rate and baseline noise
% are REPORTED for every subject. Nothing is excluded automatically on
% these numbers; the exclusion decision rests with the PI.
qc_spotlight_IDs = {'08'};   % subjects to print a detailed summary for

% Exclusions: as in the original script, except
%  - SERP27 removed from the list (PI #2: narrower window resolves it)
%  - SERP08 removed from the list pending the PI's decision on his data
%    quality numbers (PI #3). To exclude him: manual_exclude_control = {'08'};
manual_exclude_control = {'08'};
manual_exclude_pdcd    = {'26'};

% N1-amplitude outlier rule (evaluated at ref_chan, within each group)
apply_outlier_rule = true;
outlier_sd = 2;

alpha = 0.05;
plot_individual_grids = true;
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

for c = 1:numel(chans_to_analyze)
    if ~any(strcmpi({valid_erp.chanlocs.labels}, chans_to_analyze{c}))
        error('Channel %s not found. Available labels: %s', chans_to_analyze{c}, strjoin({valid_erp.chanlocs.labels}, ', '));
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
n1_amp_ref = nan(nS, nC);

bl_idx   = time_ms >= baseline_window_ms(1) & time_ms <= baseline_window_ms(2);
ref_idx  = time_ms >= ref_window_ms(1)      & time_ms <= ref_window_ms(2);
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

        % Secondary amplitude: re-referenced to -500..-200 ms (Toledo-style
        % proxy). A constant offset cannot move the peak, so latency is unchanged.
        n1_amp_ref(i, c) = find_n1_peak(tr - mean(tr(ref_idx)), time_ms, n1_window_ms, nbhd_samples);
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

%% STATS: Control vs pDCD at every electrode (Welch t, Cohen's d + 95% CI, rank-sum)
measure_names = {'Amplitude_uV', 'Latency_ms', 'AmplitudeRefBL_uV'};
measure_label = {'N1 amplitude (uV)', 'N1 latency (ms)', ...
                 sprintf('N1 amplitude re-referenced to %d..%d ms (uV)', ref_window_ms)};
measure_data  = {n1_amp, n1_lat, n1_amp_ref};

st_list = {};
for m = 1:numel(measure_names)
    for c = 1:nC
        x1 = measure_data{m}(included & is_control, c);  x1 = x1(~isnan(x1));
        x2 = measure_data{m}(included & is_pdcd,    c);  x2 = x2(~isnan(x2));
        st_list{end+1} = group_stats(measure_names{m}, chans_to_analyze{c}, x1, x2, alpha); %#ok<SAGROW>
    end
end
stats_tbl = struct2table([st_list{:}]);

% Holm correction across the electrodes, within each measure
stats_tbl.p_Holm = nan(height(stats_tbl), 1);
for m = 1:numel(measure_names)
    idx = strcmp(stats_tbl.Measure, measure_names{m});
    stats_tbl.p_Holm(idx) = holm_adjust(stats_tbl.p_Welch(idx));
end

fprintf('\n=====================================================================================================\n');
fprintf(' CONTROL vs pDCD, Bin %d, N1 window %d-%d ms  (d = Control - pDCD, pooled SD; CI via noncentral t)\n', bin_to_plot, n1_window_ms);
fprintf('=====================================================================================================\n');
for m = 1:numel(measure_names)
    fprintf('\n--- %s ---\n', measure_label{m});
    fprintf('%-5s | %-20s | %-20s | %-12s %-8s %-8s | %-22s | %-8s\n', ...
        'Elec', 'Control M (SD) n', 'pDCD M (SD) n', 'Welch t(df)', 'p', 'p_Holm', 'Cohen''s d [95% CI]', 'p_rank');
    idx = find(strcmp(stats_tbl.Measure, measure_names{m}));
    for k = idx'
        r = stats_tbl(k, :);
        fprintf('%-5s | %7.2f (%5.2f) %3d | %7.2f (%5.2f) %3d | %5.2f(%4.1f) %8.4f %8.4f | %5.2f [%5.2f, %5.2f] | %8.4f\n', ...
            r.Electrode{1}, r.M_Control, r.SD_Control, r.n_Control, r.M_pDCD, r.SD_pDCD, r.n_pDCD, ...
            r.t, r.df, r.p_Welch, r.p_Holm, r.d, r.d_CI_lo, r.d_CI_hi, r.p_ranksum);
    end
end
fprintf('\nNormality (Lilliefors) p-values are stored in stats_tbl (pNorm_Control, pNorm_pDCD).\n');
fprintf('Hedges'' g (small-sample corrected d) is stored in stats_tbl.g.\n');

%% STATS: Group x Electrode mixed ANOVA (does the group effect differ between sites?)
for m = 1:2
    fprintf('\n--- Mixed ANOVA, %s: Group (between) x Electrode (within) ---\n', measure_label{m});
    try
        wide = measure_data{m}(included, :);
        grp  = categorical(all_subj_group(included));
        ok   = all(~isnan(wide), 2);
        varnames = matlab.lang.makeValidName(chans_to_analyze);
        T = array2table(wide(ok, :), 'VariableNames', varnames);
        T.Group = grp(ok);
        within = table(categorical(chans_to_analyze(:)), 'VariableNames', {'Electrode'});
        rm = fitrm(T, sprintf('%s-%s ~ Group', varnames{1}, varnames{end}), 'WithinDesign', within);
        ranova_tbl = ranova(rm, 'WithinModel', 'Electrode');
        disp(ranova_tbl);
        fprintf('Group:Electrode row = whether the group difference depends on the site (use pValueGG).\n');
    catch ME
        fprintf('Mixed ANOVA could not be run: %s\n', ME.message);
    end
end

%% SAVE TABLES
try
    qc_tbl = table(strcat('SERP', all_subj_ids), all_subj_group, has_data, n_acc, n_rej, 100*rej_rate, ...
        included, excl_reason, 'VariableNames', {'Subject', 'Group', 'HasData', 'AcceptedTrials', ...
        'RejectedTrials', 'RejectPct', 'Included', 'ExclusionReason'});
    for c = 1:nC
        cn = matlab.lang.makeValidName(chans_to_analyze{c});
        qc_tbl.(['BaselineRMS_' cn])  = bl_rms(:, c);
        qc_tbl.(['PlusMinusRMS_' cn]) = pm_rms(:, c);
        qc_tbl.(['SNR_' cn])          = snr(:, c);
    end
    writetable(qc_tbl, fullfile(output_dir, 'QC_table.csv'));

    [ii, cc] = ndgrid(1:nS, 1:nC);
    ii = ii(:); cc = cc(:);
    long_tbl = table(strcat('SERP', all_subj_ids(ii)), all_subj_group(ii), chans_to_analyze(cc)', ...
        n1_amp(:), n1_lat(:), n1_amp_ref(:), n1_local(:), n1_edge(:), included(ii), ...
        'VariableNames', {'Subject', 'Group', 'Electrode', 'N1_Amplitude_uV', 'N1_Latency_ms', ...
        'N1_AmplitudeRefBL_uV', 'IsLocalPeak', 'AtWindowEdge', 'Included'});
    writetable(long_tbl, fullfile(output_dir, sprintf('N1_measures_%d-%dms.csv', n1_window_ms)));
    writetable(stats_tbl, fullfile(output_dir, sprintf('N1_group_stats_%d-%dms.csv', n1_window_ms)));
    fprintf('\nTables written to %s\n', output_dir);
catch ME
    warning('Could not write tables: %s', ME.message);
end

%% PLOT: 3x4 individual-subject grids for every electrode and group
if plot_individual_grids
    for c = 1:nC
        for g = 1:2
            ids = all_IDs{g};
            figure('Color', 'w', 'Position', [100, 100, 1200, 800], ...
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
                yl = ylim;
                hp = patch([n1_window_ms(1) n1_window_ms(2) n1_window_ms(2) n1_window_ms(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
                uistack(hp, 'bottom');
                line([0 0], yl, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
                line(xlim, [0 0], 'Color', 'k', 'LineWidth', 0.75);
                ylim(yl);
                ttl = sprintf('SERP%s', ids{s}); tcol = 'k';
                if ~included(i), ttl = sprintf('%s (%s)', ttl, excl_short{i}); tcol = [0.6 0.6 0.6]; end
                title(ttl, 'FontSize', 9, 'Color', tcol);
                hold off;
            end
            sgtitle(sprintf('%s - %s - Bin %d (shaded = N1 window %d-%d ms, o = N1)', ...
                groups{g}, chans_to_analyze{c}, bin_to_plot, n1_window_ms), 'FontWeight', 'bold', 'FontSize', 12);
        end
    end
end

%% PLOT: Grand averages after exclusions, all electrodes
figure('Color', 'w', 'Position', [50, 50, 1500, 850], 'Name', 'Grand averages after exclusions - all electrodes');
light_c = {control_light, pdcd_light};
dark_c  = {control_dark,  pdcd_dark};
for c = 1:nC
    subplot(2, 3, c); hold on;
    hl = []; lab = {};
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
    title(chans_to_analyze{c}, 'FontWeight', 'bold');
    if c == 1 && ~isempty(hl), legend(hl, lab, 'Location', 'southwest', 'FontSize', 8); end
    hold off;
end
sgtitle(sprintf('Bin %d grand averages after exclusions (x = N1 of grand average, shaded = %d-%d ms)', ...
    bin_to_plot, n1_window_ms), 'FontWeight', 'bold');

%% PLOT: Effect sizes with 95% CI across electrodes
figure('Color', 'w', 'Position', [200, 200, 1000, 420], 'Name', 'Cohen''s d by electrode');
forest_note = {'d > 0: Control less negative (pDCD larger N1)', 'd > 0: Control later than pDCD'};
for m = 1:2
    subplot(1, 2, m); hold on;
    idx = find(strcmp(stats_tbl.Measure, measure_names{m}));
    d  = stats_tbl.d(idx);
    lo = stats_tbl.d_CI_lo(idx);
    hi = stats_tbl.d_CI_hi(idx);
    y  = (1:numel(idx))';
    errorbar(d, y, d - lo, hi - d, 'horizontal', 'o', 'Color', 'k', ...
        'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 8);
    line([0 0], [0.5 numel(idx) + 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    set(gca, 'YTick', y, 'YTickLabel', stats_tbl.Electrode(idx), 'YDir', 'reverse', 'FontSize', 10);
    ylim([0.5 numel(idx) + 0.5]); grid on;
    xlabel('Cohen''s d (Control - pDCD), 95% CI');
    title({measure_label{m}, forest_note{m}}, 'FontSize', 10);
    hold off;
end

fprintf('\nDone.\n');

%% =====================================================================
%  LOCAL FUNCTIONS (must stay at the end of the script)
%  =====================================================================
function [peak_amp, peak_lat, is_local, at_edge] = find_n1_peak(trace, t, window, nbhd)
% FIND_N1_PEAK  N1 = most negative LOCAL peak within window (ms).
%   A sample k is a local peak if trace(k) is lower than all samples in
%   k-nbhd..k-1 and not higher than all samples in k+1..k+nbhd (neighbours
%   may lie outside the window). If no local peak exists, the absolute
%   minimum in the window is returned and is_local = false.
%   nbhd = 0 reproduces the old absolute-minimum behaviour.
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

function [d, ci, g] = cohens_d_ci(x1, x2, alpha)
% COHENS_D_CI  Cohen's d (pooled SD, n-weighted) with a (1-alpha) CI from
%   the noncentral t distribution, plus Hedges' g. Falls back to the
%   large-sample normal approximation if nctcdf is unavailable.
    n1 = numel(x1); n2 = numel(x2); df = n1 + n2 - 2;
    sp = sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / df);
    d  = (mean(x1) - mean(x2)) / sp;
    g  = d * (1 - 3 / (4 * df - 1));
    k  = sqrt(n1 * n2 / (n1 + n2));
    t_obs = d * k;
    ci = [NaN NaN];
    try
        f_lo = @(ncp) nctcdf(t_obs, df, ncp) - (1 - alpha / 2);
        f_hi = @(ncp) nctcdf(t_obs, df, ncp) - alpha / 2;
        w = 10;
        while f_lo(t_obs - w) < 0 || f_hi(t_obs + w) > 0
            w = w * 2;
            if w > 1e4, error('bracket'); end
        end
        ci = [fzero(f_lo, [t_obs - w, t_obs + w]), fzero(f_hi, [t_obs - w, t_obs + w])] / k;
    catch
        se = sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)));
        z  = sqrt(2) * erfinv(1 - alpha);
        ci = [d - z * se, d + z * se];
    end
end

function p_adj = holm_adjust(p)
% HOLM_ADJUST  Holm-Bonferroni adjusted p-values (NaNs ignored).
    p = p(:);
    p_adj = nan(size(p));
    v = find(~isnan(p));
    m = numel(v);
    if m == 0, return; end
    [ps, ord] = sort(p(v));
    adj = min(1, (m - (1:m)' + 1) .* ps);
    adj = cummax(adj);
    p_adj(v(ord)) = adj;
end

function st = group_stats(measure, electrode, x1, x2, alpha)
% GROUP_STATS  Descriptives, Welch t, Cohen's d + CI, rank-sum, normality.
    st.Measure       = measure;
    st.Electrode     = electrode;
    st.n_Control     = numel(x1);
    st.M_Control     = mean(x1);
    st.SD_Control    = std(x1);
    st.n_pDCD        = numel(x2);
    st.M_pDCD        = mean(x2);
    st.SD_pDCD       = std(x2);
    st.Diff          = st.M_Control - st.M_pDCD;
    st.Diff_CI_lo    = NaN;
    st.Diff_CI_hi    = NaN;
    st.t             = NaN;
    st.df            = NaN;
    st.p_Welch       = NaN;
    st.d             = NaN;
    st.d_CI_lo       = NaN;
    st.d_CI_hi       = NaN;
    st.g             = NaN;
    st.p_ranksum     = NaN;
    st.pNorm_Control = NaN;
    st.pNorm_pDCD    = NaN;
    if numel(x1) < 2 || numel(x2) < 2, return; end

    [~, p, ci, s] = ttest2(x1, x2, 'Vartype', 'unequal', 'Alpha', alpha);
    st.t = s.tstat; st.df = s.df; st.p_Welch = p;
    st.Diff_CI_lo = ci(1); st.Diff_CI_hi = ci(2);

    [st.d, dci, st.g] = cohens_d_ci(x1, x2, alpha);
    st.d_CI_lo = dci(1); st.d_CI_hi = dci(2);

    st.p_ranksum = ranksum(x1, x2);
    try, [~, st.pNorm_Control] = lillietest(x1); catch, end
    try, [~, st.pNorm_pDCD]    = lillietest(x2); catch, end
end