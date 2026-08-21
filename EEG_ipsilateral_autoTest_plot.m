%% 1. INITIALIZATION & USER SETTINGS
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % Open EEGLAB background workspace

% === IMPORT FILES ===
filepath = 'C:\Users\cwpan\Desktop\EEG-ipsilateral position\data';
filename_control = 'SERP40_position_45-52_單_proc_convert.cdt.cnt'; % Control file name
filename_pdcd = 'SERP36_position_45-52_單_proc_convert.cdt.cnt';       % pDCD file name
% ====================

files_to_process = {filename_control, filename_pdcd};
group_names = {'Control', 'pDCD'};
ERP_results = cell(1, 2); % This will store the final averages for both subjects

%% 2. CREATE THE BINS FILE ONCE
bin_file_path = fullfile(filepath, 'autobins.txt');
fid = fopen(bin_file_path, 'w');
fprintf(fid, 'bin 1\nTrigger 19\n.{19}\n\nbin 2\nTrigger 27\n.{27}\n\nbin 3\nTrigger 39\n.{39}\n\nbin 4\nTrigger 47\n.{47}\n');
fclose(fid);

%% 3. AUTOMATED PROCESSING LOOP FOR BOTH SUBJECTS
for s = 1:2
    fprintf('\n==========================================\n');
    fprintf('Processing %s Dataset: %s\n', group_names{s}, files_to_process{s});
    fprintf('==========================================\n');

    % A. Load Raw .cnt Data
    EEG = pop_loadcnt(fullfile(filepath, files_to_process{s}), 'dataformat', 'auto', 'memmapfile', '');

    % B. Create EventList
    EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});

    % C. Assign Bins
    EEG = pop_binlister(EEG, 'BDF', bin_file_path, 'IndexEL', 1, 'SendEL2', 'EEG', 'Voutput', 'EEG');

    % D. Extract Epochs & Baseline Correct (-200 to 0 ms)
    EEG = pop_epochbin(EEG, [-500.0  1000.0], [-200 0]);

    % E. Artifact Rejection (Threshold \pm100 \muV)
    EEG = pop_artextval(EEG, 'Channel', 1:34, 'Flag', 1, 'LowPass', -1, 'Threshold', [-100 100], 'Twindow', [-500 1000]);

    % F. Compute Averaged ERP
    ERP_results{s} = pop_averager(EEG, 'Criterion', 'good', 'DQ_flag', 0); 
end

%% 4. GENERATE "FIGURE 4" (OVERLAYING GROUPS FOR THE VR CONDITION)
fprintf('\nGenerating Figure 4 (VR Condition Comparison)...\n');
figure('Color', 'w', 'Position', [100, 100, 600, 800]); % Force outer background white
time_ms = ERP_results{1}.times;
target_chans = {'Fz', 'FCz', 'Cz'};

for i = 1:length(target_chans)
    subplot(3, 1, i); 
    hold on;
    
    % Find the channel index for this specific electrode
    ch_idx = find(strcmpi({ERP_results{1}.chanlocs.labels}, target_chans{i}));
    
    if ~isempty(ch_idx)
        % Plot CONTROL subject in BLUE
        plot(time_ms, squeeze(ERP_results{1}.bindata(ch_idx, :, 1)), 'b', 'LineWidth', 2);
        
        % Plot pDCD subject in RED
        plot(time_ms, squeeze(ERP_results{2}.bindata(ch_idx, :, 1)), 'r', 'LineWidth', 2);
        
        % Force inner graph styling to scientific standard (Black & White)
        xlim([-200 800]);
        grid on;
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'FontSize', 11);
        
        % Black reference lines
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Stimulus line
        line(xlim, [0 0], 'Color', 'k', 'LineWidth', 1); % Zero line
        
        ylabel('Amplitude (\muV)', 'Color', 'k', 'FontWeight', 'bold');
        title(sprintf('Electrode: %s', target_chans{i}), 'Color', 'k', 'FontWeight', 'bold');
        
        % Add legend only to the bottom plot
        if i == 3
            xlabel('Time (ms)', 'Color', 'k', 'FontWeight', 'bold');
            legend({'Control', 'pDCD'}, 'Location', 'southwest', 'TextColor', 'k', 'Color', 'w');
        end
    end
end
hold off;