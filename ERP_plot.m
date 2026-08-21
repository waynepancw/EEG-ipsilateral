% 1. Force MATLAB to open a clean, standard figure window
figure; 

% 2. Extract the time points (-500 to 1000 ms)
time_ms = ERP.times; 

% 3. Find the index numbers for Fz, FCz, and Cz dynamically
chan_idx = [];
for ch = {'Fz', 'FCz', 'Cz'}
    idx = find(strcmpi({ERP.chanlocs.labels}, ch{1}));
    if ~isempty(idx), chan_idx(end+1) = idx; end
end

% 4. Plot the 4 bins for the first found midline channel (usually Fz or Cz)
% This extracts: bindata(channel, time_points, bin)
plot(time_ms, squeeze(ERP.bindata(chan_idx(1), :, :)), 'LineWidth', 2);

% 5. Format the graph beautifully
xlim([-200 800]);
grid on;
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); % Stimulus onset line
line(xlim, [0 0], 'Color', 'k'); % Zero baseline line
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title(['ERP Waveforms at Channel: ' ERP.chanlocs(chan_idx(1)).labels]);
legend({'Bin 1', 'Bin 2', 'Bin 3', 'Bin 4'}, 'Location', 'best');
figure;
time_ms = ERP.times;

% Target the exact three channels from the paper
target_chans = {'Fz', 'FCz', 'Cz'};

for i = 1:length(target_chans)
    % Find channel index
    ch_idx = find(strcmpi({ERP.chanlocs.labels}, target_chans{i}));

    if ~isempty(ch_idx)
        % Create a stacked subplot (3 rows, 1 column, current position i)
        subplot(3, 1, i); 

        % Plot all 4 bins
        plot(time_ms, squeeze(ERP.bindata(ch_idx, :, :)), 'LineWidth', 1.5);

        % Styling
        xlim([-200 800]);
        grid on;
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); % Stimulus onset
        line(xlim, [0 0], 'Color', 'k'); % Zero line

        ylabel('Amplitude (\muV)');
        title(['Channel: ' target_chans{i}]);

        if i == 3
            xlabel('Time (ms)');
            legend({'Bin 1 (Trig 19)', 'Bin 2 (Trig 27)', 'Bin 3 (Trig 39)', 'Bin 4 (Trig 47)'}, 'Location', 'best');
        end
    end
end