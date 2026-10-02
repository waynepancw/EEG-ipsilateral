clear all, close all, clc

%% =====================================================================
%  N1 statistics: Control vs pDCD, passive left-foot stimulation
%  Input : N1_peaks_bin1_0-200ms.xlsx (made by EEG_ipsilateral_allSubjects_N1.m)
%            sheet "N1_peaks"         -> N1 amplitude & latency per electrode
%            sheet "Inclusion_and_QC" -> which subjects are included
%  Tests :
%   1. Mixed ANOVA for each N1 attribute (amplitude, latency):
%        Group (between: Control, pDCD) x Electrode (within: 6 electrodes)
%   2. 6 x 2 matrix (rows = electrodes, columns = amplitude, latency):
%        one-way ANOVA Control vs pDCD in every cell, with Holm correction
%        over the 6 electrodes within each attribute
%   3. Same 6 x 2 matrix with the Mann-Whitney U test (non-parametric,
%        rank-based; PI's suggestion), also Holm-corrected
%  Output: tables in the Command Window + N1_stats_bin1_0-200ms.xlsx
%  Requires the Statistics and Machine Learning Toolbox.
%  =====================================================================

%% SETTINGS
data_dir = '/home/cwpan/Documents/EEG-ipsilateral/data/N1_results';
xlsx_in  = fullfile(data_dir, 'N1_peaks_bin1_0-200ms.xlsx');
xlsx_out = fullfile(data_dir, 'N1_stats_bin1_0-200ms.xlsx');

chans      = {'Fz', 'FCz', 'Cz', 'Pz', 'C3', 'C4'};   % rows of the 6x2 matrix
attrs      = {'Amplitude_uV', 'Latency_ms'};          % columns of the 6x2 matrix
attr_names = {'Amplitude (uV)', 'Latency (ms)'};
groups     = {'Control', 'pDCD'};
alpha      = 0.05;

%% LOAD DATA AND KEEP INCLUDED SUBJECTS ONLY
peaks = readtable(xlsx_in, 'Sheet', 'N1_peaks');
info  = readtable(xlsx_in, 'Sheet', 'Inclusion_and_QC');

% "Included" may come back as 1/0 or as TRUE/FALSE text depending on Excel
inc = info.Included;
if iscell(inc) || isstring(inc)
    inc = strcmpi(string(inc), 'true') | strcmp(string(inc), '1');
else
    inc = logical(inc);
end

[found, loc] = ismember(peaks.Subject, info.Subject);
if ~all(found)
    error('Subjects in "N1_peaks" and "Inclusion_and_QC" do not match.');
end
d = peaks(inc(loc), :);
d.Group = categorical(d.Group, groups);

fprintf('Subjects used in the statistics:\n');
for g = 1:2
    s = d.Subject(d.Group == groups{g});
    fprintf('  %-8s n = %2d : %s\n', groups{g}, numel(s), strjoin(s', ', '));
end
fprintf('(Excluded subjects and reasons: sheet "Inclusion_and_QC".)\n');

%% 1. MIXED ANOVA: Group (between) x Electrode (within), per attribute
within = table(categorical(chans', chans), 'VariableNames', {'Electrode'});
mixed_rows = {};

for a = 1:numel(attrs)
    vars = strcat(chans, '_', attrs{a});
    Y  = d{:, vars};
    ok = all(~isnan(Y), 2);     % repeated measures need all 6 electrodes
    if sum(~ok) > 0
        fprintf('\n%s: %d subject(s) dropped from the mixed ANOVA (missing electrode).\n', ...
            attr_names{a}, sum(~ok));
    end
    dd = d(ok, [{'Group'}, vars]);

    rm = fitrm(dd, sprintf('%s-%s ~ Group', vars{1}, vars{end}), 'WithinDesign', within);

    % Between-subject effect: Group
    bt = anova(rm);
    rG = find(string(bt.Between) == "Group", 1);
    rE = find(string(bt.Between) == "Error", 1);
    mixed_rows(end+1, :) = {attr_names{a}, 'Group', bt.F(rG), bt.DF(rG), bt.DF(rE), ...
        bt.pValue(rG), NaN, bt.SumSq(rG) / (bt.SumSq(rG) + bt.SumSq(rE))}; %#ok<SAGROW>

    % Within-subject effects: Electrode and Group x Electrode
    % (p_GG = Greenhouse-Geisser corrected, use it if sphericity is violated)
    wt  = ranova(rm);
    rErr = find(contains(wt.Properties.RowNames, 'Error'), 1);
    eff = {'(Intercept):Electrode', 'Electrode'; 'Group:Electrode', 'Group x Electrode'};
    for k = 1:size(eff, 1)
        r = find(strcmp(wt.Properties.RowNames, eff{k, 1}), 1);
        mixed_rows(end+1, :) = {attr_names{a}, eff{k, 2}, wt.F(r), wt.DF(r), wt.DF(rErr), ...
            wt.pValue(r), wt.pValueGG(r), wt.SumSq(r) / (wt.SumSq(r) + wt.SumSq(rErr))}; %#ok<SAGROW>
    end

    mt = mauchly(rm);
    fprintf('%s: Mauchly sphericity test p = %.3f\n', attr_names{a}, mt.pValue);
end

mixed_tbl = cell2table(mixed_rows, 'VariableNames', ...
    {'Attribute', 'Effect', 'F', 'df1', 'df2', 'p', 'p_GG', 'partial_eta2'});

fprintf('\n=====================================================================\n');
fprintf(' 1. MIXED ANOVA: Group (between) x Electrode (within)\n');
fprintf('=====================================================================\n');
disp(mixed_tbl);

%% 2. 6 x 2 MATRIX: Control vs pDCD at every electrode, for each attribute
P      = nan(numel(chans), numel(attrs));   % uncorrected p
P_holm = nan(numel(chans), numel(attrs));   % Holm-corrected over the 6 electrodes
det_rows = {};

for a = 1:numel(attrs)
    for c = 1:numel(chans)
        y  = d.([chans{c} '_' attrs{a}]);
        ok = ~isnan(y);
        yc = y(ok & d.Group == 'Control');
        yp = y(ok & d.Group == 'pDCD');

        [p, tbl] = anova1(y(ok), cellstr(d.Group(ok)), 'off');
        P(c, a) = p;
        F   = tbl{2, 5};  df1 = tbl{2, 3};  df2 = tbl{3, 3};
        eta2 = tbl{2, 2} / tbl{4, 2};

        % Cohen's d (pDCD - Control) with pooled SD
        sp = sqrt(((numel(yc)-1)*var(yc) + (numel(yp)-1)*var(yp)) / (numel(yc) + numel(yp) - 2));
        cohen_d = (mean(yp) - mean(yc)) / sp;

        % Levene's test: equal variances assumption of the ANOVA
        p_lev = vartestn(y(ok), cellstr(d.Group(ok)), 'TestType', 'LeveneAbsolute', 'Display', 'off');

        det_rows(end+1, :) = {chans{c}, attr_names{a}, numel(yc), mean(yc), std(yc), ...
            numel(yp), mean(yp), std(yp), F, df1, df2, p, NaN, eta2, cohen_d, p_lev}; %#ok<SAGROW>
    end
    P_holm(:, a) = holm_correct(P(:, a));
end

det_tbl = cell2table(det_rows, 'VariableNames', {'Electrode', 'Attribute', ...
    'n_Control', 'Mean_Control', 'SD_Control', 'n_pDCD', 'Mean_pDCD', 'SD_pDCD', ...
    'F', 'df1', 'df2', 'p', 'p_Holm', 'eta2', 'Cohen_d', 'p_Levene'});
det_tbl.p_Holm = P_holm(:);   % same order as the loop (attribute, then electrode)

fprintf('\n=====================================================================\n');
fprintf(' 2. CONTROL vs pDCD at each electrode (one-way ANOVA)\n');
fprintf('    * p < %.2f after Holm correction over the 6 electrodes\n', alpha);
fprintf('=====================================================================\n');
for a = 1:numel(attrs)
    fprintf('\n--- N1 %s ---\n', attr_names{a});
    fprintf('%-6s %18s %18s %8s %8s %8s\n', 'Elec', 'Control M (SD)', 'pDCD M (SD)', 'F', 'p', 'p_Holm');
    for c = 1:numel(chans)
        r = (a - 1) * numel(chans) + c;
        star = ''; if P_holm(c, a) < alpha, star = ' *'; end
        fprintf('%-6s %10.2f (%5.2f) %10.2f (%5.2f) %8.2f %8.3f %8.3f%s\n', chans{c}, ...
            det_tbl.Mean_Control(r), det_tbl.SD_Control(r), det_tbl.Mean_pDCD(r), det_tbl.SD_pDCD(r), ...
            det_tbl.F(r), P(c, a), P_holm(c, a), star);
    end
end

%% 3. MANN-WHITNEY U TEST: Control vs pDCD at every electrode (non-parametric)
% Compares the RANKS of the values instead of the means, so it does not
% assume normal distributions or equal variances and is robust to outliers.
% U = the smaller of the two U values; exact p-values (small samples).
% r_rb = rank-biserial correlation (-1..1); positive = pDCD values larger
% (for amplitude: pDCD N1 less negative, i.e. smaller), same sign as Cohen_d.
P_mw      = nan(numel(chans), numel(attrs));
P_mw_holm = nan(numel(chans), numel(attrs));
mw_rows = {};

for a = 1:numel(attrs)
    for c = 1:numel(chans)
        y  = d.([chans{c} '_' attrs{a}]);
        ok = ~isnan(y);
        yc = y(ok & d.Group == 'Control');
        yp = y(ok & d.Group == 'pDCD');
        n1 = numel(yc);  n2 = numel(yp);

        [p, ~, st] = ranksum(yc, yp, 'method', 'exact');
        U_control = st.ranksum - n1 * (n1 + 1) / 2;
        U    = min(U_control, n1 * n2 - U_control);
        r_rb = 1 - 2 * U_control / (n1 * n2);
        P_mw(c, a) = p;

        mw_rows(end+1, :) = {chans{c}, attr_names{a}, n1, median(yc), iqr(yc), ...
            n2, median(yp), iqr(yp), U, p, NaN, r_rb}; %#ok<SAGROW>
    end
    P_mw_holm(:, a) = holm_correct(P_mw(:, a));
end

mw_tbl = cell2table(mw_rows, 'VariableNames', {'Electrode', 'Attribute', ...
    'n_Control', 'Median_Control', 'IQR_Control', 'n_pDCD', 'Median_pDCD', 'IQR_pDCD', ...
    'U', 'p', 'p_Holm', 'r_rankbiserial'});
mw_tbl.p_Holm = P_mw_holm(:);

fprintf('\n=====================================================================\n');
fprintf(' 3. CONTROL vs pDCD at each electrode (Mann-Whitney U test)\n');
fprintf('    * p < %.2f after Holm correction over the 6 electrodes\n', alpha);
fprintf('=====================================================================\n');
for a = 1:numel(attrs)
    fprintf('\n--- N1 %s ---\n', attr_names{a});
    fprintf('%-6s %20s %20s %6s %8s %8s %7s\n', 'Elec', 'Control Mdn [IQR]', 'pDCD Mdn [IQR]', ...
        'U', 'p', 'p_Holm', 'r_rb');
    for c = 1:numel(chans)
        r = (a - 1) * numel(chans) + c;
        star = ''; if P_mw_holm(c, a) < alpha, star = ' *'; end
        fprintf('%-6s %11.2f [%6.2f] %11.2f [%6.2f] %6.1f %8.3f %8.3f %7.2f%s\n', chans{c}, ...
            mw_tbl.Median_Control(r), mw_tbl.IQR_Control(r), mw_tbl.Median_pDCD(r), mw_tbl.IQR_pDCD(r), ...
            mw_tbl.U(r), P_mw(c, a), P_mw_holm(c, a), mw_tbl.r_rankbiserial(r), star);
    end
end

% The 6 x 2 matrix itself (ANOVA and Mann-Whitney side by side)
matrix_tbl = table(chans', P(:, 1), P_holm(:, 1), P_mw(:, 1), P_mw_holm(:, 1), ...
    P(:, 2), P_holm(:, 2), P_mw(:, 2), P_mw_holm(:, 2), 'VariableNames', ...
    {'Electrode', 'Amp_ANOVA_p', 'Amp_ANOVA_p_Holm', 'Amp_MW_p', 'Amp_MW_p_Holm', ...
    'Lat_ANOVA_p', 'Lat_ANOVA_p_Holm', 'Lat_MW_p', 'Lat_MW_p_Holm'});
fprintf('\n--- 6 x 2 matrix of p-values (rows = electrodes) ---\n');
disp(matrix_tbl);

%% EXPORT TO EXCEL
try
    if exist(xlsx_out, 'file'), delete(xlsx_out); end   % no leftovers from an older run
    writetable(matrix_tbl, xlsx_out, 'Sheet', 'Matrix_6x2_p');
    writetable(det_tbl,    xlsx_out, 'Sheet', 'Per_electrode_details');
    writetable(mw_tbl,     xlsx_out, 'Sheet', 'MannWhitney_details');
    writetable(mixed_tbl,  xlsx_out, 'Sheet', 'Mixed_ANOVA');
    writetable(d(:, {'Subject', 'Group'}), xlsx_out, 'Sheet', 'Subjects_used');
    fprintf('\nStatistics written to %s\n', xlsx_out);
catch ME
    warning('Could not write %s (is it open in Excel?): %s', xlsx_out, ME.message);
end

fprintf('\nDone.\n');


%% =====================================================================
%  LOCAL FUNCTIONS (must stay at the end of the script)
%  =====================================================================
function p_adj = holm_correct(p)
% HOLM_CORRECT  Holm-Bonferroni adjusted p-values (NaNs are ignored).
    p_adj = nan(size(p));
    ok = find(~isnan(p));
    [ps, order] = sort(p(ok));
    m = numel(ps);
    adj = min(1, cummax((m - (1:m)' + 1) .* ps(:)));
    p_adj(ok(order)) = adj;
end