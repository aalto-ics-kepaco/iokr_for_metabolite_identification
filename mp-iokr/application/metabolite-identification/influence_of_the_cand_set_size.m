close all;

%% Load data
outputDirMP__ = ...
        '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/commit_431c5b22cce549b07374677cbd514ade2c2b63ed/';

inclExpCand__ = false;  

% ALL CANDIDATES (MP_IOKR)
ind_fold__ = load (strcat (inputDir, 'cv_ind.txt'));
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold__), ...
                     'inner', struct ('nFolds', 10));                  

[rank_perc_all_OLD__, ranks_all_OLD__] = aggregate_results_all ('UNIMKL', outputDirMP__, cv_param__, inclExpCand__);
% NEW corresponds to the results which have been calcualted after the
% bug-fix related to the centering of candidate-feature vectors in the
% inner cross-validation. The fix basically did not effect the results.
[rank_perc_all__, ranks_all__]         = aggregate_results_all ('UNIMKL', '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/', cv_param__, inclExpCand__);

% ALL CANDIDATES SEPARATE (MP_IOKR)            
[rank_perc_all_sep_OLD__, ranks_all_sep_OLD__] = aggregate_results_all_separate ('UNIMKL', outputDirMP__, cv_param__, inclExpCand__);
[rank_perc_all_sep__, ranks_all_sep__]         = aggregate_results_all_separate ('UNIMKL', '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/', cv_param__, inclExpCand__);

% RANDOM CANDIDATES (MP-IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold__));     
perc__ = 1;
% 
[rank_perc_rand__, ranks_rand__, cand_num_rand__] = aggregate_results_random ( ...
    'UNIMKL', outputDirMP__, cv_param__, perc__, inclExpCand__);

% RANDOM CANDIDATES (MP-IOKR)
[rank_perc_rand_sep__, ranks_rand_sep__, cand_num_rand_sep__] = aggregate_results_random_separate ( ...
    'UNIMKL', outputDirMP__, cv_param__, perc__, inclExpCand__);

% REFERENCE (IOKR)
tmp__ = load ( ...
    '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/iokr/a143cc0ca62fc3c75e383cb4fb679201.mat');
ranks_ref__     = tmp__.result.ranks;
rank_perc_ref__ = tmp__.result.rank_perc;

% CFM-ID
cand_num_eval__ = load ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/candidates/num_cand_eval.txt');
dt_cfm__ = readtable ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results_cfm/res_precomp/rank_for_each_MSID.txt', ...
    'Delimiter', ' ', 'ReadVariableNames', false, 'HeaderLines', 0);
% The ranks for each molecular are in the second column
ranks_cfm__ = table2array (dt_cfm__(:, 2));
% Calculate the top-k performance for the CFM method
rank_perc_cfm__ = getRankPerc (ranks_cfm__, cand_num_eval__);

% CSI-FingerID
ranks_csifi__ = load ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results_csifi/rank_csifingerid_final_version_eval.txt');
rank_perc_csifi__ = getRankPerc (ranks_csifi__, cand_num_eval__);


fprintf ('IOKR\tMP-IOKR (all)\tMP-IOKR (all) OLD\tMP-IOKR separate (all)\tMP-IOKR separate (all) OLD\tMP-IOKR (random) [mean]\tCFM\tCSI-FingerID\n');
disp ([rank_perc_ref__(1:20),               ...
    rank_perc_all__(1:20),                  ...
    rank_perc_all_OLD__(1:20),              ...
    rank_perc_all_sep__(1:20),              ...
    rank_perc_all_sep_OLD__(1:20),          ...
    mean(rank_perc_rand__(1:20, :), 2),     ...
    mean(rank_perc_rand_sep__(1:20, :), 2), ...
    rank_perc_cfm__(1:20),                  ...
    rank_perc_csifi__(1:20)]);

fprintf ('IOKR: top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_ref__(1),  2), ...
    round (rank_perc_ref__(10), 2), ...
    round (rank_perc_ref__(20), 2));
fprintf ('MP-IOKR (all, joint): top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_all__(1),  2), ...
    round (rank_perc_all__(10), 2), ...
    round (rank_perc_all__(20), 2));
fprintf ('MP-IOKR (all, joint) OLD: top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_all_OLD__(1),  2), ...
    round (rank_perc_all_OLD__(10), 2), ...
    round (rank_perc_all_OLD__(20), 2));
fprintf ('MP-IOKR (all, separate): top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_all_sep__(1),  2), ...
    round (rank_perc_all_sep__(10), 2), ...
    round (rank_perc_all_sep__(20), 2));
fprintf ('MP-IOKR (all, separate) OLD: top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_all_sep_OLD__(1),  2), ...
    round (rank_perc_all_sep_OLD__(10), 2), ...
    round (rank_perc_all_sep_OLD__(20), 2));
fprintf ('MP-IOKR (random, joint): top-1 = %.2f%% (sd = %.1f) ; top-10 = %.2f%% (sd = %.1f) ; top-20 = %.2f%% (sd = %.1f)\n', ...
    round (mean(rank_perc_rand__(1,  :), 2), 2), round (std(rank_perc_rand__(1,  :), 0, 2), 2), ...
    round (mean(rank_perc_rand__(10, :), 2), 2), round (std(rank_perc_rand__(10, :), 0, 2), 2), ...
    round (mean(rank_perc_rand__(20, :), 2), 2), round (std(rank_perc_rand__(20, :), 0, 2), 2));
fprintf ('MP-IOKR (random, separate): top-1 = %.2f%% (sd = %.1f) ; top-10 = %.2f%% (sd = %.1f) ; top-20 = %.2f%% (sd = %.1f)\n', ...
    round (mean(rank_perc_rand_sep__(1,  :), 2), 2), round (std(rank_perc_rand_sep__(1,  :), 0, 2), 2), ...
    round (mean(rank_perc_rand_sep__(10, :), 2), 2), round (std(rank_perc_rand_sep__(10, :), 0, 2), 2), ...
    round (mean(rank_perc_rand_sep__(20, :), 2), 2), round (std(rank_perc_rand_sep__(20, :), 0, 2), 2));
fprintf ('CFM: top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_cfm__(1),  2), ...
    round (rank_perc_cfm__(10), 2), ...
    round (rank_perc_cfm__(20), 2));
fprintf ('CSI-FingerID: top-1 = %.2f%% ; top-10 = %.2f%% ; top-20 = %.2f%%\n', ...
    round (rank_perc_csifi__(1),  2), ...
    round (rank_perc_csifi__(10), 2), ...
    round (rank_perc_csifi__(20), 2));

%% Rank improvement
assert (all (isnan (ranks_ref__) == isnan (ranks_all__)));
assert (all (isnan (ranks_ref__) == isnan (ranks_all_sep__)));
assert (all (isnan (ranks_ref__) == isnan (ranks_all_OLD__)));
assert (all (isnan (ranks_ref__) == isnan (ranks_all_sep_OLD__)));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_rand__, 2))));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_rand_sep__, 2))));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_maxE__, 2))));

has_rank__ = ~ isnan (ranks_ref__);

disp ('Signtest')
ranks_IOKR__            = ranks_ref__(has_rank__); 
ranks_MP_all__          = ranks_all__(has_rank__);
ranks_MP_all_sep__      = ranks_all_sep__(has_rank__);
ranks_MP_rand__         = median (ranks_rand__(has_rank__, :), 2);
ranks_MP_rand_sep__     = median (ranks_rand_sep__(has_rank__, :), 2);

signf_mat__ = [ signtest(ranks_cfm__,          ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_cfm__,          ranks_MP_rand_sep__, 'tail', 'right');
                signtest(ranks_csifi__,        ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_csifi__,        ranks_MP_rand_sep__, 'tail', 'right');
                signtest(ranks_IOKR__,         ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_IOKR__,         ranks_MP_rand_sep__, 'tail', 'right');  
                signtest(ranks_MP_all__,       ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_MP_all__,       ranks_MP_rand_sep__, 'tail', 'right');
                signtest(ranks_MP_all_sep__,   ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_MP_all_sep__,   ranks_MP_rand_sep__, 'tail', 'right');
                signtest(ranks_MP_rand__,      ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_MP_rand__,      ranks_MP_rand_sep__, 'tail', 'right');
                signtest(ranks_MP_rand_sep__,  ranks_cfm__,         'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_csifi__,       'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_IOKR__,        'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_MP_all__,      'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_MP_all_sep__,  'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_MP_rand__,     'tail', 'right'),  ...
                signtest(ranks_MP_rand_sep__,  ranks_MP_rand_sep__, 'tail', 'right') ];


[p_all__, signf_all__]           = signtest (ranks_ref__(has_rank__), ranks_all__(has_rank__), 'tail', 'right');
[p_all_sep__, signf_all_sep__]   = signtest (ranks_ref__(has_rank__), ranks_all_sep__(has_rank__), 'tail', 'right');
[p_rand__, signf_rand__]         = signtest (ranks_ref__(has_rank__), ranks_rand__(has_rank__), 'tail', 'right');
[p_rand_sep__, signf_rand_sep__] = signtest (ranks_ref__(has_rank__), ranks_rand_sep__(has_rank__), 'tail', 'right');

fprintf ('MP-IOKR (all, joint): p-value = %e, signf = %d\n', p_all__, signf_all__);
fprintf ('MP-IOKR (all, separate): p-value = %e, signf = %d\n', p_all_sep__, signf_all_sep__);
fprintf ('MP-IOKR (random, joint): p-value = %e, signf = %d\n', p_rand__, signf_rand__);
fprintf ('MP-IOKR (random, separate): p-value = %e, signf = %d\n', p_rand_sep__, signf_rand_sep__);


% [cand_num_maxE__, IX_maxE__] = sort (cand_num_maxE__);

%% Plot top-k accuracy
maxRank__ = 100;
% 
% subplot (1, 2, 1);
% stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% stairs (rank_perc_all__(1:maxRank__), 'LineWidth', 1.5);
% stairs (rank_perc_all_sep__(1:maxRank__), 'LineWidth', 1.5);
% stairs (mean (rank_perc_rand__(1:maxRank__, :), 2), 'LineWidth', 1.5);
% stairs (mean (rank_perc_rand_sep__(1:maxRank__, :), 2), 'LineWidth', 1.5);
% stairs (mean (rank_perc_maxE__(1:maxRank__, :), 2), 'LineWidth', 1.5);
% legend ('IOKR', ...
%     'MP-IOKR (all)', 'MP-IOKR separate (all)', ...
%     'MP-IOKR (random)', 'MP-IOKR separate (random)', ...
%     'MP-IOKR (maxElement)');
% grid;
% 
% subplot (1, 2, 2);
% plot (1:maxRank__, zeros (1, maxRank__), 'k--');
% hold on;
% stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), 'LineWidth', 1.5); 
% stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), 'LineWidth', 1.5);
% stairs (mean (rank_perc_rand__(1:maxRank__, :), 2) - rank_perc_ref__(1:maxRank__), 'LineWidth', 1.5); 
% stairs (mean (rank_perc_rand_sep__(1:maxRank__, :), 2) - rank_perc_ref__(1:maxRank__), 'LineWidth', 1.5); 
% stairs (mean (rank_perc_maxE__(1:maxRank__, :), 2) - rank_perc_ref__(1:maxRank__), 'LineWidth', 1.5); 
% legend ('IOKR', ...
%     'MP-IOKR (all)', 'MP-IOKR separate (all)', ...
%     'MP-IOKR (random)', 'MP-IOKR separate (random)', ...
%     'MP-IOKR (maxElement)');
% 
% grid; 


colors__ =  [ ...
    0         0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840];
% col_rand__ = colors__(1, :);
col_all__     = colors__(1, :);
col_all_sep__ = colors__(2, :);
col_csifi__   = colors__(3, :);
col_cfm__     = colors__(4, :);
col_ref__     = 'black';

x_label__ = 'k';
y_label__ = 'top-k accuracy';

figure; 
% subplot (2, 2, 1);
% p_ref__ = stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% p_rand__ = stairs (mean (rank_perc_rand__(1:maxRank__, :), 2), ...
%     'Color', col_rand__, 'LineWidth', 1.5);
% p_all__ = stairs (rank_perc_all__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% title ('Joint');
% xlabel (x_label__);
% ylabel (y_label__);
% legend ([p_ref__, p_rand__, p_all__], 'IOKR', '1% random (mean)', 'all');
% grid;
p_ref__     = stairs (rank_perc_ref__(1:maxRank__), 'Color', col_ref__, 'LineWidth', 1.5); 
hold on;
p_all__     = stairs (rank_perc_all__(1:maxRank__), 'Color', col_all__, 'LineWidth', 1.5);
p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__), 'Color', col_all_sep__, 'LineWidth', 1.5);
p_csifi__   = stairs (rank_perc_csifi__(1:maxRank__), 'Color', col_csifi__, 'LineWidth', 1.5);
p_cfm__     = stairs (rank_perc_cfm__(1:maxRank__), 'Color', col_cfm__, 'LineWidth', 1.5);
ylim ([10, 90]);
title ('Joint');
xlabel (x_label__);
ylabel (y_label__);
legend ([p_ref__, p_all__, p_all_sep__, p_csifi__, p_cfm__], 'IOKR', 'all', 'all (separate)', 'CSI-FingerID', 'CFM-ID');
grid;

%%


figure; 

title__   = 'Performance differences in the metabolite-identification';
y_label__ = 'Relative change of top-k accuracy';
p_ref__ = plot (1:maxRank__, zeros (1, maxRank__), 'Color', col_ref__, 'LineWidth', 2);
hold on; 
p_all__ = stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_all__, 'LineWidth', 2); 
p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_all_sep__, 'LineWidth', 2);
p_csifi__ = stairs (rank_perc_csifi__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_csifi__, 'LineWidth', 2);
ylim ([-4.2, 2.1]);
xlim ([1, maxRank__]);
% title (title__);
xlabel (x_label__);
ylabel (y_label__);
lgd__ = legend ([p_ref__, p_all__, p_all_sep__, p_csifi__], ...
    {'IOKR: baseline', 'MP-IOKR: all candidates, joined', 'MP-IOKR: all candidates, separate', 'CSI:FingerID'}, ...
    'Interpreter','latex', 'Location', 'SouthEast');
title (lgd__, '\textbf{Metabolite-identification approach}');
grid;


if (output_for_paper__)
    pause (2);
%     imageOutputDir = ...
%         '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/figures/';
%     set (gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 22, 9]);
    imageOutputDir = ...
        '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/figures/';
    set (gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 16.5, 12.375]);
    oFn = strcat (imageOutputDir, ...
        '/REVIEW_performance_differences_without_title');
    print ('-dsvg', '-r150', oFn);
end % if

%%
if (~ output_for_paper__) 
    subplot (1, 2, 2);
else
    figure;
end % if

p_ref__ = plot (1:maxRank__, zeros (1, maxRank__), 'Color', 'black', 'LineWidth', 1.5);
hold on; 
rank_perc_diff_rand_sep__ = rank_perc_rand_sep__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand_sep__, 2)]);
p_rand_sep__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand_sep__, 2), ...
    std (rank_perc_diff_rand_sep__, 0, 2), ...
    {'Color', col_rand__, 'LineWidth', 1.5}, 1);
p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_all__, 'LineWidth', 1.5);
ylim ([-0.2, 2.1]);
xlim ([1, maxRank__]);
title (title__);
xlabel (x_label__);
ylabel (y_label__);
lgd__ = legend ([p_ref__, p_rand__.mainLine, p_all__], '[IOKR] no candidates', '[MP-IOKR] 1% random (mean \pm sd)', '[MP-IOKR] all');
title (lgd__, 'Candidate selection strategy');
grid; 

if (output_for_paper__)
    pause (2);
    imageOutputDir = ...
        '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/figures/';
    set (gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 22, 9]);
    oFn = strcat (imageOutputDir, ...
        '/REVIEW_comparison_all-vs-random_separate.svg');
    print ('-dsvg', '-r150', oFn);
end % if


%%
%% Get the candidate set sizes
% % ALL CANDIDATES (MP-IOKR)
% selection_param_all__ = struct ('strategy', 'all', 'inclExpCand', inclExpCand__); 
% selec_all__ = getCandidateSelection (Y_C, inchis, selection_param_all__);
% cand_num_all__ = cellfun (@(c) sum (c), selec_all__(has_rank__));
% [cand_num_all__, IX_all__] = sort (cand_num_all__);
% 
% cand_num_all_sep__ = cand_num_all__; 
% IX_all_sep__       = IX_all__;
% 
% % RANDOM CANDIDATES (MP-IOKR)
% % selection_param_rand__ = struct ('strategy', 'random', 'inclExpCand', inclExpCand__, ...
% %     'perc', perc__); 
% % cand_num_rand__ = zeros (size (ranks_rand__));
% % for rep__ = 1:size (ranks_rand__, 2)
% %     selec_tmp__ = getCandidateSelection (Y_C, inchis, selection_param_rand__);
% %     cand_num_rand__(:, rep__) = cellfun (@(c) sum (c), selec_tmp__);
% % end % for
% 
% % cand_num_rand__ = median (cand_num_rand__(has_rank__, :), 2);
% % [cand_num_rand__, IX_rand__] = sort (cand_num_rand__);
% 
% cand_num_rand_sep__ = median (cand_num_rand_sep__(has_rank__, :), 2);
% [cand_num_rand_sep__, IX_rand_sep__] = sort (cand_num_rand_sep__);
% 
% % LOWER-BOUNDED RANDOM CANDIDATES (MP-IOKR)
% cand_num_maxE__ = median (cand_num_maxE__(has_rank__, :), 2);
% 
% subplot (2, 2, 2);
% p_ref__ = stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% p_rand__ = stairs (mean (rank_perc_rand_sep__(1:maxRank__, :), 2), ...
%     'Color', col_rand__, 'LineWidth', 1.5);
% p_all__ = stairs (rank_perc_all_sep__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% title ('Separate');
% xlabel (x_label__);
% ylabel (y_label__);
% legend ([p_ref__, p_rand__, p_all__], 'IOKR', '1% random (mean)', 'all');
% grid;
% 
% subplot (2, 2, 3);
% p_ref__ = stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% p_rand__ = shadedErrorBar (1:maxRank__, mean (rank_perc_rand__(1:maxRank__, :), 2), ...
%     std (rank_perc_rand__(1:maxRank__, :), 0, 2), ...
%     {'Color', col_rand__, 'LineWidth', 1.5}, 1);
% p_all__ = stairs (rank_perc_all__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% title ('Joint');
% xlabel (x_label__);
% ylabel (y_label__);
% legend ([p_ref__, p_rand__.mainLine, p_all__], 'IOKR', '1% random (mean, sd)', 'all');
% grid;
% 
% subplot (2, 2, 4);
% p_ref__ = stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% p_rand__ = shadedErrorBar (1:maxRank__, mean (rank_perc_rand_sep__(1:maxRank__, :), 2), ...
%     std (rank_perc_rand_sep__(1:maxRank__, :), 0, 2), ...
%     {'Color', col_rand__, 'LineWidth', 1.5}, 1);
% p_all__ = stairs (rank_perc_all_sep__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% title ('Separate');
% xlabel (x_label__);
% ylabel (y_label__);
% legend ([p_ref__, p_rand__.mainLine, p_all__], 'IOKR', '1% random (mean, sd)', 'all');
% grid;
% 
% output_for_paper__ = false; 
% if (output_for_paper__)
%     close all;
% end % if

% if (~ output_for_paper__) ; subplot (1, 2, 1) ; end % if

% p_ref__ = plot (1:maxRank__, zeros (1, maxRank__), 'Color', 'black', 'LineWidth', 1.5);
% hold on; 
% rank_perc_diff_rand__ = rank_perc_rand__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand__, 2)]);
% p_rand__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand__, 2), ...
%     std (rank_perc_diff_rand__, 0, 2), ...
%     {'Color', col_rand__, 'LineWidth', 1.5}, 1);
% p_all__ = stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5); 
% ylim ([-0.2, 2.1]);
% xlim ([1, maxRank__]);
% title (title__);
% xlabel (x_label__);
% ylabel (y_label__);
% lgd__ = legend ([p_ref__, p_rand__.mainLine, p_all__], '[IOKR] no candidates', '[MP-IOKR] 1% random (mean \pm sd)', '[MP-IOKR] all');
% title (lgd__, 'Candidate selection strategy');
% grid; 

% p_ref__ = plot (1:maxRank__, zeros (1, maxRank__), 'Color', col_ref__, 'LineWidth', 1.5);
% hold on; 
% p_all__ = stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5); 
% p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_all_sep__, 'LineWidth', 1.5);
% p_csifi__ = stairs (rank_perc_csifi__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_csifi__, 'LineWidth', 1.5);
% p_cfm__ = stairs (rank_perc_cfm__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_cfm__, 'LineWidth', 1.5);
% ylim ([-20, 2.1]);
% xlim ([1, maxRank__]);
% title (title__);
% xlabel (x_label__);
% ylabel (y_label__);
% lgd__ = legend ([p_ref__, p_all__, p_all_sep__, p_csifi__, p_cfm__], ...
%     'IOKR: baseline', 'MP-IOKR, joint: all candidates', 'MP-IOKR, separate: all candidates', 'CSI-FingerID', 'CFM-ID');
% title (lgd__, 'Candidate selection strategy');
% grid; 
% figure; 
% plot (1:maxRank__, zeros (1, maxRank__), 'Color', 'black', 'LineWidth', 1.5);
% hold on; 
% 
% rank_perc_diff_rand__     = rank_perc_rand__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand__, 2)]);
% rank_perc_diff_rand_sep__ = rank_perc_rand_sep__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand_sep__, 2)]);
% 
% p_rand__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand__, 2), ...
%     2 * std (rank_perc_diff_rand__, 0, 2), ...
%     {'Color', col_rand__, 'LineWidth', 1.5, 'LineStyle', '-'}, 1);
% p_all__ = stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5, 'LineStyle', '-'); 
% 
% p_rand_sep__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand_sep__, 2), ...
%     2 * std (rank_perc_diff_rand_sep__, 0, 2), ...
%     {'Color', col_rand__, 'LineWidth', 1.5, 'LineStyle', '--'}, 1);
% p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5, 'LineStyle', '--');
% 
% lgd__ = legend ([p_rand__.mainLine, p_all__], 'random selection (1%)', 'all');
% title (lgd__, 'Candidate selection strategy');
% 
% xlim ([1, maxRank__]);
% title (title__);
% xlabel (x_label__);
% ylabel (y_label__);
% grid; 


%% Plot rank improvement
% n_data__ = sum (has_rank__);
% data__ = struct ( ...
%     'ranks_imp', [ranks_imp_all__(IX_all__) ; ranks_imp_rand__(IX_rand__) ; ranks_imp_maxE__(IX_maxE__)], ...
%     'cand_num',  [cand_num_all__ ; cand_num_rand__ ; cand_num_maxE__],                                    ... 
%     'id',        {repmat((1:n_data__)', [3, 1])},                                                         ...
%     'strategy',  {[repmat({'all'}, [n_data__, 1]) ; repmat({'random'}, [n_data__, 1]) ; repmat({'maxElement'}, [n_data__, 1])]});
% 
% data__ = struct ( ...
%     'ranks_imp', [ranks_imp_all__(IX_all__) ; ranks_imp_rand__(IX_rand__)], ...
%     'cand_num',  [cand_num_all__ ; cand_num_rand__],                        ... 
%     'id',        {repmat((1:n_data__)', [2, 1])},                           ...
%     'strategy',  {[repmat({'all'}, [n_data__, 1]) ; repmat({'random'}, [n_data__, 1])]});
% 
% figure;
% g_rank_imp__(1, 1) = gramm ('x', data__.id, 'y', data__.ranks_imp, 'color', data__.strategy);
% g_rank_imp__(1, 1).set_color_options ('map', 'matlab');
% g_rank_imp__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
% g_rank_imp__(1, 1).geom_point ('alpha', 0.20);
% g_rank_imp__(1, 1).geom_hline ('yintercept', 0);
% g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'point', 'bin_in', 25, 'dodge', 0.2);
% g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'errorbar', 'bin_in', 25, 'dodge', 0.2);
% 
% g_rank_imp__(2, 1) = gramm ('x', data__.id, 'y', data__.cand_num, 'color', data__.strategy);
% g_rank_imp__(2, 1).set_color_options ('map', 'matlab');
% g_rank_imp__(2, 1).axe_property ('XGrid', 'on', 'YGrid', 'on', 'YScale', 'log');
% g_rank_imp__(2, 1).geom_point ('alpha', 0.20);
% 
% g_rank_imp__.draw();
% 
% %% Plot rank improvement 2
% data__ = struct ( ...
%     'ranks_imp_all_sep', ranks_imp_all_sep__(IX_all_sep__), ...
%     'ranks_imp_rand_sep', ranks_imp_rand_sep__(IX_rand_sep__));
% 
% figure;
% g__ = gramm ('x', data__.ranks_imp_all_sep, 'y', data__.ranks_imp_rand_sep);
% g__.set_color_options ('map', 'matlab');
% g__.geom_point();
% g__.geom_abline();
% g__.axe_property ('XGrid', 'on', 'YGrid', 'on');
% g__.set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (random)');
% 
% g__.draw();
% %% Plot scatter: ALL vs RANDOM
% data__ = struct (                       ...
%     'ranks_imp_all',  ranks_imp_all__,  ...
%     'ranks_imp_rand', ranks_imp_rand__, ...
%     'ranks_imp_maxE', ranks_imp_maxE__);
% 
% figure;
% g_vs__(1, 1) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_rand);
% g_vs__(1, 1).set_color_options ('map', 'matlab');
% g_vs__(1, 1).geom_point();
% g_vs__(1, 1).geom_abline();
% g_vs__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
% g_vs__(1, 1).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (random)');
% 
% g_vs__(1, 2) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_maxE);
% g_vs__(1, 2).set_color_options ('map', 'matlab');
% g_vs__(1, 2).geom_point();
% g_vs__(1, 2).geom_abline();
% g_vs__(1, 2).axe_property ('XGrid', 'on', 'YGrid', 'on');
% g_vs__(1, 2).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (maxElement)');
% 
% g_vs__(1, 3) = gramm ('x', data__.ranks_imp_rand, 'y', data__.ranks_imp_maxE);
% g_vs__(1, 3).set_color_options ('map', 'matlab');
% g_vs__(1, 3).geom_point();
% g_vs__(1, 3).geom_abline();
% g_vs__(1, 3).axe_property ('XGrid', 'on', 'YGrid', 'on');
% g_vs__(1, 3).set_names ('x', 'Rank improvement (random)', 'y', 'Rank improvement (maxElement)');
% 
% g_vs__.draw();
% 
% %% Plot ranks 
% data__ = struct (                 ...
%     'ranks_ref',     ranks_ref__(has_rank__), ...
%     'ranks_all_sep', ranks_all_sep__(has_rank__));
% 
% figure;
% g_rank__ = gramm ('x', data__.ranks_ref, 'y', data__.ranks_all_sep);
% g_rank__.set_color_options ('map', 'matlab');
% g_rank__.geom_point();
% g_rank__.geom_abline();
% g_rank__.axe_property ('XGrid', 'on', 'YGrid', 'on');
% g_rank__.set_names ('x', 'Rank (reference)', 'y', 'Rank (all, separate)');
% 
% g_rank__.draw();

% Leave the workspace a bit tidy
clear *__;