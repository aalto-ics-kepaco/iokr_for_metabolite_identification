close all;

%% Load data
inclExpCand__ = false;  

% ALL CANDIDATES (MP_IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                     'inner', struct ('nFolds', 10));                  

[rank_perc_all__, ranks_all__] = aggregate_results_all ('UNIMKL', cv_param__, inclExpCand__);

% ALL CANDIDATES SEPARATE (MP_IOKR)            
[rank_perc_all_sep__, ranks_all_sep__] = aggregate_results_all_separate ('UNIMKL', cv_param__, inclExpCand__);

% RANDOM CANDIDATES (MP-IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));     
perc__ = 1;
% 
% [rank_perc_rand__, ranks_rand__, cand_num_rand__] = aggregate_results_random ( ...
%     'UNIMKL', cv_param__, perc__, inclExpCand__);
[rank_perc_rand__, ranks_rand__] = aggregate_results_random ( ...
    'UNIMKL', cv_param__, perc__, inclExpCand__);

% RANDOM CANDIDATES (MP-IOKR)
[rank_perc_rand_sep__, ranks_rand_sep__, cand_num_rand_sep__] = aggregate_results_random_separate ( ...
    'UNIMKL', cv_param__, perc__, inclExpCand__);

% LOWER-BOUNDED RANDOM CANDIDATES (MP-IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));     
maxNumCand__ = 140;

[rank_perc_maxE__, ranks_maxE__, cand_num_maxE__] = aggregate_results_maxElement ( ...
    'UNIMKL', cv_param__, maxNumCand__, inclExpCand__);

% REFERENCE (IOKR)
tmp__ = load ( ...
    '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/iokr/a143cc0ca62fc3c75e383cb4fb679201.mat');
ranks_ref__     = tmp__.result.ranks;
rank_perc_ref__ = tmp__.result.rank_perc;

fprintf ('IOKR\tMP-IOKR (all)\tMP-IOKR separate (all)\tMP-IOKR (random) [median]\tMP-IOKR (maxElement) [median]\n');
disp ([rank_perc_ref__(1:20),               ...
    rank_perc_all__(1:20),                  ...
    rank_perc_all_sep__(1:20),              ...
    mean(rank_perc_rand__(1:20, :), 2),     ...
    mean(rank_perc_rand_sep__(1:20, :), 2), ...
    mean(rank_perc_maxE__(1:20, :), 2)]);

%% Rank improvement
assert (all (isnan (ranks_ref__) == isnan (ranks_all__)));
assert (all (isnan (ranks_ref__) == isnan (ranks_all_sep__)));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_rand__, 2))));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_rand_sep__, 2))));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_maxE__, 2))));

has_rank__ = ~ isnan (ranks_ref__);

ranks_imp_all__      = ranks_ref__ - ranks_all__;
ranks_imp_all_sep__  = ranks_ref__ - ranks_all_sep__;
ranks_imp_rand__     = ranks_ref__ - median (ranks_rand__, 2);
ranks_imp_rand_sep__ = ranks_ref__ - median (ranks_rand_sep__, 2);
ranks_imp_maxE__     = ranks_ref__ - median (ranks_maxE__, 2);

ranks_imp_all__      = ranks_imp_all__(has_rank__);
ranks_imp_all_sep__  = ranks_imp_all_sep__(has_rank__);
ranks_imp_rand__     = ranks_imp_rand__(has_rank__);
ranks_imp_rand_sep__ = ranks_imp_rand_sep__(has_rank__);
ranks_imp_maxE__     = ranks_imp_maxE__(has_rank__);

%% Get the candidate set sizes
% ALL CANDIDATES (MP-IOKR)
selection_param_all__ = struct ('strategy', 'all', 'inclExpCand', inclExpCand__); 
selec_all__ = getCandidateSelection (Y_C, inchis, selection_param_all__);
cand_num_all__ = cellfun (@(c) sum (c), selec_all__(has_rank__));
[cand_num_all__, IX_all__] = sort (cand_num_all__);

cand_num_all_sep__ = cand_num_all__; 
IX_all_sep__       = IX_all__;

% RANDOM CANDIDATES (MP-IOKR)
% selection_param_rand__ = struct ('strategy', 'random', 'inclExpCand', inclExpCand__, ...
%     'perc', perc__); 
% cand_num_rand__ = zeros (size (ranks_rand__));
% for rep__ = 1:size (ranks_rand__, 2)
%     selec_tmp__ = getCandidateSelection (Y_C, inchis, selection_param_rand__);
%     cand_num_rand__(:, rep__) = cellfun (@(c) sum (c), selec_tmp__);
% end % for

% cand_num_rand__ = median (cand_num_rand__(has_rank__, :), 2);
% [cand_num_rand__, IX_rand__] = sort (cand_num_rand__);

cand_num_rand_sep__ = median (cand_num_rand_sep__(has_rank__, :), 2);
[cand_num_rand_sep__, IX_rand_sep__] = sort (cand_num_rand_sep__);

% LOWER-BOUNDED RANDOM CANDIDATES (MP-IOKR)
cand_num_maxE__ = median (cand_num_maxE__(has_rank__, :), 2);
[cand_num_maxE__, IX_maxE__] = sort (cand_num_maxE__);

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

figure; 
colors__ =  [ ...
    0         0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840];
col_rand__ = colors__(1, :);
col_all__ = colors__(2, :);

title__   = 'Performance gain using MP-IOKR over IOKR';
x_label__ = 'k';
y_label__ = 'Relative change of top-k accuracy';

% subplot (1, 2, 1);
% stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% stairs (mean (rank_perc_rand__(1:maxRank__, :), 2), ...
%     'Color', col_rand__, 'LineWidth', 1.5);
% stairs (rank_perc_all__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% xlabel (x_label__);
% ylabel ('Top-k accuracy');
% 
% subplot (1, 2, 2);
% stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
% hold on;
% stairs (mean (rank_perc_rand_sep__(1:maxRank__, :), 2), ...
%     'Color', col_rand__, 'LineWidth', 1.5);
% stairs (rank_perc_all_sep__(1:maxRank__), ...
%     'Color', col_all__, 'LineWidth', 1.5);
% ylim ([30, 90]);
% xlabel (x_label__);
% ylabel ('Top-k accuracy');

figure; 
subplot (1, 2, 1);
plot (1:maxRank__, zeros (1, maxRank__), 'Color', 'black', 'LineWidth', 1.5);
hold on; 
rank_perc_diff_rand__ = rank_perc_rand__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand__, 2)]);
% stairs (mean (rank_perc_diff_rand__, 2), ...
%     'Color', colors(2, :), 'LineWidth', 1.5, 'LineStyle', '-'); 
p_rand__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand__, 2), ...
    2 * std (rank_perc_diff_rand__, 0, 2), ...
    {'Color', col_rand__, 'LineWidth', 1.5}, 1);
p_all__ = stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_all__, 'LineWidth', 1.5); 
ylim ([-0.3, 2]);
xlim ([1, maxRank__]);
title (title__);
xlabel (x_label__);
ylabel (y_label__);
lgd__ = legend ([p_rand__.mainLine, p_all__], 'random selection (1%)', 'all');
title (lgd__, 'Candidate selection strategy');
grid; 

subplot (1, 2, 2);
plot (1:maxRank__, zeros (1, maxRank__), 'Color', 'black', 'LineWidth', 1.5);
hold on; 
rank_perc_diff_rand_sep__ = rank_perc_rand_sep__(1:maxRank__, :) - repmat (rank_perc_ref__(1:maxRank__), [1, size(rank_perc_rand_sep__, 2)]);
% stairs (mean (rank_perc_diff_rand_sep__, 2), ...
%     'Color', colors(2, :), 'LineWidth', 1.5, 'LineStyle', '-.'); 
p_rand_sep__ = shadedErrorBar (1:maxRank__, mean (rank_perc_diff_rand_sep__, 2), ...
    2 * std (rank_perc_diff_rand_sep__, 0, 2), ...
    {'Color', col_rand__, 'LineWidth', 1.5}, 1);
p_all_sep__ = stairs (rank_perc_all_sep__(1:maxRank__) - rank_perc_ref__(1:maxRank__), ...
    'Color', col_all__, 'LineWidth', 1.5);
ylim ([-0.3, 2]);
xlim ([1, maxRank__]);
title (title__);
xlabel (x_label__);
ylabel (y_label__);
lgd__ = legend ([p_rand_sep__.mainLine, p_all_sep__], 'random selection (1%)', 'all');
title (lgd__, 'Candidate selection strategy');
grid; 

% imageOutputDir = ...
%     '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/figures/';
% set (gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 44, 9]);
% oFn = strcat (imageOutputDir, ...
%     '/comparison_all-vs-random_and_separate-vs-joint.svg');
% print ('-dsvg', '-r150', oFn);

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
n_data__ = sum (has_rank__);
data__ = struct ( ...
    'ranks_imp', [ranks_imp_all__(IX_all__) ; ranks_imp_rand__(IX_rand__) ; ranks_imp_maxE__(IX_maxE__)], ...
    'cand_num',  [cand_num_all__ ; cand_num_rand__ ; cand_num_maxE__],                                    ... 
    'id',        {repmat((1:n_data__)', [3, 1])},                                                         ...
    'strategy',  {[repmat({'all'}, [n_data__, 1]) ; repmat({'random'}, [n_data__, 1]) ; repmat({'maxElement'}, [n_data__, 1])]});

data__ = struct ( ...
    'ranks_imp', [ranks_imp_all__(IX_all__) ; ranks_imp_rand__(IX_rand__)], ...
    'cand_num',  [cand_num_all__ ; cand_num_rand__],                        ... 
    'id',        {repmat((1:n_data__)', [2, 1])},                           ...
    'strategy',  {[repmat({'all'}, [n_data__, 1]) ; repmat({'random'}, [n_data__, 1])]});

figure;
g_rank_imp__(1, 1) = gramm ('x', data__.id, 'y', data__.ranks_imp, 'color', data__.strategy);
g_rank_imp__(1, 1).set_color_options ('map', 'matlab');
g_rank_imp__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_rank_imp__(1, 1).geom_point ('alpha', 0.20);
g_rank_imp__(1, 1).geom_hline ('yintercept', 0);
g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'point', 'bin_in', 25, 'dodge', 0.2);
g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'errorbar', 'bin_in', 25, 'dodge', 0.2);

g_rank_imp__(2, 1) = gramm ('x', data__.id, 'y', data__.cand_num, 'color', data__.strategy);
g_rank_imp__(2, 1).set_color_options ('map', 'matlab');
g_rank_imp__(2, 1).axe_property ('XGrid', 'on', 'YGrid', 'on', 'YScale', 'log');
g_rank_imp__(2, 1).geom_point ('alpha', 0.20);

g_rank_imp__.draw();

%% Plot rank improvement 2
data__ = struct ( ...
    'ranks_imp_all_sep', ranks_imp_all_sep__(IX_all_sep__), ...
    'ranks_imp_rand_sep', ranks_imp_rand_sep__(IX_rand_sep__));

figure;
g__ = gramm ('x', data__.ranks_imp_all_sep, 'y', data__.ranks_imp_rand_sep);
g__.set_color_options ('map', 'matlab');
g__.geom_point();
g__.geom_abline();
g__.axe_property ('XGrid', 'on', 'YGrid', 'on');
g__.set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (random)');

g__.draw();
%% Plot scatter: ALL vs RANDOM
data__ = struct (                       ...
    'ranks_imp_all',  ranks_imp_all__,  ...
    'ranks_imp_rand', ranks_imp_rand__, ...
    'ranks_imp_maxE', ranks_imp_maxE__);

figure;
g_vs__(1, 1) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_rand);
g_vs__(1, 1).set_color_options ('map', 'matlab');
g_vs__(1, 1).geom_point();
g_vs__(1, 1).geom_abline();
g_vs__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 1).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (random)');

g_vs__(1, 2) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_maxE);
g_vs__(1, 2).set_color_options ('map', 'matlab');
g_vs__(1, 2).geom_point();
g_vs__(1, 2).geom_abline();
g_vs__(1, 2).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 2).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (maxElement)');

g_vs__(1, 3) = gramm ('x', data__.ranks_imp_rand, 'y', data__.ranks_imp_maxE);
g_vs__(1, 3).set_color_options ('map', 'matlab');
g_vs__(1, 3).geom_point();
g_vs__(1, 3).geom_abline();
g_vs__(1, 3).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 3).set_names ('x', 'Rank improvement (random)', 'y', 'Rank improvement (maxElement)');

g_vs__.draw();

%% Plot ranks 
data__ = struct (                 ...
    'ranks_ref',     ranks_ref__(has_rank__), ...
    'ranks_all_sep', ranks_all_sep__(has_rank__));

figure;
g_rank__ = gramm ('x', data__.ranks_ref, 'y', data__.ranks_all_sep);
g_rank__.set_color_options ('map', 'matlab');
g_rank__.geom_point();
g_rank__.geom_abline();
g_rank__.axe_property ('XGrid', 'on', 'YGrid', 'on');
g_rank__.set_names ('x', 'Rank (reference)', 'y', 'Rank (all, separate)');

g_rank__.draw();
% Leave the workspace a bit tidy
clear *__;