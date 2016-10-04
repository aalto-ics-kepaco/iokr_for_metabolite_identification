close all;

%% Load data
inclExpCand__ = false;  

% ALL CANDIDATES (MP_IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                     'inner', struct ('nFolds', 10));                  

[rank_perc_all__, ranks_all__] = aggregate_results_all ('UNIMKL', cv_param__, inclExpCand__);

% RANDOM CANDIDATES (MP-IOKR)
cv_param__ = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));     
perc__ = 1;

[rank_perc_rand__, ranks_rand__, cand_num_rand__] = aggregate_results_random ( ...
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

fprintf ('IOKR\tMP-IOKR (all)\tMP-IOKR (random) [median]\tMP-IOKR (maxElement) [median]\n');
disp ([rank_perc_ref__(1:20),             ...
    rank_perc_all__(1:20),                ...
    median(rank_perc_rand__(1:20, :), 2), ...
    median(rank_perc_maxE__(1:20, :), 2)]);

%% Rank improvement
assert (all (isnan (ranks_ref__) == isnan (ranks_all__)));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_rand__, 2))));
assert (all (isnan (ranks_ref__) == isnan (median (ranks_maxE__, 2))));

has_rank__ = ~ isnan (ranks_ref__);

ranks_imp_all__ = ranks_ref__ - ranks_all__;
has_rank_change_all__ = ranks_imp_all__ ~= 0;

ranks_imp_rand__ = ranks_ref__ - median (ranks_rand__, 2);
has_rank_change_rand__ = ranks_imp_rand__ ~= 0;

ranks_imp_maxE__ = ranks_ref__ - median (ranks_maxE__, 2);
has_rank_change_maxE__ = ranks_imp_maxE__ ~= 0;


% ranks_imp__ = ranks_imp__(has_rank__ & has_rank_change__);
ranks_imp_all__  = ranks_imp_all__(has_rank__);
ranks_imp_rand__ = ranks_imp_rand__(has_rank__);
ranks_imp_maxE__ = ranks_imp_maxE__(has_rank__);

%% Get the candidate set sizes
% ALL CANDIDATES (MP-IOKR)
selection_param_all__ = struct ('strategy', 'all', 'inclExpCand', inclExpCand__); 
selec_all__ = getCandidateSelection (Y_C, inchis, selection_param_all__);
% cand_num__ = cellfun (@(c) sum (c), selec__(has_rank__ & has_rank_change__));
cand_num_all__ = cellfun (@(c) sum (c), selec_all__(has_rank__));
[cand_num_all__, IX_all__] = sort (cand_num_all__);

% RANDOM CANDIDATES (MP-IOKR)
% selection_param_rand__ = struct ('strategy', 'random', 'inclExpCand', inclExpCand__, ...
%     'perc', perc__); 
% cand_num_rand__ = zeros (size (ranks_rand__));
% for rep__ = 1:size (ranks_rand__, 2)
%     selec_tmp__ = getCandidateSelection (Y_C, inchis, selection_param_rand__);
%     cand_num_rand__(:, rep__) = cellfun (@(c) sum (c), selec_tmp__);
% end % for
cand_num_rand__ = median (cand_num_rand__(has_rank__, :), 2);
[cand_num_rand__, IX_rand__] = sort (cand_num_rand__);

% LOWER-BOUNDED RANDOM CANDIDATES (MP-IOKR)
cand_num_maxE__ = median (cand_num_maxE__(has_rank__, :), 2);
[cand_num_maxE__, IX_maxE__] = sort (cand_num_maxE__);

%% Plot top-k accuracy
maxRank__ = 25;

subplot (1, 2, 1);
stairs (rank_perc_ref__(1:maxRank__), 'Color', 'black', 'LineWidth', 1.5); 
hold on;
stairs (rank_perc_all__(1:maxRank__), 'Color', 'red', 'LineWidth', 1.5);
stairs (mean (rank_perc_rand__(1:maxRank__, :), 2), 'Color', 'blue', 'LineWidth', 1.5);
stairs (mean (rank_perc_maxE__(1:maxRank__, :), 2), 'Color', 'gree', 'LineWidth', 1.5);
legend ('IOKR', 'MP-IOKR (all)', 'MP-IOKR (random)', 'MP-IOKR (maxElement)');
grid;

subplot (1, 2, 2);
stairs (rank_perc_all__(1:maxRank__) - rank_perc_ref__(1:maxRank__), 'Color', 'red', 'LineWidth', 1.5); 
hold on;
stairs (mean (rank_perc_rand__(1:maxRank__, :), 2) - rank_perc_ref__(1:maxRank__), 'Color', 'blue', 'LineWidth', 1.5); 
stairs (mean (rank_perc_maxE__(1:maxRank__, :), 2) - rank_perc_ref__(1:maxRank__), 'Color', 'green', 'LineWidth', 1.5); 
legend ('MP-IOKR (all)', 'MP-IOKR (random)', 'MP-IOKR (maxElement)');
plot (1:maxRank__, zeros (1, maxRank__), 'k--');
grid; 

%% Plot rank improvement
n_data__ = sum (has_rank__);
data__ = struct ( ...
    'ranks_imp', [ranks_imp_all__(IX_all__) ; ranks_imp_rand__(IX_rand__) ; ranks_imp_maxE__(IX_maxE__)], ...
    'cand_num',  [cand_num_all__ ; cand_num_rand__ ; cand_num_maxE__],                                    ... 
    'id',        {repmat((1:n_data__)', [3, 1])},                                                         ...
    'strategy',  {[repmat({'all'}, [n_data__, 1]) ; repmat({'random'}, [n_data__, 1]) ; repmat({'maxElement'}, [n_data__, 1])]});

figure;
g_rank_imp__(1, 1) = gramm ('x', data__.id, 'y', data__.ranks_imp, 'color', data__.strategy);
g_rank_imp__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_rank_imp__(1, 1).geom_point ('alpha', 0.20);
g_rank_imp__(1, 1).geom_hline ('yintercept', 0);
g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'point', 'bin_in', 25, 'dodge', 0.2);
g_rank_imp__(1, 1).stat_summary ('type', '95percentile', 'geom', 'errorbar', 'bin_in', 25, 'dodge', 0.2);

g_rank_imp__(2, 1) = gramm ('x', data__.id, 'y', data__.cand_num, 'color', data__.strategy);
g_rank_imp__(2, 1).axe_property ('XGrid', 'on', 'YGrid', 'on', 'YScale', 'log');
g_rank_imp__(2, 1).geom_point ('alpha', 0.20);

g_rank_imp__.draw();

%% Plot scatter: ALL vs RANDOM
data__ = struct (                       ...
    'ranks_imp_all',  ranks_imp_all__,  ...
    'ranks_imp_rand', ranks_imp_rand__, ...
    'ranks_imp_maxE', ranks_imp_maxE__);

figure;
g_vs__(1, 1) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_rand);
g_vs__(1, 1).geom_point();
g_vs__(1, 1).geom_abline();
g_vs__(1, 1).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 1).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (random)');

g_vs__(1, 2) = gramm ('x', data__.ranks_imp_all, 'y', data__.ranks_imp_maxE);
g_vs__(1, 2).geom_point();
g_vs__(1, 2).geom_abline();
g_vs__(1, 2).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 2).set_names ('x', 'Rank improvement (all)', 'y', 'Rank improvement (maxElement)');

g_vs__(1, 3) = gramm ('x', data__.ranks_imp_rand, 'y', data__.ranks_imp_maxE);
g_vs__(1, 3).geom_point();
g_vs__(1, 3).geom_abline();
g_vs__(1, 3).axe_property ('XGrid', 'on', 'YGrid', 'on');
g_vs__(1, 3).set_names ('x', 'Rank improvement (random)', 'y', 'Rank improvement (maxElement)');

g_vs__.draw();

% Leave the workspace a bit tidy
clear *__;