function summarizeResults ()
    inputKernel = 'UNIMKL';
    perc = 1; 
    inclExpCand = false;
            
    inputDir = ...
        '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/';
    
    ind_fold = load (strcat (inputDir, '/cv_ind.txt'));
    cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));
    
    rank_perc = aggregate_results_random (inputKernel, cv_param, perc, inclExpCand);
    
    %% Load reference
    outputDir = ...
        '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/iokr/';
    
    param = struct ();
    
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'iokr_param', 'data_param'});
    
    param.data_param.cv_param = cv_param;
    param.data_param.inputKernel = inputKernel;
    
    settingHash = DataHash (struct (                            ...
        'cv_param',        param.data_param.cv_param,           ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.iokr_param.center,             ...
        'cv_type',         param.iokr_param.cv_type));
    
    resFn = strcat (outputDir, '/', settingHash, '.mat'); 
    fprintf ('%s exists %d\n', resFn, exist (resFn, 'file'));
    
    tmp = load (resFn);
    rank_perc_ref = tmp.result.rank_perc(1:100);
    
    %% Results
    maxRankPerc = 25;
    
    subplot (1, 2, 1);
    stairs (rank_perc_ref(1:maxRankPerc), 'k', 'LineWidth', 1.5);
    title (sprintf ('Comparison using random candidates: %.2f%%, inclExpCand %d', ...
        perc, inclExpCand));
    hold on; 
    stairs (mean (rank_perc(1:maxRankPerc, :), 2), 'red', 'LineWidth', 1.5);
    legend ('baseline (IOKR)', 'MP-IOKR', 'Location', 'SouthEast');
    xlabel ('Rank');
    ylabel ('Cumulative %');
    ylim ([30, 80]);
    xlim ([1, maxRankPerc]);
    grid;
    
    subplot (1, 2, 2);
    plot ([1, maxRankPerc], [0, 0], 'k', 'LineWidth', 1.5); hold on;
    title (sprintf ('Comparison using random candidates: %.2f%%, inclExpCand %d', ...
        perc, inclExpCand));
    shadedErrorBar (1:maxRankPerc, mean (rank_perc(1:maxRankPerc, :), 2) - rank_perc_ref(1:maxRankPerc), ...
        2 * std (rank_perc(1:maxRankPerc, :), 0, 2), {'r', 'LineWidth', 1.5}, 1); 
    xlabel ('Rank');
    ylabel ('Rel. change of %-points');
        legend ('baseline (IOKR)', 'MP-IOKR (error 2\sigma)', 'MP-IOKR (mean)', ...
        'Location', 'SouthEast');
    ylim ([-0.4, 1.2]);
    xlim ([1, maxRankPerc]);
    grid;
    
    imageOutputDir = ...
        '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/figures/';
    set (gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 40, 12]);
    oFn = strcat (imageOutputDir, ...
        '/MP-IOKR_random_perc=', num2str(perc), 'inclExpCand=', num2str(inclExpCand), '.svg');
    print ('-dsvg', '-r150', oFn);
    
    disp ('hallo')
end % function 