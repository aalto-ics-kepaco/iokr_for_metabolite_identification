function add_results_to_db (resultDir)
    files = dir ([resultDir, '/*.mat']); 
    line_info = cell (1);
    ii = 1;
    
%     figure; hold on;
    for file = files'
        res = load([resultDir, '/', file.name]);
        
        if (strcmp (res.result.param.ky_param.type, 'gaussian') && ...
                strcmp (res.result.param.ky_param.base_kernel, 'linear') && ...
                strcmp (res.result.param.mp_iokr_param.mkl, 'unimkl'))
            
            if (strcmp (res.result.param.ky_param.representation, 'feature'))
                fprintf ('rff_dimension=%d: %.3f %.3f %.3f %.3f\n', ...
                    res.result.param.ky_param.rff_dimension, ...
                    res.result.rank_perc(1), ...
                    res.result.rank_perc(5), ...
                    res.result.rank_perc(10), ...
                    res.result.rank_perc(20));
                
%                  stairs (res.result.rank_perc(1:100))
%                  line_info{ii} = sprintf ('rff_dimension=%d', res.result.param.ky_param.rff_dimension);
            else
                fprintf ('kernel: %.3f %.3f %.3f %.3f\n', ...
                    res.result.rank_perc(1), ...
                    res.result.rank_perc(5), ...
                    res.result.rank_perc(10), ...
                    res.result.rank_perc(20));
                
                
%                  stairs (res.result.rank_perc(1:100))
%                  line_info{ii} = 'gaussian';
            end % if
            
%             ii = ii + 1;
        end % if
        

    end % for
    
%     legend (line_info);
end % function

