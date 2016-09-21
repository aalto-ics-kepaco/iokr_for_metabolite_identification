%% MP_IOKR_DEFAULTS
classdef MP_IOKR_Defaults
    properties (Constant, GetAccess = private)
        paramCategories = {'opt_param', 'mp_iokr_param', 'data_param'};
    end % properties
    
    properties (Constant, GetAccess = public)
        % Optimization parameter
        opt_param = struct (                                                       ...
            'val_gamma',  sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1e3 1e4]), ...
            'val_lambda', sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1e3 1e4]));
        
        % MP_IOKR parameter
        mp_iokr_param = struct (  ...
            'center',   true,     ...
            'mkl',      'unimkl', ...
            'rev_iokr', 'separate');
        
        % Data parameter
        data_param = struct (                                                           ...
            'usePreCalcStat',    false,                                                 ...
            'inputKernel',       'unimkl',                                              ...
            'availInputKernels', {{'ALIGND', 'ALIGN', 'CEC', 'CP2Plus', 'CP2', 'CPC',   ...
                 'CPI', 'CPJB', 'CPJ', 'CPK', 'CSC', 'FIPP', 'LB', 'LC', 'LIPP', 'LI',  ...
                 'LW', 'NB', 'NI', 'NSF', 'PPKR', 'RLB', 'RLI', 'WPC'}});
    end % properties
    
    methods (Static, Access = public)
        function param = setDefaultsIfNeeded (param, paramCategories)
        %% SETDEFAULTSIFNEEDED
            if (nargin < 1)
                error ('MP_IOKR_Defaults:setDefaultsIfNeeded:InvalidInput', ...
                    'Not enough input arguments.');
            end % if
            if (nargin < 2)
                paramCategories = MP_IOKR_Defaults.paramCategories;
            end % if
            
            for fnCat = paramCategories'
                if (~ ismember (fnCat, MP_IOKR_Defaults.paramCategories))
                    warning ('No defaults available: %s\n', fnCat);
                    
                    continue;
                end % if
                
                if (~ ismember (fnCat, fieldnames (param)))
                    fprintf ('Use all default values: PARAM.%s\n', upper (fnCat));

                    param.(fnCat)= MP_IOKR_Defaults.(fnCat);
                else
                    for fn = fieldnames (MP_IOKR_Defaults.(fnCat))'
                        if ((~ ismember (fn, fieldnames (param.(fnCat)))) || isempty (param.(fnCat).(fn))) 
                            fprintf ('Use default value: PARAM.%s.%s\n', upper (fnCat), upper (fn));

                            param.(fnCat).(fn) = MP_IOKR_Defaults.(fnCat).(fn);
                        end % if
                    end % for
                end % if
            end % for
        end % function 
    end % methods
end % class