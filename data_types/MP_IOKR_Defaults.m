%% MP_IOKR_DEFAULTS
classdef MP_IOKR_Defaults
    properties (Constant, GetAccess = private)
        paramCategories = {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'};
    end % properties
    
    properties (Constant, GetAccess = public)
        % Debug parameter
        debug_param = struct (    ...
            'isDebugMode', false, ...
            'verbose',     true,  ...
            'n_debug_set', 50,    ...
            'randomSeed',  123);
        
%         Optimization parameter
%         opt_param = struct (                                                        ...
%             'val_gamma',   [0.5], ...
%             'val_lambda',  [0.5], ...
%             'nOuterFolds', 10,                                                      ...
%             'nInnerFolds', 2);
        
        opt_param = struct (                                                        ...
            'val_gamma',   sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1e3 1e4]), ...
            'val_lambda',  sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1e3 1e4]), ...
            'nOuterFolds', 10,                                                      ...
            'nInnerFolds', 10);
        
        % MP_IOKR parameter
        mp_iokr_param = struct (  ...
            'center',   true,     ...
            'mkl',      'unimkl', ...
            'rev_iokr', 'joint');
        
        % Data parameter
        % FIXME: 'availInputKernels' is specific to the metabolite
        % identification.
        data_param = struct (                                                              ...
            'usePreCalcStat',       false,                                                 ...
            'inputKernel',          'PPKR',                                                ...
            'availInputKernels',    {{'ALIGND', 'ALIGN', 'CEC', 'CP2Plus', 'CP2', 'CPC',   ...
                 'CPI', 'CPJB', 'CPJ', 'CPK', 'CSC', 'FIPP', 'LB', 'LC', 'LIPP', 'LI',     ...
                 'LW', 'NB', 'NI', 'NSF', 'PPKR', 'RLB', 'RLI', 'WPC'}},                   ...
            'selection_param',      struct ('strategy', 'all', 'inclExpCand', true),       ...
            'repetition',           -1);
    end % properties
    
    methods (Static, Access = public)
        function param = setDefaultsIfNeeded (param, paramCategories)    
        %% SETDEFAULTSIFNEEDED
%             warning ('The parameter for VAL_LAMBDA and VAL_GAMMA are in debug setting.');
        
            if (nargin < 1)
                error ('MP_IOKR_Defaults:setDefaultsIfNeeded:InvalidInput', ...
                    'Not enough input arguments.');
            end % if
            if (nargin < 2)
                paramCategories = MP_IOKR_Defaults.paramCategories;
            end % if
            
            for fnCat = paramCategories
                if (~ ismember (fnCat{1}, MP_IOKR_Defaults.paramCategories))
                    warning ('No defaults available: %s\n', fnCat{1});
                    
                    continue;
                end % if
                
                if (~ ismember (fnCat{1}, fieldnames (param)))
                    fprintf ('Use all default values: PARAM.%s\n', upper (fnCat{1}));

                    param.(fnCat{1})= MP_IOKR_Defaults.(fnCat{1});
                else
                    for fn = fieldnames (MP_IOKR_Defaults.(fnCat{1}))'
                        if ((~ ismember (fn{1}, fieldnames (param.(fnCat{1})))) || isempty (param.(fnCat{1}).(fn{1}))) 
                            fprintf ('Use default value: PARAM.%s.%s\n', upper (fnCat{1}), upper (fn{1}));

                            param.(fnCat{1}).(fn{1}) = MP_IOKR_Defaults.(fnCat{1}).(fn{1});
                        end % if
                    end % for
                end % if
            end % for
        end % function 
    end % methods
end % class