%% MP_IOKR_DEFAULTS
classdef MP_IOKR_Defaults
    properties (Constant, GetAccess = private)
        paramCategories = {  ...
            'debug_param',   ...
            'opt_param',     ... 
            'mp_iokr_param', ...
            'iokr_param',    ...
            'data_param',    ...
            'ky_param'};
    end % properties
    
    properties (Constant, GetAccess = public)
        % Debug parameter
        debug_param = struct (        ...
            'verbose',     true,      ...
            'isDebugMode', false,     ...
            'n_debug_set', 50,        ...
            'randomSeed',  'shuffle', ...
            'cand',        []);
        
        % Optimization parameter        
        opt_param = struct (                                                     ...
            'val_gamma',   sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1000]), ...
            'val_lambda',  sort ([1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100 1000]), ...
            'nOuterFolds', 10,                                                   ...
            'nInnerFolds', 10,                                                   ...
            'cv_type',     'loocv');
        
        % MP_IOKR parameter
        mp_iokr_param = struct (  ...
            'center',   true,     ...
            'mkl',      'unimkl', ...
            'rev_iokr', 'separate');
        
        % IOKR parameter
        iokr_param = struct (                 ...
            'center',               true,     ...
            'mkl',                  'unimkl', ...
            'model_representation', 'only_C');
        
        % Output kernel parameters
        ky_param = struct ( ...    
            'representation',  'kernel',   ...
            'type',            'gaussian', ...
            'base_kernel',     'tanimoto', ...
            'param_selection', 'entropy',  ...
            'rff_dimension',   -1);
        
        % Data parameter
        % TODO: Parameter regarding the candidate selection should go to
        %       the 'mp_iokr_param'.
        data_param = struct (                                                              ...
            'usePreCalcStat',       false,                                                 ...
            'inputKernel',          'ALL',                                             ...
            'availInputKernels',    {{'ALIGND', 'ALIGN', 'CEC', 'CP2Plus', 'CP2', 'CPC',   ...
                 'CPI', 'CPJB', 'CPJ', 'CPK', 'CSC', 'FIPP', 'LB', 'LC', 'LIPP', 'LI',     ...
                 'LW', 'NB', 'NI', 'NSF', 'PPKR', 'RLB', 'RLI', 'WPC'}},                   ...
            'selection_param',      struct ('strategy', 'all', 'inclExpCand', true),       ...
            'repetition',           -1);
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
        
        function str = param2str (param, simple, include_info)
        %% PARAM2STR Returns a string containing the parameter setting
            if (nargin < 1)
                error ('MP_IOKR_Defaults:param2str:InvalidInput', ...
                    'Not enough input arguments.');
            end % if
            if (nargin < 2)
                simple = true;
            end % if
            if (nargin < 3)
                include_info = {'iokr_param', 'ky_param', 'data_param'};
            end % if
            
            if simple 
                % data_param
                if (strcmp (param.data_param.inputKernel, 'ALL'))
                    n_kernel = length (param.data_param.availInputKernels);
                else
                    n_kernel = 1;
                end % if
                
                str = '';
                
                if ismember('data_param', include_info)
                    str = strcat(str, sprintf ('ikernel=%s__n=%d__', ...
                        param.data_param.inputKernel, n_kernel));
                end % if
                
                % iokr_param / mp_iokr_param
                if ismember('iokr_param', include_info)
                    center               = num2str (param.iokr_param.center);
                    mkl                  = param.iokr_param.mkl;
                    model_representation = param.iokr_param.model_representation;
                end % if
                
                if ismember('mp_iokr_param', include_info)
                    center               = param.mp_iokr_param.center;
                    mkl                  = param.mp_iokr_param.mkl;
                    model_representation = param.mp_iokr_param.model_representation;
                end % if
                str = strcat (str, sprintf ('cent=%s__mkl=%s__modelrep=%s__', ...
                        center, mkl, model_representation));
                
                % ky_param
                if ismember('ky_param', include_info)
                    str = strcat (str, sprintf ('type=%s__base=%s__parsel=%s', ...
                        param.ky_param.type, param.ky_param.base_kernel, param.ky_param.param_selection));
                end % if
            else
                error ('Not implemented yet')
            end % if
        end % function
    end % methods
end % class