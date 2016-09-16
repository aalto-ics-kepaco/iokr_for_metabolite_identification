%% STATISTICSETTINGS2STR compiles a string given a parameter-structure
function str = statisticSettings2Str (varargin)
    strategyStr    = '';
    percStr        = '';
    inclExpCandStr = '';
    foldOuterStr   = '';
    foldInnerStr   = '';
    otherStr       = '';
%     kStr           = '';
    
    for i = 1:2:(nargin-1)
        if (isempty (varargin{i + 1})) ; continue ; end % if: skip empty paramter
            
        switch (varargin{i})
            case 'strategy'
                strategyStr    = sprintf ('_strategy=%s',    ...
                    varargin{i + 1});
            case 'perc'
                percStr        = sprintf ('_perc=%.4f',      ...
                    varargin{i + 1});
            case 'inclExpCand'
                inclExpCandStr = sprintf ('_inclExpCand=%d', ...
                    varargin{i + 1});
            case 'foldOuter'
                foldOuterStr   = sprintf ('_foldOuter=%d',   ...
                    varargin{i + 1});
            case 'foldInner'
                foldInnerStr   = sprintf ('_foldInner=%d',   ...
                    varargin{i + 1});
%             case 'k'
%                 kStr           = sprintf ('_k=%d',           ...
%                     varargin{i + 1});
            otherwise
                if (isnumeric (varargin{i + 1}) || islogical (varargin{i + 1}))
                    otherStr   = sprintf ('%s_%s=%s',        ...
                        otherStr, varargin{i}, num2str (varargin{i + 1}));
                elseif (ischar (varargin{i + 1})) 
                    otherStr   = sprintf ('%s_%s=%s',        ...
                        otherStr, varargin{i}, varargin{i + 1});
                else
                    assert (0, 'Upps?!: How should I compile a filename using an object of class %s?', ...
                        class (ischar (varargin{i + 1}))
                end % if
        end % switch
    end % for
    

    str = strcat (strategyStr, inclExpCandStr, percStr, foldOuterStr, foldInnerStr, otherStr);    
%     str = strcat (strategyStr, usePredGNPSStr, percStr, kStr, foldOuterStr, foldInnerStr, otherStr);
end % function 