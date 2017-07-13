names =  {'ALIGND', 'ALIGN', 'CEC', 'CP2PLUS', 'CP2', 'CPC',   ...
          'CPI', 'CPJB', 'CPJ', 'CPK', 'CSC', 'FIPP', 'LB', 'LC', 'LIPP', 'LI',     ...
          'LW', 'NB', 'NI', 'NSF', 'PPKR', 'RLB', 'RLI', 'WPC'};                   ...

for KX_name = names
    if (exist (strcat (KX_name{1}, '.mat'), 'file'))
        continue;
    end % if
    
%     [KX, header] = loadKernel(strcat (KX_name{1}, '.txt'), 0);
    KX= loadKernel(strcat (KX_name{1}, '.txt'), 0);
    saveKernel (strcat (KX_name{1}, '.mat'), KX);
end % for