for KX_name = names
    [KX, header] = loadKernel(strcat (KX_name{1}, '.txt'));
    saveKernel (strcat (KX_name{1}, '.mat'), KX, header);
end % for