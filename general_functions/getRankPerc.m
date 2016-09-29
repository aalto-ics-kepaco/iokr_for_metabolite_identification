function rankPerc = getRankPerc (ranks, maxCandNum)
%% GETRANKPERC Calculate percentage for the different rankings.
%    INPUT:
%       ranks         ranking for each molecule
%       maxCandNum    maximum number of candidates for a molecule
%
%    OUTPUT:
%       rankPerc     percantages for rank 1 - 100
    nValidRanks = sum (~ isnan (ranks));

    nel = hist (ranks, 1:max(maxCandNum));
    rankPerc = cumsum(nel)';
    rankPerc = rankPerc / nValidRanks * 100;
end % function