function rankPerc = getRankPerc (ranks, candNum)
%% GETRANKPERC Calculate percentage for the different rankings.
%    INPUT:
%       ranks         ranking for each molecule
%       candNum       number of candidates for each example
%
%    OUTPUT:
%       rankPerc     top k accuracy for each k from 1:max(candNum)
    hasRank = (~ isnan (ranks));
    nValidRanks = sum (hasRank);

    nel = hist (ranks(hasRank), 1:max(candNum(hasRank)));
    rankPerc = cumsum(nel)';
    rankPerc = rankPerc / nValidRanks * 100;
end % function