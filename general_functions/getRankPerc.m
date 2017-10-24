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

    if     numel (candNum) == numel (ranks)
        maxCandNum = max(candNum(hasRank));
    elseif numel (candNum) == 1
        maxCandNum = candNum;
    else
        error ('getRankPerc:InvalidArgument', ...
            '"candNum" must be of same length as "ranks" or one.');
    end % if
            
    
    nel = hist (ranks(hasRank), 1:maxCandNum);
    rankPerc = cumsum(nel)';
    rankPerc = rankPerc / nValidRanks * 100;
end % function