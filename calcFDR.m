function [adjP]=calcFDR(pval)
% the function 'calcFDR' performs a multiple comparison correction using the
% Benjamini-Hochbergs Linear Step-Up Procedure (FDR)
    
% Input Arguments:
%     pval  = list of p-values

% Output Arguments:
%     adjP = FDR-adjusted p-values

    nt = numel(pval);   % number of tests 
    [p_sorted, pIdx] = sort(pval(:));
    adjP_sorted =  p_sorted .* (nt./(1:nt))';
    adjP_sorted = cummin(adjP_sorted, 'reverse');
    adjP(pIdx) = adjP_sorted;    
end
