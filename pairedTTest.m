function [tval] = pairedTTest(x,y)
    
% the function pairedTTest computes a paired-sample t-test.

% Input Arguments:
%    x = vector with the samples of the first distribution
%    y = vector with the samples of the second distribution, matched to the
%    first

% Output Arguments:
%    tval = t-values


    dim = 1;

    sampleSize = size(x,dim); 

    xy_diff = x - y;

    diff_mean = nanmean(xy_diff,dim);
    diff_sd = nanstd(xy_diff,[],dim);

    diff_stdErr = diff_sd ./ sqrt(sampleSize);

    tval = diff_mean ./ diff_stdErr;
    
end