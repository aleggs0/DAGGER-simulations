function [pval] = fisher(pvals)
pval = 1-chi2cdf(-2*sum(log(pvals)),2*length(pvals));
end