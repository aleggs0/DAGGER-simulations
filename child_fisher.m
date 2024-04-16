function [pvals_1,pvals_2] = child_fisher(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3)
%SIMES Summary of this function goes here
%   Detailed explanation goes here
pvals_2=zeros(size(layer_2));
pvals_1=zeros(size(layer_1));
for i=layer_2
    pvals_family=pvals_3(layer_3_parents==i);
    pvals_2(i) = fisher(pvals_family);
end
for i=layer_1
    pvals_family=pvals_2(layer_2_parents==i);
    pvals_1(i) = fisher(pvals_family);
end
end


