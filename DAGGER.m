function [R_1,R_2,R_3] = DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
%DAGGER Summary of this function goes here
%   Detailed explanation goes here
leaves_3=ones(size(layer_3));
desc_3=ones(size(layer_3));
leaves_2=sum((layer_2==layer_3_parents').*leaves_3',1);
desc_2=ones(size(layer_2))+sum((layer_2==layer_3_parents').*desc_3',1);
leaves_1=sum((layer_1==layer_2_parents').*leaves_2',1);
desc_1=ones(size(layer_1))+sum((layer_1==layer_2_parents').*desc_2',1);
leaves=length(layer_3);
R_num=0;
rs = 0:length(layer_1);
R_1_num = find( rs' <= sum(pvals_1 <= q.*leaves_1./leaves.*(desc_1+rs'+R_num-1)./desc_1,2) , 1, "last")-1;
R_1 = (pvals_1 <= q.*leaves_1./leaves.*(desc_1+R_1_num+R_num-1)./desc_1);
R_num=R_num+R_1_num;
rs = 0:length(layer_2);
R_2_num = find( rs' <= sum(pvals_2 <= R_1(layer_2_parents) .* q.*leaves_2./leaves.*(desc_2+rs'+R_num-1)./desc_2,2) , 1, "last")-1;
R_2 = (pvals_2 <= R_1(layer_2_parents) .* q.*leaves_2./leaves.*(desc_2+R_2_num+R_num-1)./desc_2);
R_num=R_num+R_2_num;
rs = 0:length(layer_3);
R_3_num = find( rs' <= sum(pvals_3 <= R_2(layer_3_parents) .* q.*leaves_3./leaves.*(desc_3+rs'+R_num-1)./desc_3,2) , 1, "last")-1;
R_3 = (pvals_3 <= R_2(layer_3_parents) .* q.*leaves_3./leaves.*(desc_3+R_3_num+R_num-1)./desc_3);
end

