rng default; close all;
depth=3;
widths=4.^(1:3);
layer_1=1:widths(1);
layer_2=1:widths(2);
layer_3=1:widths(3);
num_false_param=8;

%the results shall be recorded as a 5D array

%the first index refers to the target q
qs=[0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2];

%the second index refers to the method:
methods = cell(1,1);
methods{1} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
    DAGGER_child_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
methods{2} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
    DAGGER_descendant_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
methods{3} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
    DAGGER_child_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
methods{4} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
    DAGGER_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
methods{5} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
    reshaped_DAGGER_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
%methods{6} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
%    focused_BH_child_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
%methods{7} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
%    focused_BH_descendant_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
%methods{8} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
%    focused_BH_child_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
%methods{9} = @(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q)...
%    focused_BH_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);

%the third refers to the way the false hypotheses are chosen
% 1. no false hypotheses
% 2. the first leaves 8 are false
% 3. 8 randomly chosen leaves are false

%the fourth index refers to the quantity:
%FDP,all-nodes power,leaf-nodes power

%the fifth refers to the repeats
numtrials=10000;

results=zeros(length(qs),length(methods),3,3,numtrials);
%in total 2700 times a procedure is used when numtrials=20

for trialno=1:numtrials
    if mod(trialno,100)==0
        disp(trialno)
    end
    %generate a random tree with depth D=3 and blocks of width 4,16,64,
    %ensuring all leaves are at the bottom level
    layer_3_parents=sort(randi(widths(2),size(layer_3)));
    family_sizes_2 = histcounts(layer_3_parents, BinMethod='integers', BinLimits = [1 widths(2)]);
    while ~all(family_sizes_2)
        layer_3_parents=sort(randi(widths(2),size(layer_3)));
        family_sizes_2 = histcounts(layer_3_parents, BinMethod='integers', BinLimits = [1 widths(2)]);
    end
    
    layer_2_parents=sort(randi(widths(1),size(layer_2)));
    family_sizes_1 = histcounts(layer_2_parents, BinMethod='integers', BinLimits = [1 widths(1)]);
    while ~all(family_sizes_1)
        layer_2_parents=sort(randi(widths(1),size(layer_2)));
        family_sizes_1 = histcounts(layer_2_parents, BinMethod='integers', BinLimits = [1 widths(1)]);
    end
    
    for patternno=1:3
        % choose false hypotheses
        if patternno==1
            layer_3_false=[];
        elseif patternno==2
            layer_3_false=1:num_false_param;
        else
            layer_3_false = datasample(layer_3,num_false_param,'Replace',false);
        end

        layer_2_false = unique(layer_3_parents(layer_3_false));
        layer_1_false = unique(layer_2_parents(layer_2_false));
        num_leaves_false = length(layer_3_false);
        total_num_false = length(layer_3_false)+length(layer_2_false)+length(layer_1_false);
        
        %%generate p-values
        pvals_3=rand(size(layer_3));
        pvals_3(layer_3_false)=normcdf(normrnd(0,1,size(layer_3_false))-3);
        
        % run each procedure
        for qno=1:length(qs)
            q=qs(qno);
            for methodno=1:length(methods)
                [R_1,R_2,R_3] = methods{methodno}(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
                R_num=sum(R_1)+sum(R_2)+sum(R_3);
                S=sum(R_1(layer_1_false))+sum(R_2(layer_2_false))+sum(R_3(layer_3_false));
                V=R_num-S;
                results(qno,methodno,patternno,1,trialno)=V/max(R_num,1);
                results(qno,methodno,patternno,2,trialno)=S/total_num_false;
                results(qno,methodno,patternno,3,trialno)=S/num_leaves_false;
            end
        end
    end
end

%=== results processing ===

means=mean(results,5);
standarderrors=std(results,0,5)./sqrt(numtrials);

figure;hold on
plot(qs,qs,'k-')
plot(qs,means(:,1,1,1),'mx-')
plot(qs,means(:,2,1,1),'bx-')
plot(qs,means(:,3,1,1),'x-','color',"#EDB120")
plot(qs,means(:,4,1,1),'gx-')
plot(qs,means(:,5,1,1),'gx--')
legend("target level","child simes","descendant simes","child fisher","descendant fisher","descendant fisher, reshaped",'Location','northwest')
xlabel("target level q"),ylabel("average FDP")

figure;hold on
plot(qs,qs,'k-')
plot(qs,means(:,1,2,1),'mx-')
plot(qs,means(:,2,2,1),'bx-')
plot(qs,means(:,3,2,1),'x-','color',"#EDB120")
plot(qs,means(:,4,2,1),'gx-')
plot(qs,means(:,5,2,1),'gx--')
legend("target level","child simes","descendant simes","child fisher","descendant fisher","descendant fisher, reshaped",'Location','northwest')
xlabel("target level q"),ylabel("average FDP")

figure;hold on
plot(qs,qs,'k-')
plot(qs,means(:,1,3,1),'mx-')
plot(qs,means(:,2,3,1),'bx-')
plot(qs,means(:,3,3,1),'x-','color',"#EDB120")
plot(qs,means(:,4,3,1),'gx-')
plot(qs,means(:,5,3,1),'gx--')
legend("target level","child simes","descendant simes","child fisher","descendant fisher","descendant fisher, reshaped",'Location','northwest')
xlabel("target level q"),ylabel("average FDP")

figure;hold on
plot(qs,means(:,1,2,2),'mx-')
plot(qs,means(:,2,2,2),'bx-')
plot(qs,means(:,3,2,2),'x-','color',"#EDB120")
plot(qs,means(:,4,2,2),'gx-')
plot(qs,means(:,5,2,2),'gx--')
legend("child simes","descendant simes","child fisher","descendant fisher","descendant fisher, reshaped",'Location','northwest')
xlabel("target level q"),ylabel("average power")

figure;hold on
plot(qs,means(:,1,3,2),'mx-')
plot(qs,means(:,2,3,2),'bx-')
plot(qs,means(:,3,3,2),'x-','color',"#EDB120")
plot(qs,means(:,4,3,2),'gx-')
plot(qs,means(:,5,3,2),'gx--')
legend("child simes","descendant simes","child fisher","descendant fisher","descendant fisher, reshaped",'Location','northwest')
xlabel("target level q"),ylabel("average power")

%=== function definitions ===

function [R_1,R_2,R_3] = DAGGER_child_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
[pvals_1,pvals_2] = child_simes(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3);
[R_1,R_2,R_3] = DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
end

function [R_1,R_2,R_3] = DAGGER_descendant_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
[pvals_1,pvals_2] = descendant_simes(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3);
[R_1,R_2,R_3] = DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
end

function [R_1,R_2,R_3] = DAGGER_child_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
[pvals_1,pvals_2] = child_fisher(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3);
[R_1,R_2,R_3] = DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
end

function [R_1,R_2,R_3] = DAGGER_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
[pvals_1,pvals_2] = descendant_fisher(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3);
[R_1,R_2,R_3] = DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
end

function [R_1,R_2,R_3] = reshaped_DAGGER_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
[pvals_1,pvals_2] = descendant_fisher(layer_1,layer_2,layer_3,layer_2_parents,layer_3_parents,pvals_3);
[R_1,R_2,R_3] = reshaped_DAGGER(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_1,pvals_2,pvals_3,q);
end

function [R_1,R_2,R_3] = focused_BH_child_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
R_1=zeros(size(layer_1));
R_2=zeros(size(layer_2));
R_3=zeros(size(layer_3));
end

function [R_1,R_2,R_3] = focused_BH_descendant_simes(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
R_1=zeros(size(layer_1));
R_2=zeros(size(layer_2));
R_3=zeros(size(layer_3));
end

function [R_1,R_2,R_3] = focused_BH_child_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
R_1=zeros(size(layer_1));
R_2=zeros(size(layer_2));
R_3=zeros(size(layer_3));
end

function [R_1,R_2,R_3] = focused_BH_descendant_fisher(layer_1,layer_2,layer_3, layer_2_parents, layer_3_parents,pvals_3,q);
R_1=zeros(size(layer_1));
R_2=zeros(size(layer_2));
R_3=zeros(size(layer_3));
end
