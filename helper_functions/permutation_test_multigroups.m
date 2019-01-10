function pval = permutation_test_multigroups(data, group, num_sims)

%%%function for running non-parametric permutation on more than two groups

%check to see if group is grouping variable and change to 
if iscategorical(group)
    group = findgroups(group);
end

[~,tbl] = anova1(data, group,'off');
test_stat = tbl{2,5};

%create null distributions
for i=1:num_sims
    
    shuffled_idx = group(randperm(size(group,1)));
    [~,tbl] = anova1(data, shuffled_idx, 'off');
    null_dist(i) = tbl{2,5};
    
end

figure; hold on
histogram(null_dist)
plot([test_stat test_stat],[ylim],'k')
pval = sum(null_dist > test_stat)/num_sims;
title(['Perm test p=',num2str(pval)])
