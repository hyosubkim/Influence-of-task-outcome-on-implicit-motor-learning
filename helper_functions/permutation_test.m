function pval = permutation_test(data, group, test_stat, num_sims)

%%%function for running non-parametric permutation on only two groups

%check to see if group is grouping variable and change to 
if iscategorical(group)
    group = findgroups(group);
end

%create null distributions
for i=1:num_sims
    
    shuffled_idx = group(randperm(length(group)));
    
    null_dist(i,:) = mean(data(shuffled_idx==1,1) - data(shuffled_idx==2,1));

end

figure; hold on
histogram(null_dist)
plot([test_stat test_stat],[ylim],'k')
pval = sum(null_dist > abs(test_stat) | null_dist < -abs(test_stat))/num_sims
title(['Perm test p=',num2str(pval)])
