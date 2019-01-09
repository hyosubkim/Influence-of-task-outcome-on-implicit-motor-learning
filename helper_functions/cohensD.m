function [d] = cohensD(x1,x2)

%This calculates Cohen's d, a measure of effect size when comparing 2
%groups. [d = (M1 - M2)/sqrt((SD1^2 + SD2^2)/2)]. 

diffGrps = mean(x1) - mean(x2);
s1 = std(x1); s2 = std(x2);
s = sqrt((s1^2 + s2^2)/2);
d = diffGrps/s;

end

