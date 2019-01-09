function [etaSq] = etaSqCalc(anovaTable)
%This fxn calculates etaSquared values--a measure of effect size in ANOVAs
%The inputs/outputs are very specific to the way Matlab outputs an ANOVA
%table

ssb = anovaTable{2,2};
sst = anovaTable{4,2};
etaSq = ssb/sst;

end

