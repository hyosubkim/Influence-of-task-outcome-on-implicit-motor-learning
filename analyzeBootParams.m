clear all; close all; clc

load e3_4paramAM_boot_params
load e3_IA_params_boot_params

%sort parameter estimates
am4_sort = sort(bestparams_4paramAM_boot);
ia_sort = sort(bestparams_IA_boot);

%establish upper and lower bounds
lb = .025*1000;
ub = .975*1000;

%CIs for parameters
am4_ci = [am4_sort(lb,:); am4_sort(ub,:)]
ia_ci = [ia_sort(lb,:); ia_sort(ub,:)]

