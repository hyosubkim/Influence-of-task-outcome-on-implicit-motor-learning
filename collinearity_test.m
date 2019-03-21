clear all; close all; clc

% decide whether to print figs or not
printfigs = 0;

% Cols: Amiss, Ahit, Umiss, Uhit
load('e3_4paramAM_boot_params.mat')

% Cols: Aspe, Atpe, Uspe, Utpe
load('e3_IA_params_boot_params.mat')

if printfigs
    cd('/Users/hyosubkim/Dropbox/Projects/inProgress/targetSize/tgtSwitch/Figures')
end

% Place params into tables
for i=1:size(bestparams_4paramAM_boot,2)
    T_am(:,i) = table(bestparams_4paramAM_boot(:,i));
    T_ia(:,i) = table(bestparams_IA_boot(:,i));
end

T_am.Properties.VariableNames = {'A_miss','A_hit','U_miss','U_hit'};
T_ia.Properties.VariableNames = {'A_spe','A_te','U_spe','U_te'};

[r_am,pvals]=corrplot(T_am,'tail','both')
[r_ia,pvals]=corrplot(T_ia,'tail','both')

figure('position', [100 100 400 300]); hold on 
imagesc(r_am);
%imagesc places origin in bottom left; must flip rows (ie, bottom-to-top)
set(gca,'ydir','reverse')
colorbar
box on
xticks([1:4])
xticklabels({'A_{miss}','A_{hit}','U_{miss}','U_{hit}'})
yticks([1:4])
yticklabels({'A_{miss}','A_{hit}','U_{miss}','U_{hit}'}) %these are flipped 
title('Adaptation Modulation','fontweight','normal','fontsize',14)
if printfigs
    print('heatmap_am','-dpdf','-r300','-painters')
end

figure('position', [100 100 400 300]); hold on 
%imagesc places origin in bottom left; must flip rows (ie, bottom-to-top)
imagesc(r_ia);
set(gca,'ydir','reverse')
colorbar
box on
xticks([1:4])
xticklabels({'A_{spe}','A_{te}','U_{spe}','U_{te}'})
yticks([1:4])
yticklabels({'A_{spe}','A_{te}','U_{spe}','U_{te}'})
title('Dual Error','fontweight','normal','fontsize',14)
if printfigs
    print('heatmap_ia','-dpdf','-r300','-painters')
end





