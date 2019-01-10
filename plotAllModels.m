clear all; close all; clc

% default figure properties
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'defaultAxesFontSize',12)

load tgtSwitch_data
load e3_4paramAM_params_groupfit
load e3_IA_params_groupfit
load e3_popnVectorMR_params_groupfit
load e3_am_mr_hybrid_params_groupfit
load e3_hybridIA_params_groupfit


%decide whether or not to print figures
printfigs = 0;

clamptrials = 11:210;
rwd_sth = [zeros(1,10) zeros(1,120) ones(1,80) zeros(1,10)];
rwd_hts = [zeros(1,10) ones(1,120) zeros(1,80) zeros(1,10)];
rwd = [rwd_sth(clamptrials); rwd_hts(clamptrials)];

group = [1; 2];
grpclrs = parula(4);
grpclrs = [grpclrs(1,:); grpclrs(3,:)];
gray = [0.75 0.75 0.75];
clamp_angle = 1.75;

rotations = 1.75;
 
for i=1:2
    
    sims_popnVectMR(i,:) = popnVectorMR(bestparams_popnVectorMR,groupavg_hand_data(i,11:210),clamp_angle,group(i),rwd(i,:));
    sims_4AM(i,:) = adaptationModulation_4params(bestparams_4paramAM,groupavg_hand_data(i,11:210),clamp_angle,group(i),rwd(i,:));
    sims_IA(i,:) = implicitAim(bestparams_IA,groupavg_hand_data(i,11:210),clamp_angle,group(i), rwd(i,:));
    sims_am_mr_hybrid(i,:) = am_mr_hybrid_simulator(bestparams_am_mr_hybrid,groupavg_hand_data(i,11:210),clamp_angle,group(i),rwd(i,:));    
    sims_hybridIA(i,:) = hybridImplicitAim_simulator(bestparams_hybridIA,groupavg_hand_data(i,11:210),clamp_angle,group(i),rwd(i,:));
    
end


figure 
set(gcf,'units','inches','pos',[5 5 5.5 12]);
set(gcf,'PaperPositionMode','auto')
subplot(3,1,1); hold on
xlim([-2 202])
ylim([-5 35])
line(xlim,[0 0],'color','k')
s2h = shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'.','markersize',10,'color',grpclrs(1,:)},1)
h2s = shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'.','markersize',10,'color',grpclrs(2,:)},1)
model = plot(sims_popnVectMR(1,:),'k','linewidth',2)
plot(sims_popnVectMR(2,:),'k','linewidth',2)
plot([121 121], ylim, '--k')
%xlabel('Movement cycle (8 reaches)')
ylabel(['Hand angle (',char(176),')'])
title('Movement Reinforcement model','fontsize',14)
legend([s2h.mainLine, h2s.mainLine, model],{'Straddle-to-Hit','Hit-to-Straddle','Model'},...
    'location','northwest','box','off')

subplot(3,1,2); hold on
%set(gcf,'units','inches','pos',[5 5 6 4],'paperpositionmode','auto')
xlim([-2 202])
ylim([-5 35])
line(xlim,[0 0],'color','k')
s2h = shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'.','markersize',10,'color',grpclrs(1,:)},1)
h2s = shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'.','markersize',10,'color',grpclrs(2,:)},1)
model = plot(sims_4AM(1,:),'k','linewidth',2)
plot(sims_4AM(2,:),'k','linewidth',2)
plot([121 121], ylim, '--k')
%xlabel('Movement cycle')
ylabel(['Hand angle (',char(176),')'])
title('Adaptation Modulation model','fontsize',14)
legend([s2h.mainLine, h2s.mainLine, model],{'Straddle-to-Hit','Hit-to-Straddle','Model'},...
    'location','northwest','box','off')

subplot(3,1,3); hold on
xlim([-2 202])
ylim([-5 35])
line(xlim,[0 0],'color','k')
s2h = shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'.','markersize',10,'color',grpclrs(1,:)},1)
h2s = shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'.','markersize',10,'color',grpclrs(2,:)},1)
model = plot(sims_IA(1,:),'k','linewidth',2)
plot(sims_IA(2,:),'k','linewidth',2)
plot([121 121], ylim, '--k')
xlabel('Movement cycle (8 reaches)')
ylabel(['Hand angle (',char(176),')'])
title('Dual Error model','fontsize',14)
legend([s2h.mainLine, h2s.mainLine, model],{'Straddle-to-Hit','Hit-to-Straddle','Model'},...
    'location','northwest','box','off')
if printfigs
    print('e3_model_fits','-painters','-dpdf','-r300')
end

figure; hold on
set(gcf,'units','inches','pos',[5 5 6 4],'paperpositionmode','auto')
xlim([-2 202])
line(xlim,[0 0],'color',gray)
shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'color',grpclrs(1,:)},1)
shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'color',grpclrs(2,:)},1)
plot(sims_am_mr_hybrid(1,:),'k','linewidth',3)
plot(sims_am_mr_hybrid(2,:),'k','linewidth',3)
xlabel('Movement cycle')
ylabel('Hand angle (deg)')
title('AM + MR hybrid','fontsize',14)
if printfigs
    print('MRandAMhybrid','-dtiff','-r300')
end

figure; hold on
set(gcf,'units','inches','pos',[5 5 6 4],'paperpositionmode','auto')
xlim([-2 202])
line(xlim,[0 0],'color',gray)
shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'color',grpclrs(1,:)},1)
shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'color',grpclrs(2,:)},1)
plot(sims_hybridIA(1,:),'k','linewidth',3)
plot(sims_hybridIA(2,:),'k','linewidth',3)
xlabel('Movement cycle')
ylabel('Hand angle (deg)')
title('IA + MR hybrid','fontsize',14)
if printfigs
    print('IAandMRhybrid','-dtiff','-r300')
end



%plot AICs
figure; hold on
set(gcf,'units','inches','pos',[5 5 6 5],'paperpositionmode','auto')
m = [aic_popnVectorMR; aic_4paramAM; aic_IA; aic_am_mr_hybrid; aic_hybridIA]
bar(m)
set(gca,'xtick',1:5,'xticklabel',{'MR', 'AM','IA','AM + MR hybrid'},'xticklabelrotation',45)
title('AIC scores','fontsize',14)
xlabel('Model','fontsize',12)
ylabel('Score','fontsize',12) 
if printfigs
    print('aics','-dtiff','-r300')
end
