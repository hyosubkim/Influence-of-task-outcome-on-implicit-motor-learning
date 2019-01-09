%%%Analysis of E1 behavioral data from intrinsic rewards project%%%

clear all; close all; clc

load('e1_data.mat');

% default figure properties
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'defaultAxesFontSize',12)

%
subjects = unique(T.SN);
targets = unique(T.ti);
bad_trial = zeros(size(subjects,1),size(targets,1));
ntrials = max(T.TN);
blks_trials = [40 120 760 800];
blks_cycles = blks_trials/8;
fbBslCycles = 7:15;
clampcycles_one_to_ten = 16:25;
earlyclamp = 18:22; %cycles 3-7 of clamp block
prelateclamp1 = 76:85;
lateclamp = 86:95;
lateclamp = 86:95;
final_ae_cycle = 96; %first no-feedback cycle following clamp

hold_t = .500;
endptfb_t = .05;

%aesthetics
tgtclrs=parula(8);
grpclrs=parula(4);
grpclrs([1 2 3],:) = [grpclrs(1,:); grpclrs(3,:); grpclrs(2,:)];
gray=[0.8 0.8 0.8];

%comment next line if analyzing with endpoint hand angle instead of hand
%angle at peak velocity
T.hand_theta = T.hand_theta_maxradv;
T.raw_hand_theta = T.hand_theta;
T.bslc_hand_theta = T.hand_theta;

T.TT = hold_t + endptfb_t + T.ST + T.RT + T.MT;

% loop through all subjects
tic
for sn = 1:length(subjects)
    
    %find individual subject indices 
    subjidx = T.SN==subjects(sn);
    idx = find(subjidx,1);
    grp(sn,1) = T.group(idx);
    ccw = T.cond(idx) > 0;
    raw_hand_angle(sn,:) = T.hand_theta(subjidx);
    
    
    % next loop is for by-target analysis
    for ti = 1:length(targets)
        
        tgt_idx = find(subjidx & T.ti==targets(ti));
        
        %next lines are only used for identifying outliers
        smoothed_hand_theta = smooth(T.hand_theta(tgt_idx),5);
        detrended = T.hand_theta(tgt_idx) - smoothed_hand_theta;
        tgt_SDs(sn,ti) = std(detrended);
        
        % outlier removal
        for i=1:length(detrended)
            if abs(detrended(i)) > 3*tgt_SDs(sn,ti) | abs(T.hand_theta(tgt_idx(i))) > 90
                [T.hand_theta(tgt_idx(i)),T.raw_hand_theta(tgt_idx(i)),T.hand_theta_maxradv(tgt_idx(i)),T.hand_theta_50(tgt_idx(i)),...
                    T.raw_ep_hand_ang(tgt_idx(i)), T.MT(tgt_idx(i)), T.RT(tgt_idx(i)),T.ST(tgt_idx(i)),...
                    T.radvelmax(tgt_idx(i)),T.maxRadDist(tgt_idx(i))] = deal(nan);
                bad_trial(sn,ti) = bad_trial(sn,ti)+1;
            end
        end
        nofb_baseline_tgtErr(sn,ti) = nanmedian(T.raw_hand_theta(tgt_idx(1:5)));
        nofb_baseline_sd(sn,ti) = nanstd(T.hand_theta(tgt_idx(1:5)));
        fb_baseline_tgtErr(sn,ti) = nanmedian(T.raw_hand_theta(tgt_idx(fbBslCycles)));
        fb_baseline_sd(sn,ti) = nanstd(T.hand_theta(tgt_idx(fbBslCycles)));
        
        st_tgt_med(sn,ti) = nanmedian(T.ST(tgt_idx));
        rt_tgt_med(sn,ti) = nanmedian(T.RT(tgt_idx));
        mt_tgt_med(sn,ti) = nanmedian(T.MT(tgt_idx));
        tt_tgt_med(sn,:) = nanmedian(T.TT(tgt_idx));
        
        baseline_bias(sn,ti) = nanmean(T.hand_theta(tgt_idx(fbBslCycles)));
        T.bslc_hand_theta(tgt_idx) = T.hand_theta(tgt_idx) - baseline_bias(sn,ti);
       
        % detrend hand data
        validearlyidx = find(~isnan(T.hand_theta(tgt_idx(clampcycles_one_to_ten))));
        validlateidx = find(~isnan(T.hand_theta(tgt_idx(lateclamp))));
        detrend_early(sn,ti,1:length(validearlyidx)) = detrend(T.hand_theta(tgt_idx(clampcycles_one_to_ten(validearlyidx))));
        detrend_late(sn,ti,1:length(validlateidx)) = detrend(T.hand_theta(tgt_idx(lateclamp(validlateidx))));
        
        earlyhandvar_tgt(sn,ti) = std(detrend_early(sn,ti,:));
        latehandvar_tgt(sn,ti) = std(detrend_late(sn,ti,:));
        
        % flip ccw hand angle
        if ccw
            T.bslc_hand_theta(tgt_idx) = T.bslc_hand_theta(tgt_idx)*-1;
        end
        
        early_plateau(sn,ti) = mean(T.bslc_hand_theta(tgt_idx(prelateclamp1)));
        late_plateau(sn,ti) = mean(T.bslc_hand_theta(tgt_idx(lateclamp)));
        
    end
    
    st(sn,:) = T.ST(subjidx);
    rt(sn,:) = T.RT(subjidx);
    mt(sn,:) = T.MT(subjidx);
    tt(sn,:) = T.TT(subjidx);
       
    % cycle analysis
    for bn = 1:max(T.move_cycle)
        idx = (subjidx & T.move_cycle == bn);
        hand_ang_m(sn,bn) = nanmean(T.bslc_hand_theta(idx));
        mt_cycle(sn,bn) = nanmedian(T.MT(idx));
        rt_cycle(sn,bn) = nanmedian(T.RT(idx));
    end
    
end

% analyze percent outliers
num_trials = max(T.TN);
perc_outlier_ind = sum(bad_trial,2)./num_trials;
perc_outlier_mean = mean(perc_outlier_ind);
perc_outlier_max = max(perc_outlier_ind);

% calculate average baseline variability and errors in reaches across 8 
% targets (to compare big versus small target reaches)
nofb_bsl_var_ind = mean(nofb_baseline_sd,2);
bsl_var_ind = mean(fb_baseline_sd,2);
earlyhandvar = mean(earlyhandvar_tgt,2);
latehandvar = mean(latehandvar_tgt,2);
nofb_bslTgtErr = nanmean(nofb_baseline_tgtErr,2);
bslTgtErr = nanmean(fb_baseline_tgtErr,2);

%point estimates of kinematics during different phases
earlyadapt = nanmean(hand_ang_m(:,earlyclamp),2)/length(earlyclamp);
asymptote = nanmean(hand_ang_m(:,lateclamp),2);
delta_plateau = nanmean(late_plateau,2)-nanmean(early_plateau,2);
final_ae = nanmean(hand_ang_m(:,final_ae_cycle),2);

%grouping variables
[tgtgrp, tgtid] = findgroups(T.ti);
[grpnum,grpid] = findgroups(T.group);
grpid = [grpid(1) grpid(3) grpid(2)];
grplabel = {'Miss','Straddle','Hit'};

%temporal variables
mt_ind_mean = mean(mt_tgt_med,2);
rt_ind_mean = mean(rt_tgt_med,2);
st_ind_mean = mean(st_tgt_med,2);
tt_ind_mean = mean(tt_tgt_med,2);

%bsl RT and MT
mt_bsl = nanmean(mt_cycle(:,fbBslCycles),2);
rt_bsl = nanmean(rt_cycle(:,fbBslCycles),2);

toc

%% Baseline variability

xtick_grp = [1 2 3];

for i=1:3
    figure(323); hold on
    set(gcf,'units','inches','pos',[5 5 2 2]);
    set(gcf,'PaperPositionMode','auto')
    bar(i,mean(bsl_var_ind(grp==grpid(i))),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bsl_var_ind(grp==grpid(i))),nanstderr(bsl_var_ind(grp==grpid(i))),'k','linewidth',2)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,bsl_var_ind(grp==grpid(i)),'.','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:),'markersize',20)
    set(gca,'xtick',[1 2 3],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['Hand angle (',char(176),')'])
    title('Baseline variability')
    
    figure(327); hold on
    set(gcf,'units','inches','pos',[5 5 2 2]);
    set(gcf,'PaperPositionMode','auto')
    bar(i,mean(bslTgtErr(grp==grpid(i))),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bslTgtErr(grp==grpid(i))),nanstderr(bslTgtErr(grp==grpid(i))),'k','linewidth',2)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,bslTgtErr(grp==grpid(i)),'.','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:),'markersize',20)
    set(gca,'xtick',[1 2 3],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['Hand angle (',char(176),')'])
    title('Baseline error')
end

% print(323,'baselineSD','-painters','-dpdf')
% print(327,'baselineErr','-painters','-dpdf')

[p,tbl_bslVar,stats]=anova1(bsl_var_ind,grp)
[p,tbl_bslTgtErr,stats]=anova1(bslTgtErr,grp)

etaSqCalc(tbl_bslVar)
etaSqCalc(tbl_bslTgtErr)

%% Temporal variables analyses

%plot MTs and RTs for entire experiment
for i=1:length(grpid)
    
    figure(919); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'units','inches','pos',[5 5 4 2.5]);
    set(gcf,'PaperPositionMode','auto')
    rectangle('position',[0 0 10 1],'facecolor',[0.8 0.8 0.8 0.5],'edgecolor',gray)
    rtGrp(i)=shadedErrorBar(1:length(rt_cycle),nanmean(rt_cycle(grp==grpid(i),:)),...
        nanstderr(rt_cycle(grp==grpid(i),:)),{'.','markersize',12,'color',grpclrs(i,:)},1);
    xlim([-2 112])
    ylim([0 1])
    plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
        [240 240], ylim, '--k',[250 250], ylim, '--k')
    xlabel('Movement cycle (8 reaches)')
    ylabel('RT (s)')
    
    figure(324); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'units','inches','pos',[5 5 4 2.5]);
    set(gcf,'PaperPositionMode','auto')
    rectangle('position',[0 0 10 0.5],'facecolor',[0.8 0.8 0.8 0.5],'edgecolor',gray)
    mtGrp(i)=shadedErrorBar(1:length(mt_cycle),nanmean(mt_cycle(grp==grpid(i),:)),...
        nanstderr(mt_cycle(grp==grpid(i),:)),{'.','markersize',12,'color',grpclrs(i,:)},1);
    xlim([-2 112])
    ylim([0 0.5])
    plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
        [240 240], ylim, '--k',[250 250], ylim, '--k')
    xlabel('Movement cycle (8 reaches)')
    ylabel('MT (s)')
end
% print(919,'RTs','-painters','-dpdf')
% print(324,'MTs','-painters','-dpdf')


%compare average of median MTs and RTs
[p_mt,tbl_mt,stats_mt] = anova1(mt_ind_mean,grp)
[p_rt,tbl_rt,stats_rt] = anova1(rt_ind_mean,grp)
[p_tt,tbl_tt,stats_tt] = anova1(tt_ind_mean,grp)

%compare baseline MTs and RTs
[p, tbl_bsl_mt,mt_bsl_stats] = anova1(mt_bsl,grp);
[p, tbl_bsl_rt,rt_bsl_stats] = anova1(rt_bsl,grp);

%posthoc test for significant ANOVA
multcompare(rt_bsl_stats)

%calculate standardized effect size
etaSqCalc(tbl_bsl_mt)
etaSqCalc(tbl_bsl_rt)

%posthoc t-test
[h p] = ttest2(rt_ind_mean(grp=='small tgt'),rt_ind_mean(grp=='big tgt'))

%run simple linear regression to see if there is association between RT and
%total learning/aftereffects
for i=1:length(grpid)
%     xreg=[ones(sum(grp==grpid(i)),1) rt_ind_mean(grp==grpid(i))];
%     [b,bint,r,rint,stats] = regress(earlyadapt(grp==grpid(i)),xreg);
%     rt_rate_regb(i,:) = b;
%     rt_rate_r2(i) = stats(1);
%     rt_rate_p(i) = stats(3);
    
    xreg=[ones(sum(grp==grpid(i)),1) rt_bsl(grp==grpid(i))];
    [b,bint,r,rint,stats] = regress(asymptote(grp==grpid(i)),xreg);
    rt_ae_regb(i,:) = b;
    rt_ae_r2(i) = stats(1);
    rt_ae_p(i) = stats(3);    
end

%% Plot individual data

for i=1:length(subjects)
    
    % trial view
    trialidx=find(T.SN==i);
    
%     figure; hold on
%     set(gcf,'paperpositionmode','auto','units','inches','pos',[5 5 7 3.5])
%     ind = gscatter(1:length(trialidx),T.raw_hand_theta(trialidx),tgtgrp(trialidx),...
%         tgtclrs,'.',25)
%     plot(xlim, [0 0],'-k',[blks_trials(1) blks_trials(1)], ylim,'--k',...
%         [blks_trials(2) blks_trials(2)], ylim, '--k',[blks_trials(3) blks_trials(3)],...
%         ylim, '--k',[blks_trials(4) blks_trials(4)], ylim, '--k')
%     legend(ind,{'0','45','90','135','180','225','270','315'})
    
%   % cycle view
    plot(hand_ang_m(i,:),'o','markerfacecolor',grpclrs(grp(i),:))
    plot(xlim, [0 0],'-k',[5 5], ylim,'--k',[15 15], ylim, '--k',[95 95], ylim, '--k',...
        [100 100], ylim, '--k')
    title(['S',num2str(subjects(i)),' group: ',char(grp(i))])
    xlabel('Movement cycle')
    ylabel('Hand angle (deg)')
            
end
   
%% Plot group analyzed hand data
zeroline = zeros(size(hand_ang_m,1));
figure; hold on
set(gcf,'PaperPositionMode','auto','units','inches','pos',[5 3 3.5 2.63]);
rectangle('position',[0 -5 5 50],'facecolor',gray,'edgecolor',gray)
rectangle('position',[95 -5 5 50],'facecolor',gray,'edgecolor',gray)
for i=[1 2 3]
   
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==grpid(i),:)),...
        nanstderr(hand_ang_m(grp==grpid(i),:)),{'.','markersize',12,'color',grpclrs(i,:)},1);
    
end
xlim([-2 112])
ylim([-5 40])
plot(xlim, [0 0],'-k',[5 5], ylim, '--k',[15 15], ylim,'--k',[95 95], ylim, '--k',[100 100], ylim, '--k')
title(['3.5',char(176),' clamp'])
xlabel('Movement cycle (8 reaches)')
ylabel(['Hand angle (',char(176),')'])
legend([grpdata(1).mainLine,grpdata(2).mainLine,grpdata(3).mainLine],...
    {'Small','Medium','Large'},'location','northwest')
legend('boxoff')
%print('e1_learningFxns','-painters','-dpdf')

%% Plot early rate and accumulated learning

for i=1:length(grpid)
    
    asy=figure(122); hold on
    set(gcf,'paperpositionmode','auto','units','inches','pos',[5 5 2 2.63])
    bar(i,mean(asymptote(grp==grpid(i))),'facecolor',gray)
    errbar(i-.1,mean(asymptote(grp==grpid(i))),nanstderr(asymptote(grp==grpid(i))),'k','linewidth',1.5)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,asymptote(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2 3],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['Hand angle (',char(176),')'])
    title('Asymptote')
 
%     aes=figure(128); hold on
%     set(gcf,'paperpositionmode','auto','units','inches','pos',[5 5 2 2.63])
%     bar(i,mean(final_ae(grp==grpid(i))),'facecolor',gray)
%     errbar(i-.1,mean(final_ae(grp==grpid(i))),nanstderr(final_ae(grp==grpid(i))),'k','linewidth',1.5)
%     plot(.2*rand(sum(grp==grpid(i)),1)+i,final_ae(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
%         'markeredgecolor',grpclrs(i,:))
%     set(gca,'xtick',[1 2 3],'xticklabel',grplabel,'xticklabelrotation',45)
%     ylabel(['Hand angle (',char(176),')'])
%     title('Aftereffects')
    
    ea=figure(235); hold on
    set(gcf,'paperpositionmode','auto','units','inches','pos',[5 5 2 2.63])
    bar(i,mean(earlyadapt(grp==grpid(i))),'facecolor',gray)
    errbar(i-.1,mean(earlyadapt(grp==grpid(i))),nanstderr(earlyadapt(grp==grpid(i))),'k','linewidth',1.5)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,earlyadapt(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2 3],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['\DeltaHand angle/trial (',char(176),')'])
    title('Early rate')
    
end

% print('asymptote_3pt5','-painters','-dpdf')
% print(aes,'-painters','-dpdf')
% print('e1_asy',asy,'-painters','-dpdf')
% print('e1_ea',ea,'-painters','-dpdf')

% ANOVAs
[p_ea,tbl_ea,stats_ea] = anova1(earlyadapt,grp)
[p_asy1,tbl_asy1,stats_asy] = anova1(asymptote,grp)
[p_ae,tbl_ae,stats_ae] = anova1(final_ae,grp)

% posthoc tests for significant ANOVAs
c_asy1 = multcompare(stats_asy)
c_ae = multcompare(stats_ae)

%Scheffe contrasts -- for comparing small and med together versus big
MSE = tbl_asy1{3,4};
critF = finv(0.95,tbl_asy1{2,3},tbl_asy1{3,3});
k=3; %num of groups
c1=1/4;
c2=1/4;
c3=-1/2;
if (c1+c2+c3) ~= 0 
    disp('coefficients must sum to 1!')
end
c=[c1;c2;c3];
m1=mean(asymptote(grp=='small tgt'))
sd1=std(asymptote(grp=='small tgt'))
m2=mean(asymptote(grp=='med tgt'))
sd2=std(asymptote(grp=='med tgt'))
m3=mean(asymptote(grp=='big tgt'))
sd3=std(asymptote(grp=='big tgt'))

Yhat=[m1;m2;m3];
[n1 n2 n3]=deal(sum(grp=='small tgt'));
n=[n1;n2;n3];

Lhat = c'*Yhat
critDiff = sqrt((k-1)*critF)*sqrt(MSE*(sum(c.^2./n)))
if Lhat>critDiff
    hyp=1;
elseif Lhat<=critDiff
    hyp=0;
end

%effect sizes: eta^2
etaEarlyAdapt = etaSqCalc(tbl_ea)
etaAsy1 = etaSqCalc(tbl_asy1)

% t-tests
[h1, p1, ci1, stats1] = ttest2(earlyadapt(grp=='small tgt'),earlyadapt(grp=='med tgt'))
[h2, p2, ci2, stats2] = ttest2(earlyadapt(grp=='med tgt'),earlyadapt(grp=='big tgt'))
[h3, p3, ci3, stats3] = ttest2(earlyadapt(grp=='small tgt'),earlyadapt(grp=='big tgt'))

[h1, p1, ci1, stats1] = ttest2(asymptote(grp=='small tgt'),asymptote(grp=='med tgt'))
[h2, p2, ci2, stats2] = ttest2(asymptote(grp=='big tgt'),asymptote(grp=='med tgt'))
[h3, p3, ci3, stats3] = ttest2(asymptote(grp=='big tgt'),asymptote(grp=='small tgt'))

[h4, p4, ci4, stats4] = ttest2(final_ae(grp=='small tgt'),final_ae(grp=='med tgt'))
[h5, p5, ci5, stats5] = ttest2(final_ae(grp=='med tgt'),final_ae(grp=='big tgt'))
[h6, p6, ci6, stats6] = ttest2(final_ae(grp=='small tgt'),final_ae(grp=='big tgt'))

cohensD(asymptote(grp=='small tgt'),asymptote(grp=='med tgt'))
cohensD(asymptote(grp=='big tgt'),asymptote(grp=='med tgt'))
cohensD(asymptote(grp=='big tgt'),asymptote(grp=='small tgt'))


