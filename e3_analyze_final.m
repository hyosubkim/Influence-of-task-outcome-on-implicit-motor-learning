%%%Analysis of E3 behavioral data from intrinsic rewards project%%% 

clear all; close all; clc

load('e3_data.mat')

% default figure properties
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.0)

%%%Group key:
%%%grp1=(SB: small tgt to big tgt); grp2=(BS: big tgt to small tgt);
%%%grp3=(LN: lined tgt to no-lined tgt)

subjects = unique(T.SN);
targets = unique(T.ti);
bad_trial = zeros(size(subjects,1),size(targets,1));
ntrials = max(T.TN);
blks_trials = [40 80 1040 1680];
blks_cycles = blks_trials/8;
bsl_cycle = 10;
fbBslCycles = 7:10;
earlyclamp = 13:17; %cycles 3-7 of clamp block
lateclamp1 = 121:130;
lateclamp2 = 201:210;

% %for movement velocity analyses
% bsltrials = 41:80; %cycles 6-10
% earlyphase1 = 81:160; %cycles 11-20
% latephase1 = 961:1040; %cycles 121-130
% latephase2 = 1601:1680; %cycles 201-210

hold_t = .500;
endptfb_t = .05;

%aesthetics
tgtclrs=parula(8);
grpclrs = parula(4);
grpclrs(2:3,:) = [grpclrs(3,:); [1 0 1]];
gray=[0.8 0.8 0.8];

%comment next line if analyzing with endpoint hand angle instead of hand
%angle at peak velocity
T.hand_theta = T.hand_theta_maxradv;
T.raw_hand_theta = T.hand_theta;
T.bslc_hand_theta = T.hand_theta;

% T.radvelmax = T.radvelmax./10;
% T.maxRadDist = T.maxRadDist./10;

T.TT = hold_t + endptfb_t + T.ST + T.RT + T.MT;

% loop through all subjects
tic
for sn = 1:length(subjects)
    
    %find individual subject indices
    subjidx = T.SN == subjects(sn);
    idx = find(subjidx,1);
    grp(sn) = T.group(idx);
    ccw = T.cond(idx) > 0;
    raw_hand_angle(sn,:) = T.hand_theta(subjidx);
    subject_code(sn,1) = sn;
   
   
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
        fb_baseline_tgtErr(sn,ti) = nanmedian(T.raw_hand_theta(tgt_idx(fbBslCycles)));
        fb_baseline_sd(sn,ti) = nanstd(T.hand_theta(tgt_idx(fbBslCycles)));
        
        st_tgt_med(sn,ti) = nanmedian(T.ST(tgt_idx));
        rt_tgt_med(sn,ti) = nanmedian(T.RT(tgt_idx));
        mt_tgt_med(sn,ti) = nanmedian(T.MT(tgt_idx));
        tt_tgt_med(sn,:) = nanmedian(T.TT(tgt_idx));
        
        baseline_bias(sn,ti) = nanmean(T.hand_theta(tgt_idx(fbBslCycles)));
        T.bslc_hand_theta(tgt_idx) = T.hand_theta(tgt_idx) - baseline_bias(sn,ti);
        
        % flip ccw hand angle
        if ccw
            T.bslc_hand_theta(tgt_idx) = T.bslc_hand_theta(tgt_idx)*-1;
        end
        
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

% calculate average baseline variability in reaches across 8 targets (to 
% compare big versus small target reaches)
bsl_var_ind = nanmean(fb_baseline_sd,2);
bslTgtErr = nanmean(fb_baseline_tgtErr,2);

%point estimates of kinematics during different phases
earlyadapt = nanmean(hand_ang_m(:,earlyclamp),2)/length(earlyclamp);
asymptote1 = nanmean(hand_ang_m(:,lateclamp1),2);
asymptote2 = nanmean(hand_ang_m(:,lateclamp2),2);
delta_asymptote = asymptote2 - asymptote1;

%grouping variables
[tgtgrp, tgtid] = findgroups(T.ti);
[~,grpid] = findgroups(T.group);
grplabel = {'Small/Big','Big/Small','Line/No line'};

%temporal variables
mt_ind_mean = mean(mt_tgt_med,2);
rt_ind_mean = mean(rt_tgt_med,2);
st_ind_mean = mean(st_tgt_med,2);
tt_ind_mean = mean(tt_tgt_med,2);

%bsl RT and MT
mt_bsl = nanmean(mt_cycle(:,fbBslCycles),2);
rt_bsl = nanmean(rt_cycle(:,fbBslCycles),2);

toc

%% optional: save data in .mat file

handAngle_s2b = hand_ang_m(grp==1,:);
handAngle_b2s = hand_ang_m(grp==2,:);

cd('../Modeling')
save('tgtSwitch_data','handAngle_s2b','handAngle_b2s','grp')

%% Baseline target error and variability

%focus comparisons on Straddle-to-Hit vs Hit-to-Straddle (ie, groups 1 and
%2)
for i=1:2
    figure(323); hold on
    bar(i,mean(bsl_var_ind(grp==i)),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bsl_var_ind(grp==i)),nanstderr(bsl_var_ind(grp==i)),'k','linewidth',2)
    plot(.2*rand(sum(grp==i),1)+i,bsl_var_ind(grp==i),'o','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2],'xticklabel',{'Small first','Big first'})
    title('Baseline variance')
    
    figure(327); hold on
    bar(i,mean(bslTgtErr(grp==i)),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bslTgtErr(grp==i)),nanstderr(bslTgtErr(grp==i)),'k','linewidth',2)
    plot(.2*rand(sum(grp==i),1)+i,bslTgtErr(grp==i),'o','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2],'xticklabel',{'Small first','Big first'})
    title('Baseline mean target error')
    
end

%perform omnibus ANOVA including all 3 groups
[p_bsl,tbl_bsl,stats_bsl] = anova1(bsl_var_ind,grp)
etaSqCalc(tbl_bsl)

[h p ci stats] = ttest2(bslTgtErr(grp==1),bslTgtErr(grp==2))
[h p ci stats] = ttest2(bsl_var_ind(grp==1),bsl_var_ind(grp==2))

mean(bslTgtErr(grp==1))
nanstderr(bslTgtErr(grp==1))

mean(bslTgtErr(grp==2))
nanstderr(bslTgtErr(grp==2))

mean(bsl_var_ind(grp==1))
nanstderr(bsl_var_ind(grp==1))

mean(bsl_var_ind(grp==2))
nanstderr(bsl_var_ind(grp==2))


%% Temporal variables analyses

%plot MTs and RTs for entire experiment
for i=1:length(grpid)
    
    figure(111); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'units','inches','pos',[5 5 4.25 2.5]);
    rtfig(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(rt_cycle(grp==i,:)),...
        nanstderr(rt_cycle(grp==i,:)),{'.','markersize',10,'color',grpclrs(i,:)},1);
    xlim([-2 222])
    plot([blks_cycles(1) blks_cycles(1)],ylim,'--k',[blks_cycles(2) blks_cycles(2)], ylim,'--k',...
        [blks_cycles(3) blks_cycles(3)], ylim, '--k',[blks_cycles(4) blks_cycles(4)], ylim, '--k')
    xlabel('Movement cycle (8 reaches)')
    ylabel('RT (s)')
    
    figure(912); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'units','inches','pos',[5 5 4.25 2.5]);
    mtfig(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(mt_cycle(grp==i,:)),...
        nanstderr(mt_cycle(grp==i,:)),{'.','markersize',10,'color',grpclrs(i,:)},1);
    xlim([-2 222])
    plot([blks_cycles(1) blks_cycles(1)],ylim,'--k',[blks_cycles(2) blks_cycles(2)], ylim,'--k',...
        [blks_cycles(3) blks_cycles(3)], ylim, '--k',[blks_cycles(4) blks_cycles(4)], ylim, '--k')
    xlabel('Movement cycle (8 reaches)')
    ylabel('MT (s)')
    
end

%compare average of median MTs and RTs -- measure is confounded bc of transfer
[h p ci stats] = ttest2(rt_ind_mean(grp==1),rt_ind_mean(grp==2))
[h p ci stats] = ttest2(mt_ind_mean(grp==1),mt_ind_mean(grp==2))

%compare baseline MTs and RTs
[h p ci stats] = ttest2(rt_bsl(grp==1),rt_bsl(grp==2))
[h p ci stats] = ttest2(mt_bsl(grp==1),mt_bsl(grp==2))


mean(rt_bsl(grp==1))
nanstderr(rt_bsl(grp==1))

mean(rt_bsl(grp==2))
nanstderr(rt_bsl(grp==2))

mean(mt_bsl(grp==1))
nanstderr(mt_bsl(grp==1))

mean(mt_bsl(grp==2))
nanstderr(mt_bsl(grp==2))

% effect size calculations
dRT = cohensD(rt_bsl(grp==1),rt_bsl(grp==2))
dMT = cohensD(mt_bsl(grp==1),mt_bsl(grp==2))

%compare baseline MTs and RTs
[p, tbl_bsl_mt,mt_bsl_stats] = anova1(mt_bsl,grp);
[p, tbl_bsl_rt,rt_bsl_stats] = anova1(rt_bsl,grp);

%calculate standardized effect size
etaSqCalc(tbl_bsl_mt)
etaSqCalc(tbl_bsl_rt)


%% Plot individual data

for i=1:length(subjects)
    
    % trial view
    trialidx=find(T.SN==i);
    
%     figure
%     hold on
%     ind = gscatter(1:length(trialidx),T.raw_hand_theta(trialidx),tgtgrp(trialidx),...
%         tgtclrs,'.',25)
%     legend(ind,{'0','45','90','135','180','225','270','315'})
    
    % cycle view
    figure; hold on
    plot(hand_ang_m(i,:),'o','markerfacecolor',grpclrs(grp(i),:))
    plot([blks_cycles(3) blks_cycles(3)],ylim,'--k')
    title(['S',num2str(i),' grp: ',num2str(grp(i))])
    xlabel('Movement cycle')
    ylabel('Hand angle (deg)')

end


%% Plot group analyzed hand data
zeroline = zeros(size(hand_ang_m,1));
figure; hold on
% set(gcf,'units','inches','pos',[5 5 4.25 2.5]);
set(gcf,'units','inches','pos',[5 5 5.5 3]);
set(gcf,'PaperPositionMode','auto')
rectangle('position',[0 -5 5 45],'facecolor',gray,'edgecolor',gray)
for i=[1 2]
   
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==i,:)),...
        nanstderr(hand_ang_m(grp==i,:)),{'.','markersize',10,'color',grpclrs(i,:)},1);
    
end
xlim([-2 222])
ylim([-5 35])
plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[130 130], ylim, '--k',[210 210], ylim, '--k')
xlabel('Movement cycle (8 reaches)')
ylabel(['Hand angle (',char(176),')'])
title(['1.75',char(176),' clamp with target switch'],'fontsize',12)
% legend([grpdata(1).mainLine,grpdata(2).mainLine,grpdata(3).mainLine],...
%     {grplabel{1},grplabel{2},grplabel{3}},'location','northwest')
% print('bigAndBisect','-painters','-dpdf','-r300')
% print('all_grps','-dtiff','-r300')
% print('e3_learning','-r300','-painters','-dpdf')


%plot control group
zeroline = zeros(size(hand_ang_m,1));
figure; hold on
set(gcf,'PaperPositionMode','auto')
% set(gcf,'units','inches','pos',[5 5 4.25 2.5]);
set(gcf,'units','inches','pos',[5 5 5.5 3]);
rectangle('position',[0 -5 5 45],'facecolor',gray,'edgecolor',gray)
for i=[3]
    
    shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==2,:)),...
        nanstderr(hand_ang_m(grp==2,:)),{'.','markersize',10,'color',grpclrs(2,:)},1);
%     shadedErrorBar(1:130,nanmean(hand_ang_m(grp==2,1:130)),...
%         nanstderr(hand_ang_m(grp==2,1:130)),{'.','markersize',10,'color',grpclrs(2,:)},1);
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==i,:)),...
        nanstderr(hand_ang_m(grp==i,:)),{'.','markersize',10,'color',grpclrs(3,:)},1);
    
end
xlim([-2 222])
ylim([-5 35])
plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[130 130], ylim, '--k',[210 210], ylim, '--k')
xlabel('Movement cycle (8 reaches)')
ylabel(['Hand angle (',char(176),')'])
title(['1.75',char(176),' clamp with target switch'],'fontsize',12)

% cd('../../Figures')
% print('bisected_tgt_group','-dpdf','-painters')

%% Plot initial rate and accumulated learning

%main comparison between two groups
for i=1:2
    
    fig_ea = figure(984); hold on
    set(gcf,'PaperPositionMode','auto','units','inches','pos',[5 5 2 2.5]);
    bar(i,mean(earlyadapt(grp==i)),'facecolor',gray)
    errbar(i-.1,mean(earlyadapt(grp==i)),nanstderr(earlyadapt(grp==i)),'k','linewidth',1.5)
    plot(.2*rand(sum(grp==i),1)+i,earlyadapt(grp==i),'o','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2 3],'xticklabel',grplabel)
    set(gca,'xtick',[1 2 ],'xticklabel',{'Straddle','Hit'},'xticklabelrotation',45)
    ylabel(['\DeltaHand angle/trial(',char(176),')'])
    title('Early rate')
    
%     figure(121); hold on
%     set(gcf,'PaperPositionMode','auto','units','inches','pos',[5 5 2 3]);
%     bar(i,mean(asymptote1(grp==i)),'facecolor',gray)
%     errbar(i-.1,mean(asymptote1(grp==i)),nanstderr(asymptote1(grp==i)),'k','linewidth',1.5)
%     plot(.2*rand(sum(grp==i),1)+i,asymptote1(grp==i),'o','markerfacecolor',grpclrs(i,:),...
%         'markeredgecolor',grpclrs(i,:))
%     set(gca,'xtick',[1 2 3],'xticklabel',grplabel)
%     set(gca,'xtick',[1 2 ],'xticklabel',grplabel,'xticklabelrotation',45)
%     ylabel(['Hand angle (',char(176),')'])
%     title('Phase one asymptote')
%     
%     figure(131); hold on
%     set(gcf,'PaperPositionMode','auto','units','inches','pos',[5 5 2 3]);
%     bar(i,mean(asymptote2(grp==i)),'facecolor',gray)
%     errbar(i-.1,mean(asymptote2(grp==i)),nanstderr(asymptote2(grp==i)),'k','linewidth',1.5)
%     plot(.2*rand(sum(grp==i),1)+i,asymptote2(grp==i),'o','markerfacecolor',grpclrs(i,:),...
%         'markeredgecolor',grpclrs(i,:))
%     set(gca,'xtick',[1 2 3],'xticklabel',grplabel)
%     set(gca,'xtick',[1 2 ],'xticklabel',grplabel,'xticklabelrotation',45)
%     ylabel(['Hand angle (',char(176),')'])    
%     title('Phase two asymptote')

end


for i=1:2
    
    % plot within subjects changes from Phase One to Phase Two
    transfer = figure(303)
    hold on
    set(gcf,'units','inches','pos',[5 3 2.5 3]);
    set(gcf,'PaperPositionMode','auto');
    set(gca,'XTick',[1: 2*i],'XTickLabel',{'Acquisition', 'Transfer'},'xticklabelrotation',45)
    title(['\DeltaAsymptote'])
    xlim([0.5 4.5])
    ylim([0 40])
    ylabel(['Hand angle (',char(176),')'])
    plot([2*i-1, 2*i],[asymptote1(grp==i) asymptote2(grp==i)],'color',grpclrs(i,:))
    plot(repmat(2*i-1,sum(grp==i),1),asymptote1(grp==i),'.','markerfacecolor',grpclrs(i,:),'markeredgecolor',grpclrs(i,:),'markersize',18)
    plot(repmat(2*i,sum(grp==i),1),asymptote2(grp==i),'.','markerfacecolor',grpclrs(i,:),'markeredgecolor',grpclrs(i,:),'markersize',18)
    plot([2*i-1.3, 2*i-.70],[mean(asymptote1(grp==i)) mean(asymptote1(grp==i))],'color',grpclrs(i,:),'linewidth',1.5)
    plot([2*i-.3, 2*i+.3],[mean(asymptote2(grp==i)) mean(asymptote2(grp==i))],'color',grpclrs(i,:),'linewidth',1.5)
    line([2.5 2.5],[0 40],'Color',[0.5,0.5,0.5],'LineStyle','--')
    %line([4.5 4.5],[0 40],'Color',[0.5,0.5,0.5],'LineStyle','--')
    
end

%control group data
for i=3
    
    % plot within subjects changes from Phase One to Phase Two
    transfer_ctl = figure(121)
    hold on
    set(gcf,'units','inches','pos',[5 3 2 2.5]);
    set(gcf,'PaperPositionMode','auto');
    set(gca,'XTick',[1: 2],'XTickLabel',{'Acquisition', 'Transfer'},'xticklabelrotation',45)
    title(['\DeltaAsymptote'])
    xlim([0.5 2.5])
    ylim([0 40])
    ylabel(['Hand angle (',char(176),')'])
    plot([1, 2],[asymptote1(grp==i) asymptote2(grp==i)],'color','m')
    plot(repmat(1,sum(grp==i),1),asymptote1(grp==i),'.','markerfacecolor','m','markeredgecolor','m','markersize',18)
    plot(repmat(2,sum(grp==i),1),asymptote2(grp==i),'.','markerfacecolor','m','markeredgecolor','m','markersize',18)
    plot([.7, 1.3],[mean(asymptote1(grp==i)) mean(asymptote1(grp==i))],'color','m','linewidth',1.5)
    plot([1.7, 2.3],[mean(asymptote2(grp==i)) mean(asymptote2(grp==i))],'color','m','linewidth',1.5)
    
end

%%% uncomment if printing figs
% cd('../../Figures')
% print(fig_ea,'fig_ea','-dpdf','-painters')
% print(transfer,'change_in_asymptote','-dpdf','-r300','-painters')
% print(transfer_ctl,'ctl_delta_asymptote','-dpdf','-r300','-painters')


%%% Mixed-model ANOVA
% X: design matrix with four columns (future versions may allow different input configurations)
%     - first column  (i.e., X(:,1)) : all dependent variable values
%     - second column (i.e., X(:,2)) : between-subjects factor (e.g., subject group) level codes (ranging from 1:L where 
%         L is the # of levels for the between-subjects factor)
%     - third column  (i.e., X(:,3)) : within-subjects factor (e.g., condition/task) level codes (ranging from 1:L where 
%         L is the # of levels for the within-subjects factor)
%     - fourth column (i.e., X(:,4)) : subject codes (ranging from 1:N where N is the total number of subjects)
X(:,1) = [asymptote1; asymptote2];
X(:,2) = repmat(grp',[2,1]);
X(:,3) = [ones(length(asymptote1),1); ones(length(asymptote2),1)*2];
X(:,4) = repmat(subject_code,[2,1]);
[SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(X,0);
% 
% % see https://stardock.cs.virginia.edu/empirical/resources/Brown28.pdf for
% % effect sizes on mixed model ANOVA
% SSQs=cell2mat(SSQs);
% TSE = sum(SSQs);
% etasq_btwn = SSQs(1)/TSE; 
% etasq_within = SSQs(3)/TSE;
partial_etasq_int = SSQs{4}/(SSQs{4}+SSQs{5})

% Significant interaction ==> posthoc one-way ANOVAs
[p_asy1,tbl_asy1,stats_asy1] = anova1(asymptote1,grp)
[p_asy2,tbl_asy2,stats_asy2] = anova1(asymptote2,grp)
[p_delta_asy,tbl_delta_asy,stats_delta_asy] = anova1(delta_asymptote,grp)

% posthoc tests for significant ANOVAs
c_asy1 = multcompare(stats_asy1)
c_delta = multcompare(stats_delta_asy)

etaSqCalc(tbl_asy1)

% within-subjects t-tests
for i=1:3
    
  [h_ws(i,1),p_ws(i,1),ci_ws(i,:),stats_ws(i,:)] = ttest(asymptote2(grp==i),asymptote1(grp==i))
    
end

%effect size calculations
[h,p,ci,stats]=ttest2(earlyadapt(grp==1),earlyadapt(grp==2))
[h,p,ci,stats]=ttest2(earlyadapt(grp==2),earlyadapt(grp==3))
[h,p,ci,stats]=ttest2(earlyadapt(grp==3),earlyadapt(grp==1))

dBigVsBisectEA = cohensD(earlyadapt(grp==2),earlyadapt(grp==3))
dBisectVsSmallEA = cohensD(earlyadapt(grp==3),earlyadapt(grp==1))
dEarlyAdapt = cohensD(earlyadapt(grp==1),earlyadapt(grp==2))

mEAsmall = mean(earlyadapt(grp==1))
sdEAsmall = std(earlyadapt(grp==1))

mEAbig = mean(earlyadapt(grp==2))
sdEAbig = std(earlyadapt(grp==2))


[h,p,ci,stats]=ttest2(asymptote1(grp==1),asymptote1(grp==2))
[h,p,ci,stats]=ttest2(asymptote1(grp==2),asymptote1(grp==3))
[h,p,ci,stats]=ttest2(asymptote1(grp==3),asymptote1(grp==1))

dBigVsBisectAsy = cohensD(asymptote1(grp==2),asymptote1(grp==3))
dSmallVsBisectAsy = cohensD(asymptote1(grp==1),asymptote1(grp==3))
dAsymptote1 = cohensD(asymptote1(grp==1),asymptote1(grp==2))

mAsy1Small = mean(asymptote1(grp==1))
sdAsy1Small = std(asymptote1(grp==1))

mAsy1Big = mean(asymptote1(grp==2))
sdAsy1Big = std(asymptote1(grp==2))


%effect size according to standardized difference scores (see:
%http://jakewestfall.org/blog/index.php/2016/03/25/five-different-cohens-d-statistics-for-within-subject-designs/)
% d = t/sqrt(n)
dzSmallToBig = stats_ws(1).tstat/sqrt(sum(grp==1))
dzBigToSmall = stats_ws(2).tstat/sqrt(sum(grp==2))
dzBisectToBig = stats_ws(3).tstat/sqrt(sum(grp==3))


%% Compare parameter estimates of individual data -- permutation test

load e3_ind_bestParams

test_stat_a = abs(mean(best_params(grp==1,1) - best_params(grp==2,1)));
test_stat_b = abs(mean(best_params(grp==1,2) - best_params(grp==2,2)));

%create null distributions
num_sims = 1e4
for i=1:num_sims
    
    shuffled_idx = grp(randperm(size(grp,1)));
    
    null_dist_a(i,:) = mean(best_params(shuffled_idx==1,1) - best_params(shuffled_idx==2,1));
    null_dist_b(i,:) = mean(best_params(shuffled_idx==1,2) - best_params(shuffled_idx==2,2));

end

figure; hold on
histogram(null_dist_a)
plot([test_stat_a test_stat_a],[ylim],'k')
pval_a = sum(null_dist_a > test_stat_a | null_dist_a < -test_stat_a)/num_sims

figure; hold on
histogram(null_dist_b)
plot([test_stat_b test_stat_b],[ylim],'k')
pval_b = sum(null_dist_b > test_stat_b | null_dist_b < -test_stat_b)/num_sims

%% Quantify quality of model fitting

load e3_boot_bestParams_test

%without concern for group
sorted_r2 = sort(r2);

lb = floor(.025*size(best_params,1));
ub = floor(.975*size(best_params,1));

ci = [sorted_r2(lb) sorted_r2(ub)]
median(sorted_r2)

%bootstrap test
hit_params = best_params(1:1000,:); 
miss_params = best_params(1001:2000,:);

avg_diff_a = mean(hit_params(:,1)) - mean(miss_params(:,1));
avg_diff_u = mean(hit_params(:,2)) - mean(miss_params(:,2));

null_a = hit_params(:,1) - miss_params(:,1) - avg_diff_a;
null_u = hit_params(:,2) - miss_params(:,2) - avg_diff_u;

figure; hold on
histogram(null_a)

figure; hold on
histogram(null_u)

pa = sum(null_a > -avg_diff_a | null_a < avg_diff_a)/1e3
pu = sum(null_u > avg_diff_u | null_u < -avg_diff_u)/1e3

%compare r2 values when fixing A or fixing B
load e3_boot_bestParams_fixedA

%without concern for group
r2_fixedA = r2;
sorted_r2_fixedA = sort(r2_fixedA);
ci_fixedA = [sorted_r2_fixedA(ceil(.025*length(r2_fixedA))) sorted_r2_fixedA(ceil(.975*length(r2_fixedA)))] 

load e3_boot_bestParams_fixedU
r2_fixedB = r2;
sorted_r2_fixedB = sort(r2_fixedB);
ci_fixedB = [sorted_r2_fixedB(ceil(.025*length(r2_fixedB))) sorted_r2_fixedB(ceil(.975*length(r2_fixedB)))] 



