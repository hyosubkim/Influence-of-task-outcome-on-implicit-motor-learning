%%%Analysis of E2 behavioral data from intrinsic rewards project%%%

clear all; close all; clc

load('e2_data.mat');
% load('extTgtExp_movefile_table.mat');

% default figure properties
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'defaultAxesFontSize',12)

%anonymous fxn for pooled SD
pooledsd = @(n,m,x,y)(sqrt(((n-1)*std(x)^2 + (m-1)*std(y)^2)/(n+m-2)));

subjects = unique(T.SN);
targets = unique(T.ti);
bad_trial = zeros(size(subjects,1),size(targets,1));
bsl_cycle = 19:20;
fbBslCycles = 12:20;
ntrials = max(T.TN);
blks_trials = [40 80 360 640 920 960 1000];
blks_cycles = blks_trials/4;
clampcycles_onetoten = 21:30;
earlyclamp1 = 23:27;
earlyclamp2 = 22:31;
lateclamp = 221:230;
zeroclamp = 241:250;

%for movement velocity analyses
bsltrials = 41:80; %cycles 11-20
earlyphase1 = 81:160; %cycles 21-40
latephase1 = 881:960; %cycles 121-130
zeroclamptrials = 961:1000; %cycles 201-210

hold_t = .500;
endptfb_t = .05;

%
tgtclrs=parula(4);
grpclrs=[0 0 0; 1 0 0];
gray=[0.8 0.8 0.8];

% uncomment next line if analyzing with hand angle at max vel
T.hand_theta = T.hand_theta_maxradv;
T.raw_hand_theta = T.hand_theta;
T.bslc_hand_theta = T.hand_theta;

T.radvelmax = T.radvelmax./10;
T.maxRadDist = T.maxRadDist./10;

T.TT = hold_t + endptfb_t + T.ST + T.RT + T.MT;

% loop through all subjects
tic
for sn = 1:length(subjects)
    
    %indices for each subject
    subjidx = T.SN==subjects(sn);
    idx = find(subjidx,1);
    grp(sn,1) = T.group(idx);
    ccw = T.cond(idx) > 0;
    subject_code(sn,1) = sn;
        
    % next loop is for by-target analysis
    for ti = 1:length(targets)
        
        tgt_idx = find(subjidx & T.ti==targets(ti));
        
        %next lines are only used for identifying outliers
        smoothed_hand_theta = smooth(T.hand_theta(tgt_idx),5);
        smoothed = T.hand_theta(tgt_idx) - smoothed_hand_theta;
        tgt_SDs(sn,ti) = std(smoothed);
        
        % outlier removal
        for i=1:length(smoothed)
            if abs(smoothed(i)) > 3*tgt_SDs(sn,ti) | abs(T.hand_theta(tgt_idx(i))) > 90
                [T.hand_theta(tgt_idx(i)),T.raw_hand_theta(tgt_idx(i)),T.hand_theta_maxradv(tgt_idx(i)),T.hand_theta_50(tgt_idx(i)),...
                    T.raw_ep_hand_ang(tgt_idx(i)), T.MT(tgt_idx(i)), T.RT(tgt_idx(i)),T.ST(tgt_idx(i)),...
                    T.radvelmax(tgt_idx(i)),T.maxRadDist(tgt_idx(i))] = deal(nan);
                bad_trial(sn,ti) = bad_trial(sn,ti)+1;
            end
%             if abs(T.hand_theta(tgt_idx(i))) > 90
%                 T.hand_theta(tgt_idx(i)) = nan;
%                 bad_trial(sn,ti) = bad_trial(sn,ti)+1;
%             end
        end
        nofb_baseline_tgtErr(sn,ti) = nanmedian(T.raw_hand_theta(tgt_idx(1:5)));
        nofb_baseline_sd(sn,ti) = nanstd(T.hand_theta(tgt_idx(1:5)));
        fb_baseline_tgtErr(sn,ti) = nanmedian(T.raw_hand_theta(tgt_idx(fbBslCycles)));
        fb_baseline_sd(sn,ti) = nanstd(T.hand_theta(tgt_idx(fbBslCycles)));
        
        st_tgt_med(sn,ti) = nanmedian(T.ST(tgt_idx));
        rt_tgt_med(sn,ti) = nanmedian(T.RT(tgt_idx));
        mt_tgt_med(sn,ti) = nanmedian(T.MT(tgt_idx));
        tt_tgt_med(sn,:) = nanmedian(T.TT(tgt_idx));
        
        radvelmax_tgt(sn,ti) = nanmedian(T.radvelmax(tgt_idx));
        maxRadDist_tgt(sn,ti) = nanmedian(T.maxRadDist(tgt_idx));
        
        baseline_bias(sn,ti) = nanmean(T.hand_theta(tgt_idx(fbBslCycles)));
        T.bslc_hand_theta(tgt_idx) = T.hand_theta(tgt_idx) - baseline_bias(sn,ti);

        % detrend hand data
        validearlyidx = find(~isnan(T.hand_theta(tgt_idx(clampcycles_onetoten))));
        validlateidx = find(~isnan(T.hand_theta(tgt_idx(lateclamp))));
        detrend_early(sn,ti,1:length(validearlyidx)) = detrend(T.hand_theta(tgt_idx(clampcycles_onetoten(validearlyidx))));
        detrend_late(sn,ti,1:length(validlateidx)) = detrend(T.hand_theta(tgt_idx(lateclamp(validlateidx))));
        
        % flip ccw hand angle
        if ccw
            T.bslc_hand_theta(tgt_idx) = T.bslc_hand_theta(tgt_idx)*-1;
        end
        
        if ti==1
            handTgt1(sn,:) = T.bslc_hand_theta(tgt_idx);
        elseif ti==2
            handTgt2(sn,:) = T.bslc_hand_theta(tgt_idx);
        elseif ti==3
            handTgt3(sn,:) = T.bslc_hand_theta(tgt_idx);
        elseif ti==4
            handTgt4(sn,:) = T.bslc_hand_theta(tgt_idx);
        end
                
    end
    
    st(sn,:) = T.ST(subjidx);
    rt(sn,:) = T.RT(subjidx);
    mt(sn,:) = T.MT(subjidx);
    tt(sn,:) = T.TT(subjidx);
    radvelmax(sn,:) = T.radvelmax(T.SN==sn);
    maxRadDist(sn,:) = T.maxRadDist(T.SN==sn);
    
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
nofb_bsl_var_ind = mean(nofb_baseline_sd,2);
bsl_var_ind = mean(fb_baseline_sd,2);

%point estimates of kinematics during different phases
earlyadapt1 = nanmean(hand_ang_m(:,earlyclamp1),2)/length(earlyclamp1);
earlyadapt2 = nanmean(hand_ang_m(:,earlyclamp2),2)/length(earlyclamp2);
asymptote = nanmean(hand_ang_m(:,lateclamp),2);
abs_decay = hand_ang_m(:,250)-hand_ang_m(:,240);
retention = hand_ang_m(:,250)./hand_ang_m(:,240); %retention as percentage of last clamp cycle

nofb_bslTgtErr = nanmean(nofb_baseline_tgtErr,2);
bslTgtErr = nanmean(fb_baseline_tgtErr,2);

%grouping variables
[tgtgrp, tgtid] = findgroups(T.ti);
[~,grpid] = findgroups(T.group);
grpid = [grpid(2) grpid(1)];
grplabel = {'Straddle','Hit'};

%temporal variables
mt_ind_mean = mean(mt_tgt_med,2);
rt_ind_mean = mean(rt_tgt_med,2);
st_ind_mean = mean(st_tgt_med,2);
tt_ind_mean = mean(tt_tgt_med,2);

radvelmax_ind = mean(radvelmax_tgt,2);
maxRadDist_ind = mean(maxRadDist_tgt,2);

%MT, RT
rt_bsl = nanmean(rt_cycle(:,fbBslCycles),2);
earlyclampRT = nanmean(rt_cycle(:,earlyclamp1),2);
lateclampRT = nanmean(rt_cycle(:,lateclamp),2);
zeroclampRT = nanmean(rt_cycle(:,zeroclamp),2);

mt_bsl = nanmean(mt_cycle(:,fbBslCycles),2);
earlyclampMT = nanmean(mt(:, earlyclamp1),2);
lateclampMT = nanmean(mt(:, lateclamp),2);
zeroclampMT = nanmean(mt_cycle(:,zeroclamp),2);


toc

%% Normality tests

for i=1:length(grpid)
% [h p] = swtest(earlyadapt1(grp==grpid(i)))
% [h p] = swtest(asymptote(grp==grpid(i)))
[h p] = swtest(retention(grp==grpid(i)))
%[h p] = swtest(abs_decay(grp==grpid(i)))
end

p = vartestn(earlyadapt1,grp,'TestType','LeveneAbsolute')
p = vartestn(asymptote,grp,'TestType','LeveneAbsolute')
p = vartestn(retention,grp,'TestType','LeveneAbsolute')
% p = vartestn(abs_decay,grp,'TestType','LeveneAbsolute')
%% Baseline target error and variability

% load e1_baselineRTs
% grp_e1e2 = [grp_e1(grp==grpid(1) | grp==grpid(2)); grp];
% rt_ind_mean_e1e2 = [rt_e1(grp==grpid(1) | grp==grpid(2)); rt_ind_mean];
% permutation_test(rt_ind_mean_e1e2,grp_e1e2,mean(rt_ind)


for i=1:2
    
    figure(323); hold on
    set(gcf,'units','inches','pos',[5 5 2 2]);
    set(gcf,'PaperPositionMode','auto')
    bar(i,mean(bsl_var_ind(grp==grpid(i))),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bsl_var_ind(grp==grpid(i))),nanstderr(bsl_var_ind(grp==grpid(i))),'k','linewidth',2)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,bsl_var_ind(grp==grpid(i)),'.','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:),'markersize',12)
    set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['Hand angle (',char(176),')'])
    title('Baseline variability')
    
    figure(327); hold on
    set(gcf,'units','inches','pos',[5 5 2 2]);
    set(gcf,'PaperPositionMode','auto')
    bar(i,mean(bslTgtErr(grp==grpid(i))),'facecolor',gray,'edgecolor',gray)
    errbar(i-.1,mean(bslTgtErr(grp==grpid(i))),nanstderr(bslTgtErr(grp==grpid(i))),'k','linewidth',2)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,bslTgtErr(grp==grpid(i)),'.','markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:),'markersize',12)
    set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
    ylabel(['Hand angle (',char(176),')'])
    title('Baseline error')
    
end

% print(323,'baselineSD','-painters','-dpdf')
% print(327,'baselineErr','-painters','-dpdf')

%t-tests
[h,p,ci,stats]=ttest2(bsl_var_ind(grp==grpid(1)),bsl_var_ind(grp==grpid(2)))
[h,p,ci,stats]=ttest2(bslTgtErr(grp==grpid(1)),bslTgtErr(grp==grpid(2)))

%report summary statistics
mean(bsl_var_ind(grp==grpid(1)))
nanstderr(bsl_var_ind(grp==grpid(1)))

mean(bsl_var_ind(grp==grpid(2)))
nanstderr(bsl_var_ind(grp==grpid(2)))

mean(bslTgtErr(grp==grpid(1)))
nanstderr(bslTgtErr(grp==grpid(1)))

mean(bslTgtErr(grp==grpid(2)))
nanstderr(bslTgtErr(grp==grpid(2)))

%effect size calculation
dVar = cohensD(bsl_var_ind(grp=='small tgt'),bsl_var_ind(grp=='big tgt'))


%% Temporal variables analyses
for i=1:2
    
    %plot MTs and RTs for entire experiment
    figure(919); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'units','inches','pos',[5 5 4 2.5]);
    set(gcf,'PaperPositionMode','auto')
    rectangle('position',[0 0 10 1],'facecolor',[0.8 0.8 0.8 0.5],'edgecolor',gray)
    rtGrp(i)=shadedErrorBar(1:length(rt_cycle),nanmean(rt_cycle(grp==grpid(i),:)),...
        nanstderr(rt_cycle(grp==grpid(i),:)),{'.','markersize',8,'color',grpclrs(i,:)},1);
    xlim([-2 272])
    ylim([0.2 1])
    plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
        [240 240], ylim, '--k',[250 250], ylim, '--k')
    title('Reaction Times')
    xlabel('Movement cycle (4 reaches)')
    ylabel('RT (s)')
    
    figure(324); hold on
    zeroline = zeros(size(hand_ang_m,1));
    set(gcf,'units','inches','pos',[5 5 4 2.5]);
    set(gcf,'PaperPositionMode','auto')
    rectangle('position',[0 0 10 0.5],'facecolor',[0.8 0.8 0.8 0.5],'edgecolor',gray)
%     rect2=rectangle('position',[240 0 10 0.5],'facecolor',gray,'edgecolor',gray)
%     alpha(rect2,0)
    mtGrp(i)=shadedErrorBar(1:length(mt_cycle),nanmean(mt_cycle(grp==grpid(i),:)),...
        nanstderr(mt_cycle(grp==grpid(i),:)),{'.','markersize',8,'color',grpclrs(i,:)},1);
    xlim([-2 272])
    ylim([0.1 0.5])
    plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
        [240 240], ylim, '--k',[250 250], ylim, '--k')
    title('Movement times')
    xlabel('Movement cycle (4 reaches)')
    ylabel('Time (s)')
end
% legend([grpdata(1).mainLine,grpdata(2).mainLine],{grplabel{1},grplabel{2}})
% print(919,'RTs','-painters','-dpdf')
% print(324,'MTs','-painters','-dpdf')


%compare average of median MTs and RTs
[h p ci stats] = ttest2(rt_ind_mean(grp=='big tgt'),rt_ind_mean(grp=='small tgt'))
[h p ci stats] = ttest2(mt_ind_mean(grp=='big tgt'),mt_ind_mean(grp=='small tgt'))

%compare baseline MTs and RTs
[h p ci stats] = ttest2(rt_bsl(grp=='big tgt'),rt_bsl(grp=='small tgt'))
[h p ci stats] = ttest2(mt_bsl(grp=='big tgt'),mt_bsl(grp=='small tgt'))

%compare zero clamp MTs and RTs
[h p ci stats] = ttest2(zeroclampRT(grp=='big tgt'),zeroclampRT(grp=='small tgt'))
[h p ci stats] = ttest2(zeroclampMT(grp=='big tgt'),zeroclampMT(grp=='small tgt'))


%summary statistics
mean(rt_bsl(grp=='big tgt'))
nanstderr(rt_bsl(grp=='big tgt'))

mean(rt_bsl(grp=='small tgt'))
nanstderr(rt_bsl(grp=='small tgt'))

mean(mt_bsl(grp=='big tgt'))
nanstderr(mt_bsl(grp=='big tgt'))

mean(mt_bsl(grp=='small tgt'))
nanstderr(mt_bsl(grp=='small tgt'))

% effect size calculations
dRT = cohensD(rt_bsl(grp=='big tgt'),rt_bsl(grp=='small tgt'))
dMT = cohensD(mt_bsl(grp=='big tgt'),mt_bsl(grp=='small tgt'))


%%% Mixed-model ANOVA
% X: design matrix with four columns (future versions may allow different input configurations)
%     - first column  (i.e., X(:,1)) : all dependent variable values
%     - second column (i.e., X(:,2)) : between-subjects factor (e.g., subject group) level codes (ranging from 1:L where 
%         L is the # of levels for the between-subjects factor)
%     - third column  (i.e., X(:,3)) : within-subjects factor (e.g., condition/task) level codes (ranging from 1:L where 
%         L is the # of levels for the within-subjects factor)
%     - fourth column (i.e., X(:,4)) : subject codes (ranging from 1:N where N is the total number of subjects)
%X(:,1) = [rt_bsl; earlyclampRT; lateclampRT];
X(:,1) = [mt_bsl; earlyclampMT; lateclampMT];
X(:,2) = repmat(grp,[3,1]);
X(:,3) = [ones(length(rt_bsl),1); ones(length(earlyclampRT),1)*2; ones(length(lateclampRT),1)*3];
X(:,4) = repmat(subject_code,[3,1]);
[SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(X,0);
% 
% % see https://stardock.cs.virginia.edu/empirical/resources/Brown28.pdf for
% % effect sizes on mixed model ANOVA
% SSQs=cell2mat(SSQs);
% TSE = sum(SSQs);
% etasq_btwn = SSQs(1)/TSE; 
% etasq_within = SSQs(3)/TSE;
partial_etasq_int = SSQs{4}/(SSQs{4}+SSQs{5})


%Report mean +/- SEM for RTs and MTs
for ii=1:length(grpid)
    
%     %RTs
%     mean(rt_bsl(grp==grpid(ii)))
%     nanstderr(rt_bsl(grp==grpid(ii)))
%     
%     mean(earlyclampRT(grp==grpid(ii)))
%     nanstderr(earlyclampRT(grp==grpid(ii)))
%     
    mean(lateclampRT(grp==grpid(ii)))
    nanstderr(lateclampRT(grp==grpid(ii)))
% 
%     mean(zeroclampRT(grp==grpid(ii)))
%     nanstderr(zeroclampRT(grp==grpid(ii)))
    
%     %MTs
%     mean(mt_bsl(grp==grpid(ii)))
%     nanstderr(mt_bsl(grp==grpid(ii)))
%     
%     mean(earlyclampMT(grp==grpid(ii)))
%     nanstderr(earlyclampMT(grp==grpid(ii)))
%     
%     mean(lateclampMT(grp==grpid(ii)))
%     nanstderr(lateclampMT(grp==grpid(ii)))
%     
%     mean(zeroclampMT(grp==grpid(ii)))
%     nanstderr(zeroclampMT(grp==grpid(ii)))
end


%% Plot individual data

for i=1:length(subjects)
    
    %trial view
%     trialidx=find(T.SN==i);
    
%     figure; hold on
%     ind = gscatter(1:length(trialidx),T.raw_hand_theta(trialidx),tgtgrp(trialidx),...
%         tgtclrs,'.',25)
%     legend(ind,{'45','135','225','315'})
    
    %cycle view
    figure
    plot(hand_ang_m(i,:),'o','markerfacecolor',grpclrs(grp(i),:))
    title(['S',num2str(subjects(i)),' group: ',char(grp(i))])
    xlabel('Movement cycle')
    ylabel('Hand angle (deg)')

end


%% Plot group analyzed hand data
zeroline = zeros(size(hand_ang_m,1));
figure
hold on
set(gcf,'units','inches','pos',[5 5 4.75 2.63]);
set(gcf,'PaperPositionMode','auto')
rectangle('position',[0 -10 10 40],'facecolor',gray,'edgecolor',gray)
rectangle('position',[240 -5 10 40],'facecolor',gray,'edgecolor',gray)
for i=1:2
   
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==grpid(i),:)),...
        nanstderr(hand_ang_m(grp==grpid(i),:)),{'.','markersize',8,'color',grpclrs(i,:)},1);
    
end
xlim([-2 272])
ylim([-5 30])
plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
    [240 240], ylim, '--k',[250 250], ylim, '--k')
title(['Extended 1.75',char(176),' clamp'])
xlabel('Movement cycle (4 reaches)')
ylabel(['Hand angle (',char(176),')'])
% legend([grpdata(1).mainLine,grpdata(2).mainLine],{grplabel{1},grplabel{2}})
% print('smallTgtOnly','-painters','-dpdf')
print('e2_extTgt','-painters','-dpdf')


%% Initial learning rate and asymptote

% Early adaptation
for i=1:length(grpid)
    
    % early rate: calculated using clamp cycles 3-7
    ea1 = figure(111); hold on
    set(gcf,'units','inches','pos',[5 5 1.5 2.63]);
    set(gcf,'PaperPositionMode','auto');
    bar(i,mean(earlyadapt1(grp==grpid(i))),'facecolor',gray)
    errbar(i-.1,mean(earlyadapt1(grp==grpid(i))),nanstderr(earlyadapt1(grp==grpid(i))),'k')
    plot(.2*rand(sum(grp==grpid(i)),1)+i,earlyadapt1(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2 ],'xticklabel',grplabel,'xticklabelrotation',45)
    %     set(gca,'xtick',[1 2],'xticklabel',grplabel)
    xlim([.2 2.8])
    ylim([-.5 1.5])
    ylabel('\DeltaHand angle/trial (deg)')
    title('Early rate','fontsize',12)
    
%     % early rate: calculated using clamp cycles 2-11
%     ea2 = figure(115); hold on
%     set(gcf,'units','inches','pos',[5 5 1.5 2.63]);
%     set(gcf,'PaperPositionMode','auto');
%     bar(i,mean(earlyadapt2(grp==grpid(i))),'facecolor',gray)
%     errbar(i-.1,mean(earlyadapt2(grp==grpid(i))),nanstderr(earlyadapt2(grp==grpid(i))),'k','linewidth',1.5)
%     plot(.2*rand(sum(grp==grpid(i)),1)+i,earlyadapt2(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
%         'markeredgecolor',grpclrs(i,:))
%     set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
%     %     set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
%     xlim([.2 2.8])
%     ylim([-.1 .8])
%     ylabel(['\DeltaHand angle/trial (',char(176),')'])
%     title('Early rate','fontsize',12)
    
    % Asymptote
    asy = figure(121); hold on
    set(gcf,'units','inches','pos',[5 5 1.5 2.63]);
    set(gcf,'PaperPositionMode','auto');
    bar(i,mean(asymptote(grp==grpid(i))),'facecolor',gray)
    errbar(i-.1,mean(asymptote(grp==grpid(i))),nanstderr(asymptote(grp==grpid(i))),'k','linewidth',1.5)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,asymptote(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
    xlim([.2 2.8])
    ylabel(['Hand angle (',char(176),')'])
    title('Asymptote','fontsize',12)
    
    % Retention
    retfig = figure(122); hold on
    set(gcf,'units','inches','pos',[5 5 1.5 2.63]);
    set(gcf,'PaperPositionMode','auto');
    bar(i,mean(retention(grp==grpid(i))),'facecolor',gray)
    errbar(i-.1,mean(retention(grp==grpid(i))),nanstderr(retention(grp==grpid(i))),'k','linewidth',1.5)
    plot(.2*rand(sum(grp==grpid(i)),1)+i,retention(grp==grpid(i)),'.','markersize',12,'markerfacecolor',grpclrs(i,:),...
        'markeredgecolor',grpclrs(i,:))
    set(gca,'xtick',[1 2],'xticklabel',grplabel,'xticklabelrotation',45)
    xlim([.2 2.8])
    ylabel(['Retention (%)'])
    title('Retention','fontsize',12)
 
end

% print('e2_ret',retfig,'-dpdf','-painters')
% t-tests for rate and asymptote
[~,p_rate1,ci_rate1,stats_rate1] = ttest2(earlyadapt1(grp=='big tgt'),earlyadapt1(grp=='small tgt'))
% [~,p_rate2,ci_rate2,stats_rate2] = ttest2(earlyadapt2(grp=='big tgt'),earlyadapt2(grp=='small tgt'))
mEAbig = mean(earlyadapt2(grp=='big tgt'));
sdEAbig = std(earlyadapt2(grp=='big tgt'));
mEAsmall = mean(earlyadapt2(grp=='small tgt')); 
sdEAbig = std(earlyadapt2(grp=='small tgt')); 

[~,p_asy,ci_asym,stats_asy] = ttest2(asymptote(grp=='big tgt'),asymptote(grp=='small tgt'))
mean(asymptote(grp=='big tgt'))
std(asymptote(grp=='big tgt'))
mean(asymptote(grp=='small tgt'))
std(asymptote(grp=='small tgt'))

permutation_test(asymptote,grp,mean(asymptote(grp=='big tgt'))-mean(asymptote(grp=='small tgt')),1e5)

%effect size calculations
dRate1 = cohensD(earlyadapt1(grp=='small tgt'),earlyadapt1(grp=='big tgt'))
% dRate2 = cohensD(earlyadapt2(grp=='small tgt'),earlyadapt2(grp=='big tgt'))
dAsymptote = cohensD(asymptote(grp=='small tgt'),asymptote(grp=='big tgt'))

%%% uncomment if printing figures
print('ea1',ea1,'-dpdf','-painters')
% cd('../../Figures')
% print('e2_ea',ea2,'-dpdf','-painters')
% print('e2_asy',asy,'-dpdf','-painters')

%% Compare retention from end of 1.75 deg clamp to end of 0 deg clamp

for i=1:length(grpid)
    [h p] = swtest(retention(grp==grpid(i)))
end

p_levene = vartestn(retention,grp,'TestType','LeveneAbsolute')

[~,p_retention,ci_retention,stats_retention] = ttest2(retention(grp=='big tgt'),retention(grp=='small tgt'))

%permutation test
diffScore = mean(retention(grp=='big tgt') - retention(grp=='small tgt'));
permutation_test(retention,grp,diffScore,1e5)

%compare absolute decay in hand angle
[~,p_retention,ci_retention,stats_retention] = ttest2(abs_decay(grp=='small tgt'),abs_decay(grp=='big tgt'))

%effect sizes
dRetention = cohensD(retention(grp=='small tgt'),retention(grp=='big tgt'))
cohensD(abs_decay(grp=='small tgt'),abs_decay(grp=='big tgt'))
