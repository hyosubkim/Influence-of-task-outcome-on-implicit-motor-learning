%%% Note: This is the Dual Error model in the manuscript. Originally, it
%%% was referred to as the Implicit Aim ("i.e., IA model"); however, we
%%% have since went with the name "Dual Error" as we are agnostic as to
%%% whether the second model-based process is a form of aiming or not.


function tgtsize_SSM_fitting

clear all; close all; clc

% Load hand data
% Block cycles:
% 1:5 baseline nofb, 
% 6:10 baseline fb, 
% 11:130 end of first clamp, 
% 131:210 end of second clamp
% 211:220 washout

load tgtSwitch_data

groupavg_hand_data = [nanmean(handAngle_s2b); nanmean(handAngle_b2s)];
clampData = groupavg_hand_data(:,11:210);

clamp_angle = 1.75;
rotations = [ones(size(groupavg_hand_data,1),1)*clamp_angle]; % Clamp sizes for each group

grp = [1; 2];
grpclrs = parula(3);
grpclrs = grpclrs(1:2,:);

rwd_sth = [zeros(1,10) zeros(1,120) ones(1,80) zeros(1,10)];
rwd_hts = [zeros(1,10) ones(1,120) zeros(1,80) zeros(1,10)];
rwd = [rwd_sth(11:210); rwd_hts(11:210)];

% Options to enter the minimizer
options=optimset('MaxFunEvals',1e16,'TolFun',1e-16,'TolX',1e-16);

%%%%%% initialize params 
lb = [0 0 0 0];
ub = [1 1 clamp_angle clamp_angle];
num_params = length(lb);

% Number of randomly initialized starting points for the minimizer
num_inits = 10;


for iteration = 1:num_inits
    
    iteration
    
    rot = rotations;
    group = [1; 2];
    hand_data = clampData;  
    
    %Initial values
    initials = rand(1, num_params) .* ub;
    
    [params, err, simulations] = fmincon(@squared_err, initials,[],[], [], [], lb, ub, [], options, hand_data, rot(1), group, rwd);
    
    best_err = Inf;
    % Save the best fitting parameters
    if err < best_err
        best_err = err;
        bestparams = params;
        sq_err = err;
        r2 = 1 - (sq_err)/sum(sum((hand_data - mean(hand_data)).^2));
        
        %%%add AIC%%%
        aic = 2*num_params + length(hand_data(:))*log(sq_err/length(hand_data(:)));
    end
    
    temp_best_err(iteration,1) = best_err;
    temp_bestparams(iteration,:) = bestparams;
    temp_r2(iteration,1) = r2;
    temp_aic(iteration,1) = aic; 
    
end

temp_best_err
temp_bestparams
[~,bestidx] = min(temp_best_err);
bestparams_IA = temp_bestparams(bestidx,:);
r2_IA = temp_r2(bestidx);
aic_IA = temp_aic(bestidx);

save e3_IA_params_groupfit bestparams_IA groupavg_hand_data r2_IA aic_IA


for i=1:2
    
    sims(i,:) = implicitAim(bestparams_IA,groupavg_hand_data(i,11:210),clamp_angle,group(i), rwd(i,:));
    
    figure; hold on
    plot(sims(i,:),'k','linewidth',3)
    plot(groupavg_hand_data(i,11:210),'color',grpclrs(i,:),'linewidth',3)
    xlabel('Movement cycles')
    ylabel('Hand angle (deg)')
    
end
figure; hold on
plot(sims(1,:),'k','linewidth',3)
plot(sims(2,:),'k','linewidth',3)
shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'color',grpclrs(1,:)},1)
shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'color',grpclrs(2,:)},1)
xlabel('Movement cycle')
ylabel('Hand angle (deg)')

end


%Cost fxn
function [sq_err, simulations] = squared_err(params, hand_data, rot, group, rwd)

% Simulate behavior from initial parameters
simulations = implicitAim(params, hand_data, rot, group, rwd);

% Calculate Squared Error - between data and simulation
sq_err = nansum(nansum((hand_data - simulations).^2)) ;

end


%Simulator
function simulations = implicitAim(params, hand_data, rot, group, rwd)

Aspe = params(1);
Atpe = params(2);
Uspe = params(3);
Utpe = params(4);

num_trials = size(hand_data, 2);

% Load Initial Values
error = rot;

x(1) = 0; 
x_spe(1) = 0;
x_tpe(1) = 0;

for n = 1:size(hand_data,1)
 
    grp=group(n);
    reward=rwd(n,:);
    
    for i = 1:num_trials-1
        if reward(i)==1
            x_tpe(i+1) = Atpe*x_tpe(i);
        elseif reward(i)==0
            x_tpe(i+1) = Atpe*x_tpe(i) + Utpe;
        end
        
        x_spe(i+1) = Aspe*x_spe(i) + Uspe;
        
        x(i+1) = x_tpe(i+1) + x_spe(i+1);
        
    end
    
    simulations(n,:) = x;
    
end

end



