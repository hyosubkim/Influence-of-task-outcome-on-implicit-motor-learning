
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
num_params = 6;
lb = [0 0 0 0 0 0];
ub = [1 1 clamp_angle clamp_angle 1 1];

% Number of randomly initialized starting points for the minimizer
num_inits = 10;

for iteration = 1:num_inits
    
    iteration

    rot = rotations;
    group = [1; 2];
    hand_data = clampData;
    %Initial values
    initials = rand(1, num_params) .* ub;
    
    %Optimizer
    [params, err] = fmincon(@squared_err, initials,[],[], [], [], lb, ub, [], options, hand_data, rot(1), group, rwd);
    
    
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
bestparams_am_mr_hybrid = temp_bestparams(bestidx,:);
r2_am_mr_hybrid = temp_r2(bestidx);
aic_am_mr_hybrid = temp_aic(bestidx);

%save e3_am_mr_hybrid_params_groupfit bestparams_am_mr_hybrid r2_am_mr_hybrid  aic_am_mr_hybrid

for i=1:2
    
    sims(i,:) = am_mr_hybrid_simulator(bestparams_am_mr_hybrid,groupavg_hand_data(i,11:210),clamp_angle,group(i),rwd(i,:));
    
    figure; hold on
    plot(sims(i,:),'k','linewidth',3)
    plot(groupavg_hand_data(i,11:210),'color',grpclrs(i,:),'linewidth',3)
    
end
figure; hold on
plot(sims(1,:),'k','linewidth',3)
plot(sims(2,:),'k','linewidth',3)
shadedErrorBar(1:200,mean(handAngle_s2b(:,11:210)),sem(handAngle_s2b(:,11:210)),{'color',grpclrs(1,:)},1)
shadedErrorBar(1:200,mean(handAngle_b2s(:,11:210)),sem(handAngle_b2s(:,11:210)),{'color',grpclrs(2,:)},1)
xlabel('Movement cycle')
ylabel('Hand angle (deg)')

end


function [sq_err, simulations] = squared_err(params, hand_data, rot, group,rwd)

% Simulate behavior from initial parameters
simulations = am_mr_hybrid_simulator(params, hand_data, rot, group,rwd);

% Calculate Squared Error - between data and simulation
sq_err = nansum(nansum ( (hand_data - simulations).^2) ) ;

end


function simulations = am_mr_hybrid_simulator(params, hand_data, rot, group,rwd)

Amiss = params(1);
Ahit = params(2);
Umiss = params(3);
Uhit = params(4);
Ar = params(5);
s = params(6);

num_trials = size(hand_data, 2);

%initializing
nAngles = 36000;
angles = linspace(0,359.99,nAngles)'; %we want representations at every hundredth of a degree
resolution = 360/nAngles;
deg2idx = 1/resolution; %bc 0 deg is 1st element
unitvec = [cosd(angles), sind(angles)];
r = zeros(nAngles,num_trials); %initialize weights for every direction
k = zeros(nAngles,1); %this is selection vector determining which angle has its weight increased

PA(1)=0; %the population vector (P) will have a direction and magnitude
PL=0; %the length of P will also determine relative weights on adaptation and reward learning
wLmax = 1;

for n = 1:size(hand_data,1)
    
    grp=group(n);
    reward=rwd(n,:);
    
    %Load Initial Values
    [v(1) y(1) PA(1) PL(1) r(:,num_trials)] = deal(0);
    
    for i = 1:num_trials-1
        
        rwdAng = round(y(i),2); %hand angle when you were hitting target; round to 2nd decimal
        rwdAng = rwdAng*deg2idx; %to get accurate index
        if rwdAng<0
            kidx=round(rwdAng+nAngles+1); %bc 0 deg is 1st element
        else
            kidx=round(rwdAng+1);
        end
        
        %increase weight on directionally-tuned "neuron" representing rewarded angle        if reward(i)==1
        if reward(i)==1
            k(kidx)=k(kidx)+s;       
            v(i+1) = Ahit*v(i) + Uhit;
        elseif reward(i)==0
            k=k;
            v(i+1) = Amiss*v(i) + Umiss;
        end
        
        r(:,i+1) = Ar*r(:,i) + k;
        
        if sum(r(:,i)) > 1e-4
            Px(i,1) = dot(r(:,i),unitvec(:,1)); %x-component of Popn vector
            Py(i,1) = dot(r(:,i),unitvec(:,2)); %y-component of Popn vector
        else
            Px(i,1) = 0;
            Py(i,1) = 0;
        end
        PA(i+1) = atan2d(Py(i),Px(i)); %Popn vector ANGLE
        PL(i+1) = sqrt(Px(i).^2 + Py(i).^2); %Popn vector LENGTH
        %PV length max value must be 1
        if PL(i+1) > wLmax
            PL(i+1) = wLmax; 
        end
        
        %actual movement direction
        y(i+1) = (1-PL(i+1))*(v(i+1)) + PL(i+1)*PA(i+1);
                
        k = zeros(36000,1);
        
    end
    
    simulations(n,:) = y;
    
end

end



