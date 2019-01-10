function simulations = adaptationModulation_4params(params, hand_data, rot, group, rwd)

%Note: Umiss and Uhit could also be interpreted as gmiss*U or ghit*U (i.e.,
%as having separate gain factors on update (same for Amiss and Ahit); this 
%is how it's described in paper)

num_trials = size(hand_data, 2);

% Load Initial Values
error = rot;

x(1) = 0; 
Amiss = params(1);
Ahit = params(2);
Umiss = params(3);
Uhit = params(4);

for n = 1:size(hand_data,1)
    
    grp=group(n);
    reward=rwd(n,:);
    
    for i = 1:num_trials-1
        if reward(i)==1
            x(i+1) = Ahit*x(i) + Uhit;
        elseif reward(i)==0
            x(i+1) = Amiss*x(i) + Umiss;
        end
    end

    simulations(n,:) = x;
    
end

end
