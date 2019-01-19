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