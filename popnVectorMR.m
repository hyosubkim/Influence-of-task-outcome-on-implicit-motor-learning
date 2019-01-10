function simulations = popnVectorMR(params, hand_data, rot, group,rwd)

A = params(1);
U = params(2);
Ar = params(3);
s = params(4);

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
        
        %state estimate update
        v(i+1) = A*v(i) + U;
        
        rwdAng = round(y(i),2); %hand angle when you were hitting target; round to 2nd decimal
        rwdAng = rwdAng*deg2idx; %to get accurate index
        if rwdAng<0
            kidx=round(rwdAng+nAngles+1); %bc 0 deg is 1st element
        else
            kidx=round(rwdAng+1);
        end
        
        %increase weight on directionally-tuned "neuron" representing rewarded angle
        if reward(i)==1
            k(kidx)=k(kidx)+s;
        elseif reward(i)==0
            k=k;
        end
        
        %bump up weight on rewarded angle
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
        
        %actual movement direction: Popn vector LENGTH denotes weight on
        %reinforcement system
        y(i+1) = (1-PL(i+1))*(v(i+1)) + PL(i+1)*PA(i+1);
                
        k = zeros(nAngles,1);
        
    end
    
    simulations(n,:) = y;
    
end

end



