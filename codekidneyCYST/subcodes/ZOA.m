function[Best_score,Best_pos,ZOA_curve]=ZOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)

lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%% INITIALIZATION
for i=1:dimension
    X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
end

for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end
%%

for t=1:Max_iterations
    %% update the global best (fbest)
    [best , location]=min(fit);
    if t==1
        PZ=X(location,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        PZ=X(location,:);
    end
    
    %% PHASE1: Foraging Behaviour
    for i=1:SearchAgents
        
        I=round(1+rand);
        X_newP1=X(i,:)+ rand(1,dimension).*(PZ-I.* X(i,:)); %Eq(3)
        X_newP1= max(X_newP1,lowerbound);X_newP1 = min(X_newP1,upperbound);
        
        
        % Updating X_i using (5)
        f_newP1 = fitness(X_newP1);
        if f_newP1 <= fit (i)
            X(i,:) = X_newP1;
            fit (i)=f_newP1;
        end

    end
    Ps=rand;
    k=randperm(SearchAgents,1);
    AZ=X(k,:);% attacked zebra
    
    for i=1:SearchAgents
        
        if Ps<0.5
   
            R=0.1;
            X_newP2= X(i,:)+ R*(2*rand(1,dimension)-1)*(1-t/Max_iterations).*X(i,:);% Eq.(5) S1
            X_newP2= max(X_newP2,lowerbound);X_newP2 = min(X_newP2,upperbound);
      
        else
      
            
            I=round(1+rand(1,1));
            X_newP2=X(i,:)+ rand(1,dimension).*(AZ-I.* X(i,:)); 
            X_newP2= max(X_newP2,lowerbound);X_newP2 = min(X_newP2,upperbound);
             
        end
        
        f_newP2 = fitness(X_newP2); 
        if f_newP2 <= fit (i)
            X(i,:) = X_newP2;
            fit (i)=f_newP2;
        end

    end 
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end 
Best_score=fbest;
Best_pos=PZ;
ZOA_curve=best_so_far;
end

