function [Positions,Best_Score, BestFit]=updates(popsize,Max_iter,Lmax,dim, thresholds, probR)
n=popsize;
t=0;% Loop counter
p=[1:360];

Convergence_curve=zeros(1,Max_iter);
best_score = -inf;
 ub = ones(1,1) * Lmax;
 lb = ones(1,1);

throughoutgbest_pos = zeros(1,dim);
Positions = rand(popsize,dim).*(ub-lb)+lb;
Positions = floor(Positions);

 for si=1:length(Positions)
       Positions(si,:)=sort(Positions(si,:)); 
 end

BestFit=zeros(1,dim);  % bestfit
Best_Score=-inf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ABC.M   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ObjVal = otsuobj(Positions,thresholds,probR);
 fitness   = calculateFitness(ObjVal);


for i=1:size(Positions,1)          
        if fitness(i,1) > Best_Score % Change this to > for maximization problem
            Best_Score=fitness(i,1); % Update alpha
            BestFit=Positions(i,:);
        end        
    end


while t<Max_iter
  for i=1:size(Positions,1)
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
 
        
        ObjVal = otsuobj(Positions(i,:),thresholds,probR);
        fitness(i,1)   = calculateFitness(ObjVal);

        %fitness(i,1)=otsuobj(Positions(i,:),thresholds,probR);
        if fitness(i,1)>Best_Score
            Best_Score=fitness(i,1);
            BestFit=Positions(i,:);
        end
    end
    
    S=2;                                    %%% S is maximum Sensitivity range 
    rg=S-((S)*t/(Max_iter));                %%%% guides R
    
    for i=1:size(Positions,1)
        r=rand*rg;
        R=((2*rg)*rand)-rg;                 %%%%   controls to transtion phases  
        for j=1:size(Positions,2)
        teta=RouletteWheelSelection(p);
           if((-1<=R)&&(R<=1))              %%%% R value is between -1 and 1 EXPLOITATION
                if(rg>=0 && rg<=1) % vary rg range 1. 
                    y=rand();
                    randomwalk=tan(pi*(y-0.5));  
                    Rand_position=abs(rand*BestFit(j)-Positions(i,j)+randomwalk*rg);
                    Positions(i,j)=BestFit(j)-r*Rand_position*cos(teta);
                
                else
                    Rand_position=abs(rand*BestFit(j)-Positions(i,j));
                    Positions(i,j)=BestFit(j)-r*Rand_position*cos(teta);
                end
           else 
                  cp=floor(popsize*rand()+1);
                  CandidatePosition =Positions(cp,:);
                  Positions(i,j)=r*(CandidatePosition(j)-rand*Positions(i,j));
            end
        end
    end

 for i=1:size(Positions,1)
     
         for j=1:dim                 
             Positions(i,j) = abs(ceil(Positions(i,j)));
             
            if Positions(i,j) < 1 || Positions(i,j) > 255
                rndp = randperm(255);                 
               Positions(i,j)= rndp(i);
            end
         
         end
  end



 for si=1:length(Positions)
       Positions(si,:)=sort(Positions(si,:)); 
 end
  
ObjVal = otsuobj(Positions,thresholds,probR);
 fitness   = calculateFitness(ObjVal);
 
  for i=1:size(Positions,1)         
        if fitness(i,1) > Best_Score % Change this to > for maximization problem
            Best_Score = fitness(i,1); % Update alpha
            BestFit = Positions(i,:);
        end        
end
  
     t=t+1;
    Convergence_curve(t) = Best_Score;
   
end