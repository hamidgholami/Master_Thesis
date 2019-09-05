%IWO for Seru
%Model Inputs
    NumOp =8;
    NumCell =6;
    NumBatch =4;
    NumTask =5;
    NumRot =5;
 Bs = [55	50	54	49]; %Bm In Modelling
 EpsiOp = [0.18	0.19	0.2	0.21	0.2	0.2	0.2	0.22];
 NooOp = [1	1	1	4	1	3	2	1];
 B_BatchOp =[1.02	1.05	1.1	1.05
1.09	1.15	1.16	1.24
0.96	0.98	1.06	1.16
0.94	0.99	1.1	1.09
0.96	1.1	1.08	1.07
0.92	0.97	1.12	0.99
1.1	1.13	1.13	1.22
0.98	1.08	1.06	1.3];


   
 %Bmi In Modelling
                
 ProTime =[0.10	0.50	0.20	0.10	0.60	0.30	0.60	0.40	0.40	0.20	0.50	0.30	0.10	0.10	0.10	0.40	0.10	0.10	0.60	0.30	0.20	0.10	0.20	0.60	0.50	0.60	0.70	0.10	0.40	0.40	0.30	0.30	0.40	0.20	0.20	0.10	0.40	0.70	0.50	0.10	0.10	0.50	0.20	0.10	0.70	0.30	0.40	0.60	0.70	0.10	0.50	0.50	0.40	0.70	0.50	0.30	0.60	0.20	0.50	0.10	0.70	0.10	0.50	0.70	0.10	0.20	0.70	0.30	0.20	0.30	0.60	0.60	0.50	0.50	0.30	0.10	0.50	0.70	0.50	0.70	0.70	0.70	0.40	0.40	0.50	0.40	0.10	0.20	0.20	0.70	0.50	0.70	0.60	0.20	0.20	0.60	0.70	0.40	0.70	0.40	0.30	0.10	0.50	0.30	0.20	0.40	0.70	0.60	0.10	0.20	0.60	0.60	0.10	0.70	0.10	0.20	0.50	0.60	0.60	0.60
0.50	0.30	0.20	0.70	0.50	0.30	0.20	0.10	0.10	0.10	0.40	0.50	0.20	0.60	0.50	0.10	0.60	0.40	0.30	0.60	0.30	0.20	0.20	0.70	0.70	0.10	0.40	0.60	0.30	0.10	0.30	0.30	0.60	0.10	0.70	0.50	0.50	0.10	0.70	0.10	0.50	0.50	0.40	0.60	0.70	0.40	0.20	0.50	0.30	0.70	0.20	0.10	0.40	0.10	0.30	0.40	0.40	0.50	0.70	0.20	0.50	0.60	0.50	0.10	0.70	0.50	0.20	0.30	0.60	0.60	0.70	0.20	0.60	0.50	0.10	0.50	0.50	0.30	0.50	0.70	0.30	0.10	0.10	0.50	0.10	0.40	0.60	0.50	0.20	0.70	0.50	0.40	0.10	0.50	0.10	0.40	0.60	0.10	0.50	0.10	0.70	0.70	0.60	0.50	0.10	0.50	0.30	0.30	0.20	0.70	0.50	0.20	0.70	0.40	0.30	0.30	0.20	0.70	0.30	0.30
0.10	0.70	0.70	0.30	0.60	0.70	0.20	0.10	0.40	0.30	0.20	0.70	0.70	0.50	0.20	0.70	0.60	0.10	0.30	0.70	0.20	0.30	0.20	0.70	0.40	0.30	0.60	0.30	0.50	0.50	0.10	0.70	0.20	0.70	0.40	0.30	0.60	0.20	0.30	0.60	0.30	0.60	0.30	0.10	0.30	0.30	0.50	0.10	0.60	0.20	0.30	0.40	0.60	0.70	0.70	0.60	0.30	0.40	0.20	0.40	0.60	0.30	0.30	0.70	0.40	0.70	0.30	0.30	0.50	0.30	0.40	0.20	0.50	0.50	0.20	0.70	0.30	0.70	0.10	0.10	0.60	0.30	0.10	0.10	0.50	0.20	0.50	0.60	0.70	0.30	0.20	0.40	0.20	0.30	0.60	0.30	0.20	0.60	0.70	0.70	0.10	0.40	0.10	0.60	0.40	0.70	0.40	0.60	0.50	0.40	0.50	0.30	0.20	0.40	0.30	0.70	0.40	0.20	0.20	0.70
0.10	0.40	0.10	0.20	0.70	0.40	0.20	0.60	0.70	0.50	0.60	0.60	0.70	0.60	0.40	0.10	0.40	0.50	0.70	0.70	0.70	0.40	0.30	0.20	0.70	0.40	0.70	0.30	0.20	0.60	0.10	0.10	0.10	0.10	0.40	0.60	0.20	0.60	0.70	0.40	0.10	0.50	0.60	0.60	0.70	0.60	0.30	0.40	0.70	0.50	0.50	0.60	0.70	0.30	0.10	0.20	0.50	0.20	0.20	0.60	0.10	0.50	0.20	0.60	0.10	0.50	0.20	0.40	0.60	0.70	0.40	0.30	0.20	0.70	0.60	0.50	0.10	0.10	0.40	0.70	0.70	0.70	0.30	0.30	0.10	0.20	0.10	0.60	0.60	0.30	0.20	0.40	0.50	0.50	0.20	0.60	0.20	0.20	0.60	0.60	0.40	0.20	0.40	0.50	0.30	0.70	0.40	0.40	0.50	0.70	0.30	0.70	0.30	0.50	0.10	0.50	0.20	0.50	0.70	0.40
0.10	0.10	0.10	0.30	0.40	0.20	0.50	0.50	0.40	0.50	0.40	0.40	0.20	0.20	0.30	0.40	0.10	0.20	0.20	0.50	0.20	0.30	0.40	0.20	0.70	0.70	0.30	0.40	0.10	0.40	0.40	0.70	0.30	0.20	0.40	0.10	0.70	0.30	0.70	0.60	0.10	0.60	0.30	0.10	0.30	0.70	0.60	0.30	0.10	0.60	0.30	0.70	0.50	0.40	0.10	0.30	0.70	0.30	0.60	0.70	0.50	0.50	0.20	0.60	0.10	0.60	0.60	0.40	0.20	0.30	0.50	0.70	0.40	0.30	0.40	0.10	0.70	0.30	0.70	0.10	0.20	0.40	0.60	0.30	0.10	0.70	0.30	0.30	0.40	0.70	0.70	0.30	0.60	0.60	0.30	0.20	0.30	0.40	0.10	0.70	0.40	0.10	0.70	0.60	0.10	0.20	0.40	0.10	0.50	0.70	0.70	0.30	0.50	0.60	0.40	0.40	0.60	0.40	0.50	0.40
0.60	0.70	0.20	0.30	0.70	0.40	0.50	0.10	0.30	0.40	0.70	0.70	0.70	0.50	0.40	0.50	0.10	0.20	0.70	0.30	0.50	0.20	0.70	0.70	0.30	0.60	0.70	0.70	0.40	0.10	0.10	0.70	0.60	0.20	0.20	0.40	0.10	0.40	0.30	0.50	0.70	0.20	0.50	0.60	0.20	0.10	0.70	0.70	0.60	0.60	0.30	0.50	0.20	0.60	0.40	0.50	0.40	0.70	0.10	0.50	0.30	0.10	0.50	0.20	0.50	0.30	0.40	0.50	0.30	0.60	0.60	0.40	0.70	0.10	0.30	0.60	0.10	0.30	0.50	0.30	0.40	0.70	0.50	0.50	0.10	0.70	0.50	0.60	0.70	0.40	0.40	0.30	0.10	0.50	0.50	0.20	0.10	0.20	0.20	0.40	0.50	0.10	0.10	0.10	0.30	0.10	0.60	0.40	0.20	0.20	0.20	0.10	0.60	0.10	0.50	0.50	0.20	0.70	0.70	0.60
0.30	0.40	0.60	0.40	0.10	0.30	0.60	0.20	0.30	0.30	0.50	0.70	0.20	0.70	0.40	0.50	0.40	0.70	0.50	0.70	0.40	0.70	0.30	0.70	0.10	0.30	0.10	0.30	0.70	0.10	0.30	0.20	0.70	0.30	0.30	0.40	0.50	0.30	0.50	0.40	0.50	0.50	0.60	0.10	0.10	0.10	0.70	0.10	0.50	0.70	0.30	0.20	0.20	0.40	0.40	0.10	0.60	0.60	0.10	0.10	0.50	0.50	0.50	0.50	0.50	0.10	0.30	0.70	0.20	0.40	0.70	0.40	0.60	0.50	0.60	0.60	0.30	0.50	0.60	0.60	0.40	0.30	0.30	0.10	0.50	0.40	0.40	0.20	0.40	0.50	0.40	0.20	0.60	0.20	0.50	0.50	0.10	0.50	0.70	0.40	0.10	0.40	0.70	0.70	0.70	0.40	0.40	0.60	0.60	0.40	0.10	0.40	0.50	0.40	0.20	0.30	0.30	0.20	0.30	0.50
0.60	0.20	0.20	0.40	0.50	0.10	0.40	0.20	0.20	0.50	0.60	0.10	0.50	0.10	0.20	0.50	0.10	0.20	0.30	0.30	0.10	0.20	0.30	0.60	0.50	0.20	0.30	0.10	0.70	0.40	0.20	0.50	0.60	0.60	0.20	0.60	0.60	0.30	0.10	0.10	0.20	0.60	0.70	0.30	0.20	0.70	0.40	0.50	0.40	0.30	0.70	0.40	0.10	0.20	0.50	0.10	0.10	0.60	0.50	0.40	0.50	0.40	0.10	0.70	0.40	0.30	0.40	0.50	0.50	0.70	0.70	0.40	0.30	0.60	0.20	0.20	0.10	0.60	0.20	0.50	0.30	0.40	0.40	0.30	0.20	0.30	0.10	0.20	0.60	0.70	0.10	0.60	0.20	0.60	0.70	0.30	0.30	0.50	0.60	0.30	0.50	0.20	0.20	0.30	0.70	0.20	0.10	0.20	0.40	0.50	0.30	0.20	0.20	0.70	0.50	0.20	0.20	0.60	0.70	0.30
];


 cbar=1.8; 
 
 f2=floor(NumBatch/NumCell);

 
 time=[];
 %IWO Parameters
 nPop0=20;
 nPop =500 ;
 
 MaxIt=60;
 Exponent=0.5;
 Smin = 1;       % Minimum Number of Seeds
 Smax = 10;       % Maximum Number of Seeds
 NumOpx=4;
 
 %Initialization
empty_plant.AssignedOp= [];
empty_plant.AssignmentXIWO = [];
empty_plant.AssignmentZIWO = [];
empty_plant.AssignmentXMODEL = [];
empty_plant.AssignmentZMODEL = [];
empty_plant.Gijt = [];

empty_plant.Flow = [];
pop = repmat(empty_plant, nPop0, 1);
for i=1:numel(pop)
 w=randi([NumCell,NumCell],1);
 AssignedSet=randperm(NumOp,w);
 XIWO=zeros(NumRot,w);
 ZIWO=zeros(NumRot,NumBatch);
        for t=1:NumRot
     
    f1=floor(w/NumCell);
    for f=1:f1
    XIWO(t,1+(f-1)*NumCell:f*NumCell)=randperm(NumCell);
    end
    XIWO(t,f1*NumCell+1:w)=randperm(NumCell,w-f1*NumCell);
    for f=1:f2
    ZIWO(t,1+(f-1)*NumCell:f*NumCell)=randperm(NumCell);
    end
    ZIWO(t,f2*NumCell+1:NumBatch)=randperm(NumCell,NumBatch-f2*NumCell);
        end
    
 %XIWO TO XMODEL
     XMODEL=zeros(NumOp,NumCell*NumRot);
    for ii=1:w
        for j=1:NumCell
            for t=1:NumRot
                if XIWO(t,ii)==j
                    XMODEL(AssignedSet(ii),(j-1)*NumRot+t)=1;
                end
            end
        end
 
    end
  %ZIWO TO ZMODEL
    ZMODEL=zeros(NumBatch,NumCell*NumRot);
        for j=1:NumCell
            for t=1:NumRot
                for m=1:NumBatch
                    if ZIWO(t,m)==j
                        ZMODEL(m,(j-1)*NumRot+t)=1;
                    end
                end
            end
        end
        
   % Gijt Calculation
      Gijt=zeros(NumOp,NumCell*NumRot);
     % Wjt=zeros(1,NumCell*NumRot);
      Wjt=sum(XMODEL,1);
    for ii=1:w
        for j=1:NumCell
            for t=1:NumRot
                if Wjt((j-1)*NumRot+t)>NooOp(AssignedSet(ii))
                    Gijt(AssignedSet(ii),(j-1)*NumRot+t)=1+EpsiOp(AssignedSet(ii))*(Wjt((j-1)*NumRot+t)-NooOp(AssignedSet(ii)));
                else
                    Gijt(AssignedSet(ii),(j-1)*NumRot+t)=1;
                end
            end
        end
 
    end   
    
    %Flow calculations
    a1=zeros(NumBatch,1);
    a2=zeros(NumBatch,1);
    TCm=zeros(NumBatch,1);
    FCm=zeros(NumBatch,1);
    for ii=1:w
        for j=1:NumCell
            for t=1:NumRot
                for m=1:NumBatch
                   a1(m)=a1(m)+ Gijt(AssignedSet(ii),(j-1)*NumRot+t)*XMODEL(AssignedSet(ii),(j-1)*NumRot+t)*ZMODEL(m,(j-1)*NumRot+t)*B_BatchOp(AssignedSet(ii),m);
                   a2(m)=a2(m)+XMODEL(AssignedSet(ii),(j-1)*NumRot+t)*ZMODEL(m,(j-1)*NumRot+t);
                end
                
                   
                
            end
        end
 
    end 
    for m=1:NumBatch
        TCm(m)=a1(m)*cbar/a2(m);
        FCm(m)=TCm(m)*w*NumRot*Bs(m)/a2(m);
    end
   Flow=max(FCm);
 pop(i).AssignedOp=AssignedSet;
 pop(i).AssignmentXIWO=XIWO;
 pop(i).AssignmentZIWO=ZIWO;
 pop(i).AssignmentXMODEL=XMODEL;
 pop(i).AssignmentZMODEL=ZMODEL;
 pop(i).Gijt=Gijt;
 pop(i).Flow=Flow;
end 

 BestFlows = zeros(MaxIt, 1);
 time=zeros(MaxIt,1);
 t0=clock;
%IWO MAIN LOOP

for it=1:MaxIt
    
    % Update Standard Deviation
    sigma = floor(((MaxIt - it)/(MaxIt - 1))^Exponent * (NumRot - 1) + 1);
    sigma2 = floor(((MaxIt - it)/(MaxIt - 1))^Exponent * (NumOpx - 1) + 1);    
    % Get Best and Worst Flow Values
    Flows = [pop.Flow];
    BestFlow = min(Flows);
    WorstFlow = max(Flows);
    
    % Initialize Offsprings Population
    newpop = [];
    
        % Reproduction
    for i = 1:numel(pop)
        w=size(pop(i).AssignmentXIWO,2);
        ratio = (pop(i).Flow - WorstFlow)/(BestFlow - WorstFlow);
        S = floor(Smin + (Smax - Smin)*ratio);
        
    
    
        for j=1:S
            
            % Initialize Offspring
            
            newsol = empty_plant;
            
            % Generate Random Location
            rrr=rand;
           if rrr>((MaxIt - it)/(MaxIt - 1))^Exponent 
            newsol.AssignmentXIWO = pop(i).AssignmentXIWO + repmat(randi(NumCell-1,[NumRot ,1]),1,size(pop(i).AssignmentXIWO,2));
            newsol.AssignmentZIWO = pop(i).AssignmentZIWO + repmat(randi(NumCell-1,[NumRot ,1]),1,size(pop(i).AssignmentZIWO,2));
            for r=1:NumRot-sigma
                temp=randi(NumRot);
                newsol.AssignmentXIWO(temp,:)=pop(i).AssignmentXIWO(temp,:);
             
            end
            for r=1:NumRot-sigma
                temp=randi(NumRot);
                newsol.AssignmentZIWO(temp,:)=pop(i).AssignmentZIWO(temp,:);
             
            end            
            
            for r=1:NumRot
               for ii=1:size(pop(i).AssignmentXIWO,2)
                  if newsol.AssignmentXIWO(r,ii)>NumCell
                      newsol.AssignmentXIWO(r,ii)=newsol.AssignmentXIWO(r,ii)-NumCell;
                  end
               end
                for ii=1:size(pop(i).AssignmentZIWO,2)
                  if newsol.AssignmentZIWO(r,ii)>NumCell
                      newsol.AssignmentZIWO(r,ii)=newsol.AssignmentZIWO(r,ii)-NumCell;
                  end
                end              
            end
           newsol.AssignedOp = pop(i).AssignedOp;
           else
               if rrr>0.85*((MaxIt - it)/(MaxIt - 1))^Exponent 
            newsol.AssignmentXIWO = pop(i).AssignmentXIWO ;
            newsol.AssignmentZIWO = pop(i).AssignmentZIWO ;
            newsol.AssignedOp = pop(i).AssignedOp+randi(NumOp,1);
            for iii=1:size(pop(i).AssignedOp,2)
                if newsol.AssignedOp(1,iii)>NumOp
                    newsol.AssignedOp(1,iii)=newsol.AssignedOp(1,iii)-NumOp;
                end
            end
               else
           if size(pop(i).AssignedOp,2)<NumOp
           newsol.AssignmentXIWO=myCol(pop(i).AssignmentXIWO,randi(NumCell,[NumRot,1]));
           NonAssigned=[];
           gg=0;
           for iiii=1:NumOp
               if ismember(iiii,pop(i).AssignedOp)==0
                   %myCol(iiii,NonAssigned);
                   gg=gg+1;
                   NonAssigned(1,gg)=iiii;
               end
           end
           newsol.AssignedOp = myCol(NonAssigned(1,1),pop(i).AssignedOp);
           newsol.AssignmentZIWO = pop(i).AssignmentZIWO ;
           else
            newsol.AssignmentXIWO = pop(i).AssignmentXIWO ;
            newsol.AssignmentZIWO = pop(i).AssignmentZIWO ;
            newsol.AssignedOp = pop(i).AssignedOp;
           end
               end
           end
 %-----------------------------------------------------------------------------------------------------
  %XIWO TO XMODEL
  newsol.AssignmentXMODEL=zeros(NumOp,NumRot*NumCell);
     for ii=1:w
        for j1=1:NumCell
            for t=1:NumRot
                if newsol.AssignmentXIWO(t,ii)==j1
                    newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)=1;
                end
            end
        end
 
    end
  %ZIWO TO ZMODEL
  newsol.AssignmentZMODEL=zeros(NumBatch,NumRot*NumCell);
         for j1=1:NumCell
            for t=1:NumRot
                for m=1:NumBatch
                    if newsol.AssignmentZIWO(t,m)==j1
                        newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t)=1;
                    end
                end
            end
        end
        
   % Gijt Calculation
      %Gijt=zeros(NumOp,NumCell*NumRot);
     % Wjt=zeros(1,NumCell*NumRot);
      Wjt=sum(pop(i).AssignmentXMODEL,1);
    for ii=1:w
        for j1=1:NumCell
            for t=1:NumRot
                if Wjt((j1-1)*NumRot+t)>NooOp(pop(i).AssignedOp(ii))
                    newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)=1+EpsiOp(pop(i).AssignedOp(ii))*(Wjt((j1-1)*NumRot+t)-NooOp(pop(i).AssignedOp(ii)));
                else
                    newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)=1;
                end
            end
        end
 
    end 
     
 %-----------------------------------------------------------------------------------------------------
            % Evaluate Offsring
            % newsol.Flow = ?
    a1=zeros(NumBatch,1);
    a2=zeros(NumBatch,1);
    TCm=zeros(NumBatch,1);
    FCm=zeros(NumBatch,1);
    for ii=1:w
        for j1=1:NumCell
            for t=1:NumRot
                for m=1:NumBatch
                   a1(m)=a1(m)+ newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t)*B_BatchOp(pop(i).AssignedOp(ii),m);
                   a2(m)=a2(m)+newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t);
                end
                
                   
                
            end
        end
 
    end 
    for m=1:NumBatch
        TCm(m)=a1(m)*cbar/a2(m);
        FCm(m)=TCm(m)*w*NumRot*Bs(m)/a2(m);
    end
   Flow=max(FCm);    

   newsol.Flow=Flow;        
           
%-----------------------------------------------------------------------------------------------------
            % Add Offpsring to the Population
            newpop = [newpop
                      newsol];  %#ok 
           
           
        end
       
            
    end 
    % Merge Populations
    pop = [pop
           newpop];
    
    % Sort Population
    [~, SortOrder]=sort([pop.Flow]);
    pop = pop(SortOrder);

    % Competitive Exclusion (Delete Extra Members)
    if numel(pop)>nPop
        pop = pop(1:nPop);
    end
    
    % Store Best Solution Ever Found
    BestSol = pop(1);
    
    % Store Best Cost History
    BestFlows(it) = BestSol.Flow;
    time(it) =etime(clock,t0);
    
       
end
 
%% Results

figure;
%plot(BestFlows,'LineWidth',2);
 
%xlabel('Iteration');
%ylabel('Best Flow');

plot(time,BestFlows);
xlabel('Computation time(Sec)');
ylabel('Best Flow');
