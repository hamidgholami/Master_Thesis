clear all
close all
clc

%% init
%   i=3;          % *Number of operators*
%   j=5;          % *Number of jobs*
%   r=6;          % *Number of days*
%   n=2;          % *Number of assigned period*
%Cin=[1 1;1 0];    % *Operator Preferences*
%Cin=randi([0 1], 5,5); %or
Cin=[1 1 1 0 1;0 0 1 1 1;1 1 1 1 1;1 1 1 1 1;0 0 0 0 1];
%Cin=[0     1     0     1     1;0     1     1     0     0];
Pj1j2=[1 0.88 0.99 0.46 0.45;0.88 1 0.71 0.55 0.44;0.99 0.71 1 0.90 0.98;0.46 0.55 0.90 1 0.58;0.45 0.44 0.98 0.58 1];

nPop=200;
empty_individual.Assigned= [];
empty_individual.Sresults= [];
empty_individual.Xresults= [];
empty_individual.BoredomResults= [];
%empty_individual.BoredomResults.results= [];
pop=repmat(empty_individual,nPop,1);
%ga;
for i=1:nPop
   
    pop(i).Assigned=init();
    
    disp(['Iteration ' num2str(i) ': Result = ' mat2str(pop(i).Assigned.result)]);
  
end

%% Loop
s=0;
for i=1:nPop
    for ii=1:2
        for j1=1:5
            for j2=1:5
                if pop(i).Assigned.result(ii,j1) == 1 & pop(i).Assigned.result(ii,j2) == 1
                    s=s+Pj1j2;
                    pop(i).Sresults=s;
                    X=(1+s*Pj1j2)/(s+1);
                    pop(i).Xresults=X;
                    Boredom=((1-Cin)*X)+(Cin*(1-X));
                    pop(i).BoredomResults=Boredom;
                    %pop(i).BoredomResults.results=Boredom;
                    disp(['Iteration ' num2str(i) sprintf('    \n s is: ') mat2str(pop(i).Sresults) sprintf('    \n X is: ') mat2str(pop(i).Xresults) sprintf('    \n Boredom is: ')  mat2str(pop(i).BoredomResults)]);
                    %disp(['Iteration ' num2str(i) mat2str([pop(i).Assigned.result(ii,j1),pop(i).Assigned.result(ii,j2)])]);
                    %disp(['Iteration ' num2str(i) ' pop(' num2str(i) ').Assigned.result(' num2str(ii) ',' num2str(j1) 'and' num2str(j2) ')' ]);
                    %disp(s);
                end
            end
        end
    end
end

%% Min Function

% Sort Population
BR=[pop.BoredomResults];
[BR, SortOrder]=sort(BR,'descend');
pop=pop(SortOrder);

% Store Best Solution
BestSolution=pop(1);

% Array to Hold Best solution Values
BestBR=zeros(nPop,1);

disp('Best Solution is: ');disp([struct2table(BestSolution.Assigned)]);
%struct2table(BestSolution.Assigned)

%% Model Inputs
    NumOp =8;
    NumCell =6;
    NumPeriod =4;
    NumTask =5;
    NumRot =5;
 Bs = [55	50	54	49]; %Bm In Modelling
 EpsiOp = [0.18	0.19	0.2	0.21	0.2	0.2	0.2	0.22];
 %NooOp = [1	1	1	4	1	3	2	1];
 NooOp = BestSolution.Assigned.result;          
%ProTime =[1 0.88 0.99 0.46 0.45;0.88 1 0.71 0.55 0.44;0.99 0.71 1 0.90 0.98;0.46 0.55 0.90 1 0.58;0.45 0.44 0.98 0.58 1];
ProTime = pop(1).BoredomResults;

 cbar=1.8; 
 
 f2=floor(NumPeriod/NumCell);

 
 time=[];
 %% IWO Parameters
 nPop0=20;
 nPop =500;
 
 MaxIt=100;
 Exponent=0.5;
 Smin = 1;       % Minimum Number of Seeds
 Smax = 20;       % Maximum Number of Seeds
 NumOpx=4;
 %% Initialization
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
 ZIWO=zeros(NumRot,NumPeriod);
        for t=1:NumRot
     
    f1=floor(w/NumCell);
    for f=1:f1
    XIWO(t,1+(f-1)*NumCell:f*NumCell)=randperm(NumCell);
    end
    XIWO(t,f1*NumCell+1:w)=randperm(NumCell,w-f1*NumCell);
    for f=1:f2
    ZIWO(t,1+(f-1)*NumCell:f*NumCell)=randperm(NumCell);
    end
    ZIWO(t,f2*NumCell+1:NumPeriod)=randperm(NumCell,NumPeriod-f2*NumCell);
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
    ZMODEL=zeros(NumPeriod,NumCell*NumRot);
        for j=1:NumCell
            for t=1:NumRot
                for m=1:NumPeriod
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
    a1=zeros(NumPeriod,1);
    a2=zeros(NumPeriod,1);
    TCm=zeros(NumPeriod,1);
    FCm=zeros(NumPeriod,1);
    for ii=1:w
        for j=1:NumCell
            for t=1:NumRot
                for m=1:NumPeriod
                   a1(m)=a1(m)+ Gijt(AssignedSet(ii),(j-1)*NumRot+t)*XMODEL(AssignedSet(ii),(j-1)*NumRot+t)*ZMODEL(m,(j-1)*NumRot+t);%*B_BatchOp(AssignedSet(ii),m);
                   a2(m)=a2(m)+XMODEL(AssignedSet(ii),(j-1)*NumRot+t)*ZMODEL(m,(j-1)*NumRot+t);
                end
            end
        end
 
    end 
    for m=1:NumPeriod
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
  
  newsol.AssignmentZMODEL=zeros(NumPeriod,NumRot*NumCell);
         for j1=1:NumCell
            for t=1:NumRot
                for m=1:NumPeriod
                    if newsol.AssignmentZIWO(t,m)==j1
                        newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t)=1;
                    end
                end
            end
        end
        
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
    a1=zeros(NumPeriod,1);
    a2=zeros(NumPeriod,1);
    TCm=zeros(NumPeriod,1);
    FCm=zeros(NumPeriod,1);
    for ii=1:w
        for j1=1:NumCell
            for t=1:NumRot
                for m=1:NumPeriod
                   a1(m)=a1(m)+ newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t);%*B_BatchOp(pop(i).AssignedOp(ii),m);
                   a2(m)=a2(m)+newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*NumRot+t)*newsol.AssignmentZMODEL(m,(j1-1)*NumRot+t);
                end
            end
        end
 
    end 
    for m=1:NumPeriod
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
    BestFlows(it) = (BestSol.Flow/150)-2;
    time(it) =etime(clock,t0);
    
       
end
 
%% Results

figure;
plot(BestFlows,'LineWidth',2);
 xlabel('Iteration');
 ylabel('Best Flow');

figure;
plot(time,BestFlows,'LineWidth',2);
xlabel('Computation time(Sec)');
ylabel('Best Flow');


%ga;