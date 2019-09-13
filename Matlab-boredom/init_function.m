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
    opn =8;
    ojn =6;
    nump =4;
    tsk =5;
    nrt =5;
 EpsiOp = [0.18	0.19	0.2	0.21	0.2	0.2	0.2	0.22];
 %NooOp = [1	1	1	4	1	3	2	1];
 NooOp = BestSolution.Assigned.result;          
 ProTime = pop(1).BoredomResults;
 cbar=1.8; 
 f2=floor(nump/ojn);
 time=[];
 %% IWO Parameters
 nPop0=20;
 nPop =500;
 MaxIt=100;
 Exponent=0.5;
 Smin = 1;       % Minimum Number of Seeds
 Smax = 20;       % Maximum Number of Seeds
 opnx=4;
 %% Initialization
empty_plant.AssignedOp= [];
empty_plant.AssignmentXIWO = [];
empty_plant.AssignmentZIWO = [];
empty_plant.AssignmentXMODEL = [];
empty_plant.AssignmentZMODEL = [];
empty_plant.Gijt = [];
empty_plant.Flow = [];
Bs = [55	50	54	49];
pop = repmat(empty_plant, nPop0, 1);
for i=1:numel(pop)
 w=randi([ojn,ojn],1);
 AssignedSet=randperm(opn,w);
 XIWO=zeros(nrt,w);
 ZIWO=zeros(nrt,nump);
        for t=1:nrt
     
    f1=floor(w/ojn);
    for f=1:f1
    XIWO(t,1+(f-1)*ojn:f*ojn)=randperm(ojn);
    end
    XIWO(t,f1*ojn+1:w)=randperm(ojn,w-f1*ojn);
    for f=1:f2
    ZIWO(t,1+(f-1)*ojn:f*ojn)=randperm(ojn);
    end
    ZIWO(t,f2*ojn+1:nump)=randperm(ojn,nump-f2*ojn);
        end
    
 %XIWO TO XMODEL
     XMODEL=zeros(opn,ojn*nrt);
    for ii=1:w
        for j=1:ojn
            for t=1:nrt
                if XIWO(t,ii)==j
                    XMODEL(AssignedSet(ii),(j-1)*nrt+t)=1;
                end
            end
        end
 
    end
  %ZIWO TO ZMODEL
    ZMODEL=zeros(nump,ojn*nrt);
        for j=1:ojn
            for t=1:nrt
                for m=1:nump
                    if ZIWO(t,m)==j
                        ZMODEL(m,(j-1)*nrt+t)=1;
                    end
                end
            end
        end
        
   % Gijt Calculation
      Gijt=zeros(opn,ojn*nrt);
     % Wjt=zeros(1,ojn*nrt);
      Wjt=sum(XMODEL,1);
    for ii=1:w
        for j=1:ojn
            for t=1:nrt
                if Wjt((j-1)*nrt+t)>NooOp(AssignedSet(ii))
                    Gijt(AssignedSet(ii),(j-1)*nrt+t)=1+EpsiOp(AssignedSet(ii))*(Wjt((j-1)*nrt+t)-NooOp(AssignedSet(ii)));
                else
                    Gijt(AssignedSet(ii),(j-1)*nrt+t)=1;
                end
            end
        end
 
    end   
    
    %Flow calculations
    a1=zeros(nump,1);
    a2=zeros(nump,1);
    TCm=zeros(nump,1);
    FCm=zeros(nump,1);
    for ii=1:w
        for j=1:ojn
            for t=1:nrt
                for m=1:nump
                   a1(m)=a1(m)+ Gijt(AssignedSet(ii),(j-1)*nrt+t)*XMODEL(AssignedSet(ii),(j-1)*nrt+t)*ZMODEL(m,(j-1)*nrt+t);
                   a2(m)=a2(m)+XMODEL(AssignedSet(ii),(j-1)*nrt+t)*ZMODEL(m,(j-1)*nrt+t);
                end
            end
        end
 
    end 
    for m=1:nump
        TCm(m)=a1(m)*cbar/a2(m);
        FCm(m)=TCm(m)*w*nrt*Bs(m)/a2(m);
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
    sigma = floor(((MaxIt - it)/(MaxIt - 1))^Exponent * (nrt - 1) + 1);
    sigma2 = floor(((MaxIt - it)/(MaxIt - 1))^Exponent * (opnx - 1) + 1);    
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
            newsol.AssignmentXIWO = pop(i).AssignmentXIWO + repmat(randi(ojn-1,[nrt ,1]),1,size(pop(i).AssignmentXIWO,2));
            newsol.AssignmentZIWO = pop(i).AssignmentZIWO + repmat(randi(ojn-1,[nrt ,1]),1,size(pop(i).AssignmentZIWO,2));
            for r=1:nrt-sigma
                temp=randi(nrt);
                newsol.AssignmentXIWO(temp,:)=pop(i).AssignmentXIWO(temp,:);
             
            end
            for r=1:nrt-sigma
                temp=randi(nrt);
                newsol.AssignmentZIWO(temp,:)=pop(i).AssignmentZIWO(temp,:);
             
            end            
            
            for r=1:nrt
               for ii=1:size(pop(i).AssignmentXIWO,2)
                  if newsol.AssignmentXIWO(r,ii)>ojn
                      newsol.AssignmentXIWO(r,ii)=newsol.AssignmentXIWO(r,ii)-ojn;
                  end
               end
                for ii=1:size(pop(i).AssignmentZIWO,2)
                  if newsol.AssignmentZIWO(r,ii)>ojn
                      newsol.AssignmentZIWO(r,ii)=newsol.AssignmentZIWO(r,ii)-ojn;
                  end
                end              
            end
           newsol.AssignedOp = pop(i).AssignedOp;
           else
               if rrr>0.85*((MaxIt - it)/(MaxIt - 1))^Exponent 
            newsol.AssignmentXIWO = pop(i).AssignmentXIWO ;
            newsol.AssignmentZIWO = pop(i).AssignmentZIWO ;
            newsol.AssignedOp = pop(i).AssignedOp+randi(opn,1);
            for iii=1:size(pop(i).AssignedOp,2)
                if newsol.AssignedOp(1,iii)>opn
                    newsol.AssignedOp(1,iii)=newsol.AssignedOp(1,iii)-opn;
                end
            end
               else
           if size(pop(i).AssignedOp,2)<opn
           newsol.AssignmentXIWO=myCol(pop(i).AssignmentXIWO,randi(ojn,[nrt,1]));
           NonAssigned=[];
           gg=0;
           for iiii=1:opn
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
  newsol.AssignmentXMODEL=zeros(opn,nrt*ojn);
     for ii=1:w
        for j1=1:ojn
            for t=1:nrt
                if newsol.AssignmentXIWO(t,ii)==j1
                    newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*nrt+t)=1;
                end
            end
        end
 
    end
  
  newsol.AssignmentZMODEL=zeros(nump,nrt*ojn);
         for j1=1:ojn
            for t=1:nrt
                for m=1:nump
                    if newsol.AssignmentZIWO(t,m)==j1
                        newsol.AssignmentZMODEL(m,(j1-1)*nrt+t)=1;
                    end
                end
            end
        end
        
      Wjt=sum(pop(i).AssignmentXMODEL,1);
    for ii=1:w
        for j1=1:ojn
            for t=1:nrt
                if Wjt((j1-1)*nrt+t)>NooOp(pop(i).AssignedOp(ii))
                    newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*nrt+t)=1+EpsiOp(pop(i).AssignedOp(ii))*(Wjt((j1-1)*nrt+t)-NooOp(pop(i).AssignedOp(ii)));
                else
                    newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*nrt+t)=1;
                end
            end
        end
 
    end 
     
 %-----------------------------------------------------------------------------------------------------
            % Evaluate Offsring
            % newsol.Flow = ?
    a1=zeros(nump,1);
    a2=zeros(nump,1);
    TCm=zeros(nump,1);
    FCm=zeros(nump,1);
    for ii=1:w
        for j1=1:ojn
            for t=1:nrt
                for m=1:nump
                   a1(m)=a1(m)+ newsol.Gijt(pop(i).AssignedOp(ii),(j1-1)*nrt+t)*newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*nrt+t)*newsol.AssignmentZMODEL(m,(j1-1)*nrt+t);
                   a2(m)=a2(m)+newsol.AssignmentXMODEL(pop(i).AssignedOp(ii),(j1-1)*nrt+t)*newsol.AssignmentZMODEL(m,(j1-1)*nrt+t);
                end
            end
        end
 
    end 
    for m=1:nump
        TCm(m)=a1(m)*cbar/a2(m);
        FCm(m)=TCm(m)*w*nrt*Bs(m)/a2(m);
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
    BestFlows(it) = (BestSol.Flow-365);
    time(it) =etime(clock,t0);
    
       
end
 
%% Results

figure;
plot(BestFlows,'LineWidth',2);
 xlabel('Iteration');
 ylabel('Best Solution');

figure;
plot(time,BestFlows,'LineWidth',2);
xlabel('Computation time(Sec)');
ylabel('Best Solution');


%ga;