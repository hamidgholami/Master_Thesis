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

nPop=100;
empty_individual.Assigned= [];
empty_individual.Sresults= [];
empty_individual.Xresults= [];
empty_individual.BoredomResults= [];
%empty_individual.BoredomResults.results= [];
pop=repmat(empty_individual,nPop,1);

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

%%
% Store Best Solution
%WorstSolution=pop(100).BR;

% clc;
% disp('salam');
% if pop(100).Assigned.result(2,1) == 1
% disp('a bozorg ast');
% end
% 
% clc;
% for i=1:3
%     for ii=1:2
%         for jj=1:5
%             if pop(i).Assigned.result(ii,jj) == 1 & pop(i).Assigned.result(ii,jj) >= 1
%             disp('a bozorg ast');
%             end
%         end
%     end
% end

%n=input('yek adad vared konid');
%init = [randperm(3),randperm(3,2);randperm(3),randperm(3,2)]
% for i=1:numPop
%     init=[randperm(3),randperm(3,2);randperm(3),randperm(3,2)];
% end
%[~, out] = sort(rand(j,n),2);
% i=[1 2];
% j=[1 2 3 4];
% r=[1 2 3 4 5];
% n=[1 2];
% Nj=4;
% Nb=2;
% RI=5;
% %Cin=[1 1;1 0];
% Cin=randi([0 100],2,20);
% Pjk=[1 0.20 1 0.46;0.20 1 0.71 0.08;1 0.71 1 0.90;0.46 0.08 0.90 1];
% X=[1 4];
% %%  WTT  
% %Cin=[1 2;1 2];
%     Wtt=0;
%  for ii=1:20
%     for nn=1:20
%         Wtt=Wtt+Cin(ii);
%     end
%  end
%     
% %% Wi
% %Cin=[1 2;1 2];
% Wi=0;
% for nn=1:20;
%     Wi=Wi+Cin(nn);
% end
% %% WTi
% 
% WTi=Wi/Wtt;
% 
% fprintf('Wtt = %f\n',Wtt);fprintf('Wi = %f\n',Wi);fprintf('WTi = %f\n',WTi);
% %% airj
% %syms i r j
% %airj=symsum(i*r*j,i,1,5)==1
% %Airj=symsum(i*r*j,j,1,2)>=1
% %% Xn
% xn= x()
% %syms n
% %RI=symsum(x,n,1,2)
% %% l(n)
%     l=0;
%    for m=1:n-1
%        l=l+x(n); 
%    end
% n = 2;
% m = 6;
% a = rand(n,1);
% x = round(a/sum(a)*m);
% %out(1) = x(1) - sum(out) + m;
% %%
% %---------------------------------------------------------
% %A_i_j=[6 5 2 1 9;2 4 8 4 8;7 3 6 1 4;3 2 2 3 7;9 1 5 8 9];
% %B_m=[9 3 6 8];
% %C=0;
% %for mm=1:4
%  %   for ii=1:5
%   %      for jj=1:5
%    %         C=C+B_m(mm)*A_i_j(jj,ii);
%     %    end
%    % end
% %end