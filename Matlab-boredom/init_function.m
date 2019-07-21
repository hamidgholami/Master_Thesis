clear all
close all
clc

%%

i=[1 2];
j=[1 2 3 4];
r=[1 2 3 4 5];
n=[1 2];

Nj=4;
Nb=2;
RI=5;
%Cin=[1 1;1 0];
Cin=randi([0 100],2,20);
Pjk=[1 0.20 1 0.46;0.20 1 0.71 0.08;1 0.71 1 0.90;0.46 0.08 0.90 1];
X=[1 4];

%%  WTT  
%Cin=[1 2;1 2];
    Wtt=0;
 for ii=1:20
    for nn=1:20
        Wtt=Wtt+Cin(ii);
    end
 end
    
%% Wi
%Cin=[1 2;1 2];
Wi=0;
for nn=1:20;
    Wi=Wi+Cin(nn);
end
%% WTi

WTi=Wi/Wtt;

fprintf('Wtt = %f\n',Wtt);fprintf('Wi = %f\n',Wi);fprintf('WTi = %f\n',WTi);
%% airj
%syms i r j
%airj=symsum(i*r*j,i,1,5)==1
%Airj=symsum(i*r*j,j,1,2)>=1

%% Xn

xn= x()
%syms n
%RI=symsum(x,n,1,2)
%% l(n)

    l=0;
   for m=1:n-1
       l=l+x(n); 
   end

   
 
n = 2;
m = 6;
a = rand(n,1);
x = round(a/sum(a)*m);
%out(1) = x(1) - sum(out) + m;






%%
%---------------------------------------------------------
%A_i_j=[6 5 2 1 9;2 4 8 4 8;7 3 6 1 4;3 2 2 3 7;9 1 5 8 9];
%B_m=[9 3 6 8];

%C=0;
%for mm=1:4
 %   for ii=1:5
  %      for jj=1:5
   %         C=C+B_m(mm)*A_i_j(jj,ii);
    %    end
   % end
%end