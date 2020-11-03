function [y1, y2]=Crossover(x1,x2)
    
        %mat2str(randi([0 100],1,20))
    %x1=[82 91 12 92 63 9 28 55 96 97 15 98 96 49 80 14 42 92 80 96];
    %x2=[75 25 51 70 89 96 55 14 15 26 84 25 82 24 93 35 19 25 62 47];
    % Determine Length
    n=numel(x1);
    
    % Determine Crossover Point(Cut point)
    c=randi([1 n-1]);
    
    % Single-point Crossover
    y1=[x1(1:c) x2(c+1:end)];
    y2=[x2(1:c) x1(c+1:end)];
    
    R1=find(ismember(x2(c+1:end),x1(1:c)))+c;
    R2=find(ismember(x1(c+1:end),x2(1:c)))+c;

    for k=1:numel(R1)
        
        i=R1(k);
        j=R2(k);
        
        temp=y1(i);
        y1(i)=y2(j);
        y2(j)=temp;
        
   end
    
end