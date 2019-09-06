
%Difine function
function model=TSPModel()

    %mat2str(randi([0 100],1,20))
    x = [66 3 85 94 68 76 75 39 66 17 71 3 27 4 9 83 70 32 95 3]';
    y = [44 38 77 80 18 49 45 65 71 76 27 68 66 16 12 50 96 34 59 22]';
    
    n = numel(x);
    
    position = [x y];
       %Distance matrix
    D = pdist2(position,position);
    
    model.n=n;
    model.x=x;
    model.y=y;
    model.D=D;
    
end