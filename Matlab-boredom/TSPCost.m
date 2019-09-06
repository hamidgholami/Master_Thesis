%Difine cost function
function cost=TSPCost(solution,model)

    D=model.D;
    
    cost=0;
    
    if isempty(solution);
        return
    end
    
    solution=[solution solution(1)];
    for i=1:numel(solution)-1
         cost=cost+D(solution(i),solution(i+1));
    end
        
 end