function model=init()

    result = [randperm(3),randperm(3,2);randperm(3),randperm(3,2)];
   % result = [randi([1 3],1,5);randi([1 3],1,5)];
    model.result=result;

end