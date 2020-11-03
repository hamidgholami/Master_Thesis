%Difine mutate function
function y=Mutate(x)

    pSwap=0.3;
    pReversion=0.45;
    pInsertion=1-pSwap-pReversion;
    
    %»—œ«— «Õ „«·? «?Ã«œ „?ò‰?„ òÂ »« «Õ „«·0.2Œ—ÊÃ?1 Ê»« «Õ „«·0.35 Œ—ÊÃ?2
     METHOD=RouletteWheelSelection([pSwap pReversion pInsertion]);
    
     switch METHOD
         
         case 1
             y=DoSwap(x);
         case 2
             y=DoReversion(x);
         case 3
             y=DoInsertion(x);
     end
end

%Do Swap
function y=DoSwap(x)

    n=numel(x);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    y=x;
    y([i1 i2])=x([i2 i1]);
    
end

% Do Reversion
function y=DoReversion(x)
    
    n=numel(x);
    
    i=randsample(n,2);
    i1=min(i);
    i2=max(i);
    
    y=x;
    y(i1:i2)=x(i2:-1:i1);

end

%Do Insertion
function y=DoInsertion(x)

    n=numel(x);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    if i1<i2
        %1 ... i1-1 i1 i1+1 ...... i2-1 i2 i2+1 .... n
        y=x([1:i1-1 i1+1:i2 i1 i2+1:end]);
    else
        %1 ... i2-1 i2 i2+1 ...... i1-1 i1 i1+1 .... n
        y=x([1:i2 i1 i2+1:i1-1 i1+1:end]);
    end
    
end