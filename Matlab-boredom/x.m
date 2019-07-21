function x = x(n)

    n = 2;
    m = 6;
    a = rand(n,1);
    x = round(a/sum(a)*m);
    x(1) = x(1) - sum(x) + m;
end








%    n=[1 2];
 %  for i=1:n
 %      if n==1 | n==2
 %           x=3;
 %      else
 %          x=3;
 %      end
 %  end
% end
    
  %  for ii=1:2
 %       for nn=1:2
 %            Wtt=Wtt+Cin(ii);
 %       end
       
 %   end

%end