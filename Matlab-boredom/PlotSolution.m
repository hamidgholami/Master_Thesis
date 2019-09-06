function PlotSolution(sol,model)

    x=model.x;
    y=model.y;
    n=model.n;
    sol=[sol sol(1)];
    
    plot(x(sol),y(sol),'b-o','MarkerSize',12,'MarkerFaceColor','r','LineWidth',2);
    axis equal;
%      hold off;
%      grid on;
%     
%      for i=1:n
%          text(x(sol(i)),y(sol(i)),num2str(i),'FontSize',20);
%      end
end