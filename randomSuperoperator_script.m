%% Script for Figure 9
%Generate random Choi operators and plot the bounds as a function of the
%number of iterations it takes to converge.
n = 10;
neg = cell(1,n);
h = waitbar(0,'Please wait...');
ctr = 1;
ctrmax = n;
    for i = 1:n
        X = RandomSuperoperator([4 4],1,0,0,8)/4;
        X = PermuteSystems(X,[3 4 1 2],[2 2 2 2]);     %CDAB
        [neg{i},~] = improved_bound1(X,[2 2 2 2],[0,0,1]);
        waitbar(ctr/(ctrmax),h,sprintf('Progress: %5.2f %% , ctr = %d',ctr/(ctrmax)*100,ctr))
        ctr = ctr+1;      
    end
    delete(h)
hold off
for i = 1:n
    plot(0:length(neg{i})-1,neg{i},'LineWidth',1)
    xlim([0,4]);
    xlabel('Iterations','FontSize',20)
    ylabel('Negativity','FontSize',20)
    set(gca,'FontSize',20)
    hold on
end
