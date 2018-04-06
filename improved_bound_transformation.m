%% Script for figures 7 and 8.
% First compute the Negativity of the Choi-operator for the transformations
% |psi> --> |phi>, then use it to find an improved bound, which can be seen
% in the plot.
p = 0:0.1:0.5;
n = length(p);
ctr = 1;
neg = zeros(n,n);
neg1 = zeros(n,n);
neg2 = cell(n,n);
W = cell(n,n);
rho_best = W;
h = waitbar(0,'Please wait...');
for i = 1:n
    for j = 1:n        
    psi = sqrt(1-p(i)) * nn(0,0,2) + sqrt(p(i)) * nn(1,1,2); 
    phi = sqrt(1-p(j)) * nn(0,0,2) + sqrt(p(j)) * nn(1,1,2);
    rho{1} = psi*psi'; 
    sigma{1} = phi*phi'; 
    [neg(i,j),W{i,j}] = NegMapEns2(rho,2,sigma,2); 
    [neg2{i,j},rho_best{i,j}]= improved_bound1(W{i,j},[2 2 2 2],[0,0,1]);
    waitbar(ctr/n^2,h,sprintf('Progress: %5.2f %% , ctr = %d',ctr/n^2*100,ctr))
    ctr = ctr+1;
    end
end
delete(h)

for i = 1:n
    for j = 1:n
        negativity = neg2{i,j};
        neg1(i,j) = negativity(length(negativity));     
    end
end
figure
surf(p,p,neg)
alpha 0.5
hold on
surf(p,p,neg1)
xlabel('p','FontSize',20)
ylabel('q','FontSize',20) 
zlabel('Negativity','FontSize',20)
set(gca, 'XTick',0:0.1:0.5,'FontSize',17)
set(gca, 'YTick',0:0.1:0.5,'FontSize',17)
set(gca, 'ZTick',0:0.1:0.5,'FontSize',17)
zlim([0,0.5])

for i = 1:n^2
    plot(0:length(neg2{i})-1,neg2{i},'LineWidth',2)
    xlim([0,2])
    ylim([0,0.5])
    set(gca, 'XTick',0:2,'FontSize',17)
    xlabel('Iterations','FontSize',20)
    ylabel('Negativity','FontSize',20)
    set(gca,'FontSize',20)
    hold on
end
