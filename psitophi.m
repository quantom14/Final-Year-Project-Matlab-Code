%% General Case to plot FIG.3
p = 0:0.05:0.5;
n = length(p);
neg = zeros(n,n);
h = waitbar(0,'Please wait...');
ctr = 0;
for i = 1:n
    for j = 1:n
        psi = sqrt(1-p(i)) * nn(0,0,2) + sqrt(p(i)) * nn(1,1,2); 
        phi = sqrt(1-p(j)) * nn(0,0,2) + sqrt(p(j)) * nn(1,1,2); 
        rho{1} = psi*psi'; 
        sigma{1} = phi*phi'; 
   
        [neg(i,j)] = NegMapEns(rho,2,sigma,2); 
        ctr = ctr+1;
        waitbar(ctr/n^2,h,sprintf('Progress: %5.3f %% , ctr = %d' ,ctr/n^2*100,ctr))
    end
end
close(h)
figure
surf(p,p,neg')
xlabel('p','FontSize',20)
ylabel('q','FontSize',20)
zlabel('Negativity','FontSize',20)
set(gca, 'XTick',0:0.1:0.5,'YTick',0:0.1:0.5,'ZTick',0:0.025:0.125,'FontSize',18)

%% Case 1: product state --> varying state (FIG.4)
p = 0:0.001:0.5;
n = length(p);
neg = zeros(1,n);
h = waitbar(0,'Please wait...');
for i = 1:n
    psi = sqrt(1-p(i)) * nn(0,0,2) + sqrt(p(i)) * nn(1,1,2); 
    rho{1} = nn(0,0,2)*nn(0,0,2)'; 
    sigma{1} = psi*psi'; 

    [neg(i)] = NegMapEns(rho,2,sigma,2); 
    waitbar(i/n,h,sprintf('Progress: %5.2f %% , ctr = %d' ,i/n*100,i))
end
plot(p,neg(1,:),'LineWidth',2)
xlabel('q','FontSize',20)
ylabel('Negativity','FontSize',20)
ax = gca;
ax.XAxis.TickLabelFormat = '%,.1f';
set(gca, 'YTick',0:0.025:0.125,'FontSize',20)
ylim([0,0.125])

%% Case 2: varying state --> maximally entangled state (FIG.5)
p = 0:0.025:0.5;
n = length(p);
neg = zeros(1,n);
h = waitbar(0,'Please wait...');
for i = 1:n
    for j = 1:n
        psi = sqrt(1-p(i)) * nn(0,0,2) + sqrt(p(i)) * nn(1,1,2); 
        rho{1} = psi*psi'; 
        sigma{1} = MaxEntangled(2) * MaxEntangled(2)'; 
   
        [neg(i)] = NegMapEns(rho,2,sigma,2); 
        waitbar(i/n,h,sprintf('Progress: %5.3f %% , ctr = %d' ,i/n*100,i))
    end
end
close(h)
plot(p,neg(:,21),'LineWidth',2)
xlabel('p','FontSize',20)
ylabel('Negativity','FontSize',20)
ax = gca;
ax.XAxis.TickLabelFormat = '%,.1f';
set(gca, 'YTick',0:0.025:0.125,'FontSize',20)
ylim([0,0.125])
