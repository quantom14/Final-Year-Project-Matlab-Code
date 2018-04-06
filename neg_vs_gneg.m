%% 3-qubit case
p = 0:0.02:1;
m = length(p);
dim = [2 2 2];
GHZ = GHZState(2,3)*GHZState(2,3)';
zero = [1 0]';
psi = Tensor(zero*zero' , MaxEntangled(2)*MaxEntangled(2)');
neg1 = zeros(m,2);
for i = 1:m
    phi = (1-p(i)) * psi + p(i) *GHZ;
    neg1(i,1) = Multi_Negativity(phi,dim);
    neg1(i,2) = Genuine_Negativity(phi,dim);
end
%% Plot
plot(p,neg1(:,1),'LineWidth',1)
hold on
plot(p,neg1(:,2),'LineWidth',1)
title('Negativities of 3 qubits','FontSize',20)
xlabel('p','FontSize',20)
ylabel('Negativity','FontSize',20) 
set(gca, 'XTick',0:0.1:1,'FontSize',20)
set(gca, 'YTick',0:0.1:0.5,'FontSize',20)
ylim([0,0.5])