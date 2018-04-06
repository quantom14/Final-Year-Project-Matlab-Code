%% Script for Figure 12.

p = 0:0.025:0.5;
m = length(p);
neg = zeros(5,m);
n = 1:5;

for i = n
    for j = 1:m
        neg(i,j) = 1/2^i * sqrt(p(j)-p(j)^2);
    end
end

surf(p,n,neg)
xlabel('p')
ylabel('n')
zlabel('Negativity')
set(gca, 'XTick',0:0.1:0.5,'FontSize',20)
set(gca, 'YTick',n,'FontSize',20)