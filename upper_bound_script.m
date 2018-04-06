%% Script to generate Figure 10 using previous results

p = 0:0.025:0.5;
n = length(p);
neg2 = zeros(n,n);
for i = 1:n
    for j = 1:n        
        
        r0 = min((1-p(j))/(1-p(i)),1);
        chi = sqrt(r0) * nn(0,0,2) + sqrt(1-r0) * nn(1,1,2);        
        neg2(i,j) = Negativity(chi*chi');      
    end
end
[X,Y] = meshgrid(0:0.025:0.5);
Z = X+Y;
%neg1 is the matrix of negativities obtained in the script psitophi.
s1 = surf(p,p,neg1,X);
hold on
s2 = surf(p,p,neg2,-Z);
xlabel('p','FontSize',20)
ylabel('q','FontSize',20)
zlabel('Negativity','FontSize',20)
set(gca, 'XTick',0:0.1:0.5,'YTick',0:0.1:0.5,'FontSize',18)
hold on
plot3(p,0*ones(1,n),neg2(1,:),'r','LineWidth',2)
