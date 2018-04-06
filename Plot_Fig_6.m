% Code that generates plot in FIG.6
p = 0:0.025:0.5;
n = length(p);
I = eye(4);
neg = zeros(n,n);
M = cell(n,n);
ctr = 1;
h = waitbar(0,'Please wait...');

for i = 1:n
    for j = 1:n
        psi = sqrt(1-p(i)) * nn(0,0,2) + sqrt(p(i)) * nn(1,1,2);
        phi = sqrt(1-p(j)) * nn(0,0,2) + sqrt(p(j)) * nn(1,1,2);
        
        W = 1/4 *Tensor((phi*phi'), psi*psi' )  + 1/12 * Tensor(I-phi*phi',I-psi*psi' );
        neg(i,j) = Negativity(PermuteSystems(W,[1 3 2 4],[2 2 2 2]),4);

        waitbar(ctr/n^2,h,sprintf('Progress: %5.3f %% , ctr = %d' ,ctr/n^2*100,ctr))
        ctr = ctr+1;
    end
end
 close(h);  
 
 neg1 = zeros(1,n);
 neg2 = neg1;
 for i = 1:n
 neg1(i) = 0.25 * sqrt(p(i)*(1-p(i)));
 neg2(i) = 1/8 - 1/4*p(i);
 end

 hold off
 surf(p,p,neg);
 xlabel('p','FontSize',20)
 ylabel('q','FontSize',20) 
 zlabel('Negativity','FontSize',20)
 set(gca, 'XTick',0:0.1:0.5,'FontSize',17)
 set(gca, 'YTick',0:0.1:0.5,'FontSize',17)
 set(gca, 'ZTick',0:0.025:0.125,'FontSize',17)
 zlim([0,0.125])
 alpha 0.5;
 hold on
 plot3(p,0.5*ones(1,n),neg2,'r','LineWidth',2)
 plot3(0*ones(1,n),p,neg1,'r','LineWidth',2)