close all,clear,clc;
%linear equation
%very large number for G for iteration 
G=zeros(1,4);
%first index is denoting variable listing and second index the max number
%of iteration possible
for i=1:1:4
    G(1,i)=input(sprintf('Enter the Initial guess of the variable X(%d) \n',i));
end
tol=input('Enter the tolerance limit \n');
l=1;
error=1;
itr=0;
res=zeros(1,4);
while error>tol
    l=l+1;
    for i=l
    %gauss seidel equations
    G(i,1)=(9-G(i-1,2)*G(i-1,3)-G(i-1,4))/2;
    G(i,2)=(14-G(i,1)^2-G(i-1,3))/4;
    G(i,3)=(8-G(i,2)^2-G(i-1,4))/5;
    G(i,4)=(9-G(i,2)*G(i,3))/4;
    %residue finding equations (normalised)
    res(i-1,1)=(9-2*G(i,1)-G(i,2)*G(i,3)-G(i,4))/9;
    res(i-1,2)=(14-G(i,1)^2-4*G(i,2)-G(i,3))/14;
    res(i-1,3)=(8-G(i,2)^2-5*G(i,3)-G(i,4))/8;
    res(i-1,4)=(9-G(i,2)*G(i,3)-4*G(i,4))/9;
    [x,y]=max(G(i,:));
    error = 0;
    end
    error = error+abs((G(i,y)-G(i-1,y))/(G(i-1,y)));
    itr=itr+1;
    L(i-1)=itr;
end
%Graph Essentials
xlabel('Iteration Number');
ylabel('Residue');
title('Gauss Siedel Method');
for i=1:1:4
fprintf('The converged solution of X(%d) is %f \n',i,G(itr+1,i));
end
fprintf('The number of iterations required for convergence is %d \n',itr);
hold on;
grid on;
%residue plot for each and every equation
j=1;
for j=1:1:4
    plot(L,res(:,j));    
end
legend('1st equation','2nd equation','3rd equation','4th equation','Location','southeast','orientation','vertical');
hold off;