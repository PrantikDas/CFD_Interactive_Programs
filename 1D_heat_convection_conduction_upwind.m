%using diary function for printing the output in a text file
diary dynamictable.txt
%command for clearing the workspace and output and also any previous
%command window
clc,clear;
close;
%user input of grid point G and peclet number P
G=input('Enter the number of grid points \n');
P=input('Enter the number of Peclet Number \n');
%the value of interface diffusivity coefficient and length
T=.1;
L=1;
%calculation of deltaX
deltaX=(L/(G));
disp(deltaX);
%defining the counter variable
m=0;
%defining diffusion conductance and here P is the peclet number and F is
%the convective mass flux per unit area
D=(T/(deltaX));
F=P*D;
%defining the boundary conditions
J1=1;
J2=0;
%running the loop for calculation the matrix a , i is denoting the row and
%j is denoting the column
for i=2:1:G-1
    j=2+m;
    a(i,j-1)=-(D+F);
    a(i,j+1)=-(D);
    a(i,j)=-(a(i,j-1)+a(i,j+1));
    m=m+1;
end
a(1,1)=((3*D)+F);
a(1,2)=-(D);
a(G,G-1)=-(D+F);
a(G,G)=((3*D)+F);
%displaying the matrix
disp(a);
%calculation the value of x at various node points
x(1)=0;
x(2)=x(1)+deltaX/2;
for i=3:1:G+1
    x(i)=x(i-1)+deltaX;
    fprintf('%f \n',x(i));
end
x(G+2)=x(G+1)+deltaX/2;
%calculating the B matrix
B(1,1)=(2*D+F)*J1;
B(G,1)=(2*D)*J2;
%calculating the B matrix
disp(B);
%finding the value of the numerical solution
U=mldivide(a,B);
%displaying the value of the numerical solution and putting the boundary
%condition
disp(U);
W(1)=J1;
W(G+2)=J2;
fprintf('The numerical solution at x(%d) is %f \n',1,W(1));
for i=2:1:G+1
    W(i)=U(i-1);
    fprintf('The numerical solution at x(%d) is %f \n',i,W(i));
end
fprintf('The numerical solution at x(%d) is %f \n',G+2,W(G+2));
fprintf('\n');
disp(x);
%displaying the value of the analytical solution and calculating the node
%number
for i=1:1:G+2
    I(i)=1-(((exp((P/deltaX)*x(i)))-1)/((exp(P/deltaX)*L)-1));
    fprintf('The analytical solution at x(%d) is %f \n',i,I(i));
    N(i)=i;
end
fprintf('\n');
hold on;
%plotting the value of x vs numerical solution and analytical solution on
%the same curve
plot (x,W);
plot (x,I);
%putting the limit of the x curve
xlim([0 1]);
%labelling the x and y axis and also the graph itself
xlabel('deltaX (GRID SPACING)');
ylabel('Value of the function');
title('Transport Property Profile');
grid on;
%calculation of the percentage error by using O variable
O=((I-W)./I).*100;
legend ('Numerical Transport Property vs deltaX','Analytical Transport Property vs deltaX','Location','southeast','orientation','vertical');
hold off;
%constructing a dynamic table
T = table(N.',x.',W.',I.',(I-W).',O.','VariableNames',{'Node' 'Distance' 'Finite Volume Solution' 'Analytical Solution' 'Difference' 'Percentage Error'});
disp(T);
diary off;
