clc,clear; 
close all;
% Taken the initial guess of x and maximum allowable tolerance of Newton
% Raphson Method from user input
x0=input('Enter the initial guess of x :'); 
t=input('Enter the maximum allowable tolerance of Newton Raphson Method :'); 
% I am defining the function as given in the assignment
f = @(x) x^3-x-1; 
% I am defining the derivative of the function
f1 = @(x) 3*x^2-1; 
% for loop
for i=1:1:inf 
        %Calculation of new approximate root by Newton Raphson method
        x = x0-(f(x0)/f1(x0)); 
        %if condition for tolerance limit
        if(abs(x-x0)<=t) 
        break; 
        end 
        %relative error calculation
        e=abs((x-x0)/(x0)); 
        fprintf('The iteration number is %i \t The relative error is %f \t The value of x is %f \n' ,i,e,x); 
        %updating the value of x
        x0=x;   
        %storing iteration number inside an arrays
        D(i)=i; 
        %storing relative error inside an arrays
        E(i)=e; 
end 
%plotting of relative error (y axis) vs iteration step (x axis)
plot(D,E); 
%labeling the axis and title of the graph plot
xlabel('Iteration Steps'); 
ylabel('Relative Error'); 
title('Convergence History'); 
grid on ; 
%printing the root of the equation
fprintf('The root of the equation is %f \n',x); 
