clear all
close all
clc
%Entering the number of grid points in x and y direction
Nx =  input('Number of grid points in x/y direction');
% code written for equal number of grid points for both sides
Ny = Nx;
%defining the boundary
x = 1:Nx;
dx = abs(x(1)-x(2));
nx = length(x);
y = 1:Ny;
dy = abs(y(1)-y(2));
ny = length(y);

%defining the boundary conditions
t_top = 600;
t_bottom = 900;
t_left = 400;
t_right = 800;

%defining the initial conditions and the 1st guess value of temperature.
t = 300*ones(nx,ny);
t(ny,:) = t_top;
t(1,:) = t_bottom;
t(:,1) = t_left;
t(:,ny) = t_right;
t(1,1) = (t_left + t_bottom)/2;
t(1,ny) = (t_top + t_left)/2;
t(nx,ny) = (t_top + t_right)/2;
t(nx,1) = (t_right + t_bottom)/2;

%initialising the t_initial and t_old for future computation.
t_initial = t;
t_old = t;

%defining the k term for future simplification.
k =  (dx^2*dy^2)/(2*(dx^2+dy^2));

%defining the tolerance and error values
tolerance = 1e-4;
error = 1;

%starting the iteration count.
iteration = 0;

%starting the while loop that uses Jacobi method and it runs until the error value is acceptable.
while error>tolerance

    %starting the spatial loops
    for i = 2:(length(x)-1)
        for j = 2:(length(y)-1)

            %calculating the temperature of a point using Jacobi method
            term1 = (t_old(i-1,j)+t_old(i+1,j))/(dx^2);
            term2 = (t_old(i,j-1)+t_old(i,j+1))/(dy^2);
            t(i,j) = k*(term1 + term2);
        end
    end

    %calculating the error value
    error = max(max(abs(t_old - t)));

    %updating the t_old
    t_old = t;

    %updating the iteration count.
    iteration = iteration + 1;

end
hold on
%ploting for the solution obtained using the Jacobi method.
figure(1)
[c,h] = contourf(x,y,t);
clabel(c,h);
colorbar
caxis([400 900])
colormap(jet(256));
%set(gca,'YDIR','reverse');
xlabel('X Axis');
ylabel('Y Axis');
title({['Plot for Jacobi method'],['Number of iterations = ',num2str(iteration)]});
figure(2)
surface(x,y,t)
surf(x,y,t)
hold off
