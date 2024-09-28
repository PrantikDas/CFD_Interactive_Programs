close all, clear ,clc;
Nx=input('Enter the number of points along x and y axis \n');
Ny=Nx;
Lx=1;
Ly=1;
L=(Lx*Ly)^.5;
u=2;
v=2;
p_w=100;
p_n=100;
p_e=0;
p_s=0;
A=zeros(Nx*Nx,Nx*Nx);
B=zeros(Nx*Ny,1);
WE = [0 0 -3/8 (-2/8-1) (2/8+1)*p_w
      0 (7/8+1/8) (-3/8) (1/4) -(1/4)*p_w
      -1/8 6/8 0 1 -p_e
      -1/8 (6/8+1/8) -3/8 0 0];
SN = [0 0 -3/8 (-2/8-1) (2/8+1)*p_s
      0 (7/8+1/8) (-3/8) 1/4 -(1/4)*p_s
      -1/8 6/8 0 1 -p_n
      -1/8 (6/8+1/8) -3/8 0 0]; 
for i=1:Ny
    for j=1:Nx
        a_n=(i-1)*Nx+j;
        switch(j)
            case(1)
                aWW=WE(1,1);aW=WE(1,2);aE=WE(1,3);Spwe=WE(1,4);Suwe=WE(1,5);
                A(a_n,a_n+1)=-aE;
            case(2)
                aWW=WE(2,1);aW=WE(2,2);aE=WE(2,3);Spwe=WE(2,4);Suwe=WE(2,5);
                A(a_n,a_n-1)=-aW;A(a_n,a_n+1)=-aE;
            case (Nx)
                aWW=WE(3,1);aW=WE(3,2);aE=WE(3,3);Spwe=WE(3,4);Suwe=WE(3,5);
                A(a_n,a_n-2)=-aWW;A(a_n,a_n-1)=-aW;
            otherwise
                aWW=WE(4,1);aW=WE(4,2);aE=WE(4,3);Spwe=WE(4,4);Suwe=WE(4,5);
                A(a_n,a_n-2)=-aWW;A(a_n,a_n-1)=-aW;A(a_n,a_n+1)=-aE;
        end
        switch(i)
            case(Ny)
                aSS=SN(1,1);aS=SN(1,2);aN=SN(1,3);Spsn=SN(1,4);Susn=SN(1,5);
                A(a_n,a_n-Nx)=-aN;
            case(Ny-1)
                aSS=SN(2,1);aS=SN(2,2);aN=SN(2,3);Spsn=SN(2,4);Susn=SN(2,5);
                A(a_n,a_n+Nx)=-aS;A(a_n,a_n-Nx)=-aN;
            case(1)
                aSS=SN(3,1);aS=SN(3,2);aN=SN(3,3);Spsn=SN(3,4);Susn=SN(3,5);
                A(a_n,a_n+2*Nx)=-aSS;A(a_n,a_n+Nx)=-aS;
            otherwise
                aSS=SN(4,1);aS=SN(4,2);aN=SN(4,3);Spsn=SN(4,4);Susn=SN(4,5);
                A(a_n,a_n+2*Nx)=-aSS;A(a_n,a_n+Nx)=-aS;A(a_n,a_n-Nx)=-aN;
        end
        B(a_n,1)=Suwe+Susn;
        aP=aWW+aW+aE-Spwe+aSS+aS+aN-Spsn;
        A(a_n,a_n)=aP;
    end
end
T=mldivide(A,B);
hold on;
%Graph Essentials
xlabel('Distance along diagonal');
ylabel('Value of Transport Property');
title('QUICK False Diffusion Analysis');
%Analytical solution
z = linspace (0,L*1.4,1000);
for k=1:1:max(size(z))
if z(k)<(L*.7075)
    TA(k)=100;
else
    TA(k)=0;
end
end
xticks([0 .2 .4 .6 .8 1.0 1.2 1.4]);
%diagonalisation of temperature
p = linspace (0,L*1.4,Nx+2);
m=1;
plot(z,TA);
for i=1:Nx:Nx*Ny-Nx+1
Td(m+1)=T(i+m-1);
m=m+1;
end
Td(1)=100;
Td(Nx+2)=0;
plot (p,Td);
hold off
legend('Analytical Solution','Numerical Solution','Location','southeast','orientation','vertical');
%Graph Essentials
xlim([0 1.4]);
ylim([-10 110]);
xlabel('Distance along diagonal');
ylabel('Value of Transport Property');
title('QUICK False Diffusion Analysis');    