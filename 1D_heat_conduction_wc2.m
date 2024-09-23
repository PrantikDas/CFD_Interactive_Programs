clc, clear;
close;
G = input ('Enter the number of grid points');
deltaX = (0.5)/G; # in m
A = 0.01; # in m2
k = 0.5; #W/mK
V = (k)/(deltaX);
q = 0; %volumetric basis
Su = q*deltaX*A;
U = 0.1;
rho = 1;
F = rho * U;
m = 0; # a counter
# Boundary conditions
T1 = 1;
T2 = 0;
for i = 2:1:G-1
  j = 2 + m;
  a(i,j-1) = -(V+F);
  a(i,j+1) = -V;
  a(i,j) = 2*V+F;
  m = m + 1;
end
a(1,1) = 3V+F;
a(1,2) = -(F);
a(G,G-1) = -(V+F);
a(G,G) = 3V+F;
disp(a);
x(1) = 0;
x(2) = x(1) + deltaX*0.5;
for i = 3:1:G+1
  x(i) = x(i-1) + deltaX;
end
x(G+2) = x(G+1) + deltaX*0.5;
B = zeros(G,1);
B(1,1) = (2*T1+F)*V+Su;
for i=2:1:G-1
B(i,1) = Su;
end
B(G,1) = 2*T2*V+Su;
T = mldivide(a,B);
disp(T);
disp(x');
Tmod = zeros(G+2,1)
Tmod(1,1) = T1;
Tmod(G+2,1) = T2;
for i=2:1:G+1
   Tmod(i,1) = T(i-1,1);
end
plot(x',Tmod)
