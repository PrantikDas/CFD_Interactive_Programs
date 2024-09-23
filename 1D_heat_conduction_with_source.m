clc, clear;
close;
G = input ('Enter the number of grid points');
L=input('Enter the length of pipe');
deltaX = (L)/G; # in m
A = 1; # in m2
k = 0.5; #W/mK
q= 1000;
V = (k*A)/(deltaX);
S=q*A*deltaX;
m = 0; # a counter
# Boundary conditions
T1 = 100;  # left boundary temperature
T2 = 500; # right boundary temperature
for i = 2:1:G-1
  j = 2 + m;
  a(i,j-1) = -V;
  a(i,j+1) = -V;
  a(i,j) = 2*V;
  m = m + 1;
end
a(1,1) = 3*V;
a(1,2) = -V;
a(G,G-1) = -V;
a(G,G) = 3*V;
disp(a);
x(1) = 0;
x(2) = x(1) + deltaX*0.5;
for i = 3:1:G+1
  x(i) = x(i-1) + deltaX;
end
x(G+2) = x(G+1) + deltaX*0.5;
B = zeros(G,1);
B(1,1) = 2*T1*V+S;
for i=2:1:G-1
B(i,1) = S;
end
B(G,1) = 2*T2*V+S;
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
