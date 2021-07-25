/*Compiled by using GCC compiler in Linux*/
#include <stdio.h>
#include <math.h>
#include <time.h>
int main()
{
 clock_t t; 
 t = clock(); 
 int n,i,j,m,k;
 float L,k_c,A;
 float dx;
 printf("Enter the number of discretised points on the rod excluding boundary nodes\n");
 scanf ("%d",&n);
 float x[n];
 float a[1000][1000]={0},T[1000]={0};
 float Ta,Tb;
 float pivot = 0.0;
 float factor = 0.0;
 float sum = 0.0;
 printf("Enter the length of the rod \n");
 scanf ("%f",&L);
 printf("Enter the conductivity of the rod \n");
 scanf ("%f",&k_c);
 printf("Enter the cross sectonal area of the rod \n");
 scanf ("%f",&A);
 printf("Enter the temperature of left hand side of the rod \n");
 scanf ("%f",&Ta);
 printf("Enter the temperature of left hand side of the rod \n");
 scanf ("%f",&Tb);
 dx = L*pow(n,-1);
 x[0]=0;
 x[1]=x[0]+dx*0.5;
 for(i=2;i<=n;i++)
 {
 x[i]=x[i-1]+dx;
 }
 x[n+1]=x[n]+dx*0.5;
 a[1][1]=3*k_c*A*pow(dx,-1);
 a[1][2]=-1*k_c*A*pow(dx,-1);
 a[1][n+1]=2*Ta;
 for(i=2;i<=n-1;i++)
 {
    for(j=2;j<=n+1;j++)
    {
    a[j][j]=2*k_c*A*pow(dx,-1);
    a[j][j-1]=-k_c*A*pow(dx,-1);
    a[j][j+1]=-k_c*A*pow(dx,-1);
    }
 }
 a[n][n]=3*k_c*A*pow(dx,-1);
 a[n][n-1]=-1*k_c*A*pow(dx,-1);
 a[1][n+1]=2*Ta*k_c*A*pow(dx,-1);
 a[n][n+1]=2*Tb*k_c*A*pow(dx,-1);
 for(i=1;i<=n;i++)
{
    for(j=1;j<=n+1;j++)
    {
    printf("%f\t",a[i][j]);
    }
    printf("\n\n");
}
printf("Solution by Simple Gauss Elimination \n\n");
for(k=1;k<=n-1;k++)
{
if(a[k][k]==0.0)
{
 printf("error");

}
 else
{
  pivot = a[k][k];
  for(j=k;j<=n+1;j++)
    {
    a[k][j]= a[k][j]/pivot;
    } 
    for(i=k+1;i<=n;i++)
    {
    factor = a[i][k];
    for(j = k;j<=n+1;j++)
    {
     a[i][j] = a[i][j] - factor * a[k][j];
    }
    } 
} 
if(a[n][n]==0)
{
printf("error");
}
else
{
 T[n] = a[n][n+1]/a[n][n];
 for(i=n-1;i>=1;i--)
    {
    sum = 0.0;
    for(j=i+1;j<=n;j++)
    {
    sum = sum + a[i][j] * T[j];
    T[i]= (a[i][n+1]-sum)/a[i][i];
    }
    }
}
} 
for(i=1;i<=n;i++)
{
    printf("\n\tT[%1d]=%10.4f",i,T[i]); 
} 
    T[0]=Ta;
    T[n+1]=Tb;
    FILE *gnuplot = fopen("gnuplot.txt", "w");
    for (i = 0; i <= n+1; i++)
    fprintf(gnuplot, "%g %g\n", x[i], T[i]);
    fflush(gnuplot);
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
    printf("\n The program took %f seconds to execute \n", time_taken); 
return 0;
}









