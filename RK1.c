// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

#define N 5001

//declare functions outside main
double Fy(double x), f1(double a[N]), f2(double b[N]), f3(double c[N]), f4(double d[N]);

void main()	//main code

{	
	printf("\nRUNGE-KUTTA METHOD\n\nxvalues		yvalues\n");	//printing titles for values displayed
	
	double h = 0.001;	//define step size
	int n;		//the integer steps from x0 to x5
	
	//declare arrays
	double y[N];
	double x[N];
	
	//declare initial conditons for arrays
	y[0]=0;
	x[0]=0;

	//for loop for n=0,1,2,3,4
	for(n=0;n<N-1;n++)
	{
		x[n+1] = x[n]+h;			//declared how x(n+1) relates to x(n)

							//declared how y(n+1) relates to y(n)
		y[n+1] = y[n]+(h/6)*(Fy(x[n])+2*(Fy(x[n]+h/2)+(h/2)*Fy(x[n]))+2*(Fy(x[n]+h/2)+(h/2)*(Fy(x[n]+h/2)+(h/2)*Fy(x[n])))+(Fy(x[n]+h/2)+h*(Fy(x[n]+h/2)+(h/2)*(Fy(x[n]+h/2)+(h/2)*Fy(x[n])))));
		printf("x=%f	y=%f\n",x[n+1],y[n+1]);	//printed values for x and y respectively
	}
}

//declared how y (as Fy) is a function of x
double Fy(double x)
{
	double t = x*x;		//input function y=x^2=t 
	
	return t;		//return the input function t
}