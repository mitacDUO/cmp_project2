// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

double Fy(double x);	//declare functions outside main

void main()	//main code

{	
	printf("\nEULER METHOD\n\nx values	y values\n");	//printing titles for values displayed
	
	double h = 1;	//define step size
	int N = 5;	//define x(subscript n) 
	int n;		//the integer steps from x0 to x5
	
	//declare arrays
	double y[N];
	double x[N];
	
	//declare initial conditons for arrays
	y[0]=0;
	x[0]=0;

	//for loop for n=0,1,2,3,4
	for(n=0;n<N;n++)
	{
		x[n+1] = x[n]+h;			//declared how x(n+1) relates to x(n)
		y[n+1] = y[n]+h*(Fy(x[n]));		//declared how y(n+1) relates to y(n)
	printf("x=%f	y=%f\n",x[n+1],y[n+1]);		//printed values for x and y respectively
	}
}

//declared how y (as Fy) is a function of x
double Fy(double r)
{
	double row;
	
	double m = r*r*row;		//input function y=x^2=t 
	
	return m;		//return the input function t
}