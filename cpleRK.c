// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

#define N 5001

//declare functions outside main
double F(double x, double y, double z), G(double x, double y, double z);
double F2(double x, double y, double z), F3(double x, double y, double z), F4(double x, double y, double z);
double G2(double x, double y, double z), G3(double x, double y, double z), G4(double x, double y, double z);

double h = 0.001;	//define step size

void main()	//main code

{	
	printf("\nRUNGE-KUTTA METHOD\n\nxvalues		yvalues		zvalues\n");	//printing titles for values displayed
	
	int n;		//the integer steps from x0 to x5
	
	//declare arrays
	double x[N];
	double y[N];
	double z[N];
	
	//declare initial conditons for arrays
	x[0] = 0;
	y[0] = 0;
	z[0] = 0;
	
	//for loop for n=0,1,2,3,4
	for(n=0;n<N-1;n++)
	{
		//declared how x(n+1) relates to x(n), y(n+1) relates to y(n), z(n+1) relates to z(n)
		x[n+1] = x[n]+h;
		y[n+1] = y[n]+(h/6)*(F(x[n],y[n],z[n])+2*F2(x[n],y[n],z[n])+2*F3(x[n],y[n],z[n])+F4(x[n],y[n],z[n]));
		z[n+1] = z[n]+(h/6)*(G(x[n],y[n],z[n])+2*G2(x[n],y[n],z[n])+2*G3(x[n],y[n],z[n])+G4(x[n],y[n],z[n]));
		
		printf("x=%f	y=%f	z=%f\n",x[n+1],y[n+1],z[n+1]);		//printed values for x and y respectively
	}
}

//declared how y (as Fy) is a function of x
double F(double x, double y, double z)
{
	double t = x*x;		//input function y=x^2=t 
	
	return t;		//return the input function t
}

double G(double x, double y, double z)
{
	double b = pow(x,4.);
	
	return b;
}

//declared the function for M2 from M
double F2(double x, double y, double z)
{
	double Fn = F(x+h/2,y+(h/2)*F(x,y,z),z+(h/2)*G(x,y,z));
	
	return Fn;
}

//declared the function for rho2 from rho
double G2(double x, double y, double z)
{	
	double Fn = G(x+h/2,y+(h/2)*F(x,y,z),z+(h/2)*G(z,y,z));
	
	return Fn;
}

//declared the function for M3 from M2
double F3(double x, double y, double z)
{
	double Fn = F(x+h/2,y+(h/2)*F2(x,y,z),z+(h/2)*G2(x,y,z));
	
	return Fn;
}

//declared the function for rho3 from rho2
double G3(double x, double y, double z)
{	
	double Fn = G(x+h/2,y+(h/2)*F2(x,y,z),z+(h/2)*G2(z,y,z));
	
	return Fn;
}

//declared the function for M4 from M3
double F4(double x, double y, double z)
{
	double Fn = F(x+h,y+h*F3(x,y,z),z+h*G3(x,y,z));
	
	return Fn;
}

//declared the function for rho4 from rho3
double G4(double x, double y, double z)
{	
	double Fn = G(x+h/2,y+h*F3(x,y,z),z+h*G3(z,y,z));
	
	return Fn;
}