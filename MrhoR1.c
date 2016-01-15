// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

#define N 200

//declare functions outside main
double rho(double x), M(double a), gamma(double n);

void main()	//main code

{	
	printf("\nRUNGE-KUTTA METHOD\n\nxvalues		yvalues		zvalues\n");	//printing titles for values displayed
	
	double s = 0.01;	//define step size
	double h = 0.01;	//the initail parameter h
	double c = 10.0;	//the parameter row(0)
	int n;			//the integer steps from x0 to x5
	
	//declare arrays
	double r[N];
	double y[N];
	double z[N];
	
	//declare initial conditons for arrays
	r[0] = h;				//first array is for r=h
	y[0] = (h*h*h/3)*c;			//initial conditon for scaled mass (m)
	z[0] = c*(1-((h*h*c)/(6*gamma(c))));	//initial conditon for rho
	
	//for loop for n=0,1,...,200
	for(n=0;n<N-1;n++)
	{
		double f1 = M(r[n]);				//added sustitutes for f1,f2,f3
		double f2 = M(r[n]+s/2)+(s/2)*f1;
		double f3 = M(r[n]+s/2)+(s/2)*f2;
		double f4 = M(r[n]+s/2)+s*f3;
		
		double g1 = rho(r[n]);				//added sustitutes for g1,g2,g3
		double g2 = rho(r[n]+s/2)+(s/2)*g1;
		double g3 = rho(r[n]+s/2)+(s/2)*g2;
		double g4 = rho(r[n]+s/2)+s*g3;
								//declared how y(n+1) relates to y(n)
		r[n+1] = r[n]+s;				//declared how x(n+1) relates to x(n)
		y[n+1] = y[n]+(s/6)*(f1+2*f2+2*f3+f4);
		z[n+1] = z[n]+(s/6)*(g1+2*g2+2*g3+g4);
		
		printf("r=%.2f		M=%.2f		rho=%.2f\n",r[n+1],y[n+1],z[n+1]);		//printed values for x and y respectively
	}
}

//declared the function for dm/dr
double M(double a)
{
	double s = 0.01;			//define step size
	double h = 0.01;			//the initail parameter h
	double c = 10.0;			//the parameter rho(c)

	double b,p;
	p = c*(1-(h*h*c/(6*gamma(c))));		//initial conditon for rho
	
	b = a*a*p;				//coupled equation for dm/dr
	
	return b;				//return dm/dr
}

//declared the function for d(rho)/dr
double rho(double x)
{
	double s = 0.01;			//define step size
	double h = 0.01;			//the initail parameter h
	double c = 10.0;			//the parameter rho(c)
	double m,p,t;
	
	m = ((h*h*h)/3)*c;			//initial conditon for scaled mass (m)
	p = c*(1-(h*h*c/(6*gamma(c))));		//initial condition for rho
	
	t = -(m*p)/(gamma(p)*x*x);		//coupled equation for d(rho)/dr
	
	return t;				//return d(rho)/dr
}

//declared the function for gamma(x) = x^2/(3*sqrt(1+x^2))
double gamma(double d)
{
	double x = ((d*d)/(3*(sqrt(1+d*d))));
	
	return x;
}