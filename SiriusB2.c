// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

#define N 20000

//declared functions outside main
double gamma(double p), M(double r, double m, double p), rho(double r, double m, double p);
double M2(double r, double m, double p), M3(double r, double m, double p), M4(double r, double m, double p);
double rho2(double r, double m, double p), rho3(double r, double m, double p), rho4(double r, double m, double p);

double h = 0.01;	//the initail parameter h, the step size
double Y = 0.464;		//Y(e) the number of electrons per nucleon

void main()	//main code

{	
	double Rs = Y*((7.72*pow(10,8))/(6.95*pow(10,10)));
	double Ms = Y*Y*((5.67*pow(10,33))/(1.98*pow(10,33)));
	
	double c = 82;	//the parameter rho(c)
	
	printf("\nSIRIUS B\n");
	printf("\nRunge-Kutta Method\nMs (solar masses)	Rs (solar radii)	Rho(c) (solar densities)	Y(e)\n");	//printing titles for values displayed
	
	int n;		//the integer steps, n
	
	//declare arrays
	float x[N];
	float y[N];
	float z[N];
	float x2[N];
	float y2[N];
	float z2[N];
	
	//r,m,p as the radius, mass and density
	double r,m,p,r2,m2,p2;
		
		//declare initial conditons for arrays
		x[0] = h;				//first array is for r=h
		y[0] = (h*h*h/3)*c;			//initial conditon for scaled mass (m)
		z[0] = c*(1-((h*h*c)/(6*gamma(c))));	//initial conditon for rho
	
		//for loop for n=0,1,...,N
		for(n=0;n<N;n++)
		{
			//declared how x(n+1) relates to x(n), y(n+1) relates to y(n), z(n+1) relates to z(n)
			x[n+1] = x[n]+h;
			y[n+1] = y[n]+(h/6)*(M(x[n],y[n],z[n])+2*M2(x[n],y[n],z[n])+2*M3(x[n],y[n],z[n])+M4(x[n],y[n],z[n]));
			z[n+1] = z[n]+(h/6)*(rho(x[n],y[n],z[n])+2*rho2(x[n],y[n],z[n])+2*rho3(x[n],y[n],z[n])+rho4(x[n],y[n],z[n]));

			if(isnan(z[n+1]))
			{
				break;
			}
		
			//r,m,p will be declared in pg-plot
			r = x[n+1];
			m = y[n+1];
			p = z[n+1];
				
		}
	
		printf("%.3e		%.3e		%.2f				%.3f\n",y[n]*Ms,x[n]*Rs,c,Y);	//printed values for x and y respectively
	
	printf("\nEuler Method\nMs (solar masses)	Rs (solar radii)	Rho(c) (solar densities)	Y(e)\n");	//printing titles for values displayed
		
		//declare initial conditons for arrays
		x2[0] = h;				//first array is for r=h
		y2[0] = (h*h*h/3)*c;			//initial conditon for scaled mass (m)
		z2[0] = c*(1-((h*h*c)/(6*gamma(c))));	//initial conditon for rho
	
		//for loop for n=0,1,...,200
		for(n=0;n<N;n++)
		{
			//declared how x2(n+1) relates to x2(n), y2(n+1) relates to y2(n), z2(n+1) relates to z2(n)
			x2[n+1] = x2[n]+h;
			y2[n+1] = y2[n]+h*(M(x2[n],y2[n],z2[n]));
			z2[n+1] = z2[n]+h*(rho(x2[n],y2[n],z2[n]));
		
			if(z2[n+1]<0)
			{
				break;
			}
		
			//r2,m2,p2 will be declared in pg-plot
			r2 = x2[n+1];
			m2 = y2[n+1];
			p2 = z2[n+1];
		
		}
	
	//printed values for x and y respectively
	printf("%.3e		%.3e		%.2f				%.3f\n",y2[n]*Ms,x2[n]*Rs,c,Y);
}

//declared the function for gamma(x) = x^2/(3*sqrt(1+x^2))
double gamma(double p)
{
	double e = cbrt(p);
	double x = ((e*e)/(3*(sqrt(1+e*e))));
	
	return x;
}

//declared the function for dm/dr
double M(double r, double m, double p)
{
	double b = r*r*p;		//coupled equation for dm/dr
	
	return b;			//return dm/dr
}

//declared the function for d(rho)/dr
double rho(double r, double m, double p)
{	
	double t = -((m*p)/(gamma(p)*r*r));	//coupled equation for d(rho)/dr
	
	return t;				//return d(rho)/dr
}

//declared the function for M2 from M
double M2(double r, double m, double p)
{
	double Fn = M(r+h/2,m+(h/2)*M(r,m,p),p+(h/2)*rho(r,m,p));
	
	return Fn;
}

//declared the function for rho2 from rho
double rho2(double r, double m, double p)
{	
	double Fn = rho(r+h/2,m+(h/2)*M(r,m,p),p+(h/2)*rho(r,m,p));
	
	return Fn;
}

//declared the function for M3 from M2
double M3(double r, double m, double p)
{
	double Fn = M(r+h/2,m+(h/2)*M2(r,m,p),p+(h/2)*rho2(r,m,p));
	
	return Fn;
}

//declared the function for rho3 from rho2
double rho3(double r, double m, double p)
{	
	double Fn = rho(r+h/2,m+(h/2)*M2(r,m,p),p+(h/2)*rho2(r,m,p));
	
	return Fn;
}

//declared the function for M4 from M3
double M4(double r, double m, double p)
{
	double Fn = M(r+h,m+h*M3(r,m,p),p+h*rho3(r,m,p));
	
	return Fn;
}

//declared the function for rho4 from rho3
double rho4(double r, double m, double p)
{	
	double Fn = rho(r+h,m+h*M3(r,m,p),p+h*rho3(r,m,p));
	
	return Fn;
}