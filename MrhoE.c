// declare libraries
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>

#define N 200

double gamma(double p), M(double r, double m, double p), rho(double r, double m, double p);	//declare functions outside main

double h = 0.01;	//define step size

void main()	//main code

{	
	printf("\nEULER METHOD\n\nx values	y values\n");	//printing titles for values displayed
	
	double c = 10.0; 	//the parameter rho(c)
	int n;			//the integer steps from x0 to x2
	
	//declare arrays
	float y[N];
	float x[N];
	float z[N];
	
	//r,m,p as the radius, mass and density
	double r,m,p;
	
	//declare initial conditons for arrays
	x[0] = h;				//first array is for r=h
	y[0] = (h*h*h/3)*c;			//initial conditon for scaled mass (m)
	z[0] = c*(1-((h*h*c)/(6*gamma(c))));	//initial conditon for rho

	//for loop for n=0,1,...,200
	for(n=0;n<N;n++)
	{
		x[n+1] = x[n]+h;			//declared how x(n+1) relates to x(n)
		y[n+1] = y[n]+h*(M(x[n],y[n],z[n]));		//declared how y(n+1) relates to y(n)
		z[n+1] = z[n]+h*(rho(x[n],y[n],z[n]));
		
		if(z[n+1]<0)
		{
			break;
		}
		
		r = x[n+1];
		m = y[n+1];
		p = z[n+1];
		
	printf("r=%.2e	M=%.2e	rho=%.2e\n",x[n+1],y[n+1],z[n+1]);		//printed values for x and y respectively
	}
	
	/**** Use pgplot to plot this function ****/

  // cpgbeg starts a plotting page, in this case with 2x1 panels
  cpgbeg(0,"?",2,1);

  // sets colour: 1-black, 2-red, 3-green, 4-blue
  cpgsci(1);

  // sets line style 1-solid, 2-dashed, 3-dot-dashed, 4-dotted
  cpgsls(1);

  // sets charachter height, larger number = bigger
  cpgsch(1.);

  // cpgpage() moves to the next panel
  cpgpage();

  // sets the axes limits in the panel
  cpgswin(0,r,0,m);

  // draw the axes
  cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

  // label the bottom axis
  cpgmtxt("B",2.,.5,.5,"radius");
  
  // label the left axis
  cpgmtxt("L",2.5,.5,.5,"scaled mass");

  // connect N points in ax and ay with a line
  cpgline(n,x,y);

  // cpgpage() moves to the next panel
  cpgpage();

  // sets the axes limits in the panel
  cpgswin(0,r,0,c);

  // draw the axes
  cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

  // label the bottom axis
  cpgmtxt("B",2.,.5,.5,"radius");
  
  // label the left axis
  cpgmtxt("L",2.5,.5,.5,"density");
  
  // connect N points in ax and ay with a line
  cpgline(n,x,z);

  // close all pgplot applications
  cpgend();
  
  // end program
  return;
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
	double t = -(m*p)/(gamma(p)*r*r);	//coupled equation for d(rho)/dr
	
	return t;				//return d(rho)/dr
}