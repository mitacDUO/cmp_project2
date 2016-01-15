/*
 *  Author:  A. Phan, 2010
 *    modified from S.Wyithe
 *
 *  This code demonstrates some simple comands,
 *  functions and plotting in pgplot
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cpgplot.h>


// Defined constants
#define N 100  // number of loop iterations


// Declare functions
float fy(float x);


// Main code
int main()
{
  /**** Count to 10 in integers ****/

  // Declare integer for loop counting
  int i;

  // Loop from 0 to 10, printing at each step
  for(i=0; i<10; i++)
  {
    printf("i= %d\n", i);
  }

  /**** Plot a function y=f(x) ****/

  // Declare arrays of N real numbers
  float ax[N]; // x
  float ay[N]; // y

  float aylow[N];  // lower error in y
  float ayhigh[N]; // upper error in y

  // Set minimum and maximum for x
  float xmin = 0.0;
  float xmax = 10.0;

  // Assigning ax with N values for x between xmin and xmax
  for(i=0;i<N;i++)
  {
    ax[i] = xmin + (xmax-xmin)*(float)i/(float)(N-1);
  }

  // Fill ay using function fy
  for(i=0;i<N;i++)
  {
    ay[i] = fy(ax[i]);
  }

  // Fill aylow and ayhigh using sqrt(y) as the error
  for(i=0;i<N;i++)
  {
    aylow[i]  = ay[i] - sqrt(ay[i]);
    ayhigh[i] = ay[i] + sqrt(ay[i]);
  }

  /**** Use pgplot to plot this function ****/

  // cpgbeg starts a plotting page, in this case with 2x1 panels
  cpgbeg(0,"?",2,1);

  // sets colour: 1-black, 2-red, 3-green, 4-blue
  cpgsci(1);

  // sets line style: 1-solid, 2-dashed, 3-dot-dashed, 4-dotted
  cpgsls(1);

  // sets charachter height, larger number = bigger
  cpgsch(2.);

  // cpgpage() moves to the next panel
  cpgpage();

  // sets the axes limits in the panel
  cpgswin(xmin,xmax,0.,100.);

  // draw the axes
  cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

  // label the bottom axis
  cpgmtxt("B",2.,.5,.5,"x");
  // label the left axis
  cpgmtxt("L",2.5,.5,.5,"f");

  // connect N points in ax and ay with a line
  cpgline(N,ax,ay);

  // cpgpage() moves to the next panel
  cpgpage();

  // sets the axes limits in the panel
  cpgswin(xmin,xmax,0.,100.);

  // draw the axes
  cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

  // label the bottom axis
  cpgmtxt("B",2.,.5,.5,"x");
  // label the left axis
  cpgmtxt("L",2.5,.5,.5,"f");

  // draw N points in ax and ay
  //   17 - filled circles, 16 - filled squares, 13 - filled triangles
  cpgpt(N,ax,ay,17);

  // draw y error bars on the points
  cpgerry(N,ax,aylow,ayhigh,1.0);

  // close all pgplot applications
  cpgend();

  // end program
  return 0;
}


// Returns the function y=x^2
float fy(float x)
{
  float val;

  val = x*x;

  return val;
}
