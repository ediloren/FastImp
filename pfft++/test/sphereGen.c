/*
Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA.
All rights reserved.

This Agreement gives you, the LICENSEE, certain rights and obligations.
By using the software, you indicate that you have read, understood, and
will comply with the terms.

Permission to use, copy and modify for internal, noncommercial purposes
is hereby granted.  Any distribution of this program or any part thereof
is strictly prohibited without prior written consent of M.I.T.

Title to copyright to this software and to any associated documentation
shall at all times remain with M.I.T. and LICENSEE agrees to preserve
same.  LICENSEE agrees not to make any copies except for LICENSEE'S
internal noncommercial use, or to use separately any portion of this
software without prior written consent of M.I.T.  LICENSEE agrees to
place the appropriate copyright notice on any such copies.

Nothing in this Agreement shall be construed as conferring rights to use
in advertising, publicity or otherwise any trademark or the name of
"Massachusetts Institute of Technology" or "M.I.T."

M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By
way of example, but not limitation, M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS OR DOCUMENTATION WILL
NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
M.I.T. shall not be held liable for any liability nor for any direct,
indirect or consequential damages with respect to any claim by LICENSEE
or any third party on account of or arising from this Agreement or use
of this software.
*/
/* This program wil generate input files with constant-type 
 * elements (quads and triangles) for a simple case.
 * At the moment most of the stuff is hardcoded. 
 * Eventually, most should be reached through command line arguments */

# include <stdlib.h>

# include <stdio.h>
# include <float.h>
# include <math.h>
# include <string.h>
# include <assert.h>

#include "sphereGen.h"

#define PI 3.141592653589793
#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )

#define XPROD(V,V1,V2){ \
V[0] = V1[1]*V2[2]-V1[2]*V2[1]; \
V[1] = V1[2]*V2[0]-V1[0]*V2[2]; \
V[2] = V1[0]*V2[1]-V1[1]*V2[0];}

/**********************************************************************
 * sphereGen --
 **********************************************************************/

void sphereGen(
	       const char* filename,
	       double radius,
	       int nDense)
{
  int i,j, iside, ic, iplot, swapNormal, gradeMesh[6], igrade,
    *icc[2], icc0[3], icc1[3], ie, dostats, nn;
  double sideLength, center[3], thisrad;
  double dx,du, u[4],v[4], xout[4],yout[4],zout[4],
    xu,xv,x0, yu,yv,y0, zu,zv,z0, *uvals, *ugvals, *uuse,
    d02, d13,
    s0[3], s1[3], normal[3], area, 
    minarea, maxarea, meanarea;
  int neles, writeToFile;
  char inStr[]="inwards", outStr[]="outwards", *normalStr;
  FILE *outFile;

  /* Write to file "sphere.qui" or to stdout?*/
  writeToFile=1;

  /* Generate a cube with equal size mesh and outward pointing normal */
  center[0]=0.0;
  center[1]=0.0;
  center[2]=0.0;
  /* nDense by nDense elements on each side. */
  /*  nDense = 8; */
  /* Possibly swap normal to point inwards (for "exterior problems"): */
  swapNormal=1;

  /* Set to 0 if no statistics are needed 
   * (written as comments at end of file) */
  dostats = 1;

  /* ---- Below here everything should be OK ---*/

  /*if (writeToFile) outFile = fopen("sphere.qui","w");
  else outFile = stdout;
*/

  if (writeToFile) outFile = fopen(filename, "w");
  else outFile = stdout;

  /* Possibly grade the mesh using a cosine scaling towards the edges.
   * (This is a VERY bad idea for the sphere. - KEEP igrade=0)  */
  igrade = 0;
  gradeMesh[0] = igrade; /* normal = ( 1, 0, 0);*/
  gradeMesh[1] = igrade; /* normal = (-1, 0, 0);*/
  gradeMesh[2] = igrade; /* normal = ( 0, 1, 0);*/
  gradeMesh[3] = igrade; /* normal = ( 0,-1, 0);*/
  gradeMesh[4] = igrade; /* normal = ( 0, 0, 1); "free surface"*/
  gradeMesh[5] = igrade; /* normal = ( 0, 0,-1);*/

  /* DO NOT change this length! (sideLength SHOULD be TWO!) */
  sideLength = 2.0;

  /* Mesh size (if uniform) */
  dx = sideLength/(double) nDense;


  /* Initiate statistics settings */
  minarea=DBL_MAX;
  maxarea=0.0;
  meanarea=0.0;
  neles=0;
  /* Setup number of different coordinate values to be used, 
   * and set those values */
  uvals = (double *) calloc(nDense+1,sizeof(double));
  ugvals = (double *) calloc(nDense+1,sizeof(double));
  du = 2.0/(double) nDense;
  /* Points located from [-1:1] */
  uvals[0] = -1.0;
  for (i=1; i<nDense; i++)uvals[i] = i*du - 1.0;
  uvals[nDense] = 1.0;

  /* Cosine scaling of the vector */
  for (i=0; i<nDense+1; i++)
    ugvals[i] = cos(0.5*PI*(uvals[i]-1.0));
  /* Scale mesh with sidelength */
  for (i=0; i<nDense+1; i++) {
    uvals[i] *= sideLength/2.0;
    ugvals[i] *= sideLength/2.0;
  }

  /* Scale mesh to make approximately uniform grid on a sphere         */
  for (i=0; i<nDense+1; i++) {
    uvals[i]  = tan(PI*uvals[i]/4.0);
    ugvals[i] = tan(PI*ugvals[i]/4.0);
  }

  /* Expected number of elements: */
  nn = 2*6*nDense*nDense;

  normalStr = outStr;
  if (swapNormal) normalStr = inStr;
  /* Write header: */
  fprintf(outFile,
	  "0 Sphere w/ %s pointing normal. %d triangular elements.\n",
	  normalStr,nn);
  fprintf(outFile,
	  "* This file was automatically generated by \"spheregen\"\n");

  /* Plot # elements to be generated: */
  fprintf(stderr,"* Generating sphere with %d elements\n",nn);

  /* if I can set up ONE side, then the rest can be found from 
   * rotations (and possible translations), but no mirroring(!) of 
   * the first side. */
  for (iside=0; iside<6; iside++){
    switch (iside){
      /* Note: (u,v,n) is a right-hand-system. (outward normal) */
    case 0: 
      yu=1.0; yv=0.0; y0=0.0;
      zu=0.0; zv=1.0; z0=0.0;
      xu=0.0; xv=0.0; x0=sideLength/2.0;
      break;
    case 1: 
      yu=1.0; yv=0.0; y0=0.0;
      zu=0.0; zv=-1.0; z0=0.0;
      xu=0.0; xv=0.0; x0=-sideLength/2.0;
      break;
    case 2: 
      xu=1.0; xv=0.0; x0=0.0;
      zu=0.0; zv=-1.0; z0=0.0;
      yu=0.0; yv=0.0; y0=sideLength/2.0;
      break;
    case 3: 
      xu=1.0; xv=0.0; x0=0.0;
      zu=0.0; zv=1.0; z0=0.0;
      yu=0.0; yv=0.0; y0=-sideLength/2.0;
      break;
    case 4: 
      xu=1.0; xv=0.0; x0=0.0;
      yu=0.0; yv=1.0; y0=0.0;
      zu=0.0; zv=0.0; z0=sideLength/2.0;
      break;
    case 5: 
      xu=1.0; xv=0.0; x0=0.0;
      yu=0.0; yv=-1.0; y0=0.0;
      zu=0.0; zv=0.0; z0=-sideLength/2.0;
      break;
    }
    /* Check if we should grade mesh on this side: */
    if (gradeMesh[iside])uuse = ugvals; /* Graded mesh      */
    else uuse = uvals;                  /* Equidistant mesh */
    for (i=0; i<nDense; i++){
      for (j=0; j<nDense; j++){
	/* Coordinate for this element in u,v;
	 * normal upwards, clockwise ordering  */
	u[0] = uuse[i]; 
	v[0] = uuse[j]; 
	u[1] = uuse[i]; 
	v[1] = uuse[j+1]; 
	u[2] = uuse[i+1]; 
	v[2] = uuse[j+1]; 
	u[3] = uuse[i+1]; 
	v[3] = uuse[j]; 

	/* For each corner:*/ 
	for (ic=0; ic<4; ic++){
	/* Correlate the (u,v) to (x,y,z) for the given side */
	  xout[ic] = xu*u[ic] +xv*v[ic] +x0;
	  yout[ic] = yu*u[ic] +yv*v[ic] +y0;
	  zout[ic] = zu*u[ic] +zv*v[ic] +z0;
	  /* Transform from cube to sphere by changing the radius of 
	   * the present point, i.e. performing a normalization 
	   * based on the center at the origin                          */
	  thisrad = sqrt(xout[ic]*xout[ic]+yout[ic]*yout[ic]+zout[ic]*zout[ic]);
	  xout[ic] *= radius/thisrad;
	  yout[ic] *= radius/thisrad;
	  zout[ic] *= radius/thisrad;
	  /* Add possible translation of center  */
	  xout[ic] += center[0];
	  yout[ic] += center[1];
	  zout[ic] += center[2];
	}
	/* The quad is (most likely) skew. Thus it needs to be 
	 * "cracked" into two triangular elements.
	 * Crack element along shortest diagonal.                      */
	/* First diagonal: */
	d02 =  (xout[0]-xout[2])*(xout[0]-xout[2])
	      +(yout[0]-yout[2])*(yout[0]-yout[2])
	      +(zout[0]-zout[2])*(zout[0]-zout[2]);
	d13 =  (xout[1]-xout[3])*(xout[1]-xout[3])
	      +(yout[1]-yout[3])*(yout[1]-yout[3])
	      +(zout[1]-zout[3])*(zout[1]-zout[3]);
	/* Corners of the two triangular elements: */
	if (d02<=d13){
	  icc0[0] = 0; icc0[1] = 1; icc0[2] = 2;
	  icc1[0] = 0; icc1[1] = 2; icc1[2] = 3;
	}
	else {
	  icc0[0] = 0; icc0[1] = 1; icc0[2] = 3;
	  icc1[0] = 1; icc1[1] = 2; icc1[2] = 3;
	}
	icc[0] = icc0;
	icc[1] = icc1;
	for (ie=0; ie<2; ie++){
	  /* Print line */
	  fprintf(outFile,"T"); /* Triangular element */
	  fprintf(outFile," %d",1);  /*Surface number ("conducter number") */
	  for (ic=0; ic<3; ic++){
	    iplot = icc[ie][ic];
	    if (swapNormal) /* Change the order of printing */
	      iplot = icc[ie][2-ic]; 
	    /* Print this corner to file */
	    fprintf(outFile," %16.8e %16.8e %16.8e",
		    xout[iplot],yout[iplot],zout[iplot]);
	  }
	  fprintf(outFile,"\n");
	} /* Plot next triangle */
	if (dostats){ /* Add two triangles to statistics */
	  for (ie=0; ie<2; ie++){
	    /* First side: */
	    s0[0] = xout[icc[ie][1]]-xout[icc[ie][0]];
	    s0[1] = yout[icc[ie][1]]-yout[icc[ie][0]];
	    s0[2] = zout[icc[ie][1]]-zout[icc[ie][0]];
	    /* Second side: */
	    s1[0] = xout[icc[ie][2]]-xout[icc[ie][0]];
	    s1[1] = yout[icc[ie][2]]-yout[icc[ie][0]];
	    s1[2] = zout[icc[ie][2]]-zout[icc[ie][0]];
	    /* Normal vector (sign may be wrong) */
	    XPROD(normal,s0,s1);
	    /* Area is half of normal length: */
	    area = 0.5 * sqrt(normal[0]*normal[0]+
			      normal[1]*normal[1]+
			      normal[2]*normal[2]);
	    if (area<minarea) minarea = area;
	    if (area>maxarea) maxarea = area;
	    meanarea += area;
	    neles++;
	  }
	}
      }
    }
  }
  /* For statistics complete the calculations arnd write output */
  if (dostats){
    meanarea /= (double) neles;
    fprintf(outFile,"* Statistics:\n");
    fprintf(outFile,"*  %d elements\n",neles);
    fprintf(outFile,"*  Mean size %16.8e\n",meanarea);
    fprintf(outFile,"*  Minimum size %16.8e\n",minarea);
    fprintf(outFile,"*  Maximum size %16.8e\n",maxarea);
    fprintf(outFile,"*  Ratio between largest and smallest area %16.8e\n",
	    maxarea/minarea);
    area = 4.0*PI*radius*radius;
    fprintf(outFile,"*  Relative error on sphere area %16.8e\n",
	    (neles*meanarea-area)/area );
  }

  if (writeToFile) fclose(outFile);
  printf(" finished !\n");
}
