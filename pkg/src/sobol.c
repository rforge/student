/* C function for computing a Sobol sequence **********************************/

#include "sobol.h"

// for the binray to gray function

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// function to convert binary string
// to gray string

/**
 * @title Generate n Points of a d-dimensional Sobol Sequence
 * @param n number of points
 * @param d dimension
 * @param randomize logical (here: integer) indicating whether a digital
 *        shift should be included
 * @param res pointer to the result matrix
 * @skip  number of initial terms in the sequence to be skipped (skip = 0 means
 *        that the sequence starts at the origin)
 * @return void
 * @author Marius Hofert based on C. Lemieux's RandQMC
 */

// C++ program for Binary To Gray
// and Gray to Binary conversion

void bintogray(int n1, int size, int *b)
{

  int a[size],k=1,i=0;

	 while(n1!=0)   /* converting number to its binary equivalent */
	 {  a[i]=n1 % 2;
		 n1/=2;
		 i++;
	 }
	/* printing binary equivalent */
	 //printf("\nThe binary code of the given number is :");
	 //for(i=size-1;i>=0;i--)
	 //printf("%d",a[i]);
	 /* gray code conversion */
	 b[0]=a[size-1];

	 for(i=size-1;i>=0;i--)
	 {  if(a[i]==0 && a[i-1]==0)
		 b[k]=0;
		 if(a[i]==1 && a[i-1]==1)
		 b[k]=0;
		 if(a[i]==0 && a[i-1]==1)
		 b[k]=1;
		 if(a[i]==1 && a[i-1]==0)
		 b[k]=1;
		 k++;

	 }

}

void sobol(int n, int d, int randomize, double *res, int skip)
{
  int i, count, numcols, j, k, initcount;
  int newv, temp, degree, sizeskip;
	int column; /* the column to use */
	unsigned int maxn;
	double U; /* temporal storage for a single element of point */
	double recipd = 0.0;
	unsigned int *lastpoint = NULL;
	unsigned long long int randint, point; /* the int representation of U */
	int rmaxcoeff = 52;
	double rmaxint = pow (2, rmaxcoeff);
	double *rvector;
	int *gskip = NULL;
        unsigned int v[sobolMaxDim][sobolMaxCol];
        lastpoint = (unsigned int *) R_alloc (d, sizeof(unsigned int));

	sizeskip=0;
	if(skip>0)
	  sizeskip = ceil(log(skip+1)/log(2));
	gskip = (int *) R_alloc(sizeskip, sizeof(int));


        memset (lastpoint, 0, sizeof(unsigned int) * d);
	memset (gskip,0, sizeof(int) * sizeskip);

        /* Initialize the V array */

	/* First the first dimension */
	numcols = ceil(log(n + skip)/log(2));
	for (i = 0; i < numcols; i++)
		v[0][i] = 1;

	/* Now the rest of the dimensions */
	for (i = 1; i < d; i++)
	{
		/* Find the degree of polynomial i */
		degree = sobolMaxDegree;
		while (!((poly[i] >> degree) & 0x1))
			degree--;

		/* Copy known values from minit */

                for (j = 0; j < degree; j++){
		  v[i][j]=minit[i-1][j];
		}
		/* Find the rest of the values */
		for (j = degree; j < numcols; j++)
		{
			newv = v[i][j-degree];
			temp = 1;
                        for (k = degree-1; k >= 0; k--) {
				if ((poly[i] >> k) & 0x1)
					newv = newv ^ (v[i][j-(degree-k)] << temp);
				temp++;
			}
			v[i][j] = newv;
		}
	}

	/* Calculate corresponding v<i,j>, v<i,j> = m<i,j>/2^j */
	temp = 1;
	for (j = numcols-2; j >= 0; j--){
		for (i = 0; i < d; i++)
			v[i][j] = v[i][j] << temp;
		temp++;
	}
        /* => v is now updated */

	/* Generate a shift */
	maxn =  ((unsigned int)(1 << numcols)); /* 2^numcols */
	double dblmaxn = (double) maxn;
        if(randomize){
		double u, bigu;
		int iszero;
		rvector = (double *) R_alloc (d, sizeof(double));
		GetRNGstate();
		for (i = 0; i < d; i++) {
			iszero = 1;
			while(iszero == 1){
				u = unif_rand();
				bigu = u*dblmaxn;
				if (floor(bigu)!=bigu) iszero = 0;
			}
			*(rvector+i) = u;
		}
		PutRNGstate();
	}

        /* Compute the recipd */
        recipd = ((double) (1.0 / dblmaxn));


	/* Init result (and randomization) */
	//initcount=skip-1;

	//if(skip == 0){

	if(skip>0){
	  bintogray(skip,sizeskip, gskip);
	}
	initcount=skip;



          for(i=0; i<d; i++){
	        if(skip>0){
	           for(j=0;j<sizeskip;j++){
	              if (gskip[sizeskip-j-1]==1)
	                   *(lastpoint+i) = *(lastpoint+i) ^ v[i][j];
		   }
		}
		//res[i*n] = 0.0;
		if(randomize){
			point = *(lastpoint+i);
                        point = point << (rmaxcoeff - numcols);
			U = *(rvector+i);
			U *= rmaxint;
			randint = (unsigned long long int) U;
			point = point ^ randint;
			res[i*n] = ((double) point)/rmaxint;
		}
		else res[i*n] =  ((double) *(lastpoint+i)) * recipd;
           }
	  //}

	/* Main loop */
        for(count=initcount; count<n + skip-1; count++){
		column=0;
		while ((count >> column) & 0x1) /* '>>' = bitwise right shift by column-many bits and fill with 0s from the left; hexadecimal value of 1 */
			column++;
		/* Calculate the next point using Gray Code */
		for (i = 0; i < d; i++)
		{
			*(lastpoint+i) = *(lastpoint+i) ^ v[i][column];
			if(randomize){
				point = *(lastpoint+i);
				/* point = point * (rmaxint*recipd); */
				point = point << (rmaxcoeff - numcols);
				U = *(rvector+i);
				U *= rmaxint;
				randint = (unsigned long long int) U;
				point = point ^ randint;
				res[i*n+count+1-skip] = ((double) point)/rmaxint;
			}
			else res[i*n+count+1-skip] = ((double) *(lastpoint+i)) * recipd;
		}
	}
	}







/**
 * @title R Interface to C for Generating a Sobol Sequence
 * @param n number of points
 * @param d dimension
 * @param randomize logical indicating whether a digital shift should be
          included
 * @param number of initial terms in the sequence to be skipped (skip = 0 means
 *        that the sequence starts at the origin)
 * @return (n, d)-matrix
 * @author Marius Hofert
 */
SEXP sobol_(SEXP n, SEXP d, SEXP randomize, SEXP skip)
{
	/* Input parameters */
	int n_ = asInteger(n); /* numeric(1) */
	int d_ = asInteger(d); /* numeric(1) */
	int randomize_ = asLogical(randomize); /* 0 (= FALSE), 1 (= TRUE) */
	int skip_ = asInteger(skip); /* numeric(1) */

	/* Create result object */
	SEXP res = PROTECT(allocMatrix(REALSXP, n_, d_)); /* (n,d)-matrix */
	double *res_ = REAL(res); /* pointer to the values of res */

	/* Main */
	sobol(n_, d_, randomize_, res_, skip_);

	/* Return */
	UNPROTECT(1); /* clean-up */
	return res;
}