//
//  evalfbonly.c
//  
//
//  Created by Erik Hintz on 6/7/18.
//
//

#include "evalfbonly.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/* this function is used to evaluate f if all a_i are -Inf. */
double evalfbonly(int n, int q, double *U, double *b,
                  double *C, double nu, double ONE, double ZERO)
{
    double finit;
    double y[q-1], sumysq, f, cone, ctwo, diff, tmp;
    /* the following variables (ending in -a) are used to store the corresponding antithetic values */
    double ya[q-1], sumysqa, fa, conea, ctwoa, diffa, tmpa;
    double mean=0;
    int i, j, l;
    
    /* now that .Call is being used, U and C are vectors. To access their elements, we use the rule
     U[s,k] = U[k*numrows+s] */
    
    /* note that C starts initializing at 0 */
    
    /* initialize f. These values can be reused */
    finit = pt( b[0] / C[0], nu, 1, 0);
    
    /* for each row of U */
    for(j=0; j < n; j++){
        
        /* initialize y */
        
        /* U[j] corresponds to U[j,0] in the orginal matrix */
        y[0]  = qt(U[j]*finit, nu, 1, 0);
        ya[0] = qt((1-U[j])*finit, nu, 1, 0);
        
        sumysq  = y[0] * y[0];
        sumysqa = ya[0] * ya[0];
        
        /* C[1] corresponds to C[1,0] in the original matrix */
        cone  = y[0] * C[1];
        conea = ya[0] * C[1];
        
        /* C[q+1] corresponds to C[1,1] in the original matrix */
        ctwo  = sqrt( (nu+1)/(nu+sumysq) ) / C[q+1];
        ctwoa = sqrt( (nu+1)/(nu+sumysqa) ) / C[q+1];
        
        diff  = pt( ( b[1] - cone) * ctwo, nu+1, 1, 0);
        diffa = pt( ( b[1] - conea) * ctwoa, nu+1, 1, 0);
        
        f  = diff * finit;
        fa = diffa * finit;
        
        /* and then go through the columns */
        /* first column has already been taken care of */
        
        for(i=1; i< (q-1); i++){
            
            /* U[i*n+j] corresponds to U[j,i] in the orginal matrix */
            tmp = U[i*n+j]*diff;
            
            /* check if too close to 1 or 0 */
            if(tmp > ONE){
                tmp = ONE;
            }
            if(tmp < ZERO){
                tmp = ZERO;
            }
            
            y[i] = qt(tmp, nu+i, 1, 0) * sqrt( (nu + sumysq) / (nu+i));
            sumysq += y[i] * y[i];
            
            
            /* now the same for the antithetic value */
            tmpa = (1-U[i*n+j])*diffa;
            
            if(tmpa > ONE){
                tmpa = ONE;
            }
            if(tmpa < ZERO){
                tmpa = ZERO;
            }
            
            ya[i] = qt(tmpa, nu+i, 1, 0) * sqrt( (nu + sumysqa) / (nu+i));
            sumysqa += ya[i] * ya[i];
            
            /* calculate the scalar product */
            cone  = 0;
            conea = 0;
            
            for(l=0; l<(i+1); l++){
                /* C[l*q+i+1] corresponds to C[i+1,l] in the original matrix */
                cone  += y[l]* C[l*q+i+1];
                conea += ya[l]* C[l*q+i+1];
            }
            
            /* C[(i+1)*(q+1)] corresponds to C[i+1,i+1] in the original matrix */
            ctwo  = sqrt( (nu+i+1) / (nu + sumysq)) / C[(i+1)*(q+1)];
            ctwoa = sqrt( (nu+i+1) / (nu + sumysqa)) / C[(i+1)*(q+1)];
            
            diff  = pt( (b[i+1] - cone)*ctwo, nu+i+1, 1, 0);
            diffa = pt( (b[i+1] - conea)*ctwoa, nu+i+1, 1, 0);
            
            f  *= diff;
            fa *= diffa;
        }
        mean += (f+fa)/2;
    }
    mean = mean/ n;
    return(mean);
}



SEXP evalfbonly_(SEXP n, SEXP q, SEXP U, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO)
{
    double res = evalfbonly(INTEGER(n)[0], INTEGER(q)[0], REAL(U), REAL(b), REAL(C), REAL(nu)[0], REAL(ONE)[0], REAL(ZERO)[0]);
    
    return ScalarReal(res);
    
}



