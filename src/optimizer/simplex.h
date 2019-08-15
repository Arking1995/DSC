#ifndef _SIMPLEX_H_H_
#define _SIMPLEX_H_H_

void amoeba(float **p, float y[], int ndim, float ftol, 
            float (*funk)(float []), int *nfunk);

#endif