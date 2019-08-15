#ifndef _POWEL_H_H_
#define _POWEL_H_H_

bool powell(float p[], float **xi, int n, float retval, float ftol, int *iter,
	float *fret, float(*func)(float[]));
static void linmin(float p[], float xi[], int n, float *fret, float(*func)(float[]));
static float f1dim(float x);
static void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float(*func)(float));
static float brent(float ax, float bx, float cx, float(*f)(float), float tol, float *xmin);

#endif