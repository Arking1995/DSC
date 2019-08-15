#ifndef _nrutil_already_included
#define _nrutil_already_included
#define __cpluscplus 1
#ifdef __cplusplus
extern "C" {
#endif

static float sqrarg;
#define TSQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)  //used to be SOR

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#define TSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))


float *vector(int,int);
float **matrix(int,int,int,int);
float **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
float **submatrix();
void free_vector(float *,int,int);
void free_dvector();
void free_ivector();
void free_matrix(float **,int,int,int,int);
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror(char []);

#ifdef __cplusplus
};
#endif

#endif
