#ifndef CDSOLVER_H
#define CDSOLVER_H
    #include "CRSMat_types.h"
    void solve_cholesky(CRSMatrix *A, double *b, double *x);
    void solve_cholesky_wlds(CRSMatrix *A, double *b, double *x);
#endif