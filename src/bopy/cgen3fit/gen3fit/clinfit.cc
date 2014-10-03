#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include "clinfit.h"

#ifdef __cplusplus
extern "C" {
#endif

int getLinAPhiDC(const double* pinv,
            const double* data,
            const unsigned int* axes,
            const unsigned int nOfData,
            double* work_arr,
            double* A,
            double* phi,
            double* DC,
            double* delDC,
            double* delA,
            double* delPhi,
            double* chi2ByDOF) {

    size_t isize=axes[0];
    size_t jsize=axes[1];
    size_t ksize=axes[2];

    const unsigned int npar = 3;
    double par[npar];

    // Initialize output data
    for (unsigned int i=0; i<isize*jsize; ++i) {
        A[i]            = static_cast<double>(NAN);
        DC[i]           = static_cast<double>(NAN);
        phi[i]          = static_cast<double>(NAN);
        delA[i]         = static_cast<double>(NAN);
        delPhi[i]       = static_cast<double>(NAN);
        delDC[i]        = static_cast<double>(NAN);
        chi2ByDOF[i]    = static_cast<double>(NAN);
    }

    size_t p;
    size_t i, j, k;
    for (i=0; i<isize; ++i) {
        for (j=0; j<jsize; ++j) {
            for (p=0; p<npar; ++p) {
                // Kahan summation
                double s = 0.0;
                double c = 0.0;
                for (k=0; k<nOfData; ++k) {
                    double y = data[i+isize*(j+jsize*k)]*pinv[k+nOfData*p] - c;     // So far, so good: c is zero.
                    double t = s + y;                               // Alas, sum is big, y small, so low-order digits of y are lost.
                    c = (t - s) - y;                                // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
                    s = t;                                          // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
                    // Next time around, the lost low part will be added to y in a fresh attempt.
                }
                par[p] = s;
            }
            DC[i+isize*j] = par[0];
            A[i+isize*j] = sqrt(par[1]*par[1] + par[2]*par[2]);
            phi[i+isize*j] = atan2(par[2], par[1]);
        }
    }
    return 0;
}

#ifdef __cplusplus
}
#endif

