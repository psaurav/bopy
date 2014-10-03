#ifndef _CLINFIT_H_
#define _CLINFIT_H_
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
            double* chi2ByDOF);


#ifdef __cplusplus
}
#endif
#endif
