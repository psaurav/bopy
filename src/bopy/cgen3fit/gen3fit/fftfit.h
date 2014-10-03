#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>

#ifndef _GEN3FFT_H_
#define _GEN3FFT_H_

class Gen3FFT {
    public:
        Gen3FFT(int n, double freq, double dt);

        ~Gen3FFT();

        void fftfit(const std::vector<double>& data, std::vector<double>& dap);

    private:
        int n;
        double freq;
        double dt;
        double *in;
        double *out;
        int imin;
        fftw_plan p;

        Gen3FFT(void);

        void ap_r2hc_index(void);
};

#endif

