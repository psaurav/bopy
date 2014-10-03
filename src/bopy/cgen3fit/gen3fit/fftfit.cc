#include "fftfit.h"

Gen3FFT::Gen3FFT(int n, double freq, double dt): n(n), freq(freq), dt(dt) {
    in = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n));
    out = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n));
    p = fftw_plan_r2r_1d(n, in, out, FFTW_R2HC, FFTW_MEASURE);
    ap_r2hc_index();
    std::cerr << "------------------------" << std::endl;
    std::cerr << "N      : " << n << std::endl;
    std::cerr << "freq   : " << freq << std::endl;
    std::cerr << "dt     : " << dt << std::endl;
    std::cerr << "imin   : " << imin << std::endl;
    std::cerr << "------------------------" << std::endl;
}

Gen3FFT::~Gen3FFT() {
    fftw_free(in);
    fftw_free(out);
}

void Gen3FFT::fftfit(const std::vector<double>& data, std::vector<double>& dap) {
    for (int i=0; i<n; ++i)
        in[i] = data[i];
    fftw_execute(p);
    dap.clear();
    dap.resize(3);
    dap[0] = out[0]/n;
    dap[1] = std::sqrt(out[imin]*out[imin] + out[n-imin]*out[n-imin])/(n/2);
    dap[2] = std::atan2(out[n-imin], out[imin]);
}

void Gen3FFT::ap_r2hc_index() {
    double fmax = 1.0/dt/2.0;               // the maximum frequency
    int nf = n/2 + 1;                       // the number of frequencies
    double df = fmax/(nf-1.0);              // the frequency resolution
    double dffmin = fmax;
    imin = 0;
    for (int i=1; i<nf; ++i) {
        double dff = std::abs(i*df - freq);
        if (dff < dffmin) {
            dffmin = dff;
            imin = i;
        }
    }
    return;
}
