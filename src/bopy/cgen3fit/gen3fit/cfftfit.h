#include<iostream>
#include<cmath>
#include<fftw3.h>

class Gen3FFTfit {
    public:
        Gen3FFit(int n, double freq, double delt) : n(n), freq(freq), delt(delt) {

            
            // Get n
            if (n < 10) {
                std::cerr << "Cannot perform FFT fit" << std::endl;
            }
            n = (n >= 20) ? 20 : 10;
            n_fft = n/2 + 1;                // Number of frequencies

            // get freq index
            double f_nyq = 1.0/delt/2.0;           // Nyquist frequency
            double dfreq = f_nyq/(n_fft - 1.0);
            int freq_index = 0;
            double min = -1.0;
            for (i=0; i<n_fft; ++i)
                double dummy = std::abs(i*dfreq - freq);
                if (dummy < min) {
                    min = dummy;
                    freq_index = i;
                }
            }
            std::cerr << "Freq index: " << freq_index << std::endl;
            in = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n));
            out = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n));
            p = fftw_plan_r2r_1d(n, in, out, FFTW_R2HC, FFTW_MEASURE);
            
        }

        void getDAP(double* d, int ii, int jj, int isize, int jsize, int ksize) { // NONONONO
            for (unsigned int k=0; k<n; ++k) {
                in[k] = data.data[i+isize*(j+jsize*k)];
            }
            fftw_execute(p);
            DC.data[idx] = std::abs(out[0])/n;
            A.data[idx] =  std::sqrt(out[si]*out[si] + out[n-si]*out[n-si])/(n/2);
            phi.data[idx] = std::atan2(out[n-si], out[si]);
            

            
    private:
        int n;
        int n_fft;
        double freq;
        double delt;
        double* in;
        double* out;
        int freq_index;
        fftw_plan p;
}
