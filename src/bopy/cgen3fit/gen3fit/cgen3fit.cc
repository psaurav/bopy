/*

Changes:
    * gslnonlinfit-maxminphase-centroid.h
        * use centroid to identify the max DC pixel.
    * gslnonlinfit-maxminphase.h
        * use minimum to calculate phase (when max is at the cusp)
    * gslnonlinfit-maxphase.h
        * Get the initial guesses for A, DC,
        * Fit for the maximum DC
        * Use that phase for all pixels
        * Do no phase unwrapping in the end.
    * gslnonlinfit-ndelw-estimatedsigma-prevsol.h
        * use previous sol as init value
    * gslnonlinfit-ndelw-estimatedsigma.h
    * outputs init values

*/
#ifndef _GSLNONLINFIT_H_
#define _GSLNONLINFIT_H_
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_statistics.h>
#include "cgen3fit.h"
#include "fftfit.h"

#ifdef __cplusplus
extern "C" {
#endif

int centroid(const double* data, const unsigned int* axes, double topfraction, unsigned int& row_centroid, unsigned int& col_centroid) {

    //if (data.axes.size() != 2) {
    //    std::cerr << "centroid: dim=" << data.axes.size() << ". This works only for 2D data" << std::endl;
    //    return 1;
    //}

    // create an indexed vector of the data
    //std::cerr << "centroid: axes: " << axes[0] << " " << axes[1] << " " << axes[2] << std::endl;
    std::vector<std::pair<double, unsigned int> > vDIndx;
    for (unsigned int i=0; i<axes[0]*axes[1]; ++i) {
        vDIndx.push_back(std::pair<double, unsigned int>(data[i], i));
    }

    //std::sort(vDIndx.begin(), vDIndx.end(), CmpDIndx);
    //
    // sort it in decreasing order
    std::sort(vDIndx.rbegin(), vDIndx.rend());

    unsigned int csize = axes[0];
    unsigned int rsize = axes[1];
    unsigned int tsize = csize*rsize;
    unsigned int top   = static_cast<unsigned int>(static_cast<double>(tsize)*topfraction); // the topfraction elements

    double rmean=0, cmean=0;
    for (unsigned int i=0; i<top; ++i) {
        unsigned int r = vDIndx[i].second / csize;
        unsigned int c = vDIndx[i].second % csize;
        rmean += (static_cast<double>(r) - rmean)/(i+1);
        cmean += (static_cast<double>(c) - cmean)/(i+1);
    }
    row_centroid = static_cast<unsigned int>(rmean+0.5);
    col_centroid = static_cast<unsigned int>(cmean+0.5);

    return 0;
}

struct APhiDCData {
    const double* dat;
    unsigned int nOfData;
    unsigned int isize;
    unsigned int jsize;
    unsigned int ksize;
    unsigned int i;
    unsigned int j;
    //std::vector<T> sigma;
    double omega; // is generally 1Hz*2*pi
    double del_t; // inder takes this as 0.1, which is 1/10th of 1Hz.
    std::vector<double> wt;
    std::vector<double> t;
    APhiDCData(const double* data, const unsigned int* axes): dat(data) {
        isize = axes[0];
        jsize = axes[1];
        ksize = axes[2];
    }
    inline double getData(unsigned k) {
        return dat[i+isize*(j+jsize*k)];
    } 
    inline double getSigma() {
        return 1;
    } 
    inline void setWT() {
        wt.resize(nOfData);
        t.resize(nOfData);
        for (size_t i=0; i<nOfData; ++i) {
            wt[i] = i*omega*del_t;
            t[i] = i*del_t;
        }
        return;
    }
    
    private: 
        APhiDCData() {}
};

inline double fn_trend(double wt, double d, double c, double p, double alpha, double t) {
    return d + c*c*std::cos(wt + p) + alpha*t;
}

inline double dfna_trend(double wt, double c, double p) {
    return 2*c*std::cos(wt + p);
}

inline double dfnp_trend(double wt, double c, double p) {
    return -c*c*std::sin(wt + p);
}  

int func_cos_trend(const gsl_vector *x, void *data, gsl_vector *f) {
    APhiDCData *d = static_cast<APhiDCData*>(data);
     
    double DC    = gsl_vector_get (x, 0);
    double C     = gsl_vector_get (x, 1);
    double phi   = gsl_vector_get (x, 2);
    double alpha = gsl_vector_get (x, 3);

    for (size_t i=0; i<d->nOfData; ++i) {
        //double t    = d->t[i];
        //double Yi   = C*C*std::cos(d->omega*t+phi) + DC;
        double Yi = fn_trend(d->wt[i], DC, C, phi, alpha, d->t[i]);
        double diff = (Yi - d->getData(i));
        gsl_vector_set (f, i, diff);
    }
    return GSL_SUCCESS;
}

int dfunc_cos_trend(const gsl_vector *x, void *data, gsl_matrix *J) {
    APhiDCData *d = static_cast<APhiDCData*>(data);

    double C    = gsl_vector_get (x, 1);
    double phi  = gsl_vector_get (x, 2);

    for (size_t i=0; i<d->nOfData; ++i) {
        //double t    = d->del_t*i;
        //double arg  = d->omega*t + phi;
        gsl_matrix_set (J, i, 0, 1);
        gsl_matrix_set (J, i, 1, dfna_trend(d->wt[i], C, phi));
        gsl_matrix_set (J, i, 2, dfnp_trend(d->wt[i], C, phi));
        gsl_matrix_set (J, i, 3, d->t[i]);
    }

    return GSL_SUCCESS;
}

int fdfunc_cos_trend(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {

    func_cos_trend(x, data, f);
    dfunc_cos_trend(x, data, J);

    return GSL_SUCCESS;
}



inline double fn(double wt, double d, double c, double p) {
    return d + c*c*std::cos(wt + p);
}

inline double dfna(double wt, double c, double p) {
    return 2*c*std::cos(wt + p);
}

inline double dfnp(double wt, double c, double p) {
    return -c*c*std::sin(wt + p);
}  

// the fitted function is cos(wt+phi) for w < w_g.
int func_cos(const gsl_vector *x, void *data, gsl_vector *f) {
    APhiDCData *d = static_cast<APhiDCData*>(data);
     
    double DC  = gsl_vector_get (x, 0);
    double C   = gsl_vector_get (x, 1);
    double phi = gsl_vector_get (x, 2);

    for (size_t i=0; i<d->nOfData; ++i) {
        double t    = d->del_t*i;
        //double Yi   = C*C*std::cos(d->omega*t+phi) + DC;
        double Yi = fn(d->wt[i], DC, C, phi);
        double diff = (Yi - d->getData(i));
        gsl_vector_set (f, i, diff);
    }
    return GSL_SUCCESS;
}

//int GSLSinDF (const gsl_vector * x, void * data, gsl_matrix * J)
int dfunc_cos(const gsl_vector *x, void *data, gsl_matrix *J) {
    APhiDCData *d = static_cast<APhiDCData*>(data);

    double C    = gsl_vector_get (x, 1);
    double phi  = gsl_vector_get (x, 2);

    for (size_t i=0; i<d->nOfData; ++i) {
        double t    = d->del_t*i;
        double arg  = d->omega*t + phi;
        gsl_matrix_set (J, i, 0, 1);
        //gsl_matrix_set (J, i, 1, 2*C*std::cos(arg));
        //gsl_matrix_set (J, i, 2, -C*C*std::sin(arg));
        gsl_matrix_set (J, i, 1, dfna(d->wt[i], C, phi));
        gsl_matrix_set (J, i, 2, dfnp(d->wt[i], C, phi));
    }

    return GSL_SUCCESS;
}

//int GSLSinFDF (const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
int fdfunc_cos(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {

    func_cos(x, data, f);
    dfunc_cos(x, data, J);

    return GSL_SUCCESS;
}

double get_pcts(APhiDCData& apdData, unsigned int nOfData, unsigned int p, double DC, double C, double phi) {
    double pcts = 0.0;     // Pearson's Cummulative Test Statistic
    for (size_t i=0; i<nOfData; ++i) {
        double dummy = apdData.getData(i)-fn(apdData.wt[i], DC, C, phi);
        dummy = dummy*dummy/fn(apdData.wt[i], DC, C, phi);
        pcts += dummy;
        //std::cerr << apdData.t[i] << " " << apdData.getData(i) << " " << fn(apdData.wt[i], DC, C, phi) << std::endl;
    }
    pcts /= static_cast<double>(nOfData-p);
    return pcts;
}

int get_stats(APhiDCData& apdData, unsigned int nOfData, unsigned int p, double DC, double C, double phi, 
            double& pcts, double& rse, double& nrse) {
    
    pcts = 0.;
    rse = 0.;
    nrse = 0.;
    double A = C*C;
    for (size_t i=0; i<nOfData; ++i) {
        double fnorm = (apdData.getData(i)-DC)/A;
        double dnorm = (fn(apdData.wt[i], DC, C, phi)-DC)/A;
        double dummy = apdData.getData(i)-fn(apdData.wt[i], DC, C, phi);
        rse += dummy*dummy;
        nrse += (fnorm-dnorm)*(fnorm-dnorm);
        //dummy = dummy*dummy/fn(apdData.wt[i], DC, C, phi);
        pcts += (fnorm-dnorm)*(fnorm-dnorm)/fnorm;
    }
    pcts /= static_cast<double>(nOfData-p);
    rse = std::sqrt(rse/static_cast<double>(nOfData-p));
    nrse = std::sqrt(nrse/static_cast<double>(nOfData-p));
    
    return 0;
}

int get_stats_trend(APhiDCData& apdData, unsigned int nOfData, unsigned int p, double DC, double C, double phi, double alpha,
            double& pcts, double& rse, double& nrse) {
    
    pcts = 0.;
    rse = 0.;
    nrse = 0.;
    double A = C*C;
    for (size_t i=0; i<nOfData; ++i) {
        double fnorm = (apdData.getData(i)-DC)/A;
        double dnorm = (fn_trend(apdData.wt[i], DC, C, phi, alpha, apdData.t[i])-DC)/A;
        double dummy = apdData.getData(i)-fn_trend(apdData.wt[i], DC, C, phi, alpha, apdData.t[i]);
        rse += dummy*dummy;
        nrse += (fnorm-dnorm)*(fnorm-dnorm);
        pcts += (fnorm-dnorm)*(fnorm-dnorm)/fnorm;
    }
    pcts /= static_cast<double>(nOfData-p);
    rse = std::sqrt(rse/static_cast<double>(nOfData-p));
    nrse = std::sqrt(nrse/static_cast<double>(nOfData-p));
    
    return 0;
}


#ifdef __cplusplus
extern "C"
#endif
int getAPhiDC(const double* data,
            const unsigned int* axes,
            const double del_t, 
            const double freq, 
            const unsigned int nOfData, 
            const double absErr, 
            const double relErr, 
            const unsigned int maxIter, 
            double* A, 
            double* phi, 
            double* DC,
            unsigned long int* iters,
            long int* Status,
            double* delDC,
            double* delA,
            double* delPhi,
            double* chi2ByDOF) {
            //double* initA,
            //double* initPhi,
            //double* initDC) {

    size_t isize=axes[0];
    size_t jsize=axes[1];
    size_t ksize=axes[2];

    std::vector<double> initA(isize*jsize, static_cast<double>(NAN));
    std::vector<double> initPhi(isize*jsize, static_cast<double>(NAN));
    std::vector<double> initDC(isize*jsize, static_cast<double>(NAN));

    // Initialize output data
    for (unsigned int i=0; i<isize*jsize; ++i) {
        A[i]            = static_cast<double>(NAN);
        DC[i]           = static_cast<double>(NAN);
        phi[i]          = static_cast<double>(NAN);
        iters[i]        = static_cast<unsigned long int>(0);
        Status[i]       = 100;
        delA[i]         = static_cast<double>(NAN);
        delPhi[i]       = static_cast<double>(NAN);
        delDC[i]        = static_cast<double>(NAN);
        chi2ByDOF[i]    = static_cast<double>(NAN);
    }

    APhiDCData apdData(data, axes);
    double omega    = 2*M_PI*freq;   // * freq which = 1 here..
    apdData.omega   = omega;
    apdData.del_t   = del_t;
    apdData.nOfData = nOfData;
    apdData.setWT();
    //apdData.sigma.resize(nOfData, 100.0);

    // init some GSL
    const size_t p      = 3;
    const gsl_multifit_fdfsolver_type *Tp;
    gsl_multifit_fdfsolver *s;
    gsl_multifit_function_fdf funcs;
    funcs.f     = func_cos;
    funcs.df    = dfunc_cos;
    funcs.fdf   = fdfunc_cos;
    funcs.n     = nOfData;
    funcs.p     = p;
    funcs.params= &apdData;

    double est_C = 0.0;
    double est_DC = 0.0;
    double est_phi = 0.0;

    // FFT stuff
    std::vector<double> vec_phi_fft;        // a vector to hold the phi from ffts
    std::vector<double> fft_data(nOfData);  // vector to hold data for fft
    std::vector<double> dap(3);             // vector for dc, ac, phi
    //Gen3FFT::Gen3FFT g3fft(nOfData, freq, del_t);
    Gen3FFT g3fft(nOfData, freq, del_t);

    double common_est_phi;
    
    // initial values
    double dc_max = 0;
    unsigned int i_dc_max, j_dc_max;
    for (unsigned int i=0; i<axes[0]; ++i) { 
        for (unsigned int j=0; j<axes[1]; ++j) { 
            unsigned int idx    = i+isize*j;
            double dmax         = gsl_stats_max(&data[idx], isize*jsize, nOfData); 
            double dmin         = gsl_stats_min(&data[idx], isize*jsize, nOfData); 
            est_C               = std::sqrt((dmax-dmin)/2.0);
            est_DC              = (dmax+dmin)/2.0;
            //size_t dmax_index   = gsl_stats_max_index(&data.data[idx], isize*jsize, nOfData);
            //double est_phi      = omega*del_t*dmax_index;
            // use the minima of the signal to calculate the phase
            size_t dmin_index   = gsl_stats_min_index(&data[idx], isize*jsize, nOfData);
            double est_phi      = omega*del_t*dmin_index;

            // Do phase unwrapping.
            while (est_phi > 2*M_PI) { 
                est_phi -= 2*M_PI;
            }
            est_phi -= M_PI;    // because the minima is being used
            // Point it to the right direction: cos(wt + phi)
            est_phi = 2*M_PI - est_phi;
            
            initA[idx]     = est_C*est_C;
            initPhi[idx]   = est_phi;
            initDC[idx]    = est_DC;
            if (dc_max < est_DC) {
                dc_max = est_DC;
                i_dc_max = i;
                j_dc_max = j;
                common_est_phi = est_phi;
            }
        }
    }

    // Check the row/column conundrum
    unsigned int r, c;
    centroid(&initDC[0], axes, 0.1, r, c);
    std::vector<double> vec_phi;
    std::cerr << "Centroid c r: " << c << " " << r << std::endl;
    //unsigned int cmax = (c+1>=isize-1) ? c+1 : isize-1;
    unsigned int rmin = r - 1;
    unsigned int rmax = r + 1;
    unsigned int cmin = c - 1;
    unsigned int cmax = c + 1;
    rmin = (rmin <= 0) ? 0 : rmin;
    cmin = (cmin <= 0) ? 0 : cmin;
    rmax = (rmax >= jsize - 1) ? jsize - 1 : rmax;
    cmax = (cmax >= isize - 1) ? isize - 1 : cmax;

    std::cerr << "rmin: " << rmin << std::endl;
    std::cerr << "rmax: " << rmax << std::endl;
    std::cerr << "cmin: " << cmin << std::endl;
    std::cerr << "cmax: " << cmax << std::endl;

 
    for (unsigned int i=cmin; i<=cmax; ++i) {
        for (unsigned int j=rmin; j<=rmax; ++j) {
            unsigned int idx = i+isize*j;
            size_t dmin_index   = gsl_stats_min_index(&data[idx], isize*jsize, nOfData);
            double est_phi      = omega*del_t*dmin_index;
            // Do phase unwrapping.
            while (est_phi > 2*M_PI) {
                est_phi -= 2*M_PI;
            }
            est_phi -= M_PI;
            //est_phi = 2*M_PI - est_phi;
            est_phi = - est_phi;
            vec_phi.push_back(est_phi);

            // FFT stuff
            for (int k=0; k<nOfData; ++k) {
                fft_data[k] = data[idx +k*isize*jsize];
            }
            g3fft.fftfit(fft_data, dap);
            vec_phi_fft.push_back(dap[2]);
            std::cerr << dap[2] << " ";
        }
    }
    std::cerr << std::endl;
    
    std::sort(vec_phi.begin(), vec_phi.end());
    // use the "median" value of the phi estimates of the centroid and its nearest neighbors
    common_est_phi = vec_phi[static_cast<unsigned int>(vec_phi.size()/2)];
    std::cerr << "Common estimated phi     = " << common_est_phi << std::endl;

    // FFT stuff
    std::sort(vec_phi_fft.begin(), vec_phi_fft.end());
    double common_est_phi_fft = vec_phi_fft[static_cast<unsigned int>(vec_phi_fft.size()/2)];
    std::cerr << "Common estimated phi FFT = " << common_est_phi_fft << std::endl;
    
    for (unsigned int i=0; i<axes[0]; ++i) { 
        for (unsigned int j=0; j<axes[1]; ++j) { 

            apdData.i = i;
            apdData.j = j;

            unsigned int idx    = i+isize*j;
            est_C = std::sqrt(initA[idx]);
            est_DC = initDC[idx];
            est_phi = common_est_phi_fft;

            gsl_matrix *covar   = gsl_matrix_alloc(p, p);
            double x_init[p]    = {est_DC, est_C, est_phi};
            gsl_vector_view vv  = gsl_vector_view_array(x_init, p);
            Tp                  = gsl_multifit_fdfsolver_lmsder;
            s                   = gsl_multifit_fdfsolver_alloc (Tp, nOfData, p);
            gsl_multifit_fdfsolver_set (s, &funcs, &vv.vector);

            size_t iter = 0;
            int status; 
            do {
                ++iter;
                status          = gsl_multifit_fdfsolver_iterate(s);
                if (status)
                    break;
                status          = gsl_multifit_test_delta (s->dx, s->x, absErr, relErr);
            } while (status == GSL_CONTINUE && iter < maxIter);

            iters[idx]     = static_cast<unsigned long int>(iter);
            Status[idx]    = static_cast<long int>(status);
            DC[idx]        = gsl_vector_get(s->x, 0);
            double C       = gsl_vector_get(s->x, 1);
            A[idx]         = C*C;                          // A = C^2
            phi[idx]       = gsl_vector_get(s->x, 2); 

            if (i == c && j == r) { // centroid
                std::cerr << "Centroid fitted phi: " << phi[idx] << std::endl;
            }

            // save the solutions for the next pixel
            /*
            while (phi.data[idx] > (2*M_PI)) { 
                phi.data[idx]  -= (2*M_PI);
            }
            while (phi.data[idx] < 0.0) { 
                phi.data[idx]  += (2*M_PI);
            }
            //phi.data[idx] = -phi.data[idx];
            //phi.data[idx]  = (2*M_PI) - phi.data[idx];
            */

            // errors etc
            //gsl_multifit_covar(s->J, relErr, covar);  // Inderpreet
            gsl_multifit_covar(s->J, 0.0, covar);       // GSL example

            double chi     = gsl_blas_dnrm2(s->f);
            double dof     = nOfData - p;
            //double c     = GSL_MAX_DBL(1, chi/std::sqrt(dof));
            double c       = chi/std::sqrt(dof);
            delDC[idx]     = c*std::sqrt(gsl_matrix_get(covar,0,0));
            double delC    = c*std::sqrt(gsl_matrix_get(covar,1,1));
            delA[idx]      = 2*delC*std::abs(C);                                 // delA/A = 2*(delC/C)
            delPhi[idx]    = c*std::sqrt(gsl_matrix_get(covar,2,2));
            chi2ByDOF[idx] = chi*chi/dof;

            // free
            gsl_multifit_fdfsolver_free (s);
            gsl_matrix_free (covar);
        }
    }
    return 0;
}

int getAPhiDC2(const double* data,
            const unsigned int* axes,
            const double del_t, 
            const double freq, 
            const double phi_init,
            const unsigned int nOfData, 
            const double absErr, 
            const double relErr, 
            const unsigned int maxIter, 
            double* A, 
            double* phi, 
            double* DC,
            unsigned long int* iters,
            long int* Status,
            double* delDC,
            double* delA,
            double* delPhi,
            double* chi2ByDOF,
            double* pcts,
            double* rse,
            double* nrse) {

    
    size_t isize=axes[0];
    size_t jsize=axes[1];
    size_t ksize=axes[2];

    APhiDCData apdData(data, axes);
    double omega    = 2*M_PI*freq;   // * freq which = 1 here..
    apdData.omega   = omega;
    apdData.del_t   = del_t;
    apdData.nOfData = nOfData;
    apdData.setWT();

    // Initialize output data
    for (unsigned int i=0; i<isize*jsize; ++i) {
        A[i]            = static_cast<double>(NAN);
        DC[i]           = static_cast<double>(NAN);
        phi[i]          = static_cast<double>(NAN);
        iters[i]        = static_cast<unsigned long int>(0);
        Status[i]       = 100;
        delA[i]         = static_cast<double>(NAN);
        delPhi[i]       = static_cast<double>(NAN);
        delDC[i]        = static_cast<double>(NAN);
        chi2ByDOF[i]    = static_cast<double>(NAN);
    }

    // init some GSL
    const size_t p      = 3;
    const gsl_multifit_fdfsolver_type *Tp;
    gsl_multifit_fdfsolver *s;
    Tp              = gsl_multifit_fdfsolver_lmsder;
    s               = gsl_multifit_fdfsolver_alloc (Tp, nOfData, p);

    gsl_multifit_function_fdf funcs;
    funcs.f         = func_cos;
    funcs.df        = dfunc_cos;
    funcs.fdf       = fdfunc_cos;
    funcs.n         = nOfData;
    funcs.p         = p;
    funcs.params    = &apdData;

    for (unsigned int i=0; i<axes[0]; ++i) { 
        for (unsigned int j=0; j<axes[1]; ++j) { 

            apdData.i = i;
            apdData.j = j;
            unsigned int idx    = i+isize*j;

            // Initial values
            double dmax         = gsl_stats_max(&data[idx], isize*jsize, nOfData); 
            double dmin         = gsl_stats_min(&data[idx], isize*jsize, nOfData); 
            double est_C        = std::sqrt((dmax-dmin)/2.0);
            double est_DC       = (dmax+dmin)/2.0;
            double x_init[p]    = {est_DC, est_C, phi_init};

            gsl_vector_view vv  = gsl_vector_view_array(x_init, p);
            gsl_multifit_fdfsolver_set (s, &funcs, &vv.vector);

            size_t iter = 0;
            int status; 
            do {
                ++iter;
                status          = gsl_multifit_fdfsolver_iterate(s);
                if (status)
                    break;
                status          = gsl_multifit_test_delta (s->dx, s->x, absErr, relErr);
            } while (status == GSL_CONTINUE && iter < maxIter);

            iters[idx]     = static_cast<unsigned long int>(iter);
            Status[idx]    = static_cast<long int>(status);
            DC[idx]        = gsl_vector_get(s->x, 0);
            double C       = gsl_vector_get(s->x, 1);
            A[idx]         = C*C;                          // A = C^2
            phi[idx]       = gsl_vector_get(s->x, 2); 

            // errors etc
            //gsl_multifit_covar(s->J, relErr, covar);  // Inderpreet
            gsl_matrix *covar   = gsl_matrix_alloc (p, p);
            gsl_multifit_covar(s->J, 0.0, covar);       // GSL example

            double chi     = gsl_blas_dnrm2(s->f);
            double dof     = nOfData - p;
            double c       = chi/std::sqrt(dof);
            delDC[idx]     = c*std::sqrt(gsl_matrix_get(covar,0,0));
            delA[idx]      = 2.*c*std::sqrt(gsl_matrix_get(covar,1,1))*C;
            delPhi[idx]    = c*std::sqrt(gsl_matrix_get(covar,2,2));
            chi2ByDOF[idx] = chi*chi/dof;
            double _pcts, _rse, _nrse;
            
            get_stats(apdData, nOfData, p, DC[idx], C, phi[idx], _pcts, _rse, _nrse);
            //std::cerr << _rse << " ";
            pcts[idx] = _pcts;
            rse[idx] = _rse; 
            nrse[idx] = _nrse; 

            // free
            gsl_matrix_free (covar);
        }
    }
    std::cerr << std::endl;
    gsl_multifit_fdfsolver_free (s);

    return 0;
}

int getPixelAPhiDC2(const double* data,
            const unsigned int* axes,
            const unsigned long int row,
            const unsigned long int col,
            const double del_t,
            const double freq,
            const double phi_init,
            const unsigned int nOfData,
            const double absErr,
            const double relErr,
            const unsigned int maxIter,
            double& A,
            double& phi,
            double& DC,
            unsigned long int& iters,
            long int& Status,
            double& delDC,
            double& delA,
            double& delPhi,
            double& chi2ByDOF,
            double& pcts,
            double& rse,
            double& nrse,
            double& ssr) {
    
    size_t isize=axes[0];
    size_t jsize=axes[1];
    size_t ksize=axes[2];

    APhiDCData apdData(data, axes);
    double omega    = 2*M_PI*freq;   // * freq which = 1 here..
    apdData.omega   = omega;
    apdData.del_t   = del_t;
    apdData.nOfData = nOfData;
    apdData.i = col;
    apdData.j = row;
    apdData.setWT();

    // init some GSL
    const size_t p  = 3;
    const gsl_multifit_fdfsolver_type *Tp;
    gsl_multifit_fdfsolver *s;
    Tp              = gsl_multifit_fdfsolver_lmsder;
    s               = gsl_multifit_fdfsolver_alloc (Tp, nOfData, p);

    gsl_multifit_function_fdf funcs;
    funcs.f         = func_cos;
    funcs.df        = dfunc_cos;
    funcs.fdf       = fdfunc_cos;
    funcs.n         = nOfData;
    funcs.p         = p;
    funcs.params    = &apdData;
    unsigned int idx    = col+isize*row;
    
    // Initial values
    double dmax         = gsl_stats_max(&data[idx], isize*jsize, nOfData); 
    double dmin         = gsl_stats_min(&data[idx], isize*jsize, nOfData); 
    double est_C        = std::sqrt((dmax-dmin)/2.0);
    double est_DC       = (dmax+dmin)/2.0;
    double x_init[p]    = {est_DC, est_C, phi_init};

    gsl_vector_view vv  = gsl_vector_view_array(x_init, p);
    gsl_multifit_fdfsolver_set (s, &funcs, &vv.vector);

    size_t iter = 0;
    int status; 
    do {
        ++iter;
        status          = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        status          = gsl_multifit_test_delta (s->dx, s->x, absErr, relErr);
    } while (status == GSL_CONTINUE && iter < maxIter);

    iters     = static_cast<unsigned long int>(iter);
    Status    = static_cast<long int>(status);
    DC        = gsl_vector_get(s->x, 0);
    double C  = gsl_vector_get(s->x, 1);
    A         = C*C;                          // A = C^2
    phi       = gsl_vector_get(s->x, 2); 
    std::cerr << A << " " << phi << " " << DC << std::endl;
            
    // errors etc
    //gsl_multifit_covar(s->J, relErr, covar);  // Inderpreet
    gsl_matrix *covar   = gsl_matrix_alloc (p, p);
    gsl_multifit_covar(s->J, 0.0, covar);       // GSL example

    double chi     = gsl_blas_dnrm2(s->f);
    double dof     = nOfData - p;
    double c       = chi/std::sqrt(dof);
    delDC     = c*std::sqrt(gsl_matrix_get(covar,0,0));
    delA      = 2.*c*std::sqrt(gsl_matrix_get(covar,1,1))*C;
    delPhi    = c*std::sqrt(gsl_matrix_get(covar,2,2));
    chi2ByDOF = chi*chi/dof;

    get_stats(apdData, nOfData, p, DC, C, phi, pcts, rse, nrse); 
    ssr = nrse*std::sqrt(static_cast<double>(nOfData-p));
    std::cerr << pcts;

    gsl_multifit_fdfsolver_free (s);

    return 0;
}


int getPixelAPhiDCTr(const double* data,
            const unsigned int* axes,
            const unsigned long int row,
            const unsigned long int col,
            const double del_t,
            const double freq,
            const double phi_init,
            const unsigned int nOfData,
            const double absErr,
            const double relErr,
            const unsigned int maxIter,
            double& A,
            double& phi,
            double& DC,
            double& alpha,
            unsigned long int& iters,
            long int& Status,
            double& delDC,
            double& delA,
            double& delPhi,
            double& delAlpha,
            double& chi2ByDOF,
            double& pcts,
            double& rse,
            double& nrse,
            double& ssr) {
    
    size_t isize=axes[0];
    size_t jsize=axes[1];
    size_t ksize=axes[2];

    APhiDCData apdData(data, axes);
    double omega    = 2*M_PI*freq;   // * freq which = 1 here..
    apdData.omega   = omega;
    apdData.del_t   = del_t;
    apdData.nOfData = nOfData;
    apdData.i = col;
    apdData.j = row;
    apdData.setWT();

    // init some GSL
    const size_t p  = 4;
    const gsl_multifit_fdfsolver_type *Tp;
    gsl_multifit_fdfsolver *s;
    Tp              = gsl_multifit_fdfsolver_lmsder;
    s               = gsl_multifit_fdfsolver_alloc (Tp, nOfData, p);

    gsl_multifit_function_fdf funcs;
    funcs.f         = func_cos_trend;
    funcs.df        = dfunc_cos_trend;
    funcs.fdf       = fdfunc_cos_trend;
    funcs.n         = nOfData;
    funcs.p         = p;
    funcs.params    = &apdData;
    unsigned int idx    = col+isize*row;
    
    // Initial values
    double dmax         = gsl_stats_max(&data[idx], isize*jsize, nOfData); 
    double dmin         = gsl_stats_min(&data[idx], isize*jsize, nOfData); 
    double est_C        = std::sqrt((dmax-dmin)/2.0);
    double est_DC       = (dmax+dmin)/2.0;
    double alpha_init   = 0.;                                       // the trend coefficient, init to zero
    double x_init[p]    = {est_DC, est_C, phi_init, alpha_init};

    gsl_vector_view vv  = gsl_vector_view_array(x_init, p);
    gsl_multifit_fdfsolver_set (s, &funcs, &vv.vector);

    size_t iter = 0;
    int status; 
    do {
        ++iter;
        status          = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        status          = gsl_multifit_test_delta (s->dx, s->x, absErr, relErr);
    } while (status == GSL_CONTINUE && iter < maxIter);

    iters     = static_cast<unsigned long int>(iter);
    Status    = static_cast<long int>(status);
    DC        = gsl_vector_get(s->x, 0);
    double C  = gsl_vector_get(s->x, 1);
    A         = C*C;                          // A = C^2
    phi       = gsl_vector_get(s->x, 2); 
    alpha     = gsl_vector_get(s->x, 3); 
    std::cerr << A << " " << phi << " " << DC << " " << alpha << std::endl;
            
    // errors etc
    //gsl_multifit_covar(s->J, relErr, covar);  // Inderpreet
    gsl_matrix *covar   = gsl_matrix_alloc (p, p);
    gsl_multifit_covar(s->J, 0.0, covar);       // GSL example

    double chi     = gsl_blas_dnrm2(s->f);
    double dof     = nOfData - p;
    double c       = chi/std::sqrt(dof);
    delDC     = c*std::sqrt(gsl_matrix_get(covar,0,0));
    delA      = 2.*c*std::sqrt(gsl_matrix_get(covar,1,1))*C;
    delPhi    = c*std::sqrt(gsl_matrix_get(covar,2,2));
    delAlpha  = c*std::sqrt(gsl_matrix_get(covar,3,3));
    chi2ByDOF = chi*chi/dof;

    get_stats_trend(apdData, nOfData, p, DC, C, phi, alpha, pcts, rse, nrse); 
    ssr = nrse*std::sqrt(static_cast<double>(nOfData-p));

    gsl_multifit_fdfsolver_free (s);

    return 0;
}



#ifdef __cplusplus
}
#endif

#endif
