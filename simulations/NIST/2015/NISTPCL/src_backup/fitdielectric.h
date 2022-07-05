//
//  fitdielectric.h
//  cppWork
//
//  Created by SJH on 11/7/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#ifndef FITDIELECTRIC_H
#define FITDIELECTRIC_H

#include "functions.h"
#include "structureclass.h"
#include "alglibfittingkernel.h"
#include <algorithm>

namespace autoWork {
    
    class FitDielectric: public AlglibFittingKernel
    {
        
    private:
        
        std::string dielectric_format; // "NIST" or "NPIC"
        std::string source;            // "Kremer" or "wiki"
        
        int n_fit_HN;
        int n_used_HN;
        int n_check;
        int n_each_side;
        
        bool is_use_dc_free;
        
        double corrFacforT;
        double precision;
        double r2_threshold;
        double frac_threshold;
        
        double afreq,afreqi; // angular frequecy 2*pi()*sfreq
        double sfreq,sfreqi; // typical frequency 1/sec
        double storage;      // eps'
        double loss;         // eps''
        
        double afreq_MIN,afreq_MAX;
        double sfreq_MIN,sfreq_MAX;
        double log_afreq_MIN,log_afreq_MAX;
        double log_sfreq_MIN,log_sfreq_MAX;

        std::vector<double> afreq_max;
        std::vector<double> sfreq_max;
        std::vector<double> loss_max;        
        
        std::vector<double> cutoff_lo_afreq;
        std::vector<double> cutoff_lo_sfreq;
        std::vector<double> cutoff_hi_afreq;
        std::vector<double> cutoff_hi_sfreq;
        
        alglib::spline1dinterpolant losscurve;
        
        /* useful data containers */
        std::vector<std::vector<std::vector<double>>> dielectric_freq_sortincreasing;
        std::vector<std::vector<double>> dielectric_freq_overall;
        std::vector<std::vector<double>> normal_loss_derivatives;
        std::vector<std::vector<double>> dc_free_derivatives;
        
        
        /* ALGLIB kernel */
        //----------------------------------------------------------------------
        /* create alglib solver object */
        void alglib_solverobj(const alglib::real_2d_array& x,
                              const alglib::real_1d_array& y,
                              const alglib::real_1d_array& c,
                              alglib::lsfitstate& state);
        /* lsfit form */
        void alglib_lsfit_form(const std::string&,
                               alglib::lsfitstate& state);
        /* fitting parameters */
        void set_fitParams(const std::string& model);
        /* correct fitParams on the fly if needed */
        void fitParams_correction(const std::string& model);
        /* alglib nonlinear least-square fit routine */
        void alglib_lsfit(const std::string& model);
        //----------------------------------------------------------------------
        
        
        // The loss part of the Havriliak-Negami relation
        // Source: wiki
        static void HN_loss_wiki(const alglib::real_1d_array &c,
                                 const alglib::real_1d_array &x,
                                 double &func,
                                 void *ptr);
        // Source: Kremer
        static void HN_loss_Kremer(const alglib::real_1d_array &c,
                                   const alglib::real_1d_array &x,
                                   double &func,
                                   void *ptr);
        
        // Calculate HN function value when given afreq
        double HN_func_value(const alglib::real_1d_array&,
                             const double afreq,
                             const std::string& tcalc);
        
        // tau_alpha definitions for Havriliak-Negami relation
        double HN_tauFit_calc(const alglib::real_1d_array&,
                              const std::string& tcalc);
        
        void read_all_dielectric_data(const StructureClass& sysVar,
                                      const int process,
                                      const double Temp_d,
                                      const std::string& form);
        
        void write_fit_HN(std::ofstream& outputFile,
                          const std::string& tcalc,
                          const int whichpeak);
        
        void write_HN_params(std::ofstream& outputFile,
                             const double Temp_d);
        
        void write_peak_info(std::ofstream& outputFile,
                             const int n_peaks,
                             const int whichpeak);
        
        void peak_sampling_byValue();
        void peak_sampling_bySlope();
        void calc_loss_derivatives();
        void clear_containers();
        
    public:
        
        /* Constructors */
        FitDielectric()=delete;
        FitDielectric(const StructureClass& sysVar);
        FitDielectric(const FitDielectric&)=default;
        /* Assignment */
        FitDielectric& operator= (const FitDielectric&)=default;
        /* Destructor */
        ~FitDielectric()=default;
        
        
        
        /* preprocess data to a format readable by the program */
        std::vector<std::vector<double>> dielectricData_preprocess(const StructureClass& sysVar,
                                                                   const std::string& form);
        
        /* fit relaxation profile to the Havriliak-Negami relation */
        void fit_HN(StructureClass&,
                    const int process,
                    const double T,
                    const std::string& tcalc,
                    const std::string& format);
        
        
        
        /* public setters */
        void set_is_use_dc_free(const bool);
        
        /* public getters */
        bool get_is_use_dc_free(){return is_use_dc_free;}
    };
    
}
#endif /* FITDIELECTRIC_H */





