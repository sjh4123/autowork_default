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
        
        std::string form;
        std::string func;
        std::string model;
        std::string func_HN; // "wiki" or "Kremer"
        
        bool is_use_dc_free;
        bool is_manual_cutlowT;
        bool is_manualFitSingleDS;
        
        int n_fit_HN;
        int n_points;
        int n_peaks;
        int n_used_HN;
        int n_each_side;
        int n_check,n_total,n_count;
        
        double slope_deviation;
        double cutlowT;
        double cutlowT_manual;
        double tau_sectops;
        
        double r2_threshold;
        double r2_peak_threshold;
        double frac_threshold;
        double df_threshold;
        double dc_threshold;
        
        double afreq,afreqi; // angular frequecy 2*pi()*sfreq
        double sfreq,sfreqi; // typical frequency 1/sec
        double storage;      // eps'
        double loss;         // eps''
        
        double afreq_MIN,afreq_MAX;
        double sfreq_MIN,sfreq_MAX;
        double log_afreq_MIN,log_afreq_MAX;
        double log_sfreq_MIN,log_sfreq_MAX;
        double Theq,Tleq;
        double T_actual,T_nominal;
        
        std::vector<double> loss_max,loss_base;
        std::vector<double> afreq_max,afreq_base;
        std::vector<double> sfreq_max,sfreq_base;
        std::vector<double> cutoff_lo_afreq,cutoff_lo_sfreq;
        std::vector<double> cutoff_hi_afreq,cutoff_hi_sfreq;
        
        std::vector<int> peak_indices_sampled;
        
        std::vector<bool> is_log_data;
        
        alglib::spline1dinterpolant interpolant;
        alglib::spline1dinterpolant interpolantLossCurve;
        alglib::spline1dinterpolant interpolantReducedLoss;
        alglib::spline1dinterpolant interpolantDCLoss;
        alglib::spline1dinterpolant interpolantMasterCurve;
        alglib::spline1dinterpolant interpolanttaueqCurve;
        
        alglib::real_1d_array c_dc_loss;
        alglib::lsfitreport rep_dc_loss;
        std::vector<alglib::real_1d_array> c_losspeaks;
        
        /* useful data containers */
        std::vector<std::vector<std::vector<double>>> dielectric_freq_sortincreasing;
        std::vector<std::vector<double>> dielectric_freq;
        std::vector<std::vector<double>> dielectric_freq_peak;
        std::vector<std::vector<double>> dielectric_freq_overall;
        std::vector<std::vector<double>> dielectric_freq_overall_reduced;
        std::vector<std::vector<double>> dc_curve_freq_sortincreasing;
        std::vector<std::vector<double>> normal_loss_derivatives;
        std::vector<std::vector<double>> dc_free_derivatives;
        std::vector<std::vector<double>> masterloss;
        
        /* dummy variables */
        int    indexii;
        int    peak_indexii,base_indexii;
        int    whichpeak;
        double afreq_now=0,afreq_pre=0,afreq_prr;
        double f_now=0,df_now=0,d2f_now=0;
        double f_pre=0,df_pre=0,d2f_pre=0;
        double f_prr=0,df_prr=0,d2f_prr=0;
        
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
        
        
        /* fit dielectric relaxation */
        //----------------------------------------------------------------------
        /* HN from wikipedia */
        static void HN_loss_wiki(const alglib::real_1d_array &c,
                                 const alglib::real_1d_array &x,
                                 double &func,
                                 void *ptr);
        /* HN from Kremer's book */
        static void HN_loss_Kremer(const alglib::real_1d_array &c,
                                   const alglib::real_1d_array &x,
                                   double &func,
                                   void *ptr);
        /* dc effect on the loss curve */
        static void dc_loss(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr);
        /* calculate HN function value when given afreq */
        double HN_func_value(const alglib::real_1d_array&,
                             const double afreq,
                             const std::string& func);
        /* calculate HN function value when given afreq */
        double dc_loss_value(const alglib::real_1d_array&,
                             const double afreq);
        /* tau_alpha definitions for Havriliak-Negami relation */
        std::vector<double> HN_tauFit_calc(const alglib::real_1d_array&,
                                           const std::string& func);
        /* conductivity curve detection */
        void dc_curve_sampling(const std::vector<std::vector<double>>&);
        /* Peak sampling methods */
        void peak_sampling_byValue(const std::vector<std::vector<double>>&);
        void peak_sampling_bySlope(const std::vector<std::vector<double>>&);
        bool is_peaksRepeat(const int, std::vector<int>&);
        bool is_peaksTooClose(const int, const std::vector<int>&);
        void peak_validation(const std::vector<std::vector<double>>&, const int);
        void base_validation(const std::vector<std::vector<double>>&, const int);
        void base_sampling(const std::vector<std::vector<double>>&, const int);
        void backward_sampling(const std::vector<std::vector<double>>&, const int);
        void forward_sampling(const std::vector<std::vector<double>>&, const int);
        
        /* curve fitting */
        void fit_dc_loss();
        void fit_loss_peak();
        void mastercurve_superposition(const std::vector<std::vector<double>>&);
        //----------------------------------------------------------------------
        
        
        /* fit taueq(T): T dependence of relaxation time */
        //----------------------------------------------------------------------
        /* VFT */
        static void VFT_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr);
        static void VFT_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             alglib::real_1d_array &grad,
                             void *ptr);
        /* COOP */
        static void COOP_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr);
        static void COOP_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr);
        double calc_time_given_temp(const alglib::real_1d_array& c,
                                    const double temp,
                                    const std::string& model);
        std::vector<double> compute_Tg_fragility(StructureClass& sysVar,
                                                 const alglib::real_1d_array& c,
                                                 const std::string& model);
        std::vector<double> calc_temp_given_time(const alglib::real_1d_array& c,
                                                 const double time,
                                                 const std::string& model);
        //----------------------------------------------------------------------
        
        
        /* Useful functions */
        //----------------------------------------------------------------------
        void clear_containers();
        void cout_peaksinfo();
        void write_peaksinfo(const StructureClass&);
        void calc_taueq_derivatives();
        void get_Theq(const StructureClass&);
        void decide_cutlowT(const std::vector<std::vector<double>>&);
        void calc_derivatives(const std::vector<std::vector<double>>&,
                              const int index_x,
                              const int index_y);
        void calc_loglog_derivatives(const std::vector<std::vector<double>>&,
                                     const int index_x,
                                     const int index_y);
        void lossCurveReduction(const std::vector<std::vector<double>>& A,
                                const std::vector<std::vector<double>>& B);
        //----------------------------------------------------------------------
        
        
        /* read functions */
        //----------------------------------------------------------------------
        void read_all_dielectric_data(const StructureClass& sysVar,
                                      int process=2,
                                      double Temp_d=0);
        void read_individual_taueq_data(const StructureClass& sysVar);
        //----------------------------------------------------------------------
        
        
        /* write functions */
        //----------------------------------------------------------------------
        void write_fit_HN(std::ofstream& outputFile);
        void write_dc_loss(std::ofstream& outputFile);
        void write_master_curve(std::ofstream& outputFile);
        void write_HN_params(std::ofstream& outputFile);
        void write_peak_info(std::ofstream& outputFile);
        void write_taueqToFile(std::ofstream& taueq,
                               std::ofstream& teqinv,
                               std::ofstream& alltaueq,
                               std::ofstream& alltaueqinv);
        void write_fitcutoff(std::ofstream& outputFile);
        void write_Tg_fragility(std::ofstream& outputFile,
                                const alglib::real_1d_array& c,
                                const std::string& model);
        void write_tauFitfit_vs_T(std::ofstream& outputFile,
                                  const alglib::real_1d_array& c,
                                  const std::string& model);
        void write_tauFitfit_vs_invT(std::ofstream& outputFile,
                                     const alglib::real_1d_array& c,
                                     const std::string& model);
        //----------------------------------------------------------------------
        
        
    public:
        
        /* Constructors */
        FitDielectric()=delete;
        FitDielectric(const StructureClass& sysVar);
        FitDielectric(const FitDielectric&)=default;
        /* Assignment */
        FitDielectric& operator= (const FitDielectric&)=default;
        /* Destructor */
        ~FitDielectric()=default;
        
        
        
        /* intialize internal values */
        void initialize();
        
        /* make folder for source dielectric data */
        void make_dielectric_folder(const StructureClass&);
        
        /* preprocess data to a format readable by the program */
        std::vector<std::vector<double>>
        dielectricData_preprocess(const StructureClass& sysVar);
        
        /* fit relaxation profile to the Havriliak-Negami relation */
        void fit_HN(StructureClass&,int process=2,double T=0);
        
        /* fit temperature denpendence of relaxation time */
        void fit_tauFit(StructureClass&);
        
        
        
        /* public setters */
        void set_form(const std::string&);
        void set_func(const std::string&);
        void set_model(const std::string&);
        void set_is_use_dc_free(const bool);
        void set_is_manual_cutlowT(const bool);
        void set_is_manualFitSingleDS(const bool);
        void set_cutlowT_manual(const double);
        void set_is_log_data(const std::vector<bool>&);
        
        /* public getters */
        bool get_is_use_dc_free(){return is_use_dc_free;}
        bool get_is_manual_cutlowT(){return is_manual_cutlowT;}
        bool get_is_manualFitSingleDS(){return is_manualFitSingleDS;}
    };
    
}
#endif /* FITDIELECTRIC_H */





