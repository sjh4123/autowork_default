/*******************************************************************************
 PreSQ & AUTOWORK
 Copyright (c) Jui-Hsiang (Sean) Hung. All rights reserved.
 
 >>> SOURCE LICENSE >>>
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation (www.fsf.org); either version 2 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 A copy of the GNU General Public License is available at
 http://www.fsf.org/licensing/licenses
 >>> END OF LICENSE >>>
 ******************************************************************************/

#ifndef FITDIELECTRIC_H
#define FITDIELECTRIC_H

#include "alglibfittingkernel.h"
#include "structureclass.h"

namespace autoWork
{
    class FitDielectric: public AlglibFittingKernel
    {
        bool is_use_dc_free;
        bool is_manual_cutlowT;
        bool is_manualFitSingleDS;
        
        /* counting paramters */
        int n_fit_HN;
        int n_check;
        int n_used_HN;
        int n_each_side;
        int n_points;
        int n_dcpoints;
        int n_peaks;
        int n_total;
        int n_count;
        int n_reduction;
        int counter;
        int process_id;
        int n_points_min;
        
        double r2_threshold;
        double r2_peak_threshold;
        double r2_combined;
        double r2_global_pre;
        double r2_global_now;
        double n_peaks_pre;
        double n_peaks_now;
        double lossfrac_cut;
        double npoints_frac;
        double df_threshold;
        double dc_threshold;
        double slope_deviation;
        double cutlowT;
        double cutlowT_manual;
        
        /* range parameters */
        double afreq_MIN;
        double afreq_MAX;
        double sfreq_MIN;
        double sfreq_MAX;
        double log_afreq_MIN;
        double log_afreq_MAX;
        double log_sfreq_MIN;
        double log_sfreq_MAX;
        double Theq;
        double Tleq;
        double T_actual;
        double T_nominal;
        
        /* others */
        double tau_sectops;
        
        std::string form;
        std::string func;
        std::string model;
        std::string func_HN;
        std::string process;
        
        /** STL containers **/
        //----------------------------------------------------------------------
        std::vector<bool>   is_log_data;
        std::vector<int>    spc_id,T_id;
        std::vector<double> loss_max,loss_base;
        std::vector<double> afreq_max,afreq_base;
        std::vector<double> sfreq_max,sfreq_base;
        std::vector<double> cutoff_lo_afreq,cutoff_lo_sfreq;
        std::vector<double> cutoff_hi_afreq,cutoff_hi_sfreq;
        std::vector<double> sfreq_range;        
        
        /** internal dummy variables **/
        //----------------------------------------------------------------------
        int    indexii;
        int    peak_indexii;
        int    base_indexii;
        int    whichpeak;
        double afreq;
        double afreqi;    // angular frequecy 2*pi()*sfreq
        double sfreq;
        double sfreqi;    // typical frequency 1/sec
        double storage;   // eps'
        double loss;      // eps''
        double afreq_now;
        double afreq_pre;
        double afreq_prr;
        double f_now;
        double df_now;
        double d2f_now;
        double f_pre;
        double df_pre;
        double d2f_pre;
        double f_prr;
        double df_prr;
        double d2f_prr;
        
        /** ALGLIB fitting DS **/
        //----------------------------------------------------------------------
        alglib::real_1d_array c_dc_loss;
        alglib::lsfitreport   rep_dc_loss;
        alglib::spline1dinterpolant interpolantRawLoss;
        alglib::spline1dinterpolant interpolantDCLoss;
        alglib::spline1dinterpolant interpolanttaueqCurve;
        std::vector<alglib::real_1d_array> c_combined;
        std::vector<alglib::real_1d_array> c_global;
        std::vector<alglib::real_1d_array> c_losspeaks;
        std::vector<alglib::lsfitreport> rep_combined;
        std::vector<alglib::lsfitreport> rep_global;
        std::vector<alglib::lsfitreport> rep_losspeaks;
        std::vector<alglib::ae_int_t> info_losspeaks;
        std::vector<std::vector<std::time_t>> usedt_losspeaks;
        std::vector<int> cyclecount_losspeaks;
        
        /** Dielectric data containers **/
        //----------------------------------------------------------------------
        std::vector<int> peak_indices_original;
        std::vector<int> peak_indices_sampled;
        std::vector<std::vector<std::vector<double>>> dielectric_freq_sortincreasing;
        std::vector<std::vector<std::vector<double>>> dielectric_freq_domains;
        std::vector<std::vector<std::vector<double>>> dielectric_freq_reduced_all;
        std::vector<std::vector<std::vector<double>>> dielectric_freq_global;
        std::vector<std::vector<double>> dc_curve_freq_sortincreasing;
        std::vector<std::vector<double>> dielectric_freq_overall;
        std::vector<std::vector<double>> dielectric_freq_dcfree;
        std::vector<std::vector<double>> dielectric_freq_reduced;
        std::vector<std::vector<double>> dielectric_freq;
        std::vector<std::vector<double>> dielectric_freq_peak;
        std::vector<std::vector<double>> masterloss;
        std::vector<std::vector<double>> spc_Temps;
        std::vector<std::vector<std::vector<std::vector<double>>>> dielectric_freq_NPIC;
        
        
        /** internal functions **/
        
        
        /* ALGLIB kernel */
        //----------------------------------------------------------------------
        /* create alglib solver object */
        void alglib_solverobj(const alglib::real_2d_array&,const alglib::real_1d_array&,const alglib::real_1d_array&,alglib::lsfitstate&);
        /* lsfit form */
        void alglib_lsfit_form(const std::string&,alglib::lsfitstate&);
        /* fitting parameters */
        void set_fitParams(const std::string&);
        /* correct fitParams on the fly if needed */
        void fitParams_correction(const std::string&);
        /* alglib nonlinear least-square fit routine */
        void alglib_lsfit(const std::string&);
        
        /* Virtuals */
        //----------------------------------------------------------------------
        double calc_y_given_x(const alglib::real_1d_array&,const double,const std::string&);
        std::vector<double> calc_x_given_y(const alglib::real_1d_array&,const double,const std::string&);
        void write_tauFitfit_vs_T(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&);
        void write_tauFitfit_vs_invT(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&);
        void build_function_interpolant(const std::string&,const std::vector<std::vector<double>>&,const int,const int);
        
        /* fit dielectric relaxation: functional forms */
        //----------------------------------------------------------------------
        /* HN from wikipedia */
        static void HN_loss_wiki(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        /* HN from Kremer's book */
        static void HN_loss_Kremer(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        /* dc effect on the loss curve */
        static void dc_loss(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        /* combined functional forms */
        static void static_combined_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        
        /* fit dielectric relaxation: function values */
        //----------------------------------------------------------------------
        /* calculate HN function value when given afreq */
        double HN_func_value(const alglib::real_1d_array&,const double);
        /* calculate HN function value when given afreq */
        double dc_loss_value(const alglib::real_1d_array&,const double);
        /* calculate combined function value when given afreq */
        double combined_func_value(const alglib::real_1d_array&,const double);
        /* calculate tau_alpha based on Havriliak-Negami relation */
        std::vector<double> HN_tauFit_calc(const alglib::real_1d_array&,const std::string&);
        
        /* sampling methods for dielectric relaxation */
        //----------------------------------------------------------------------
        void dc_curve_sampling(const std::vector<std::vector<double>>&);
        void peak_sampling_byValue(const std::vector<std::vector<double>>&);
        void peak_sampling_bySlope(const std::vector<std::vector<double>>&);
        bool is_afreqRepeat(const double, std::vector<double>&);
        bool is_peaksRepeat(const int, std::vector<int>&);
        bool is_peaksTooClose(const int, const std::vector<int>&);
        void peak_validation(const std::vector<std::vector<double>>&, const int);
        void base_validation(const std::vector<std::vector<double>>&, const int);
        void base_sampling(const std::vector<std::vector<double>>&, const int);
        void backward_sampling(const std::vector<std::vector<double>>&, const int);
        void forward_sampling(const std::vector<std::vector<double>>&, const int);
        
        /* curve fitting processes */
        //----------------------------------------------------------------------
        void fit_dc_loss();
        void fit_loss_peak();
        void fit_local_domains();
        void fit_combined_local_domains();
        void fit_global_domain();
        
        /* utility functions */
        //----------------------------------------------------------------------
        void combine_domain_data_to_fit();
        void combine_global_data_to_fit();
        void reduce_raw_loss();
        void update_fitinfo(const int);
        double find_loss_given_logafreq(const double);
        
        /* fit T dependence of relaxation time */
        //----------------------------------------------------------------------
        /* VFT */
        static void VFT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void VFT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /* COOP */
        static void COOP_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void COOP_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        std::vector<double> compute_Tg_fragility(StructureClass&,const alglib::real_1d_array&,const std::string&);
        
        /* more utility functions */
        //----------------------------------------------------------------------
        void clear_containers();
        void cout_peaksinfo();
        void write_peaksinfo(const StructureClass&);
        void write_localFitinfo(const StructureClass&);
        void write_combinedFitinfo(const StructureClass&);
        void write_globalFitinfo(const StructureClass&);
        void calc_taueq_derivatives();
        void get_Theq(const StructureClass&);//cf.get_Theq() of FitData
        void decide_cutlowT(const std::vector<std::vector<double>>&);
        void lossCurveReduction(const std::vector<std::vector<double>>&,const std::vector<std::vector<double>>&);
        
        /* read functions */
        //----------------------------------------------------------------------
        void read_all_dielectric_data(const StructureClass&,int process=2,double Temp_d=0);
        void read_all_taueq_data(const StructureClass&);
        
        /* write functions */
        //----------------------------------------------------------------------
        void write_fit_HN(std::ofstream&);
        void write_dc_loss(std::ofstream&);
        void write_combined_curve(std::ofstream&);
        void write_global_curve(std::ofstream&);
        void write_reducedLoss(std::ofstream&);
        void write_HN_params(std::ofstream&);
        void write_peak_info(std::ofstream&);
        void write_taueqToFile(std::ofstream&,const double);
        void write_blankToFile(std::ofstream&);
        void write_fitcutoff(std::ofstream&);
        void write_Tg_fragility(std::ofstream&,const alglib::real_1d_array&,const std::string&);
        
    public:
        
        FitDielectric()=default;
        FitDielectric(const StructureClass&);
        FitDielectric(const FitDielectric&)=default;
        FitDielectric& operator= (const FitDielectric&)=default;
        ~FitDielectric()=default;
        
        /** transfer raw dielectric data and rename files **/
        void transfer_raw_data(const StructureClass&);
        
        /** intialize internal values **/
        void initialize();
        
        /** preprocess data to a format readable by the program **/
        std::vector<std::vector<double>> dielectricData_preprocess(const StructureClass&);
        
        /** fit relaxation data to the Havriliak-Negami relation **/
        void fit_dielectric_loss(const StructureClass&,int process=2,double T=0);
        
        /** fit temperature denpendence of relaxation time **/
        void fit_tauFit(const StructureClass&);
        
        /** public setters **/
        //----------------------------------------------------------------------
        void set_form(const std::string&);
        void set_func(const std::string&);
        void set_model(const std::string&);
        void set_process(const std::string&);
        void set_is_use_dc_free(const bool);
        void set_is_manual_cutlowT(const bool);
        void set_is_manualFitSingleDS(const bool);
        void set_process_id(const int);
        void set_cutlowT_manual(const double);
        void set_is_log_data(const std::vector<bool>&);
        void set_spc_id(const std::vector<int>&);
        void set_T_id(const std::vector<int>&);
        void set_sfreq_range(const std::vector<double>&);        
        
        /** public getters **/
        //----------------------------------------------------------------------
        bool get_is_use_dc_free() const {return is_use_dc_free;}
        bool get_is_manual_cutlowT() const {return is_manual_cutlowT;}
        bool get_is_manualFitSingleDS() const {return is_manualFitSingleDS;}
        int  get_process_id() const {return process_id;}
        std::string get_form() const {return form;}
        std::string get_func() const {return func;}
        std::string get_model() const {return model;}
        std::string get_process() const {return process;}
        std::vector<int> get_spc_id() const {return spc_id;}
        std::vector<int> get_T_id() const {return T_id;}
        std::vector<double> get_sfreq_range() const {return sfreq_range;}
    };
}
#endif /* FITDIELECTRIC_H */




