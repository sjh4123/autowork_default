//
//  Created by SJH on 8/11/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef FITDATA_H
#define FITDATA_H

#include "stdafx.h"
#include "interpolation.h"

#include <algorithm>
#include "alglibfittingkernel.h"
#include "structureclass.h"
#include "workscripts.h"
#include "amdatanalysis.h"

namespace autoWork {
    
    class FitData: public AlglibFittingKernel
    {
        
    private: // the unique part of the fitData class
        
        /* Algorithm Kernel (Assistive Extrapolation of New Temperatures) */
        void shootForNewTemperatures(StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const alglib::real_1d_array& c_org,
                                     const std::string& model);
        
    protected: // members that are also shared in derived classes
        
        std::string amdat_version;
        std::string analysispart;
        std::string relaxation_target;
        
        bool is_fit_sExp;
        bool is_find_DWF;
        bool is_fit_Arrhenius;
        bool is_fit_tauFit;
        bool is_shootForNext;
        bool is_avgFitParams;
        bool is_fit_by_TA;
        bool is_fitByEveryPoint;
        bool is_fit_by_Tc;
        bool is_use_viscosity;
        bool is_applytauFitcut;
        bool is_fit_full_alpha;
        bool is_use_sExp_cut_hi;
        bool is_use_plateautime;
        bool is_use_gammafunc;
        bool is_imposeStrictTeq;
        bool is_calc_thermoData;
        
        int current_regime;
        int n_regime;
        int n_trial;
        int n_system;
        int n_sys_beg;
        int n_sys_end;
        int index;
        int index_largest_Tg;
        int index_largest_m;
        int sExp_cut_index_hi;
        int sExp_cut_index_lo;
        int n_fit_sExp;
        int n_fit_chopping;
        int n_fit_sliding;
        
        double DWF_time;
        double DWF;
        double wavenumber;
        double n_equ_blocks;
        double n_prd_blocks;
        double sExp_tau_lo;      // fit KWW full alpha part: time > 1 tau
        double sExp_cut_hi;      // fit KWW only portion < cutoff value
        double sExp_cut_hi_org;
        double sExp_cut_lo;      // fit KWW only portion < cutoff value
        double sExp_tauFit;      // at this value tau_alpha is defined
        double tauFit_calc;
        double relativeStdev;
        double tauFit_cut;
        double shootratio;
        double largest_Tg,largest_m;
        double Theq,Tleq;
        double T_actual,tau_T,tau_hiT;
        double slope_fit;
        double r2_threshold;
        
        /* Model fit */
        double tau0_model,tau0_Arrhenius,tau0_fit;
        
        /* Arrhenius fit */
        double deviateTA; // percent deviation of VFT from Arrhenius
        double TA_avg,tauA_avg,tauA_Arrhenius;
        double tauA,uA2,Ea,TA,LA;
        double cutTforArrhenius;
        double res_Arrhenius;
        
        /* MCT fit */
        double Tc_MCT,Tc_SOU,Tc_percent;
        
        
        /* useful data containers to class FitData */
        //----------------------------------------------------------------------
        std::vector<double> c_avgVec;
        std::vector<double> r2_Arrhenius;
        std::vector<double> r2_Model;
        std::vector<std::vector<double>> ArrheniusCoeffs;
        std::vector<std::vector<double>> ModelCoeffs;
        std::vector<std::vector<double>> ExtrpCoeffs;
        std::vector<std::vector<double>> correlation_original;
        std::vector<std::vector<double>> correlation_sortincreasing; // (time,correlation)
        std::vector<std::vector<double>> dwf_sortdecreasing;         // (Teq,DWF)
        std::vector<std::vector<double>> dwf_invT_sortdecreasing;    // (1/Teq,DWF)
        std::vector<std::vector<double>> dwf_tau_data;               // (DWF,tau)
        std::vector<std::vector<double>> param_sortdecreasing;
        std::vector<std::vector<double>> thermo_sortdecreasing;
        //----------------------------------------------------------------------
        
        
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
        
        
        /* Fitting Functions */
        //----------------------------------------------------------------------
        
        // KWW form
        static void KWW_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr);
        
        static void KWW_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             alglib::real_1d_array &grad,
                             void *ptr);
        
        // KWW_pwr form
        static void KWW_pwr_func(const alglib::real_1d_array &c,
                                 const alglib::real_1d_array &x,
                                 double &func,
                                 void *ptr);
        
        static void KWW_pwr_grad(const alglib::real_1d_array &c,
                                 const alglib::real_1d_array &x,
                                 double &func,
                                 alglib::real_1d_array &grad,
                                 void *ptr);
        
        // mKWW form
        static void mKWW_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr);
        
        static void mKWW_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr);
        
        // mKWW_pwr form
        static void mKWW_pwr_func(const alglib::real_1d_array &c,
                                  const alglib::real_1d_array &x,
                                  double &func,
                                  void *ptr);
        
        static void mKWW_pwr_grad(const alglib::real_1d_array &c,
                                  const alglib::real_1d_array &x,
                                  double &func,
                                  alglib::real_1d_array &grad,
                                  void *ptr);
        
        //  tau_alpha definitions for stretched exponential function
        double sExp_tauFit_interpolateFvalue(const alglib::real_1d_array&,
                                             const std::string& tcalc);
        double sExp_tauFit_gammafunction(const alglib::real_1d_array&,
                                         const std::string& tcalc);
        
        // Calculate sExp function value when given time
        double sExp_func_value(const alglib::real_1d_array&,
                               const double time,
                               const std::string& tcalc);
        
        /* linear model */
        static void linear_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr);
        static void linear_grad(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                alglib::real_1d_array &grad,
                                void *ptr);
        
        /* exponential model for string length distribution */
        static void exp_string_func(const alglib::real_1d_array &c,
                                    const alglib::real_1d_array &x,
                                    double &func,
                                    void *ptr);
        static void exp_string_grad(const alglib::real_1d_array &c,
                                    const alglib::real_1d_array &x,
                                    double &func,
                                    alglib::real_1d_array &grad,
                                    void *ptr);
        
        /* Power-exponential model for string length distribution */
        static void pwr_exp_string_func(const alglib::real_1d_array &c,
                                        const alglib::real_1d_array &x,
                                        double &func,
                                        void *ptr);
        static void pwr_exp_string_grad(const alglib::real_1d_array &c,
                                        const alglib::real_1d_array &x,
                                        double &func,
                                        alglib::real_1d_array &grad,
                                        void *ptr);
        
        /* Arrhenius */
        static void Arrhenius_func(const alglib::real_1d_array &c,
                                   const alglib::real_1d_array &x,
                                   double &func,
                                   void *ptr);
        static void Arrhenius_grad(const alglib::real_1d_array &c,
                                   const alglib::real_1d_array &x,
                                   double &func,
                                   alglib::real_1d_array &grad,
                                   void *ptr);
        
        /* 3-PARAMETER MODELS */
        
        /* MCT */
        static void MCT_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr);
        static void MCT_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             alglib::real_1d_array &grad,
                             void *ptr);
        
        /* SOU */
        static void SOU_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr);
        static void SOU_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             alglib::real_1d_array &grad,
                             void *ptr);
        
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
        
        /* Double Exponential */
        static void Mauro_func(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               void *ptr);
        static void Mauro_grad(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               alglib::real_1d_array &grad,
                               void *ptr);
        
        /* Avramov–Milchev (AM) */
        static void AM_func(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr);
        static void AM_grad(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            alglib::real_1d_array &grad,
                            void *ptr);
        
        /* Dyre–Granato */
        static void DG_func(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr);
        static void DG_grad(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            alglib::real_1d_array &grad,
                            void *ptr);
        
        /* ArrheniusII */
        static void ArrheniusII_func(const alglib::real_1d_array &c,
                                     const alglib::real_1d_array &x,
                                     double &func,
                                     void *ptr);
        
        /* Generalied Localization Model (GLM) -- 3 free parameters */
        static void GLM_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr);
        static void GLM_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             alglib::real_1d_array &grad,
                             void *ptr);
        
        /* Generalied Localization Model (GLM) -- 1 free parameter */
        static void GLM1_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr);
        static void GLM1_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr);
        
        /* 4-PARAMETER MODELS */
        
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
        
        
        /* DEAG */
        static void DEAG_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr);
        static void DEAG_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr);
        
        /* Cohen-Grest (CG) */
        static void CG_func(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr);
        static void CG_grad(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            alglib::real_1d_array &grad,
                            void *ptr);
        
        /* Four-Param VFT */
        static void FourParamVFT_func(const alglib::real_1d_array &c,
                                      const alglib::real_1d_array &x,
                                      double &func,
                                      void *ptr);
        static void FourParamVFT_grad(const alglib::real_1d_array &c,
                                      const alglib::real_1d_array &x,
                                      double &func,
                                      alglib::real_1d_array &grad,
                                      void *ptr);
        
        /* ArrheniusIII */
        static void ArrheniusIII_func(const alglib::real_1d_array &c,
                                      const alglib::real_1d_array &x,
                                      double &func,
                                      void *ptr);
        //----------------------------------------------------------------------
        
        
        /* Data Retrieving Functions */
        //----------------------------------------------------------------------
        /* Debye-Waller Factor */
        void read_individual_DWF(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys);
        void read_individual_equ_DWF(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys);
        void read_all_DWF(const StructureClass& sysVar,
                          const int n_sys);
        void read_all_equ_DWF(const StructureClass& sysVar,
                              const int n_sys);
        void read_all_equ_DWF(const StructureClass& sysVar,
                              const int n_sys,
                              const double TA_d);
        
        /* Density Correlation Function */
        void skim_individual_correlation_data(const StructureClass& sysVar,
                                              const int n_trl,
                                              const int n_sys,
                                              const double Temp_d);
        void read_individual_correlation_data(const StructureClass& sysVar,
                                              const int n_trl,
                                              const int n_sys,
                                              const double Temp_d);
        
        /* KWW fit paramters */
        void read_all_sExp_params(const StructureClass& sysVar,
                                  const int n_sys);
        
        /* taueq */
        void read_individual_taueq_data(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys);
        void read_individual_taueq_data(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys,
                                        const std::string& model);
        void read_all_taueq_data(const StructureClass& sysVar,
                                 const int n_sys);
        void read_all_taueq_data(const StructureClass& sysVar,
                                 const int n_sys,
                                 const double TA_d);
        void read_all_taueq_data(const StructureClass& sysVar,
                                 const int n_sys,
                                 const std::string& cond);
        
        /* Teq */
        void read_individual_temps_data(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys);
        void read_all_temps_data(const StructureClass& sysVar,
                                 const int n_sys);
        //----------------------------------------------------------------------
        
        
        /* Data Output Functions */
        //----------------------------------------------------------------------
        void write_fitcutoff(std::ofstream& outputFile,
                             const StructureClass& sysVar,
                             const std::string& model);
        
        void write_errcurve_actual(std::ofstream& outputFile,
                                   const std::string& model,
                                   const std::vector<std::vector<double>>&);
        void write_errcurve_log(std::ofstream& outputFile,
                                const std::string& model,
                                const std::vector<std::vector<double>>&);
        void write_fitavg(StructureClass& sysVar,
                          const int n_sys,
                          const std::string& model);
        
        /* Glass formation properties */
        void write_Tg_fragility(std::ofstream& outputFile,
                                const alglib::real_1d_array& c,
                                const std::string& model);
        
        void write_fit_correlation(std::ofstream& outputFile,
                                   const std::string& tcalc);
        
        void write_sExp_params(std::ofstream& outputFile,
                               const double Temp_d);
        
        void write_sExp_params_avg(const StructureClass& sysVar,
                                   const int n_sys,
                                   const std::string& tcalc);
        
        void write_tauFitfit_vs_T(std::ofstream& outputFile,
                                  const alglib::real_1d_array& c,
                                  const std::string& model);
        void write_tauFitfit_vs_invT(std::ofstream& outputFile,
                                     const alglib::real_1d_array& c,
                                     const std::string& model);
        
        void write_TA(const StructureClass& sysVar,
                      const int n_sys);
        
        void write_tau0(StructureClass& sysVar,
                        const int n_sys,
                        const std::string& model);
        
        double write_Tc(const StructureClass& sysVar,
                        const int n_sys);
        //----------------------------------------------------------------------
        
        
        /* Additional Fiiting Functions */
        //----------------------------------------------------------------------
        /* Data Processing */
        void fit_every_postprocessing();
        void fit_every_preprocessing(const StructureClass& sysVar,
                                     const int n_sys,
                                     const std::string& cond);
        void fit_every_processing(const int,const int);
        
        /* Fitting Kernel */
        void fit_continuousRange(const StructureClass& sysVar,
                                 const int n_sys,
                                 const std::string& model);
        void fit_neighboringNpoints(const StructureClass& sysVar,
                                    const int n_sys,
                                    const std::string& model);
        void fit_partial(const StructureClass& sysVar,
                         const int n_sys,
                         const std::string& model);
        //----------------------------------------------------------------------
        
        
        /* Utilities Functions */
        //----------------------------------------------------------------------
        void save_fit_coeffs(StructureClass&);
        void load_fit_coeffs(const StructureClass&);
        void Tc_correction(const std::string& model);
        
        void useDynamicsRange(StructureClass& sysVar,
                              const int n_sys);
        
        double calc_time_given_temp(const alglib::real_1d_array& c,
                                    const double temp,
                                    const std::string& model);
        
        std::vector<double> calc_temp_given_time(const alglib::real_1d_array& c,
                                                 const double time,
                                                 const std::string& model);
        
        void compute_Tg_fragility(StructureClass& sysVar,
                                  const alglib::real_1d_array& c,
                                  const std::string& model);
        
        void find_largest_Tg_fragility(StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const alglib::real_1d_array& c,
                                       const std::string& model);
        
        double error_Tg(const std::string& model);
        double error_fragility(const std::string& model,
                               const double Tg_value,
                               const double Tg_error);
        
        double get_2pt_slope(const std::vector<double>& xy0,
                             const std::vector<double>& xy1)
        {return (xy1[1]-xy0[1])/(xy1[0]-xy0[0]);}
        
        std::vector<double> get_interpolated_Ts(StructureClass& sysVar,
                                                const int n_trl,
                                                const int n_sys,
                                                const int n_temps,
                                                const double Tl,
                                                const double Th);
        
        void avgfitParams(const StructureClass& sysVar,
                          const int n_trl,
                          const int n_sys,
                          const std::string& model);
        
        void avgfitParams_extrp(const StructureClass& sysVar,
                                const int n_trl,
                                const int n_sys,
                                const std::string& extrp);
        
        void find_cutTforArrhenius(const StructureClass& sysVar,
                                   const int n_sys);
        
        void find_TA(const StructureClass& sysVar,
                     const int n_sys,
                     const std::vector<std::vector<double>> tauVec);
        
        void find_tau0(StructureClass& sysVar,
                       const int n_sys,
                       const std::string& model);
        
        double get_Theq(StructureClass& sysVar,
                        const int n_trl,
                        const int n_sys);
        double get_Theq(StructureClass& sysVar,
                        const int n_sys);
        
        double get_Tleq(StructureClass& sysVar,
                        const int n_trl,
                        const int n_sys);
        double get_Tleq(StructureClass& sysVar,
                        const int n_sys);
        //----------------------------------------------------------------------
        
        
    public:
        
        /* Constructors */
        FitData()=delete;
        FitData(const StructureClass&,
                const WorkScripts&,
                const AmdatAnalysis&);
        FitData(const FitData&)=default;
        /* Assignment */
        FitData& operator= (const FitData&)=default;
        /* Destructor */
        ~FitData()=default;
        
        
        /* write out quenching temperatures */
        void write_qchTs(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys);
        
        /* write out in-equilibrium temperatures */
        void write_equTs(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys);
        
        /* check # of in-equilibrium T's to see if more simulations are needed */
        bool check_is_retry(StructureClass&);
        
        /* fit relaxation profile to stretched exponetial function */
        void fit_sExp(StructureClass&,
                      const int n_trl,
                      const int n_sys,
                      const double T,
                      const std::string& tcalc);
        
        /* find the Debye-Waller factor at each temperature */
        void find_individual_msd_data(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d);
        void write_DWF(const StructureClass&,
                       const int n_trl,
                       const int n_sys);
        void write_DWF_equ(const StructureClass& sysVar);
        void write_DWF_equ_avg(const StructureClass& sysVar);
        
        /* fit T dependence of relaxation time to the Arrhenius relation */
        void fit_Arrhenius(StructureClass&,
                           const int n_trl,
                           const int n_sys);
        
        /* fit T dependence of relaxation time to a chosen model */
        void fit_tauFit(StructureClass&,
                        const int n_trl,
                        const int n_sys,
                        const std::vector<std::string>& relmodel);
        
        
        /*==( public setters )==*/
        
        
        /* string */
        void set_analysispart(const std::string&);
        void set_relaxation_target(const std::string&);
        
        /* bool */
        void set_is_fit_sExp(const bool);
        void set_is_find_DWF(const bool);
        void set_is_fit_Arrhenius(const bool);
        void set_is_fit_tauFit(const bool);
        void set_is_shootForNext(const bool);
        void set_is_avgFitParams(const bool);
        void set_is_fit_by_TA(const bool);
        void set_is_fitByEveryPoint(const bool);
        void set_is_fit_by_Tc(const bool);
        void set_is_applytauFitcut(const bool);
        void set_is_fit_full_alpha(const bool);
        void set_is_use_gammafunc(const bool);
        void set_is_imposeStrictTeq(const bool);
        void set_is_calc_thermoData(const bool);
        
        /* int */
        void set_index_largest_Tg(const int);
        void set_index_largest_m(const int);
        
        /* double */
        void set_sExp_tauFit(const double);
        void set_shootratio(const double);
        void set_cutTforArrhenius(const double);
        void set_largest_Tg(const double);
        void set_largest_m(const double);
        void set_Theq(const double);
        void set_Tleq(const double);
        void set_TA_avg(const double);
        void set_deviateTA(const double);
        void set_Tc_MCT(const double);
        void set_Tc_SOU(const double);
        void set_tauFit_cut(const double);
        
        /*==( public getters )==*/
        
        /* string */
        const std::string get_analysispart() const {return analysispart;}
        const std::string get_relaxation_target() const {return relaxation_target;}
        
        /* bool */
        bool get_is_fit_sExp() const {return is_fit_sExp;}
        bool get_is_find_DWF() const {return is_find_DWF;}
        bool get_is_fit_Arrhenius() const {return is_fit_Arrhenius;}
        bool get_is_fit_tauFit() const {return is_fit_tauFit;}
        bool get_is_shootForNext() const {return is_shootForNext;}
        bool get_is_avgFitParams() const {return is_avgFitParams;}
        bool get_is_fit_by_TA() const {return is_fit_by_TA;}
        bool get_is_fitByEveryPoint() const {return is_fitByEveryPoint;}
        bool get_is_fit_by_Tc() const {return is_fit_by_Tc;}
        bool get_is_applytauFitcut() const {return is_applytauFitcut;}
        bool get_is_fit_full_alpha() const {return is_fit_full_alpha;}
        bool get_is_calc_thermoData() const {return is_calc_thermoData;}
        
        /* int */
        int get_index_largest_Tg() const {return index_largest_Tg;}
        int get_index_largest_m() const {return index_largest_m;}
        
        /* double */
        double get_n_equ_blocks() const {return n_equ_blocks;}
        double get_n_prd_blocks() const {return n_prd_blocks;}
        double get_sExp_tauFit() const {return sExp_tauFit;}
        double get_tauFit_cut() const {return tauFit_cut;}
        double get_shootratio() const {return shootratio;}
        double get_cutTforArrhenius() const {return cutTforArrhenius;}
        double get_largest_Tg() const {return largest_Tg;}
        double get_largest_m() const {return largest_m;}
        double get_TA_avg() const {return TA_avg;}
        double get_tauA_avg() const {return tauA_avg;}
        double get_deviateTA() const {return deviateTA;}
        double get_Tc_MCT() const {return Tc_MCT;}
        double get_Tc_SOU() const {return Tc_SOU;}
    };
    
}
#endif
