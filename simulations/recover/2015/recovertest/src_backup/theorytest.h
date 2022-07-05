//
//  theorytest.h
//  cppWork
//
//  Created by SJH on 9/30/15.
//  Copyright (c) 2015 Sean Jh H. All rights reserved.
//

#ifndef THEORYTEST_H
#define THEORYTEST_H

#include "functions.h"
#include "structureclass.h"
#include "workscripts.h"
#include "amdatanalysis.h"
#include "fitdata.h"

namespace autoWork {
    
    class TheoryTest: public FitData
    {
        
    private:
        
        bool is_use_ngp_peak_frame;
        bool is_use_counting1;
        bool is_data_smoothing;
        
        int ngp_peak_frame;
        int ngp_block_frame;
        int peak_frame;
        
        double frame_time;
        double mean_strings;
        double mean_length;
        double mean_length_count1;
        double order_parameter;
        double stringlen;
        double peak_time;
        double fast_threshold;
        double fit_tau_threshold;
        double stringlen_threshold;
        double n_distr_cutoff;
        double p1_extrp;
        
        double delG,delmiu,delHa,delSa,Ea_COOP;
        
        alglib::real_1d_array c_COOP;
        
        std::string model_stringlen_distr; // exp_string, pwr_exp_string
        
        std::string model;
        std::string amdat_version;
        std::string analysispart;
        std::vector<std::vector<double>> tinfo;
        std::vector<std::vector<double>> dataxy_sortdecreasing;    // (x,y)
        std::vector<std::vector<double>> stringlen_sortdecreasing; // (T,stringlen)
        std::vector<std::vector<double>> ngp_data_sortdecreasing;  // (time,ngp)
        std::vector<std::vector<double>> fastfrac_sortdecreasing;  // (T,fast_fraction)        
        std::vector<std::vector<double>> taueq_sortdecreasing_fit;
        std::vector<std::vector<double>> stringlen_sortdecreasing_fit;
        std::vector<std::vector<double>> dwf_sortdecreasing_fit;
        
    public:
        
        /* Constructors */
        TheoryTest()=delete;
        TheoryTest(const StructureClass&,
                   const WorkScripts&,
                   const AmdatAnalysis&,
                   const std::vector<std::vector<double>>&,
                   const std::string& model);
        TheoryTest(const TheoryTest&)=default;
        /* Assignment */
        TheoryTest& operator= (const TheoryTest&)=default;
        /* Destructor */
        ~TheoryTest()=default;
        
        
        /* string analysis by AMDAT */
        void amdatstringanalysis(const StructureClass& sysVar,
                                 const WorkScripts& ws,
                                 AmdatAnalysis& aa);
        
        /* Test the Adam-Gibbs (AG) theory */
        void AGtest(const StructureClass& sysVar);
        
        /* Test the Generalized Localization Model (GLM) */
        void GLMtest(const StructureClass& sysVar);

        
        /* find peak frame of ngp */
        void find_ngp_peak_frame(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d);
        
        
        /* find frame index for a block */
        void find_ngp_block_frame(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d);
        
        
        /* fit string length distribution P(n) to expoenential function */
        void find_mean_stringlen(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d);
        
        
        /* find L(T) data; single point */
        void find_string_length(const StructureClass& sysVar,
                                const int n_trl,
                                const int n_sys,
                                const double Temp_d);
        
        
        /* find L(T) data; f(time) */
        void find_time_stringlen(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d,
                                 const int frame);
        
        
        /* find fast_frac; single point */
        void find_fastfrac(const StructureClass& sysVar,
                           const int n_trl,
                           const int n_sys,
                           const double Temp_d);
        
        
        /* find fast_frac; f(time) */
        void find_time_fastfrac(const StructureClass& sysVar,
                                const int n_trl,
                                const int n_sys,
                                const double Temp_d,
                                const int frame);
        
        
        /* fit L(T) vs T by COOP */
        void fit_stringlen(const StructureClass& sysVar,
                           const int n_sys);
        void fit_stringlen(const StructureClass& sysVar,
                           const int n_sys,
                           const double TA_d);
        
        /* fit L(T) vs T by COOP */
        void fit_DWF(const StructureClass& sysVar,
                     const int n_sys);
        void fit_DWF(const StructureClass& sysVar,
                     const int n_sys,
                     const double TA_d);
        
        /* fit tau(T) vs T by model */
        void fit_taueq(const StructureClass& sysVar,
                       const int n_sys);
        void fit_taueq(const StructureClass& sysVar,
                       const int n_sys,
                       const double TA_d);
        
        
        /* find/write the peak string length from L(T,time) file */
        void write_peak_stringlen(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d);
        
        
        /* find/write the peak string length from L(T,time) file */
        void write_peak_fastfrac(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d);
        
        
        /* symmetric functions for retrieving string length information */
        void read_individual_stringlen(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys);
        void read_individual_equ_stringlen(const StructureClass& sysVar,
                                           const int n_trl,
                                           const int n_sys);
        void read_all_stringlen(const StructureClass& sysVar,
                                const int n_sys);
        void read_all_equ_stringlen(const StructureClass& sysVar,
                                    const int n_sys);
        void read_all_equ_stringlen(const StructureClass& sysVar,
                                    const int n_sys,
                                    const double TA_d);
        void write_stringlen_equ(const StructureClass& sysVar);
        void write_stringlen_equ_avg(const StructureClass& sysVar);
        
        
        /* symmetric functions for retrieving fast_frac information */
        void read_individual_fastfrac(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys);
        void read_individual_equ_fastfrac(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys);
        void read_all_fastfrac(const StructureClass& sysVar,
                               const int n_sys);
        void read_all_fastfrac(const StructureClass& sysVar,
                               const int n_sys,
                               const double TA_d);
        void read_all_equ_fastfrac(const StructureClass& sysVar,
                                   const int n_sys);
        void read_all_equ_fastfrac(const StructureClass& sysVar,
                                   const int n_sys,
                                   const double TA_d);
        
        void write_fastfrac_equ(const StructureClass& sysVar);
        void write_fastfrac_equ_avg(const StructureClass& sysVar);
        
        
        /* find string length at TA */
        void find_LA(const StructureClass& sysVar,
                     const int n_sys);
        
        
        /* find DWF at TA */
        void find_uA2(const StructureClass& sysVar,
                      const int n_sys);
        
        
        /* find results from Arrhenius fit */
        void find_individual_Arrhenius(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys);
        void find_all_Arrhenius(const StructureClass& sysVar,
                                const int n_sys);
        
        
        /* find results from model fit */
        void find_individual_tau0(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys);
        void find_all_tau0(const StructureClass& sysVar,
                           const int n_sys);
        
        
        /* find tau_inf as a fitting parameter */
        void find_AG_tau0_fit(const StructureClass& sysVar,
                              const int n_sys);
        void read_AG_raw(const StructureClass& sysVar,
                         const int n_sys);
        
        /* find tau0_fit as a fitting parameter */
        void find_GLM_tau0_fit(const StructureClass& sysVar,
                               const int n_sys);
        void read_GLM_raw(const StructureClass& sysVar,
                          const int n_sys);
        
        
        /*==( public setters )==*/
        void set_is_use_ngp_peak_frame(const bool);
        void set_is_use_counting1(const bool);
        void set_is_data_smoothing(const bool);
    };
}
#endif /* THEORYTEST_H */
