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

#ifndef THEORYTEST_H
#define THEORYTEST_H

#include "structureclass.h"
#include "workscripts.h"
#include "amdatanalysis.h"
#include "fitdata.h"

namespace autoWork
{
    class TheoryTest: public FitData
    {
        bool is_AGtest;
        bool is_RFOTtest;
        bool is_GLMtest;
        bool is_fitu2T;
        bool is_Leporinitest;
        bool is_HWtest;
        bool is_reltaur;
        bool is_check_amdat_version;
        bool is_use_stringLen_smoothed;
        bool is_use_fast_frac_smoothed;
        bool is_use_counting1;
        bool is_backwardExtrp;
        bool is_data_smoothing;
        bool is_find_fastfrac;
        bool is_use_peak_stringlen;
        bool is_stringmb;
        bool is_use_displacement_list;
        bool is_fit_univILP;
        bool is_use_Lbeta;
        
        int peak_frame;
        int n_trajs;
        int n_total_frames;
        int n_moments;
        
        double time_stringlen;
        double DWF;
        double mean_strings;
        double mean_length;
        double mean_length_count1;
        double order_parameter;
        double stringlen;
        double stringsRg;
        double stringsRg_cutlo;
        double stringsRg_cuthi;
        double peak_time;
        double fast_threshold;
        double fit_tau_threshold;
        double stringlen_threshold;
        double time_stringlen_threshold;
        double n_distr_cutoff;
        double p1_extrp;
        double delG;
        double delmiu;
        double delmiuTA;
        double Ea_COOP;
        double delHa;
        double delSa;
        double tau0_vib;
        double tau0_shift;
        double strings_threshold;
        double fixedfastpercentage;
        
        std::string analysispart;
        std::string model_stringlen_distr; // exp_string, pwr_exp_string
        std::string stringtype;
        std::string cusstrfolder;
        
        /** STL containers **/
        //----------------------------------------------------------------------
        std::vector<std::vector<double>> dataxy_sortdecreasing;   //(x,y)
        std::vector<std::vector<double>> stringsRg_sortdecreasing;//(T,stringsRg)
        std::vector<std::vector<double>> stringlen_sortdecreasing;//(T,stringlen)
        std::vector<std::vector<double>> stringsdf_sortdecreasing;//(T,stringsdf)
        std::vector<std::vector<double>> fast_frac_sortdecreasing;//(T,fast_frac)
        std::vector<std::vector<double>> taueq_sortdecreasing_fit;
        std::vector<std::vector<double>> stringlen_sortdecreasing_fit;
        std::vector<std::vector<double>> strings_rgvsn_sortincreasing;
        std::vector<std::vector<double>> strings_rgvsnloglog_sortincreasing;
        std::vector<std::vector<double>> dwf_sortdecreasing_fit;
        std::vector<std::vector<double>> sigmaMatrix;
        
        /** Penalizaed spline interpolation in ALGLIB **/
        //----------------------------------------------------------------------
        double rho_Lt;
        
        /** GLM **/
        //----------------------------------------------------------------------
        double tau0_GLM;
        double u02;
        double alpha;
        double lntau0;
        double dfit;
        
        /** Leporini **/
        //----------------------------------------------------------------------
        double logtauref;//fixed ref timescale in ps
        double Tref;
        double u2ref;
        double Xref;
        double u2g;
        double Yg;
        double mg;
        
        alglib::real_1d_array c_COOP;
        
        /* find mean string length  */
        void find_mean_stringlen
        (const StructureClass&,const int,const int,const double,const double,const int);
        void find_mean_fast_frac
        (const StructureClass&,const int,const int,const double,const double,const int);
        
        /* find n_trajectories from amdat screen file */
        void find_n_trajs(const StructureClass&,const int,const int,const double);
        
        /** find frame_time give frame **/
        double get_frame_time(const StructureClass&,const int,const int,const double,const int);
        
        /* find and write L(T) data; f(time) */
        void find_time_stringlen(const StructureClass&,const int,const int,const double,const int);
        void find_time_fast_frac(const StructureClass&,const int,const int,const double,const int);
        
        /* find and write strings Rg data; f(time) */
        void find_time_stringsRg(const StructureClass&,const int,const int,const double,const int);
        
        /* fit L(T) vs T by COOP */
        void fit_stringlen(const StructureClass&,const int,const std::string&);
        void fit_stringlen(const StructureClass&,const int,const double,const std::string&);
        
        /* fit L(T) vs T by COOP */
        void fit_DWF(const StructureClass&,const int,const std::string&);
        void fit_DWF(const StructureClass&,const int,const double,const std::string&);
        
        /* fit tau(T) vs T by model */
        void fit_taueq(const StructureClass&,const int,const std::string&);
        void fit_taueq(const StructureClass&,const int,const double,const std::string&);
        
        /* fit Rg of stings vs mass to get fractal dimension (df) */
        void fit_all_strings_rgvsm(const StructureClass&,const int,const int,const double,const int);
        
        /* fit string model with 2 parameters */
        void fit_stringmodel(const StructureClass&,const int,const double,const std::string&);
        
        /* fit GLM model */
        void fit_GLMmodel(const StructureClass&,const int,const bool,const std::string&);
        void write_GLMref(const StructureClass&,const int,const std::string&);
        
        /* fit Hall-Wolynes model */
        void fit_HWmodel(const StructureClass&,const int,const bool,const std::string&);
        
        /* fit Leporini universal scaling */
        void fit_Leporinimodel(const StructureClass&,const int,const bool,const std::string&);
        
        /* find/write the peak string length from L(T,time) file */
        void write_time_stringlen(const StructureClass&,const int,const int,const double);
        void write_peak_stringlen(const StructureClass&,const int,const int,const double);
        void write_peak_fast_frac(const StructureClass&,const int,const int,const double);
        
        /* symmetric functions for retrieving string Rg information */
        //TODO
        void read_all_equ_strings_rg(const StructureClass&,const int);
        void read_frame_strings_rgvsn(const StructureClass&,const int,const int,const double,const int);
        void write_frame_stringsRg(const StructureClass&,const int,const int,const double,const int);
        
        /* symmetric functions for retrieving string length information */
        void read_individual_stringlen(const StructureClass&,const int,const int);
        void read_individual_equ_stringlen(const StructureClass&,const int,const int);
        void read_all_stringlen(const StructureClass&,const int,const bool is_cuttime=true,const double cuttime=1000);
        void read_all_equ_stringlen(const StructureClass&,const int);
        void read_all_equ_stringlen(const StructureClass&,const int,const double,const bool is_ltTA=true);
        void write_stringlen_equ(const StructureClass&);
        void write_stringlen_equ_avg(const StructureClass&);
        
        /* find string length at TA */
        void find_LA(const StructureClass&,const int,const double);
        
        /* find results from Arrhenius fit */
        void find_individual_Arrhenius(const StructureClass&,const int,const int);
        void find_all_Arrhenius(const StructureClass&,const int);
        
        /* find results from model fit */
        //void find_individual_tau0(const StructureClass&,const int,const int);
        //void find_all_tau0(const StructureClass&,const int);
        
        /* find tau_inf as a fitting parameter */
        //void find_AG_diagonal(const StructureClass&,const int,const std::string&);
        void find_diagonal(const StructureClass&,const int,const std::string&);
        
        /* find tau0_fit as a fitting parameter */
        void fit_GLM_diagonal(const StructureClass&,const int,const std::string&);
        
        /* fit u2(T) to a straight line */
        void fit_u2_Tr(const StructureClass&,const int,std::ofstream&,const bool,const std::string&,const std::string&);
        void fit_u2T_linear(const StructureClass&,const int,const std::string&,const std::string&);
        void read_rawxydata(const StructureClass&,const int,const std::string&);
        
        /* write universal scaling relation between <u2> and T (based on TA) */
        void write_ref_scaling_u2(const StructureClass&,const double,const double,const std::vector<double>&);
        void write_universal_scaling_u2(const StructureClass&,const std::vector<std::string>&,const std::vector<alglib::real_1d_array>&,const alglib::real_1d_array&,const std::vector<std::vector<double>>&,std::ofstream&);
        void write_universal_scaling_dfit(const StructureClass&,const std::vector<std::string>&,const alglib::real_1d_array&,const alglib::real_1d_array&,const std::vector<double>&);
        
    public:
        
        TheoryTest()=default;
        TheoryTest(const StructureClass&,const WorkScripts&,const AmdatAnalysis&,const FitData&);
        TheoryTest(const TheoryTest&)=default;
        TheoryTest& operator= (const TheoryTest&)=default;
        ~TheoryTest()=default;
        
        /* string analysis by AMDAT */
        void amdatstringanalysis(StructureClass&,const WorkScripts&,AmdatAnalysis&,const std::vector<std::vector<double>>&);
        
        /* Test the RFOT theory */
        void RFOTtest(const StructureClass&);
        
        /* Test the Adam-Gibbs (AG) theory */
        void AGtest(const StructureClass&);
        
        /* Test the Generalized Localization Model (GLM) */
        void GLMtest(const StructureClass&);
        
        /* Test the Hall-Wolynes relation */
        void HWtest(const StructureClass&);
        
        /* Test the universal scaling proposed by Leporini */
        void Leporinitest(const StructureClass&);
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_AGtest(const bool);
        void set_is_RFOTtest(const bool);
        void set_is_GLMtest(const bool);
        void set_is_Leporinitest(const bool);
        void set_is_HWtest(const bool);
        void set_is_reltaur(const bool);
        void set_is_use_counting1(const bool);
        void set_is_backwardExtrp(const bool);
        void set_is_data_smoothing(const bool);
        /* string */
        void set_stringtype(const std::string&);
        
        /** public getters **/
        //----------------------------------------------------------------------
        /* bool */
        bool get_is_AGtest() const {return is_AGtest;}
        bool get_is_RFOTtest() const {return is_RFOTtest;}
        bool get_is_GLMtest() const {return is_GLMtest;}
        bool get_is_Leporinitest() const {return is_Leporinitest;}
        bool get_is_HWtest() const {return is_HWtest;}
        bool get_is_reltaur() const {return is_reltaur;}
        /* string */
        const std::string get_stringtype() const {return stringtype;}
    };
}
#endif /* THEORYTEST_H */




