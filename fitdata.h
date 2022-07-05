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

#ifndef FITDATA_H
#define FITDATA_H

//ALGLIB
#include "./alglib/stdafx.h"
#include "./alglib/interpolation.h"

#include "alglibfittingkernel.h"
#include "structureclass.h"
#include "workscripts.h"
#include "amdatanalysis.h"

namespace autoWork
{
    class FitData: public AlglibFittingKernel
    {
        /** PreSQ Kernel **/
        void shootForNewTemperatures(StructureClass&,const int,const int,
                                     const alglib::real_1d_array&,
                                     const std::string&);
    protected:
        
        bool is_check_amdat_version;
        bool is_applytauFitcut;
        bool is_imposeStrictTeq;
        bool is_avgFitParams;
        bool is_fitByEveryPoint;
        bool is_fit_by_TA;
        bool is_fit_by_Tc;
        bool is_fit_sExp;
        bool is_fit_lnFs;
        bool is_1stmoment_tau;
        bool is_find_DWF;
        bool is_find_NGP;
        bool is_fit_Arrhenius;
        bool is_fit_tauFit;
        bool is_shootForNext;
        bool is_use_viscosity;
        bool is_fit_Fs_by_spline;
        bool is_fit_full_alpha;
        bool is_use_sExp_cut_hi;
        bool is_use_plateautime;
        bool is_use_sExp_cut_hi_org;
        bool is_use_plateautime_org;
        bool is_use_gammafunc;
        bool is_use_KWWassist;
        bool is_calc_thermoData;
        bool is_1paramLeporini;
        bool is_fixed_tauA;
        bool is_write_derivativeT_files;
        bool is_normalModes;
        bool is_qspectrum;
        bool is_use_ngp_peak_frame;
        bool is_use_ngp_smoothed;
        bool is_binning;
        
        int index;
        int index_largest_Tg;
        int index_largest_m;
        int current_regime;
        int n_regime;
        int n_trial;
        int n_system;
        int n_sys_beg;
        int n_sys_end;
        int sExp_cut_index_hi;
        int sExp_cut_index_lo;
        int n_fit_sExp;
        int n_fit_chopping;
        int n_fit_sliding;
        int taualpha_frame;
        int ngp_peak_frame;
        int ngp_block_frame;
        
        double DWF_time;
        double wavenumber;
        double n_equ_blocks;
        double n_prd_blocks;
        double sExp_tau_lo;
        double sExp_cut_hi;
        double sExp_cut_lo;
        double sExp_tauFit;
        double tauFit_calc;
        double relativeStdev;
        double tauFit_cut;
        double shootratio;
        double largest_Tg;
        double largest_m;
        double Theq;
        double Tleq;
        double T_actual;
        double tau_T;
        double slope_fit;
        double intcp_fit;
        double r2_threshold;
        double tauA_fixed;
        double rho_arrNorm;
        double rho_tau;
        double rho_dwf;
        double rho_ngp;
        double xtauA;
        double xTA;
        double taualpha_time;
        double taualpha_frame_time;
        double frame_time;
        
        std::string amdat_svn;
        std::string amdat_version;
        std::string relaxation_target;
        std::string deviationscale;
        std::string definitionofTA;
        std::string fcorr_model;
        std::string extrp_model;
        std::string presq_model;
        std::string GLM1mode;
        std::string analysispart;
        std::string analysispart_bin;
        std::string segmode;
        
        /** Model fit **/
        double tau0_model;
        double tau0_Arrhenius;
        
        /** Arrhenius fit **/
        double deviate_ArrFit;
        double deviate_ArrNor;
        double cutTforArrhenius;
        double res_Arrhenius;
        double tauA;
        double u2A;
        double rhoA;//density*
        double Ea;
        double TA;
        double LA;
        double LT;
        double TA_avg;
        double tauA_avg;
        double tauA_Arrhenius;
        double tauA_deviation;
        
        /** MCT fit **/
        double Tc_MCT;
        double Tc_SOU;
        double Tc_percent;
        
        /** STL containers **/
        std::vector<int>    waveindices;
        std::vector<double> logtselect;
        std::vector<double> c_avgVec;
        std::vector<double> r2_Arrhenius;
        std::vector<double> r2_Model;
        std::vector<std::vector<double>> ArrheniusCoeffs;
        std::vector<std::vector<double>> ModelCoeffs;
        std::vector<std::vector<double>> ExtrpCoeffs;
        std::vector<std::vector<double>> correlation_original;        //(time,corr)
        std::vector<std::vector<double>> correlation_semilog;         //(log.time,log.corr)
        std::vector<std::vector<double>> correlation_loglog;          //(log.time,log.corr)
        std::vector<std::vector<double>> correlation_sortincreasing;  //(time,corr)
        std::vector<std::vector<double>> lncorrelation_sortincreasing;//(time,ln.corr)
        std::vector<std::vector<double>> msd_sortincreasing;          //(time,msd)
        std::vector<std::vector<double>> loglogmsd_sortincreasing;    //(log.time,log.msd)
        std::vector<std::vector<double>> semilogmsd_sortincreasing;   //(log.time,msd)
        std::vector<std::vector<double>> msdt_sortdecreasing;         //(T,msdt)
        std::vector<std::vector<double>> dwf_sortdecreasing;          //(Teq,DWF)
        std::vector<std::vector<double>> dwf_invT_sortdecreasing;     //(1/Teq,DWF)
        std::vector<std::vector<double>> dwf_tau_data;                //(DWF,tau)
        std::vector<std::vector<double>> param_sortdecreasing;
        std::vector<std::vector<double>> thermo_sortdecreasing;       //(T,thermo)
        std::vector<std::vector<double>> betaKWW_sortincreasing;      //(log.q2,beta)
        std::vector<std::vector<double>> ngp_data_raw;                //(time,ngp,frame)
        std::vector<std::vector<double>> ngp_data_raw_smoothed;       //(time,ngp,frame)
        std::vector<std::vector<double>> ngp_data_sortincreasing;     //(in-block-time,ngp,frame)
        std::vector<std::vector<double>> ngp_avg_data_sortDecreasing; //(T,in-block-time,ngp,frame)
        
        /** ALGLIB kernel (virtuals) **/
        //----------------------------------------------------------------------
        /** create alglib solver object **/
        void alglib_solverobj(const alglib::real_2d_array&,const alglib::real_1d_array&,const alglib::real_1d_array&,alglib::lsfitstate&);
        /** lsfit form **/
        void alglib_lsfit_form(const std::string&,alglib::lsfitstate&);
        /** fitting parameters **/
        void set_fitParams(const std::string&);
        /** correct fitParams on the fly if needed **/
        void fitParams_correction(const std::string&);
        /** alglib nonlinear least-square fit routine **/
        void alglib_lsfit(const std::string&);
        
        /* Virtuals */
        //----------------------------------------------------------------------
        double calc_y_given_x(const alglib::real_1d_array&,const double temp,const std::string&);
        std::vector<double> calc_x_given_y(const alglib::real_1d_array&,const double time,const std::string&);
        void write_tauFitfit_vs_T(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&);
        void write_tauFitfit_vs_invT(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&);
        void build_function_interpolant(const std::string&,const std::vector<std::vector<double>>&,const int,const int);
        
        /** functional forms for lsfit **/
        //----------------------------------------------------------------------
        /** KWW **/
        static void KWW_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void KWW_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void KWW_lnFs_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void KWW_lnFs_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** KWW_full **/
        static void KWW_full_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void KWW_full_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void KWW_full_lnFs_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void KWW_full_lnFs_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** KWW_pwr **/
        static void KWW_pwr_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void KWW_pwr_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** mKWW **/
        static void mKWW_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void mKWW_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** mKWW_pwr **/
        static void mKWW_pwr_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void mKWW_pwr_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** simple linear **/
        static void linear_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void linear_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void linear_func2(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void linear_grad2(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void linear_func3(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void linear_grad3(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** 2nd order polynomial **/
        static void secondpoly_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void secondpoly_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Exponential model for string length distribution **/
        static void exp_string_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void exp_string_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Power-exponential model for string length distribution **/
        static void pwr_exp_string_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void pwr_exp_string_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Arrhenius **/
        static void Arrhenius_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void Arrhenius_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** String model (delH and delS are fit parameters) **/
        static void string2param_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void string2param_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void string3param_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void string3param_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void transitState_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void transitState_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Leporini (with a reference) **/
        double get_Xr_leporini(const alglib::real_1d_array&,const double);
        static void leporiniref_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void leporiniref_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void leporini_universal_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void leporini_universal_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Hall-Wolynes **/
        static void hallwolynes_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void hallwolynes_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        static void hallwolynes1param_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void hallwolynes1param_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        
        /** 3-PARAMETER MODELS **/
        //----------------------------------------------------------------------
        /** MCT **/
        static void MCT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void MCT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** SOU **/
        static void SOU_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void SOU_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** VFT **/
        static void VFT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void VFT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
       
        /** Generic_MYEGA (6/27/2022) **/
        static void Generic_MYEGA_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void Generic_MYEGA_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Generic_COOP (6/27/2022) **/
        static void Generic_COOP_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void Generic_COOP_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        
        /** GLMVFT **/
        static void GLMVFT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void GLMVFT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Double Exponential **/
        static void Mauro_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void Mauro_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Avramov–Milchev (AM) **/
        static void AM_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void AM_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Dyre–Granato (DG) **/
        static void DG_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void DG_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** ArrheniusII **/
        static void ArrheniusII_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        /** Generalied Localization Model (GLM) -- 3 free parameters **/
        static void GLM_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void GLM_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Generalied Localization Model (GLM) -- 1 free parameter **/
        static void GLM1_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void GLM1_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Universal relation between tau and DWF based on an inferred localization potential **/
        static void univILP_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        /** Universal Relation from GLM -- 1 free parameter **/
        static void univGLM_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void univGLM_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** taur, <u2>r, Tr **/
        static void univGT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void univGT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** fitting <u2> to Tr with slope equal to <u2>A -- 1 free parameter **/
        static void univu2T_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void univu2T_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        
        /** 4-PARAMETER MODELS **/
        //----------------------------------------------------------------------
        /** COOP **/
        static void COOP_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void COOP_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** DEAG **/
        static void DEAG_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void DEAG_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Cohen-Grest (CG) **/
        static void CG_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void CG_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** Four-Param VFT **/
        static void FourParamVFT_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        static void FourParamVFT_grad(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,alglib::real_1d_array&,void*);
        /** ArrheniusIII **/
        static void ArrheniusIII_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        
        /** Heterogeneous Rouse Model (friction coeff version) **/
        static void HRM1_func(const alglib::real_1d_array&,const alglib::real_1d_array&,double&,void*);
        
        /** Data Processing/Retrieval/IO **/
        
        /** tau_alpha calculation from stretched exponential function **/
        //----------------------------------------------------------------------
        double sExp_tauFit_interpolateFvalue(const alglib::real_1d_array&,const std::string& tcalc);
        double sExp_tauFit_gammafunction(const alglib::real_1d_array&,const std::string& tcalc);
        
        /** calculate value of a sExp function at given time **/
        //----------------------------------------------------------------------
        double sExp_func_value(const alglib::real_1d_array&,const double time,const std::string& tcalc);
        
        /** MSD & DWF **/
        //----------------------------------------------------------------------
        void read_raw_msd_data(const StructureClass&,const int,const int,const double,const int frame=0);
        void read_individual_DWF(const StructureClass&,const int,const int);
        void read_individual_equ_DWF(const StructureClass&,const int,const int);
        void read_all_DWF(const StructureClass&,const int);
        void read_all_equ_DWF(const StructureClass&,const int);
        void read_all_equ_DWF(const StructureClass&,const int,const double);
        
        /** Non-Gaussian Parameter (NGP) **/
        //----------------------------------------------------------------------
        /* find peak frame of ngp */
        void find_ngp_peak_frame(const StructureClass&,const int,const int,const double);
        void find_avg_ngp_peak_frame(const StructureClass&);
        void write_ngp_peak_frame(const StructureClass&,const int,const int,const double);
        void write_avg_ngp_peak_frame(const StructureClass&);
        /* find tau_alpha frame */
        void find_taualpha_frame(const StructureClass&,const int,const int,const double);
        /* find frame index for a block */
        void find_ngp_block_frame(const StructureClass&,const int,const int,const double);
        
        /** get intermediate scattering function values **/
        //----------------------------------------------------------------------
        bool skim_individual_correlation_data(const StructureClass&,const int,const int,const double,const int frame=0,const int waveindex=-1,const std::string& target="");
        void read_individual_correlation_data(const StructureClass&,const int,const int,const double,const int frame=0,const int waveindex=-1,const std::string& target="");
        
        /** get fit paramters of sExp function **/
        //----------------------------------------------------------------------
        void read_all_sExp_params(const StructureClass&,const int);
        
        /** fit dynamic correlation by stretched exponential function **/
        //----------------------------------------------------------------------
        void fit_corr_by_sExp(StructureClass&,const int,const int,const double,const std::string&,bool&,const int frame=0);
        void fit_corr_by_spline(StructureClass&,const int,const int,const double,bool&,const int frame=0);
        
        /** get diffusion coefficient **/
        //----------------------------------------------------------------------
        double get_DiffCoeff(const StructureClass&,const int,const int,const double);
        
        /** fit q dependences of isfs relaxation **/
        //----------------------------------------------------------------------
        void fit_q_dependences(const StructureClass&,const int,const int,const double,const std::string&,const int);
        
        /** get T-dependence of normal modes decoupling relations **/
        //----------------------------------------------------------------------
        void fit_normal_modes_decoupling(const StructureClass&,const int,const int,const double,const std::string&);
        
        /** get taueq **/
        //----------------------------------------------------------------------
        void read_individual_taueq_data(const StructureClass&,const int,const int);
        void read_individual_taueq_data(const StructureClass&,const int,const int,const std::string& model);
        void read_all_taueq_data(const StructureClass&,const int);
        void read_all_r2teq_data(const StructureClass&,const int);
        void read_all_taueq_data(const StructureClass&,const int,const double,bool is_ltTA=true);
        void read_all_r2teq_data(const StructureClass&,const int,const double,bool is_ltTA=true);
        void read_all_taueq_data(const StructureClass&,const int,const std::string&);
        
        double get_individual_taueq_at_T(const StructureClass&,const int,const int,const double);
        
        /** read LAMMPS themo data **/
        //----------------------------------------------------------------------
        void read_all_thermo_data(const StructureClass&,const int,const std::string&);
        void read_all_thermo_data(const StructureClass&,const int,const std::string&,const double,bool is_ltTA=true);
        
        /** get Teq **/
        //----------------------------------------------------------------------
        void read_individual_temps_data(const StructureClass&,const int,const int);
        void read_all_temps_data(const StructureClass&,const int);
        
        /** find TA and related **/
        //----------------------------------------------------------------------
        void find_TA(const StructureClass&,const int,const std::vector<std::vector<double>> tauVec);
        void find_u2A(const StructureClass&,const int,const double);/* find DWF at TA */
        void find_rhoA(const StructureClass&,const int,const double);/* find density at TA */
        
        /** find transition time from ballistic to caging regimes **/
        //----------------------------------------------------------------------
        void find_balltocage_time(const StructureClass&,const int,const int,const double);
        
        /** Stickel derivative analysis **/
        //----------------------------------------------------------------------
        void stickelanalysis(const StructureClass&);
        
        /** internal utility functions **/
        //----------------------------------------------------------------------
        void write_TA(const StructureClass&,const int);
        void write_tau0(StructureClass&,const int,const std::string&);
        void write_fitcutoff(std::ofstream&,const StructureClass&,const std::string&);
        void write_spline_actual(std::ofstream&,const std::string&,const std::vector<std::vector<double>>&);
        void write_errcurve_actual(std::ofstream&,const std::string&,const std::vector<std::vector<double>>&);
        void write_errcurve_log(std::ofstream&,const std::string&,const std::vector<std::vector<double>>&);
        void write_errcurve2d_log(std::ofstream&,const std::string&,const std::vector<std::vector<double>>&);
        void write_fitavg(StructureClass&,const int,const std::string&);
        
        void write_fit_correlation(std::ofstream&,const std::string&);
        void write_fit_correlation_spline(std::ofstream&,const std::vector<std::vector<double>>&);
        void write_sExp_params(std::ofstream&,const double);
        void write_tauFitfit_vs_dwf(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&);
        void write_tauFitfit_vs_reduceddwf(std::ofstream&,const alglib::real_1d_array&,const std::vector<std::vector<double>>&,const std::vector<double>&,const double,const std::string&);
        
        void fit_every_postprocessing();
        void fit_every_preprocessing(const StructureClass&,const int,const std::string&);
        void fit_every_processing(const int,const int);
        void fit_continuousRange(const StructureClass&,const int,const std::string&);
        void fit_neighboringNpoints(const StructureClass&,const int,const std::string&);
        void fit_partial(const StructureClass&,const int,const std::string&);
        
        void save_fit_coeffs(StructureClass&);
        void Tc_correction(const std::string&);
        void compute_Tg_fragility(StructureClass&,const alglib::real_1d_array&,const std::string&);
        void compute_VFT_Tg_fragility(StructureClass&,const alglib::real_1d_array&);
        void compute_COOP_Tg_fragility(StructureClass&,const alglib::real_1d_array&);
        void find_largest_Tg_fragility(StructureClass&,const int,const int,const alglib::real_1d_array&,const std::string&);
        void avgfitParams(const StructureClass&,const int,const int,const std::string&);
        void avgfitParams_extrp(const StructureClass&,const int,const int,const std::string& extrp);
        void find_cutTforArrhenius(const StructureClass&,const int);
        void use_ArrFit();
        void use_ArrNor();
        void use_MsdCag();
        void use_FsqCag();
        void find_tau0(StructureClass&,const int,const std::string&);
        
        bool load_fit_coeffs(const StructureClass&);
        
        const std::vector<double>& get_Theq(StructureClass&,const int,const int);
        const std::vector<double>& get_Theq(StructureClass&,const int);
        const std::vector<double>& get_Tleq(StructureClass&,const int,const int);
        const std::vector<double>& get_Tleq(StructureClass&,const int);
        double calc_time_given_vecinput(const alglib::real_1d_array&,const std::vector<double>& v,const std::string&);
        double error_Tg(const std::string&);
        double error_fragility(const std::string&,const double Tg_value,const double Tg_error);
        double write_Tc(const StructureClass&,const int);
        
        std::vector<double> write_Tg_fragility(std::ofstream&,const alglib::real_1d_array&,const std::string&);
        std::vector<double> write_Yg_fragility(std::ofstream&,const alglib::real_1d_array&,const std::vector<double>&,const std::string&);
        std::vector<double> calc_Y_given_time(const alglib::real_1d_array&,const double time,const std::vector<double>& ref,const std::string&);
        std::vector<double> get_interpolated_Ts(StructureClass&,const int,const int,const int n_temps,const double Tl,const double Th);
        
        
        void find_MinMaxbySlopeSignChange(const StructureClass&,const std::vector<std::vector<double>>&,std::vector<double>&);
        void find_MinMaxbySlopeInterplation(const StructureClass&,const std::vector<std::vector<double>>&,std::vector<double>&);
        
        
        /** inline **/
        //----------------------------------------------------------------------
        double get_2pt_slope(const std::vector<double>& xy0,const std::vector<double>& xy1){return (xy1[1]-xy0[1])/(xy1[0]-xy0[0]);}
        
    public:
        
        FitData()=default;
        FitData(const StructureClass&,const WorkScripts&,const AmdatAnalysis&);
        FitData(const FitData&)=default;
        FitData& operator= (const FitData&)=default;
        ~FitData()=default;
        
        /** write quenched temperatures to file **/
        void write_qchTs(StructureClass&,const int,const int);
        
        /** write in-equilibrium temperatures to file **/
        void write_equTs(StructureClass&,const int,const int);
        
        /** check if current regime needs a retry based on # of in-equ. data **/
        bool check_is_retry(StructureClass&);
        
        /** fit relaxation data to a stretched exponetial function **/
        void fit_sExp(StructureClass&,const int,const int,const double T,const std::string&,const int frame=0);
        void write_sExp_params_avg(const StructureClass&,const int,const std::string&);
        
        /** MSD & DWF **/
        void write_DWF(const StructureClass&,const int,const int,const double,const int frame=0);
        void write_raw_msd_data(const StructureClass&,const int,const int,const double);
        void write_diffusion_coeff(const StructureClass&,const int,const int,const double);
        void write_MSDt_equ(const StructureClass&);
        void write_MSDt_equ_avg(const StructureClass&);
        void write_DWF_equ(const StructureClass&);
        void write_DWF_equ_avg(const StructureClass&);
        void write_loglogMSD_derivatives(const StructureClass&,const int,const int,const double);
        
        /** NGP **/
        void write_NGP(const StructureClass&,const int,const int,const double);
        void write_avg_NGP(const StructureClass&);
        
        /** Normal modes decoupling test **/
        void write_avg_normal_modes_decoupling(const StructureClass&);
        void write_avg_qmatch(const StructureClass&);
        void write_avg_qmatch_KWW(const StructureClass&);
        void write_smo_qmatch_KWW(const StructureClass&);
        void write_bKWW_vs_T(const StructureClass&);
        
        /** fit (highT) relaxation data Arrhenius relation **/
        void fit_Arrhenius(StructureClass&,const int,const int);
        
        /** fit relaxation data to a chosen model **/
        void fit_tauFit(StructureClass&,const int,const int,const std::string&,const std::string&);
        
        /** invoke independent ALGLIB fit routine **/
        //---------------------------------------------
        void alglibfit_routine(const StructureClass&);
        //---------------------------------------------
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_fit_sExp(const bool);
        void set_is_fit_lnFs(const bool);
        void set_is_1stmoment_tau(const bool);
        void set_is_find_DWF(const bool);
        void set_is_find_NGP(const bool);
        void set_is_fit_Arrhenius(const bool);
        void set_is_fit_tauFit(const bool);
        void set_is_shootForNext(const bool);
        void set_is_avgFitParams(const bool);
        void set_is_fit_by_TA(const bool);
        void set_is_fitByEveryPoint(const bool);
        void set_is_fit_by_Tc(const bool);
        void set_is_applytauFitcut(const bool);
        void set_is_fit_Fs_by_spline(const bool);
        void set_is_fit_full_alpha(const bool);
        void set_is_use_gammafunc(const bool);
        void set_is_imposeStrictTeq(const bool);
        void set_is_calc_thermoData(const bool);
        void set_is_use_KWWassist(const bool);
        void set_is_normalModes(const bool);
        void set_is_qspectrum(const bool);
        void set_is_use_ngp_peak_frame(const bool);
        void set_is_use_ngp_smoothed(const bool);
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
        void set_Tc_MCT(const double);
        void set_Tc_SOU(const double);
        void set_tauFit_cut(const double);
        void set_xtauA(const double);
        void set_xTA(const double);
        /* string */
        void set_relaxation_target(const std::string&);
        void set_definitionofTA(const std::string&);
        void set_fcorr_model(const std::string&);
        void set_extrp_model(const std::string&);
        void set_presq_model(const std::string&);
        /* STL */
        void set_waveindices(const std::vector<int>&);
        void set_logtselect(const std::vector<double>&);
        
        /** public getters **/
        //----------------------------------------------------------------------
        /* bool */
        bool get_is_fit_sExp() const {return is_fit_sExp;}
        bool get_is_find_DWF() const {return is_find_DWF;}
        bool get_is_find_NGP() const {return is_find_NGP;}
        bool get_is_fit_Arrhenius() const {return is_fit_Arrhenius;}
        bool get_is_fit_tauFit() const {return is_fit_tauFit;}
        bool get_is_shootForNext() const {return is_shootForNext;}
        bool get_is_avgFitParams() const {return is_avgFitParams;}
        bool get_is_fit_by_TA() const {return is_fit_by_TA;}
        bool get_is_fitByEveryPoint() const {return is_fitByEveryPoint;}
        bool get_is_fit_by_Tc() const {return is_fit_by_Tc;}
        bool get_is_applytauFitcut() const {return is_applytauFitcut;}
        bool get_is_fit_Fs_by_spline() const {return is_fit_Fs_by_spline;}
        bool get_is_fit_full_alpha() const {return is_fit_full_alpha;}
        bool get_is_use_gammafunc() const {return is_use_gammafunc;}
        bool get_is_calc_thermoData() const {return is_calc_thermoData;}
        bool get_is_use_KWWassist() const {return is_use_KWWassist;}
        bool get_is_normalModes() const {return is_normalModes;}
        bool get_is_qspectrum() const {return is_qspectrum;}
        bool get_is_fit_lnFs() const {return is_fit_lnFs;}
        bool get_is_1stmoment_tau() const {return is_1stmoment_tau;}
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
        double get_Tc_MCT() const {return Tc_MCT;}
        double get_Tc_SOU() const {return Tc_SOU;}
        double get_xtauA() const {return xtauA;}
        double get_xTA() const {return xTA;}
        /* string */
        const std::string get_definitionofTA() const {return definitionofTA;}
        const std::string get_relaxation_target() const {return relaxation_target;}
        const std::string get_fcorr_model() const {return fcorr_model;}
        const std::string get_extrp_model() const {return extrp_model;}
        const std::string get_presq_model() const {return presq_model;}
        /* STL */
        const std::vector<int>& get_waveindices() const {return waveindices;}
        const std::vector<double>& get_logtselect() const {return logtselect;}
    };
}
#endif /* FITDATA_H */




