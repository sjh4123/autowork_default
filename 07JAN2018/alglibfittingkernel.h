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

#ifndef ALGLIBFITTINGKERNEL_H
#define ALGLIBFITTINGKERNEL_H

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <deque>
#include <fstream>
#include <ctime>
#include <time.h>

//ALGLIB
#include "./alglib/stdafx.h"
#include "./alglib/interpolation.h"

#define inf INFINITY
#define Epsilon 1e-10

namespace autoWork
{
    class AlglibFittingKernel
    {
        /** This is an abstract base class **/
        
    protected:
        
        bool is_use_FG;
        
        int equal_pieces;
        int cyclecount;
        int indexi;
        int indexii;
        
        double NewtonGuessTemp;
        double epsf;
        double epsx;
        double diffstep;
        double corrFacforT;
        double precision;
        double convert;
        double extrp_time;
        double compu_time;
        
        std::string systemUnit;
        
        /** ALGLIB fitting DS **/
        //----------------------------------------------------------------------
        std::vector<std::string> xdata;
        std::vector<std::string> ydata;
        std::vector<std::string> xraw;
        std::vector<std::string> yraw;
        std::vector<double> xorg;
        std::vector<double> yorg;
        std::string fit_xData;
        std::string fit_yData;
        std::vector<std::time_t> usedt;
        alglib::real_2d_array x;
        alglib::real_1d_array y;
        alglib::real_1d_array c;
        alglib::real_1d_array s;
        alglib::real_1d_array bndl;
        alglib::real_1d_array bndu;
        alglib::ae_int_t maxits;
        alglib::ae_int_t info;
        alglib::lsfitstate state;
        alglib::lsfitreport rep;
        std::string coeffs;
        std::string coeffs_scale;
        std::string coeffs_bndl;
        std::string coeffs_bndu;
        std::vector<double> coeffs_vD;
        std::vector<double> coeffs_scale_vD;
        std::vector<double> coeffs_bndl_vD;
        std::vector<double> coeffs_bndu_vD;
        alglib::spline1dinterpolant interpolant;
        alglib::barycentricinterpolant interpolantpolynom;
        
        /** STL containers **/
        //----------------------------------------------------------------------
        std::vector<std::vector<double>> taueq_sortdecreasing;    //(T,  taueq)
        std::vector<std::vector<double>> taueqinvT_sortdecreasing;//(1/T,taueq)
        std::vector<std::vector<double>> taueq_sortincreasing;    //(T,  taueq)
        std::vector<std::vector<double>> temps_sortdecreasing;    //temperatures
        std::vector<std::vector<double>> temps_sortincreasing;    //temperatures
        std::vector<std::vector<double>> taueq2d_sortdecreasing;
        std::vector<std::vector<double>> standard2dcontainer;
        
        std::deque<std::deque<double>> sorted_avg,sorted_stdev;
        
        /** Pure Virtual Functions **/
        //----------------------------------------------------------------------
        /* create alglib solver object */
        virtual void alglib_solverobj(const alglib::real_2d_array&,const alglib::real_1d_array&,const alglib::real_1d_array&,alglib::lsfitstate&)=0;
        /* lsfit form */
        virtual void alglib_lsfit_form(const std::string&,alglib::lsfitstate&)=0;
        /* fitting parameters */
        virtual void set_fitParams(const std::string&)=0;
        /* correct fitParams on the fly if needed */
        virtual void fitParams_correction(const std::string&)=0;
        /* alglib nonlinear least-square fit routine */
        virtual void alglib_lsfit(const std::string&)=0;
        
        /* fitting information */
        //----------------------------------------------------------------------
        const std::string get_fit_info(const int);
        
        /** Useful Functions **/
        //----------------------------------------------------------------------
        /* Pure Virtuals */
        virtual double calc_y_given_x(const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual std::vector<double> calc_x_given_y(const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual void write_tauFitfit_vs_T(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual void write_tauFitfit_vs_invT(std::ofstream&,const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual void build_function_interpolant(const std::string&,const std::vector<std::vector<double>>&,const int,const int)=0;
        
        /* Shared Functions */
        void fitdata_processing_lsfit(const std::vector<std::vector<double>>&);
        void fitdata_processing_lsfit_2dx(const std::vector<std::vector<double>>&);
        void fitdata_processing_1d(const std::vector<std::vector<double>>&);
        void write_relaxation_target(std::ofstream&, const std::string&);
        void write_message(std::ofstream&,const std::string&);
        void write_fitModel(std::ofstream&,const std::string&);
        void write_fitarrays(std::ofstream&);
        void write_stopcond(std::ofstream&);
        void write_fitinfo(std::ofstream&);
        void write_yvsxby2dvecdub(std::ofstream&,const std::vector<std::vector<double>>&);
        void write_errcurve(std::ofstream&,const std::string&);
        void write_errcurve(std::ofstream&,const std::string&,const std::vector<std::vector<double>>&);
        void write_tauFit_vs_T(std::ofstream&);
        void write_tauFit_vs_invT(std::ofstream&);
        void write_badFit(std::ofstream&,const int code=0);
        void write_insufficientData(std::ofstream&,const int,const int);
        
        /** calculate value of t-dependent function at given time **/
        double tfunction_at_given_time(const std::vector<std::vector<double>>&,const double,bool& flag);
        double tfunction_at_given_valu(const std::vector<std::vector<double>>&,const double,bool& flag);
        
        /* sort 2D data from large to small */
        static bool sortDecreasing0(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[0]>vd2[0];} // NOTE: sort 1st column
        static bool sortDecreasing1(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[1]>vd2[1];} // NOTE: sort 2nd column
        static bool sortDecreasing2(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[2]>vd2[2];} // NOTE: sort 3rd column
        
        /* sort data from small to large */
        static bool sortIncreasing0(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[0]<vd2[0];} // NOTE: sort 1st column
        static bool sortIncreasing1(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[1]<vd2[1];} // NOTE: sort 2nd column
        static bool sortIncreasing2(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[2]<vd2[2];} // NOTE: sort 3rd column
        
        /* sort 1D data from large to small */
        static bool sortDecreasing1D(const double d1,const double d2){return d1>d2;}
        /* sort 1D data from small to large */
        static bool sortIncreasing1D(const double d1,const double d2){return d1<d2;}
        
        void averaging_sorted_data(const std::vector<std::vector<double>>&,std::vector<std::vector<double>>&);
        void cout_c_parms();
        void cout_r2();
        
        /* setters */
        void set_coeffs(const std::vector<double>&);
        void set_coeffs_scale(const std::vector<double>&);
        void set_coeffs_bndl(const std::vector<double>&);
        void set_coeffs_bndu(const std::vector<double>&);
        
        /* getters */
        const std::string get_coeffs() const {return coeffs;}
        const std::string get_coeffs_scale() const {return coeffs_scale;}
        const std::string get_coeffs_bndl() const {return coeffs_bndl;}
        const std::string get_coeffs_bndu() const {return coeffs_bndu;}
        const std::vector<double>& get_coeffs_vD() const {return coeffs_vD;}
        const std::vector<double>& get_coeffs_scale_vD() const {return coeffs_scale_vD;}
        const std::vector<double>& get_coeffs_bndl_vD() const {return coeffs_bndl_vD;}
        const std::vector<double>& get_coeffs_bndu_vD() const {return coeffs_bndu_vD;}
        
        /** make two container consistent in 1d elements, shrink larger to smaller **/
        template <typename T>
        void make1dconsistent(std::vector<std::vector<T>>& vec1,
                              std::vector<std::vector<T>>& vec2)
        {
            bool is_v1_big=false;
            int  n_large=0;
            int  n_small=0;
            size_t sizev1=vec1.size();
            size_t sizev2=vec2.size();
            std::vector<std::vector<T>> vvl,tmpl;
            std::vector<std::vector<T>> vvs,tmps;
            if (sizev1<sizev2) {
                n_small=(int)sizev1; vvs=vec1;
                n_large=(int)sizev2; vvl=vec2;
            } else {
                is_v1_big=true;
                n_small=(int)sizev2; vvs=vec2;
                n_large=(int)sizev1; vvl=vec1;
            }
            for (int i=0;i<n_small;++i) {
                for (int ii=0;ii<n_large;++ii) {
                    if (vvs.at(i).at(0)==vvl.at(ii).at(0)) {
                        tmps.push_back({vvs.at(i).at(0),vvs.at(i).at(1)});
                        tmpl.push_back({vvl.at(ii).at(0),vvl.at(ii).at(1)});
                        break;
                    }
                }
            }
            if (is_v1_big) {
                vec1=tmpl;
                vec2=tmps;
            } else {
                vec1=tmps;
                vec2=tmpl;
            }
        }
        
    public:
        
        AlglibFittingKernel();
        AlglibFittingKernel(const AlglibFittingKernel&)=default;
        AlglibFittingKernel& operator= (const AlglibFittingKernel&)=default;
        virtual ~AlglibFittingKernel()=default;
        
        /* Public Shared Functions */
        void build_cubicspline_interpolant(const std::vector<std::vector<double>>&,const int,const int);
        void build_cubicspline_loglog_interpolant(const std::vector<std::vector<double>>&,const int,const int);
        void build_polynomial_interpolant(const std::vector<std::vector<double>>&,const int,const int,const alglib::ae_int_t);
        void build_penalizedspline_interpolant(const std::vector<std::vector<double>>&,const int,const int,const alglib::ae_int_t m=50,const double rho=2.0);
        void build_penalizedspline_loglog_interpolant(const std::vector<std::vector<double>>&,const int,const int,const alglib::ae_int_t m=50,const double rho=0.5);
        
        /** public setters **/
        //----------------------------------------------------------------------
        void set_is_use_FG(const bool);
        
        /** public getters **/
        //----------------------------------------------------------------------
        bool get_is_use_FG() const {return is_use_FG;}
        alglib::spline1dinterpolant get_interpolant() const {return interpolant;}
    };
}
#endif /* ALGLIBFITTINGKERNEL_H */




