//
//  alglibfittingkernel.h
//  cppWork
//
//  Created by SJH on 11/11/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#ifndef ALGLIBFITTINGKERNEL_H
#define ALGLIBFITTINGKERNEL_H

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <time.h>

#include "stdafx.h"
#include "interpolation.h"
#include "alglibfittingkernel.h"

#define inf INFINITY
#define epsilon 1e-10

namespace autoWork {
    
    class AlglibFittingKernel
    {
        
    protected:
        
        int indexi,indexii;
        
        int cyclecount;
        std::vector<std::string> xdata,ydata;
        std::vector<std::string> xraw,yraw;
        std::vector<double> xorg,yorg;
        std::string fit_xData,fit_yData;
        std::vector<std::time_t> usedt;
        
        alglib::real_2d_array x;
        alglib::real_1d_array y;
        alglib::real_1d_array c;
        alglib::real_1d_array s;
        alglib::real_1d_array bndl;
        alglib::real_1d_array bndu;
        double epsf;
        double epsx;
        double diffstep;
        alglib::ae_int_t maxits;
        alglib::ae_int_t info;
        alglib::lsfitstate state;
        alglib::lsfitreport rep;
        
        std::string systemUnit;
        double corrFacforT;
        double precision;
        double extrp_time;
        double compu_time;
        double NewtonGuessTemp;
        
        bool is_use_FG; // true: use FG-form; false: use F-form
        
        std::string coeffs;
        std::string coeffs_scale;
        std::string coeffs_bndl;
        std::string coeffs_bndu;
        
        std::vector<double> coeffs_vD;
        std::vector<double> coeffs_scale_vD;
        std::vector<double> coeffs_bndl_vD;
        std::vector<double> coeffs_bndu_vD;
        
        std::vector<std::vector<double>> taueq_sortdecreasing; // (Teq,taueq)
        std::vector<std::vector<double>> taueq_sortincreasing; // (Teq,taueq)
        std::vector<std::vector<double>> temps_sortdecreasing; // temperatures
        std::vector<std::vector<double>> temps_sortincreasing; // temperatures
        
        /* Constructors */
        AlglibFittingKernel();
        AlglibFittingKernel(const AlglibFittingKernel&)=default;
        /* Assignment */
        AlglibFittingKernel& operator= (const AlglibFittingKernel&)=default;
        /* Destructor */
        ~AlglibFittingKernel()=default;
        
        
        /* Pure Virtuals */
        //----------------------------------------------------------------------
        /* create alglib solver object */
        virtual void alglib_solverobj(const alglib::real_2d_array& x,
                                      const alglib::real_1d_array& y,
                                      const alglib::real_1d_array& c,
                                      alglib::lsfitstate& state)=0;
        /* lsfit form */
        virtual void alglib_lsfit_form(const std::string&,
                                       alglib::lsfitstate& state)=0;
        /* fitting parameters */
        virtual void set_fitParams(const std::string& model)=0;
        /* correct fitParams on the fly if needed */
        virtual void fitParams_correction(const std::string& model)=0;
        /* alglib nonlinear least-square fit routine */
        virtual void alglib_lsfit(const std::string& model)=0;
        //----------------------------------------------------------------------
        
        /* fitting information */
        std::string get_fit_info(const int);
        
        
        /* Useful Functions */
        //----------------------------------------------------------------------
        
        /* Pure Virtuals */
        virtual double calc_time_given_temp(const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual std::vector<double> calc_temp_given_time(const alglib::real_1d_array&,const double,const std::string&)=0;
        virtual void write_tauFitfit_vs_T(std::ofstream&,const alglib::real_1d_array&,const std::string&)=0;
        virtual void write_tauFitfit_vs_invT(std::ofstream&,const alglib::real_1d_array&,const std::string&)=0;
        
        /* Shared Functions */
        void fitdata_processing_lsfit(const std::vector<std::vector<double>>&);
        void fitdata_processing_1d(const std::vector<std::vector<double>>&);
        void write_relaxation_target(std::ofstream&, const std::string&);
        void write_fitModel(std::ofstream& outputFile,const std::string& model);
        void write_fitarrays(std::ofstream& outputFile);
        void write_stopcond(std::ofstream& outputFile);
        void write_fitinfo(std::ofstream& outputFile);
        void write_errcurve(std::ofstream& outputFile,const std::string& model);
        void write_tauFit_vs_T(std::ofstream& outputFile);
        void write_tauFit_vs_invT(std::ofstream& outputFile);
        void write_badR2(std::ofstream& outputFile,const double r2_threshold);
        void write_insufficientDara(std::ofstream& outputFile,const int n_fit,const int n_threshold);        
        
        /* data sorted from large to small */
        static bool sortDecreasing(const std::vector<double>& vd1,
                                   const std::vector<double>& vd2)
        {return vd1[0]>vd2[0];} // NOTE: sort 1st column
        
        static bool sortDecreasing1(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[1]>vd2[1];} // NOTE: sort 2nd column
        
        static bool sortDecreasing2(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[2]>vd2[2];} // NOTE: sort 3rd column
        
        /* data sorted from small to large */
        static bool sortIncreasing(const std::vector<double>& vd1,
                                   const std::vector<double>& vd2)
        {return vd1[0]<vd2[0];} // NOTE: sort 1st column
        
        static bool sortIncreasing1(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[1]<vd2[1];} // NOTE: sort 2nd column
        
        static bool sortIncreasing2(const std::vector<double>& vd1,
                                    const std::vector<double>& vd2)
        {return vd1[2]<vd2[2];} // NOTE: sort 3rd column
        
        void averaging_sorted_data(const std::vector<std::vector<double>>&,
                                   std::vector<std::vector<double>>&);
        void cout_c_parms();
        void cout_r2();
        
        /* setters */
        void set_coeffs(const std::vector<double>&);
        void set_coeffs_scale(const std::vector<double>&);
        void set_coeffs_bndl(const std::vector<double>&);
        void set_coeffs_bndu(const std::vector<double>&);
        
        /* getters */
        const std::vector<double>& get_coeffs_vD() const {return coeffs_vD;}
        const std::vector<double>& get_coeffs_scale_vD() const {return coeffs_scale_vD;}
        const std::vector<double>& get_coeffs_bndl_vD() const {return coeffs_bndl_vD;}
        const std::vector<double>& get_coeffs_bndu_vD() const {return coeffs_bndu_vD;}
        const std::string get_coeffs() const {return coeffs;}
        const std::string get_coeffs_scale() const {return coeffs_scale;}
        const std::string get_coeffs_bndl() const {return coeffs_bndl;}
        const std::string get_coeffs_bndu() const {return coeffs_bndu;}
        
    public:
        
        /* public setters */
        void set_is_use_FG(const bool);
        
        /* public getters */
        bool get_is_use_FG() const {return is_use_FG;}
        
    };
}
#endif /* ALGLIBFITTINGKERNEL_H */
