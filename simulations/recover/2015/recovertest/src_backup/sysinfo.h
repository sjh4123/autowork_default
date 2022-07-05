//
//  Created by Sean Jh H. on 6/3/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef SYSINFO_H
#define SYSINFO_H

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <time.h>
#include "stdafx.h"
#include "interpolation.h"

namespace autoWork {
    
    class SysInfo {
        
    protected:
        
        std::string Path;
        std::string simType;
        std::string year;
        std::string usicID;
        std::string systemUnit;
        std::string sys_targets;
        
        bool is_makeFileSystem;
        bool is_defaultLmpData;
        bool is_moltempLmpData;
        bool is_LammpsPhase;
        bool is_Simulations;
        bool is_AMDAT;
        bool is_fitData;
        bool is_directSub;
        bool is_watch_hold;
        bool is_GPU;            // enables the GPU features
        bool is_fullquench;     // quench from highest to lowest or by regime
        bool is_singleTempTest; // single temperature test: gen-->prd-->amdat
        bool is_useDynamicRange;
        bool is_use_viscosity;
        bool is_set_CheckPoint;
        bool is_updatetinfo;
        bool is_theorytest;
        bool is_fit_dielectric;
        
        int indexi,indexii,indexiii;
        int n_digits;
        int n_trial;
        int n_system;
        int n_sys_beg;   // start regime (inclusive)
        int n_sys_end;   // final regime (inclusive)
        int n_startTemp;
        int n_Temp;
        int regime_beg;
        int regime_end;
        int current_regime;
        int n_regime;
        int quenchMode;
        
        double n_equ_blocks;
        double n_prd_blocks;
        
        double startTemp;    // default = 1000.0
        double crossTemp;    // default = 700.0
        double finalTemp;    // default = 200.0
        double hTres;        // default = 75.0
        double lTres;        // default = 50.0
        double corrFacforT;  // default = 1.0
        double precision;    // default = 100.0
        double cutTforArrhenius;
        
        double dynamicRange_tau0;
        double dynamicRange_tauA;
        
        std::vector<int> n_regime_temps;
        
        std::vector<double> ts_regime;
        std::vector<double> equilibration_times;
        std::vector<double> equilibration_times_org;
        
        
        /* temperature containers */
        //------------------------------------------------
        std::vector<std::vector<double>> initialtemps;
        std::vector<std::vector<double>> temperaturesInfo;
        std::vector<std::vector<double>> equilibratedTs;
        std::vector<std::vector<double>> quenchingTs;
        std::vector<double> qchTdiff;
        std::vector<std::vector<std::vector<double>>> tauFit; //[trial][T][tau]
        std::vector<std::vector<std::vector<double>>> tauEqu; //[trial][T][tau]
        //------------------------------------------------
        
        
        
        /* variables used in the fianl reports */
        //------------------------------------------------
        std::string startdatetime;
        std::string relaxationModel;
        std::string extrpolateModel;
        std::string message_dynRange;
        double extrp_time;
        double compu_time;
        double timestep_size;
        double quenchRate;
        double timediff_overall;
        double wavenumber;
        std::vector<double> elapsedtime;
        std::vector<std::vector<double>> Tg_extrp;
        std::vector<std::vector<double>> Tg_compu;
        std::vector<std::vector<double>> m_extrp;
        std::vector<std::vector<double>> m_compu;
        
        bool is_node0_info,is_node1_info;
        std::vector<int> node0_info,node1_info;
        std::vector<std::vector<int>> simulation_retry; //{regime,times retried}
        //------------------------------------------------
        
        
        /* variables for simple stats */
        //------------------------------------------------
        int max;
        double mean;
        std::vector<double> value;
        //------------------------------------------------
        
        
    public:
        
        /* Constructors */
        SysInfo();
        SysInfo(const SysInfo&)=default;
        /* Assignment */
        SysInfo& operator= (const SysInfo&)=default;
        /* Destructor */
        ~SysInfo()=default;
        
        /*==( public setters )==*/
        
        /* string */
        void set_Path(const std::string&);
        void set_simType(const std::string&);
        void set_year(const std::string&);
        void set_usicID(const std::string&);
        void set_systemUnit(const std::string&);
        void set_sys_targets(const std::string&);
        
        /* bool */
        void set_is_makeFileSystem(const bool);
        void set_is_defaultLmpData(const bool);
        void set_is_moltempLmpData(const bool);
        void set_is_LammpsPhase(const bool);
        void set_is_Simulations(const bool);
        void set_is_AMDAT(const bool);
        void set_is_fitData(const bool);
        void set_is_directSub(const bool);
        void set_is_watch_hold(const bool);
        void set_is_GPU(const bool);
        void set_is_fullquench(const bool);
        void set_is_singleTempTest(const bool);
        void set_is_useDynamicRange(const bool);
        void set_is_use_viscosity(const bool);
        void set_is_set_CheckPoint(const bool);
        void set_is_updatetinfo(const bool);
        void set_is_theorytest(const bool);
        void set_is_fit_dielectric(const bool);
        
        /* int */
        void set_n_digits(const int);
        void set_n_trial(const int);
        void set_n_system(const int);
        void set_n_sys_beg(const int);
        void set_n_sys_end(const int);
        void set_n_startTemp(const int);
        void set_n_Temp(const int);
        void set_regime_beg(const int);
        void set_regime_end(const int);
        void set_current_regime(const int);
        void set_n_regime(const int);
        void set_quenchMode(const int);
        
        /* double */
        void set_n_equ_blocks(const double);
        void set_n_prd_blocks(const double);
        void set_startTemp(const double);
        void set_crossTemp(const double);
        void set_finalTemp(const double);
        void set_hTres(const double);
        void set_lTres(const double);
        void set_corrFacforT(const double);
        void set_precision(const double);
        void set_cutTforArrhenius(const double);
        void set_dynamicRange_tau0(const double);
        void set_dynamicRange_tauA(const double);
        
        /* vector<int> */
        void set_n_regime_temps(const std::vector<int>&);
        
        /* vector<double> */
        void set_ts_regime(const std::vector<double>&);
        void set_equilibration_times(const std::vector<double>&);
        void set_equilibration_times_org(const std::vector<double>&);
        
        /* vector<vector<double>> */
        void set_initialtemps(const std::vector<std::vector<double>>&);
        
        // temperature containers (no use for setters, use getters)
        //----------------------------------------------------------------------
        //void set_temperaturesInfo(const std::vector<std::vector<double>>&);
        //void set_equilibratedTs(const std::vector<std::vector<double>>&);
        //void set_quenchingTs(const std::vector<std::vector<double>>&);
        //void set_qchTdiff(const std::vector<double>&);
        //----------------------------------------------------------------------
        
        /* variables used in fitData for fit KWW */
        //----------------------------------------------------------------------
        std::vector<double> coeffs_bndl_vD;
        std::vector<double> coeffs_bndu_vD;
        std::string coeffs_bndl;
        std::string coeffs_bndu;
        std::vector<double> c_vD;
        alglib::lsfitreport rep;
        double epsf;
        double epsx;
        double diffstep;
        alglib::ae_int_t maxits;
        //----------------------------------------------------------------------
        
        /* variables used in the fianl reports */
        //----------------------------------------------------------------------
        /* string */
        void set_startdatetime(const std::string&);
        void set_relaxationModel(const std::string&);
        void set_extrpolateModel(const std::string&);
        void set_message_dynRange(const std::string&);
        /* double */
        void set_extrp_time(const double);
        void set_compu_time(const double);
        void set_timestep_size(const double);
        void set_quenchRate(const double);
        void set_timediff_overall(const double);
        void set_wavenumber(const double);
        /* bool */
        void set_is_node0_info(const bool);
        void set_is_node1_info(const bool);
        /* vector<int> */
        void set_node0_info(const std::vector<int>&);
        void set_node1_info(const std::vector<int>&);
        //----------------------------------------------------------------------
        
        
        /*==( public getters )==*/
        
        /* string */
        const std::string get_Path() const {return Path;}
        const std::string get_simType() const {return simType;}
        const std::string get_usicID() const {return usicID;}
        const std::string get_year() const;
        const std::string get_usic() const;
        const std::string get_systemUnit() const {return systemUnit;}
        const std::string get_sys_targets() const {return sys_targets;}
        
        /* bool */
        bool get_is_makeFileSystem() const {return is_makeFileSystem;}
        bool get_is_defaultLmpData() const {return is_defaultLmpData;}
        bool get_is_moltempLmpData() const {return is_moltempLmpData;}
        bool get_is_LammpsPhase() const {return is_LammpsPhase;}
        bool get_is_Simulations() const {return is_Simulations;}
        bool get_is_AMDAT() const {return is_AMDAT;}
        bool get_is_fitData() const {return is_fitData;}
        bool get_is_directSub() const {return is_directSub;}
        bool get_is_watch_hold() const {return is_watch_hold;}
        bool get_is_GPU() const {return is_GPU;}
        bool get_is_fullquench() const {return is_fullquench;}
        bool get_is_singleTempTest() const {return is_singleTempTest;}
        bool get_is_useDynamicRange() const {return is_useDynamicRange;}
        bool get_is_use_viscosity() const {return is_use_viscosity;}
        bool get_is_set_CheckPoint() const {return is_set_CheckPoint;}
        bool get_is_updatetinfo() const {return is_updatetinfo;}
        bool get_is_theorytest() const {return is_theorytest;}
        bool get_is_fit_dielectric() const {return is_fit_dielectric;}
        
        /* int */
        int get_n_digits() const {return n_digits;}
        int get_n_trial() const {return n_trial;}
        int get_n_system() const;
        int get_n_sys_beg() const {return n_sys_beg;}
        int get_n_sys_end() const {return n_sys_end;}
        int get_n_startTemp() const {return n_startTemp;}
        int get_n_Temp() const {return n_Temp;}
        int get_regime_beg() const {return regime_beg;}
        int get_regime_end() const {return regime_end;}
        int get_current_regime() const {return current_regime;}
        int get_n_regime() const {return n_regime;}
        int get_quenchMode() const {return quenchMode;}
        
        /* double */
        double get_n_equ_blocks() const {return n_equ_blocks;}
        double get_n_prd_blocks() const {return n_prd_blocks;}
        double get_startTemp() const {return startTemp;}
        double get_crossTemp() const {return crossTemp;}
        double get_finalTemp() const {return finalTemp;}
        double get_hTres() const {return hTres;}
        double get_lTres() const {return lTres;}
        double get_corrFacforT() const {return corrFacforT;}
        double get_precision() const {return precision;}
        double get_cutTforArrhenius() const {return cutTforArrhenius;}
        double get_dynamicRange_tau0() const {return dynamicRange_tau0;}
        double get_dynamicRange_tauA() const {return dynamicRange_tauA;}
        
        
        /* vector<int> */
        const std::vector<int>& get_n_regime_temps() const
        {return n_regime_temps;}
        
        /* vector<double> */
        const std::vector<double>& get_ts_regime() const
        {return ts_regime;}
        const std::vector<double>& get_equilibration_times() const
        {return equilibration_times;};
        const std::vector<double>& get_equilibration_times_org() const
        {return equilibration_times_org;}
        
        
        const std::vector<std::vector<double>>& get_initialtemps() const
        {return initialtemps;}
        
        // temperature containers
        // vector<vector<double>>& -- direct reference to memory of object
        //----------------------------------------------------------------------
        std::vector<std::vector<double>>& get_temperaturesInfo()
        {return temperaturesInfo;}
        
        std::vector<std::vector<double>>& get_equilibratedTs()
        {return equilibratedTs;}
        
        std::vector<std::vector<double>>& get_quenchingTs()
        {return quenchingTs;}
        
        std::vector<double>& get_qchTdiff()
        {return qchTdiff;}
        
        std::vector<std::vector<std::vector<double>>>& get_tauFit()
        {return tauFit;}
        std::vector<std::vector<std::vector<double>>>& get_tauEqu()
        {return tauEqu;}
        //----------------------------------------------------------------------
        
        
        /* variables used in the fianl reports */
        //----------------------------------------------------------------------
        /* string */
        const std::string get_startdatetime() const {return startdatetime;}
        const std::string get_relaxationModel() const {return relaxationModel;}
        const std::string get_extrpolateModel() const {return extrpolateModel;}
        const std::string get_message_dynRange() const {return message_dynRange;}
        /* double */
        double get_extrp_time() const {return extrp_time;}
        double get_compu_time() const {return compu_time;}
        double get_timestep_size() const {return timestep_size;}
        double get_quenchRate() const {return quenchRate;}
        double get_timediff_overall() const {return timediff_overall;}
        double get_wavenumber() const {return wavenumber;}
        
        /* vector<vector<double>>& -- direct reference to memory of object */
        std::vector<double>& get_elapsedtime(){return elapsedtime;}
        
        /* vector<double>& -- direct reference to memory of object */
        std::vector<std::vector<double>>& get_Tg_extrp() {return Tg_extrp;}
        std::vector<std::vector<double>>& get_Tg_compu() {return Tg_compu;}
        std::vector<std::vector<double>>& get_m_extrp()  {return m_extrp;}
        std::vector<std::vector<double>>& get_m_compu()  {return m_compu;}
        
        bool get_is_node0_info(){return is_node0_info;}
        bool get_is_node1_info(){return is_node1_info;}
        const std::vector<int>& get_node0_info(){return node0_info;}
        const std::vector<int>& get_node1_info(){return node1_info;}
        
        int times_retry;
        int times_write_tempsinfo;
        int max_times_retry;
        std::vector<std::vector<int>>& get_simulation_retry(){return simulation_retry;}
        //----------------------------------------------------------------------
        
        
        /* simple stats */
        //----------------------------------------------------------------------
        void set_calcvector(const std::vector<double>&);
        double calc_mean();
        double calc_variance();
        double calc_sampleVariance();
        double get_stdev() {return sqrt(calc_variance());}
        double get_sample_stdev(){return sqrt(calc_sampleVariance());}
        //----------------------------------------------------------------------
    };
}

#endif