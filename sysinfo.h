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

#ifndef SYSINFO_H
#define SYSINFO_H

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <random>

//ALGLIB
#include "./alglib/stdafx.h"
#include "./alglib/interpolation.h"

namespace autoWork
{
    class SysInfo
    {
        
    protected:
        
        bool is_makeFileSystem;
        bool is_defaultLmpData;
        bool is_moltempLmpData;
        bool is_LammpsPhase;
        bool is_use_prepScripts;
        bool is_Simulations;
        bool is_AMDAT;
        bool is_fitData;
        bool is_directSub;
        bool is_watch_hold;
        bool is_GPU;
        bool is_fullquench;
        bool is_singleTempTest;
        bool is_use_viscosity;
        bool is_set_CheckPoint;
        bool is_updatetinfo;
        bool is_theorytest;
        bool is_fit_dielectric;
        bool is_specify_GPU;
        bool is_cutoff_traject;
        bool is_alglibfit;
        bool is_use_named_nameString;
        bool is_makeNewAnanlysisFolder;
        bool is_backupfitdata;
        bool is_backupstatistics;
        bool is_return;
        bool is_amdatinp;
        bool is_copy_from_to;
        bool is_aging;
        bool is_cancelRetryAlgo;
        
        int indexi;
        int indexii;
        int indexiii;
        int n_digits;
        int n_trial;
        int n_system;
        int n_sys_beg;
        int n_sys_end;
        int n_startTemp;
        int n_Temp;
        int regime_beg;
        int regime_end;
        int current_regime;
        int n_regime;
        int times_retry;
        int max_times_retry;
        int n_cum_taueq;
        int which_GPU;
        int sim_restart;
        int gen_restart;
        int qch_restart;
        int equ_restart;
        int res_restart;
        int prd_restart;
        int ana_restart;
        int fd_res;
        int tt_res;
        int max_initialized_regime;
        int n_kuhnsegs;
        
        double n_equ_blocks;
        double n_prd_blocks;
        double n_relaxation;
        double startTemp;
        double crossTemp;
        double finalTemp;
        double hTres;
        double lTres;
        double corrFacforT;
        double precision;
        double cutTforArrhenius;
        double dynamicRange_tau0;
        double dynamicRange_tauA;
        
        std::string userName;
        std::string masterFolder;
        std::string Path;
        std::string simType;
        std::string year;
        std::string usicID;
        std::string systemUnit;
        std::string sys_targets;
        std::string computeCluster;
        std::string copy_source;
        std::string copy_destin;
        
        /** STL containers **/
        //----------------------------------------------------------------------
        std::vector<int>    n_regime_temps;
        std::vector<double> ts_regime;
        std::vector<double> equilibration_times;
        std::vector<double> equilibration_times_org;
        std::vector<double> n_equ_blocks_plus;
        std::vector<double> n_prd_blocks_plus;
        std::vector<double> n_relaxation_plus;
        std::vector<double> tautargets;
        std::vector<std::vector<double>> extrpTmax;
        
        /** simple stats **/
        //----------------------------------------------------------------------
        int    max;
        double mean;
        std::vector<double> value;
        
        /** variables used in fitData for KWW fit **/
        //----------------------------------------------------------------------
        double epsf;
        double epsx;
        double diffstep;
        std::vector<double> coeffs_bndl_vD;
        std::vector<double> coeffs_bndu_vD;
        std::string coeffs_bndl;
        std::string coeffs_bndu;
        std::vector<double> c_vD;
        alglib::lsfitreport rep;
        alglib::ae_int_t maxits;
        
        /** variables used in fianl report **/
        //----------------------------------------------------------------------
        bool   is_node0_info;
        bool   is_node1_info;
        bool   is_node2_info;
        double extrp_time;
        double compu_time;
        double timestep_size;
        double quenchRate;
        double timediff_overall;
        double wavenumber;
        std::string startdatetime;
        std::string message_dynRange;
        std::vector<double> elapsedtime;
        std::vector<std::string> extrp_model;
        std::vector<std::string> presq_model;
        std::vector<std::vector<double>> Tg_extrp;
        std::vector<std::vector<double>> Tg_compu;
        std::vector<std::vector<double>> m_extrp;
        std::vector<std::vector<double>> m_compu;
        std::vector<std::vector<double>> Tg_VFT;
        std::vector<std::vector<double>> Tg_COOP;
        std::vector<std::vector<double>> m_VFT;
        std::vector<std::vector<double>> m_COOP;
        std::vector<int> node0_info;
        std::vector<int> node1_info;
        std::vector<int> node2_info;
        std::vector<std::vector<int>> simulation_retry; //{regime,times retried}
        
        /** temperature containers **/
        //----------------------------------------------------------------------
        std::vector<std::vector<double>> initialtemps;
        std::vector<std::vector<double>> temperaturesInfo;
        std::vector<std::vector<double>> equilibratedTs;
        std::vector<std::vector<double>> quenchingTs;
        std::vector<double> qchTdiff;
        std::vector<std::vector<std::vector<double>>> tauFit;  //[trial][T][tau]
        std::vector<std::vector<std::vector<double>>> tauEqu;  //[trial][T][tau]
        std::vector<std::vector<std::vector<double>>> inptinfo;//[regime][trial][T]
        
    public:
        
        SysInfo();
        SysInfo(const SysInfo&)=default;
        SysInfo& operator= (const SysInfo&)=default;
        ~SysInfo()=default;
        
        /** simple stats **/
        void   set_calcvector(const std::vector<double>&);
        double calc_mean();
        double calc_variance();
        double calc_sampleVariance();
        double get_stdev() {return sqrt(calc_variance());}
        double get_sample_stdev(){return sqrt(calc_sampleVariance());}
        
        /** variables used in fitData for KWW fit **/
        /** setters **/
        void set_coeffs_bndl_vD(const std::vector<double>&);
        void set_coeffs_bndu_vD(const std::vector<double>&);
        void set_coeffs_bndl(const std::string&);
        void set_coeffs_bndu(const std::string&);
        void set_c_vD(const std::vector<double>&);
        void set_rep(const alglib::lsfitreport&);
        void set_epsf(const double);
        void set_epsx(const double);
        void set_diffstep(const double);
        void set_maxits(const alglib::ae_int_t&);
        /** getters **/
        std::vector<double> get_coeffs_bndl_vD() const {return coeffs_bndl_vD;}
        std::vector<double> get_coeffs_bndu_vD() const {return coeffs_bndu_vD;}
        std::string get_coeffs_bndl() const {return coeffs_bndl;}
        std::string get_coeffs_bndu() const {return coeffs_bndu;}
        std::vector<double> get_c_vD() const {return c_vD;}
        alglib::lsfitreport get_rep() const {return rep;}
        double get_epsf() const {return epsf;}
        double get_epsx() const {return epsx;}
        double get_diffstep() const {return diffstep;}
        alglib::ae_int_t get_maxits() const {return maxits;}
        
        /** variables used in fianl report **/
        /* int */
        void set_times_retry(const int);
        void set_max_times_retry(const int);
        void set_n_cum_taueq(const int);
        void set_which_GPU(const int);
        void add_n_cum_taueq(const int);
        int  get_times_retry() const {return times_retry;}
        int  get_max_times_retry() const {return max_times_retry;}
        int  get_n_cum_taueq() const {return n_cum_taueq;}
        int  get_which_GPU() const {return which_GPU;}
        /* bool */
        void set_is_node0_info(const bool);
        void set_is_node1_info(const bool);
        void set_is_node2_info(const bool);
        bool get_is_node0_info() const {return is_node0_info;}
        bool get_is_node1_info() const {return is_node1_info;}
        bool get_is_node2_info() const {return is_node2_info;}
        /* double */
        void   set_extrp_time(const double);
        void   set_compu_time(const double);
        void   set_timestep_size(const double);
        void   set_quenchRate(const double);
        void   set_timediff_overall(const double);
        void   set_wavenumber(const double);
        double get_extrp_time() const {return extrp_time;}
        double get_compu_time() const {return compu_time;}
        double get_timestep_size() const {return timestep_size;}
        double get_quenchRate() const {return quenchRate;}
        double get_timediff_overall() const {return timediff_overall;}
        double get_wavenumber() const {return wavenumber;}
        /* string */
        void  set_startdatetime(const std::string&);
        void  set_extrp_model(const std::string&);
        void  set_presq_model(const std::string&);
        void  set_message_dynRange(const std::string&);
        const std::string get_startdatetime() const {return startdatetime;}
        const std::string get_message_dynRange() const {return message_dynRange;}
        /* vector<int> */
        void set_node0_info(const std::vector<int>&);
        void set_node1_info(const std::vector<int>&);
        void set_node2_info(const std::vector<int>&);
        const std::vector<int>& get_node0_info() const {return node0_info;}
        const std::vector<int>& get_node1_info() const {return node1_info;}
        const std::vector<int>& get_node2_info() const {return node2_info;}
        /** direct reference to memory of object **/
        std::vector<std::vector<int>>& get_simulation_retry()
        {return simulation_retry;}
        
        /** direct reference to memory of object (special data containers) **/
        std::vector<double>& get_qchTdiff()             {return qchTdiff;}
        std::vector<double>& get_elapsedtime()          {return elapsedtime;}
        std::vector<std::string>& get_extrp_model() {return extrp_model;}
        std::vector<std::string>& get_presq_model() {return presq_model;}
        std::vector<std::vector<double>>& get_Tg_extrp(){return Tg_extrp;}
        std::vector<std::vector<double>>& get_Tg_compu(){return Tg_compu;}
        std::vector<std::vector<double>>& get_m_extrp() {return m_extrp;}
        std::vector<std::vector<double>>& get_m_compu() {return m_compu;}
        std::vector<std::vector<double>>& get_Tg_VFT()  {return Tg_VFT;}
        std::vector<std::vector<double>>& get_Tg_COOP() {return Tg_COOP;}
        std::vector<std::vector<double>>& get_m_VFT()   {return m_VFT;}
        std::vector<std::vector<double>>& get_m_COOP()  {return m_COOP;}
        std::vector<std::vector<double>>& get_temperaturesInfo()
        {return temperaturesInfo;}
        std::vector<std::vector<double>>& get_equilibratedTs()
        {return equilibratedTs;}
        std::vector<std::vector<double>>& get_quenchingTs()
        {return quenchingTs;}
        std::vector<std::vector<double>>& get_extrpTmax()
        {return extrpTmax;}
        std::vector<std::vector<std::vector<double>>>& get_tauFit()
        {return tauFit;}
        std::vector<std::vector<std::vector<double>>>& get_tauEqu()
        {return tauEqu;}
        
        /** const reference to memory of object **/
        void clear_inptinfo(){inptinfo.clear();}
        void set_inptinfo(const std::vector<std::vector<std::vector<double>>>&);
        void set_regime_tinfo(const int,const std::vector<std::vector<double>>&);
        void show_inptinfo();
        void show_regime_tinfo(const int);
        const std::vector<std::vector<std::vector<double>>>& get_inptinfo() const
        {return inptinfo;}
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_makeFileSystem(const bool);
        void set_is_defaultLmpData(const bool);
        void set_is_moltempLmpData(const bool);
        void set_is_LammpsPhase(const bool);
        void set_is_use_prepScripts(const bool);
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
        void set_is_cutoff_traject(const bool);
        void set_is_alglibfit(const bool);
        void set_is_use_named_nameString(const bool);
        void set_is_makeNewAnanlysisFolder(const bool);
        void set_is_backupfitdata(const bool);
        void set_is_backupstatistics(const bool);
        void set_is_return(const bool);
        void set_is_amdatinp(const bool);
        void set_is_copy_from_to(const bool);
        void set_is_aging(const bool);
        void set_is_cancelRetryAlgo(const bool);
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
        void set_sim_restart(const int);
        void set_gen_restart(const int);
        void set_qch_restart(const int);
        void set_equ_restart(const int);
        void set_res_restart(const int);
        void set_prd_restart(const int);
        void set_ana_restart(const int);
        void set_fd_res(const int);
        void set_tt_res(const int);
        void set_max_initialized_regime(const int);
        void set_n_kuhnsegs(const int);
        /* double */
        void set_n_equ_blocks(const double);
        void set_n_prd_blocks(const double);
        void set_n_relaxation(const double);
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
        /* string */
        void set_userName(const std::string&);
        void set_masterFolder(const std::string&);
        void set_Path(const std::string&);
        void set_simType(const std::string&);
        void set_year(const std::string&);
        void set_usicID(const std::string&);
        void set_systemUnit(const std::string&);
        void set_sys_targets(const std::string&);
        void set_computeCluster(const std::string&);
        void set_copy_source(const std::string&);
        void set_copy_destin(const std::string&);
        /* vector<int> */
        void set_n_regime_temps(const std::vector<int>&);
        /* vector<double> */
        void set_n_equ_blocks_plus(const std::vector<double>&);
        void set_n_prd_blocks_plus(const std::vector<double>&);
        void set_n_relaxation_plus(const std::vector<double>&);
        void set_tautargets(const std::vector<double>&);
        void set_ts_regime(const std::vector<double>&);
        void set_equilibration_times(const std::vector<double>&);
        void set_equilibration_times_org(const std::vector<double>&);
        /* vector<vector<double>> */
        void set_initialtemps(const std::vector<std::vector<double>>&);
        
        
        /** public getters **/
        //----------------------------------------------------------------------
        /* bool */
        bool get_is_makeFileSystem() const {return is_makeFileSystem;}
        bool get_is_defaultLmpData() const {return is_defaultLmpData;}
        bool get_is_moltempLmpData() const {return is_moltempLmpData;}
        bool get_is_LammpsPhase() const {return is_LammpsPhase;}
        bool get_is_use_prepScripts() const {return is_use_prepScripts;}
        bool get_is_Simulations() const {return is_Simulations;}
        bool get_is_AMDAT() const {return is_AMDAT;}
        bool get_is_fitData() const {return is_fitData;}
        bool get_is_directSub() const {return is_directSub;}
        bool get_is_watch_hold() const {return is_watch_hold;}
        bool get_is_GPU() const {return is_GPU;}
        bool get_is_fullquench() const {return is_fullquench;}
        bool get_is_singleTempTest() const {return is_singleTempTest;}
        bool get_is_use_viscosity() const {return is_use_viscosity;}
        bool get_is_set_CheckPoint() const {return is_set_CheckPoint;}
        bool get_is_updatetinfo() const {return is_updatetinfo;}
        bool get_is_theorytest() const {return is_theorytest;}
        bool get_is_fit_dielectric() const {return is_fit_dielectric;}
        bool get_is_specify_GPU() const {return is_specify_GPU;}
        bool get_is_cutoff_traject() const {return is_cutoff_traject;}
        bool get_is_alglibfit() const {return is_alglibfit;}
        bool get_is_use_named_nameString() const {return is_use_named_nameString;}
        bool get_is_makeNewAnanlysisFolder() const {return is_makeNewAnanlysisFolder;}
        bool get_is_backupfitdata() const {return is_backupfitdata;}
        bool get_is_backupstatistics() const {return is_backupstatistics;}
        bool get_is_return() const {return is_return;}
        bool get_is_amdatinp() const {return is_amdatinp;}
        bool get_is_copy_from_to() const {return is_copy_from_to;}
        bool get_is_aging() const {return is_aging;}
        bool get_is_cancelRetryAlgo() const {return is_cancelRetryAlgo;}
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
        int get_sim_restart() const {return sim_restart;}
        int get_gen_restart() const {return gen_restart;}
        int get_qch_restart() const {return qch_restart;}
        int get_equ_restart() const {return equ_restart;}
        int get_res_restart() const {return res_restart;}
        int get_prd_restart() const {return prd_restart;}
        int get_ana_restart() const {return ana_restart;}
        int get_fd_res() const {return fd_res;}
        int get_tt_res() const {return tt_res;}
        int get_max_initialized_regime() const {return max_initialized_regime;}
        int get_n_kuhnsegs() const {return n_kuhnsegs;}
        /* double */
        double get_n_equ_blocks() const {return n_equ_blocks;}
        double get_n_prd_blocks() const {return n_prd_blocks;}
        double get_n_relaxation() const {return n_relaxation;}
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
        /* string */
        const std::string get_userName() const {return userName;}
        const std::string get_masterFolder() const {return masterFolder;}
        const std::string get_Path() const {return Path;}
        const std::string get_simType() const {return simType;}
        const std::string get_usicID() const {return usicID;}
        const std::string get_year() const;
        const std::string get_usic() const;
        const std::string get_systemUnit() const {return systemUnit;}
        const std::string get_sys_targets() const {return sys_targets;}
        const std::string get_computeCluster() const {return computeCluster;}
        const std::string get_copy_source() const {return copy_source;}
        const std::string get_copy_destin() const {return copy_destin;}
        /* STL */
        const std::vector<double>& get_n_equ_blocks_plus() const {return n_equ_blocks_plus;}
        const std::vector<double>& get_n_prd_blocks_plus() const {return n_prd_blocks_plus;}
        const std::vector<double>& get_n_relaxation_plus() const {return n_relaxation_plus;}
        const std::vector<double>& get_tautargets() const {return tautargets;}
        const std::vector<int>& get_n_regime_temps() const
        {return n_regime_temps;}
        const std::vector<double>& get_ts_regime() const
        {return ts_regime;}
        const std::vector<double>& get_equilibration_times() const
        {return equilibration_times;};
        const std::vector<double>& get_equilibration_times_org() const
        {return equilibration_times_org;}
        const std::vector<std::vector<double>>& get_initialtemps() const
        {return initialtemps;}
    };
}
#endif /* SYSINFO_H */




