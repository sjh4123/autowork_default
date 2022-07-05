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

#include "fitdata.h"
#include "functions.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

FitData::FitData(const StructureClass& sysVar,
                 const WorkScripts& ws,
                 const AmdatAnalysis& aa):
/* Base class */
AlglibFittingKernel(),

/* bool */
is_check_amdat_version(aa.get_is_check_amdat_version()),
is_applytauFitcut(true),
is_imposeStrictTeq(false),
is_avgFitParams(false),// true: shoot by avg; false: by trial
is_fitByEveryPoint(true),
is_fit_by_TA(false),// true: fit only T<TA in final regime
is_fit_by_Tc(false),// true: fit only T<Tc in final regime
is_fit_sExp(false),
is_fit_lnFs(false),
is_1stmoment_tau(false),
is_find_DWF(false),
is_find_NGP(false),
is_fit_Arrhenius(false),
is_fit_tauFit(false),
is_shootForNext(true),// true: extrapolate new temperatures
is_use_viscosity(sysVar.get_is_use_viscosity()),
is_fit_Fs_by_spline(false),
is_fit_full_alpha(false),
is_use_sExp_cut_hi(true),//NOTE
is_use_plateautime(false),//NOTE
is_use_sExp_cut_hi_org(is_use_sExp_cut_hi),
is_use_plateautime_org(is_use_plateautime),
is_use_gammafunc(false),
is_use_KWWassist(false),
is_calc_thermoData(false),
is_1paramLeporini(false),
is_fixed_tauA(false),
is_write_derivativeT_files(false),
is_qspectrum(false),
is_use_ngp_peak_frame(false),
is_use_ngp_smoothed(true),

/* int */
index(0),
index_largest_Tg(0),
index_largest_m(0),
current_regime(sysVar.get_current_regime()),
n_regime(sysVar.get_n_regime()),
n_trial(sysVar.get_n_trial()),
n_system(sysVar.get_n_system()),
n_sys_beg(sysVar.get_n_sys_beg()),
n_sys_end(sysVar.get_n_sys_end()),
sExp_cut_index_hi(0),
sExp_cut_index_lo(0),
n_fit_sExp(5),
n_fit_chopping(7),
n_fit_sliding(10),
taualpha_frame(0),
ngp_peak_frame(0),
ngp_block_frame(0),

/* double */
DWF_time(aa.get_DWF_time()),
wavenumber(1.0),
n_equ_blocks(ws.get_n_equ_blocks()),
n_prd_blocks(ws.get_n_prd_blocks()),
sExp_tau_lo(1.0),//NOTE
sExp_cut_hi(0.6),//NOTE
sExp_cut_lo(0.01),
sExp_tauFit(0.2), //at this value tau_alpha is defined
tauFit_calc(0),
relativeStdev(0),
tauFit_cut(sysVar.get_compu_time()),//see also set_is_applytauFitcut()
shootratio(0),
largest_Tg(0),
largest_m(0),
Theq(0),
Tleq(0),
T_actual(0),
tau_T(0),
slope_fit(0),
r2_threshold(0.95),
tauA_fixed(0),//log10-based
rho_arrNorm(2.0),
rho_tau(2.0),
rho_dwf(2.0),
rho_ngp(0.0),//NOTE
xtauA(1.0),
xTA(1.0),
taualpha_time(0),
taualpha_frame_time(0),
frame_time(0),

/* string */
amdat_svn(aa.get_amdat_svn()),
amdat_version(""),
analysispart(aa.get_analysispart()),
relaxation_target(aa.get_relaxation_target()),
deviationscale("raw"),
definitionofTA("ArrNor"),
fcorr_model("KWW"),
extrp_model("COOP"),
presq_model("VFT"),
GLM1mode("mode1"),

/** Model fit **/
//------------------------------------------------------------------------------
tau0_model(0),
tau0_Arrhenius(0),

/** Arrhenius fit **/
//------------------------------------------------------------------------------
deviate_ArrFit(0.1),
deviate_ArrNor(0.15),
cutTforArrhenius(sysVar.get_cutTforArrhenius()),//cutoff T in Arrhenius fit
res_Arrhenius(1.1),
tauA(0),
u2A(0),
rhoA(0),
Ea(0),
TA(0),
LA(0),
LT(0),
TA_avg(0),
tauA_avg(0),
tauA_Arrhenius(0),
tauA_deviation(0),

/** MCT fit **/
//------------------------------------------------------------------------------
Tc_MCT(0),
Tc_SOU(0),
Tc_percent(0.99),

/** STL containers **/
waveindices(),
logtselect(),
c_avgVec(),
r2_Arrhenius(),
r2_Model(),
ArrheniusCoeffs(),
ModelCoeffs(),
ExtrpCoeffs(),
correlation_original(),
correlation_semilog(),
correlation_loglog(),
correlation_sortincreasing(),
msd_sortincreasing(),
loglogmsd_sortincreasing(),
semilogmsd_sortincreasing(),
msdt_sortdecreasing(),
dwf_sortdecreasing(),
dwf_invT_sortdecreasing(),
dwf_tau_data(),
param_sortdecreasing(),
thermo_sortdecreasing(),
betaKWW_sortincreasing(),
ngp_data_raw(),
ngp_data_sortincreasing(),
ngp_avg_data_sortDecreasing()
{
    if (sysVar.get_systemUnit()=="real") {
        sExp_tau_lo *= pow(10,3);//fs
    } else if (sysVar.get_systemUnit()=="lj") {
        sExp_tau_lo *= pow(10,0);//tau
    } else if (sysVar.get_systemUnit()=="metal") {
        sExp_tau_lo *= pow(10,0);//ps
    }
    
    /** assign values to members in parent class (AlglibFittingKernel) **/
    systemUnit       = sysVar.get_systemUnit();
    precision        = sysVar.get_precision();
    corrFacforT      = sysVar.get_corrFacforT();
    convert          = 1000.0/(corrFacforT/precision);
    extrp_time       = sysVar.get_extrp_time();
    compu_time       = sysVar.get_compu_time();
    NewtonGuessTemp *= precision; // need precision first be assigned
    
    //cout << "Log.tb("<<sysVar.get_usic()<<") "<<log10(DWF_time)<<"\n";
}





void FitData::fit_sExp(StructureClass& sysVar,
                       const int n_trl,
                       const int n_sys,
                       const double Temp_d,
                       const string& tcalc,
                       const int frame)
{
    int index=n_trl*n_system+(n_sys-n_sys_beg);
    
    /** output file paths **/
    //--------------------------------------------------------------------------
    /** tauFit stores relaxation time of all simulated temperatures **/
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tauFit_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tauFit(o.c_str(),ofstream::app);//'append'
    //--------------------------------------------------------------------------
    
    double T_actual=Temp_d*pow(corrFacforT,-1);//NOTE: convert!
    double log_tauFit=0;
    
    bool is_badfit=false;
    if (is_fit_Fs_by_spline) {
        fit_corr_by_spline(sysVar,n_trl,n_sys,Temp_d,is_badfit,frame);
        if (is_badfit) return;
        tauFit_calc=spline1dcalc(interpolant,get_sExp_tauFit());
    } else {
        fit_corr_by_sExp(sysVar,n_trl,n_sys,Temp_d,tcalc,is_badfit,frame);
        if (is_badfit) return;
        if (is_use_gammafunc) tauFit_calc=sExp_tauFit_gammafunction(c,tcalc);
        else tauFit_calc=sExp_tauFit_interpolateFvalue(c,tcalc);
    }
    
    /** store "log(tauFit)" data into container **/
    log_tauFit=log10(tauFit_calc);
    write_tauFit << T_actual << " " << log_tauFit << "\n";
    sysVar.get_tauFit().at(index).push_back({Temp_d,log_tauFit});
    
    /** store "taueq" data into container **/
    // in-equilibrium T's are written to file in write_equTs()
    if (tauFit_calc<tauFit_cut) {
        sysVar.get_equilibratedTs().at(index).push_back(Temp_d);
        sysVar.get_tauEqu().at(index).push_back({Temp_d,log_tauFit});
    } write_tauFit.close();
}





void FitData::fit_corr_by_spline(StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d,
                                 bool& is_badfit,
                                 const int frame)
{
    /** output file paths **/
    //--------------------------------------------------------------------------
    /** tauFit stores relaxation time of all simulated temperatures **/
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tauFit_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tauFit(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/sExp_params_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_params(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_sExp_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_sExp(o.c_str());
    //--------------------------------------------------------------------------
    
    /** NOTE: only fit KWW when data is full relaxed **/
    bool is_fullyRelaxed=false;
    
    is_fullyRelaxed=
    skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
    if (is_fullyRelaxed) {
        read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
    } else {
        is_badfit=true;
        write_badFit(write_sExp,1); return;
    }
    
    /** record q value used **/
    sysVar.set_wavenumber(wavenumber);
    
    /** process data to ALGLIB readable form **/
    fitdata_processing_lsfit(correlation_sortincreasing);
    
    /** interpolation and record fit results **/
    build_penalizedspline_interpolant(correlation_sortincreasing,1,0);//NOTE: time vs F
    write_relaxation_target(write_sExp,relaxation_target);
    write_fitModel(write_sExp,"penalizedspline");
    write_fitarrays(write_sExp);
    write_fit_correlation_spline(write_sExp);
    
    write_tauFit.close();
    write_params.close();
    write_sExp.close();
}





void FitData::fit_corr_by_sExp(StructureClass& sysVar,
                               const int n_trl,
                               const int n_sys,
                               const double Temp_d,
                               const string& tcalc,
                               bool& is_badfit,
                               const int frame)
{
    int index=n_trl*n_system+(n_sys-n_sys_beg);
    
    /** output file paths **/
    //--------------------------------------------------------------------------
    /** tauFit stores relaxation time of all simulated temperatures **/
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tauFit_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tauFit(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/sExp_params_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_params(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_sExp_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_sExp(o.c_str());
    //--------------------------------------------------------------------------
    
    /** setup initial values of fitting parameters **/
    if (is_use_KWWassist)
    {
        /** if r2_pre >= r2_threshold, use coeffs from previous fit;
         ** otherwise, use default coeffs  **/
        if (sysVar.get_rep().r2<r2_threshold) {
            set_fitParams(tcalc);
        } else if (sysVar.get_rep().r2>=r2_threshold &&
                   sysVar.get_rep().r2<=1.0) {
            /** if unsuccessfully loading in coeffs, use default **/
            if(!load_fit_coeffs(sysVar)) set_fitParams(tcalc);
        } else {
            set_fitParams(tcalc);
        }
    } else {
        /** use default coeffs **/
        set_fitParams(tcalc);
    }
    
    /** read relaxation data for fit **/
    //--------------------------------------------------------------------------
    /** NOTE: only fit KWW when data is full relaxed **/
    bool is_fullyRelaxed=false;
    
    if (is_fit_full_alpha) { /* fit full alpha relaxation */
        
        is_fullyRelaxed=
        skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
        if (is_fullyRelaxed) {
            read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
        } else {
            is_badfit=true;
            write_badFit(write_sExp,1); return;
        }
    } else { /* data interpolated in sExp_cut_lo < corr < sExp_cut_hi */
        
        is_fullyRelaxed=
        skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
        if (is_fullyRelaxed) {
            /* default switch: cutFsbyValue */
            is_use_sExp_cut_hi=true;
            is_use_plateautime=false;
            read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame);
        } else {
            is_badfit=true;
            write_badFit(write_sExp,1); return;
        }
    }
    //--------------------------------------------------------------------------
    
    /** save q value used **/
    sysVar.set_wavenumber(wavenumber);
    
    /** process data to ALGLIB readable form **/
    if (is_fit_lnFs) {
        fitdata_processing_lsfit(lncorrelation_sortincreasing);
    } else {
        fitdata_processing_lsfit(correlation_sortincreasing);
    }
    
    /** curve-fitting and record fit results **/
    if (correlation_sortincreasing.size()>=n_fit_sExp)
    {
        alglib_lsfit(tcalc);
        write_relaxation_target(write_sExp,relaxation_target);
        write_fitModel(write_sExp,tcalc);
        write_fitarrays(write_sExp);
        write_fitcutoff(write_sExp,sysVar,tcalc);
        write_stopcond(write_sExp);
        write_fitinfo(write_sExp);
        write_fit_correlation(write_sExp,tcalc);
        write_sExp_params(write_params,Temp_d);
    }
    else
    {
        string i;
        i.append(relaxation_target+"_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append("_T"+to_string((long long int)Temp_d));
        cout
        << "in FitData::fit_sExp():\n"
        << i << "\n"
        << "correlation data too few ("<<correlation_sortincreasing.size()<<"), "
        << "less than the threshold number ("<<n_fit_sExp<<")\n";
        system_wait(5);
    }
    
    if (rep.r2>r2_threshold)
    {
        save_fit_coeffs(sysVar);
    }
    
    /** handling poor sExp fits **/
    //--------------------------------------------------------------------------
    if ((rep.r2<r2_threshold)||(correlation_sortincreasing.size()<n_fit_sExp))
    {
        is_badfit=true;
        double faketau=log10(tauFit_cut)+1;//NOTE
        write_tauFit << Temp_d*pow(corrFacforT,-1) <<" "<< faketau << "\n";
        sysVar.get_tauFit().at(index).push_back({Temp_d,faketau});//use expanded form
        if (rep.r2<r2_threshold)
        {
            write_badFit(write_sExp);
        }
        if (correlation_sortincreasing.size()<n_fit_sExp)
        {
            write_insufficientData
            (write_sExp,(int)correlation_sortincreasing.size(),n_fit_sExp);
        }
    }
    
    /** handling normal mode dynamics analyses **/
    //--------------------------------------------------------------------------
    if (is_normalModes)
    {
        if (relaxation_target=="isfs")
        {
            if (is_qspectrum) { /** fit q dependences of isfs **/
                fit_q_dependences(sysVar,n_trl,n_sys,Temp_d,tcalc,frame);
            } else {
                fit_normal_modes_decoupling(sysVar,n_trl,n_sys,Temp_d,tcalc);
            }
        }
    }
    //--------------------------------------------------------------------------
    write_tauFit.close();
    write_params.close();
    write_sExp.close();
}





void FitData::fit_normal_modes_decoupling(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d,
                                          const string& tcalc)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/decoupling_vs_T_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("."+analysispart+".dat");
    ofstream write_coupling(o.c_str(),ofstream::app);//'append'
    
    bool   is_fullyRelaxed=false;
    bool   is_qspectrum_org=is_qspectrum;
    double tau_baf=0;
    double tau_isfs_std=0;
    
    /** find tau_baf **/
    //--------------------------------------------------------------------------
    if (false)
    {
        is_qspectrum=false;//NOTE
        is_fullyRelaxed=
        skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"baf");
        if (is_fullyRelaxed) {
            read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"baf");
        } is_qspectrum=is_qspectrum_org;//NOTE
        
        if (is_fit_lnFs) {
            fitdata_processing_lsfit(lncorrelation_sortincreasing);
        } else {
            fitdata_processing_lsfit(correlation_sortincreasing);
        }
        set_fitParams(tcalc);
        alglib_lsfit(tcalc);
        
        if (rep.r2>r2_threshold)
        {
            tau_baf=log10(sExp_tauFit_interpolateFvalue(c,tcalc));
        }
    }
    //--------------------------------------------------------------------------
    
    //static int fp=0;
    //if (fp==0) {
    //    write_coupling << "T 1/T D(p=1) D(p) Log(t.isfs) Log(t.baf) \n";
    //} ++fp;
    
    /* find tau_isfs & write out decoupling to file */
    //--------------------------------------------------------------------------
    is_qspectrum=false;//NOTE
    is_fullyRelaxed=
    skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"isfs");
    if (is_fullyRelaxed) {
        read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"isfs");
    } is_qspectrum=is_qspectrum_org;//NOTE
    
    if (is_fit_lnFs) {
        fitdata_processing_lsfit(lncorrelation_sortincreasing);
    } else {
        fitdata_processing_lsfit(correlation_sortincreasing);
    }
    set_fitParams(tcalc);
    alglib_lsfit(tcalc);
    
    if (rep.r2>r2_threshold)
    {
        string analysispart_org=analysispart;//NOTE: backup
        /** get diffusion coeff from mode1 **/
        analysispart="mode1";
        double dc_mode1=get_DiffCoeff(sysVar,n_trl,n_sys,Temp_d);
        analysispart=analysispart_org;//NOTE: resume value
        
        tau_isfs_std=log10(sExp_tauFit_interpolateFvalue(c,tcalc));
        
        write_coupling
        << Temp_d*pow(corrFacforT,-1)         << " "
        << pow(Temp_d*pow(corrFacforT,-1),-1) << " "
        << dc_mode1 << " "
        << get_DiffCoeff(sysVar,n_trl,n_sys,Temp_d) << " "
        << tau_isfs_std << " "
        << tau_baf      << " "
        << "\n";
    }
    //--------------------------------------------------------------------------
}





void FitData::fit_q_dependences(const StructureClass& sysVar,
                                const int n_trl,
                                const int n_sys,
                                const double Temp_d,
                                const string& tcalc,
                                const int frame)
{
    double f=0,df=0,d2f=0;
    string o;//output file path
    
    /** write raw isfs at all tested q to a single file **/
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_isfs_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_isfs(o.c_str(),ofstream::app);//'append'
    
    //write_isfs << "T q q2 Log(t) Fs Fs(psp) DFs D2Fs \n";
    
    for (int waveind=waveindices.at(0);waveind<=waveindices.at(1);++waveind)
    {
        read_individual_correlation_data
        (sysVar,n_trl,n_sys,Temp_d,frame,waveind);
        
        build_penalizedspline_interpolant(correlation_semilog,0,1);
        
        for (int i=0;i<(int)correlation_semilog.size();++i)
        {
            double Lt=correlation_semilog.at(i).at(0);//log(time)
            double Fs=correlation_semilog.at(i).at(1);//raw Fs
            spline1ddiff(interpolant,Lt,f,df,d2f);
            
            if (sysVar.get_systemUnit()=="real") {
                write_isfs
                << Temp_d*pow(corrFacforT,-1) << " "
                << wavenumber << " "
                << pow(wavenumber,2) << " "
                << Lt-3 << " "
                << Fs   << " "
                << f    << " "
                << df   << " "
                << d2f  << "\n";
            } else {
                write_isfs
                << Temp_d*pow(corrFacforT,-1) << " "
                << wavenumber << " "
                << pow(wavenumber,2) << " "
                << Lt << " "
                << Fs   << " "
                << f    << " "
                << df   << " "
                << d2f  << "\n";
            }
        } write_isfs << "\n"; //write_isfs.close();
    }
    //--------------------------------------------------------------------------
    
    
    /** write isfs vs q at different instances of observation **/
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_isfs_vs_q_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_isfsq(o.c_str(),ofstream::app);//'append'
    
    //write_isfsq << "T Log(t) t q q2 Fs(psp) \n";
    
    for (int ts=0;ts<(int)logtselect.size();++ts)
    {
        for (int waveind=waveindices.at(0);waveind<=waveindices.at(1);++waveind)
        {
            read_individual_correlation_data
            (sysVar,n_trl,n_sys,Temp_d,frame,waveind);
            
            build_penalizedspline_interpolant(correlation_semilog,0,1);
            
            double logts=logtselect.at(ts);
            
            if (sysVar.get_systemUnit()=="real") {
                write_isfsq
                << Temp_d*pow(corrFacforT,-1) << " "
                << logts-3 << " "
                << pow(10,logts-3) << " "
                << wavenumber << " "
                << pow(wavenumber,2) << " "
                << spline1dcalc(interpolant,logts) << " "
                << "\n"; //write_isfsq.close();
            } else {
                write_isfsq
                << Temp_d*pow(corrFacforT,-1) << " "
                << logts << " "
                << pow(10,logts) << " "
                << wavenumber << " "
                << pow(wavenumber,2) << " "
                << spline1dcalc(interpolant,logts) << " "
                << "\n"; //write_isfsq.close();
            }
        } write_isfsq << "\n";
    }
    //--------------------------------------------------------------------------
    
    
    /** write relaxation time & fit KWW params at each q to file **/
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_tau_vs_q_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_tauq(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_sExp_vs_q_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_sExp(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_betaKWW_vs_T_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //ofstream write_beta_all(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_sExp_vs_q_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_sExp_all(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_bKWW_vs_q2_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_bKWW_all(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_q_match_vs_T_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_q_match(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_q_match_KWW_vs_T_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_q_KWW(o.c_str(),ofstream::app);//'append'
    
    betaKWW_sortincreasing.clear();//used to contain betaKWW for interpolation
    
    /** find tau_baf **/
    //--------------------------------------------------------------------------
    bool   is_fullyRelaxed=false;
    bool   is_qspectrum_org=is_qspectrum;
    double tau_baf=0;
    
    alglib::real_1d_array c_baf;
    alglib::lsfitreport rep_baf;
    
    if (analysispart!="mode12")
    {
        is_qspectrum=false;//NOTE
        
        is_fullyRelaxed=
        skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"baf");
        if (is_fullyRelaxed) {
            read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,0,-1,"baf");
        } is_qspectrum=is_qspectrum_org;//NOTE
        
        if (is_fit_lnFs) {
            fitdata_processing_lsfit(lncorrelation_sortincreasing);
        } else {
            fitdata_processing_lsfit(correlation_sortincreasing);
        }
        set_fitParams(tcalc);
        alglib_lsfit(tcalc);
        
        c_baf=c;
        rep_baf=rep;
        
        if (rep.r2>r2_threshold)
        {
            tau_baf=log10(sExp_tauFit_interpolateFvalue(c,tcalc));
        }
    }
    //--------------------------------------------------------------------------
    
    /** tau_isfs used for storing tau_isfs of all q's **/
    vector<vector<double>> bKWWvsq;
    vector<vector<double>> tau_isfs;//{(q,tau_isfs)}
    vector<double> qfound;
    vector<real_1d_array> cfound;
    vector<lsfitreport> repfound;
    
    double sExp_tau_lo_org=sExp_tau_lo;//NOTE
    double sExp_cut_hi_org=sExp_cut_hi;//NOTE
    
    /** apply data truncation **/
    //---------------------------------------------
    /** cutoff time **/
    sExp_tau_lo=1.0;
    //sExp_tau_lo=10.0;
    if (sysVar.get_systemUnit()=="real") {
        sExp_tau_lo *= pow(10,3);//fs
    } else if (sysVar.get_systemUnit()=="lj") {
        sExp_tau_lo *= pow(10,0);//tau
    } else if (sysVar.get_systemUnit()=="metal") {
        sExp_tau_lo *= pow(10,0);//ps
    }
    
    for (int waveind=waveindices.at(0);waveind<=waveindices.at(1);++waveind)
    {
        /** cutoff value **/
        if (analysispart=="mode1") {
            sExp_cut_hi=0.6;
        } else if (analysispart=="mode2") {
            sExp_cut_hi=0.55;
        } else if (analysispart=="mode3") {
            sExp_cut_hi=0.5;
        } else if (analysispart=="mode4") {
            sExp_cut_hi=0.45;
        } else if (analysispart=="mode6") {
            sExp_cut_hi=0.4;
        } else if (analysispart=="mode12") {
            sExp_cut_hi=0.3;
        } else {
            sExp_cut_hi=sExp_cut_hi_org;
        }
        int repcount=0;
        do //dynamic decision on lowering sExp_cut_hi based on r2 (if <0.99)
        {
            ++repcount;
            if (repcount>1) {
                sExp_cut_hi -= 0.05;
            }
            
            /** NOTE: only fit KWW when data is full relaxed **/
            bool is_fullyRelaxed=
            skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame,waveind);
            if (!is_fullyRelaxed) continue;
            
            read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d,frame,waveind);
            
            if (is_fit_lnFs) {
                fitdata_processing_lsfit(lncorrelation_sortincreasing);
            } else {
                fitdata_processing_lsfit(correlation_sortincreasing);
            }
            
            /** KWW assisit **/
            if (sysVar.get_rep().r2<r2_threshold) {
                set_fitParams(tcalc);
            } else if (sysVar.get_rep().r2>=r2_threshold &&
                       sysVar.get_rep().r2<=1.0) {
                /** if unsuccessfully loading in coeffs, use default **/
                if(!load_fit_coeffs(sysVar)) set_fitParams(tcalc);
            } else {
                set_fitParams(tcalc);
            }
            /** KWW fit routine **/
            alglib_lsfit(tcalc);
            
        } while((rep.r2<0.99)&&(repcount<=10));
        
        if (rep.r2>r2_threshold)
        {
            //    write_tauq     << "T q q2 tau fit_r2 \n";
            //    write_sExp     << "T q q2 KWW_Params fit_r2 \n";
            //    write_sExp_all << "T q q2 KWW_Params fit_r2 \n";
            
            write_tauq
            << Temp_d*pow(corrFacforT,-1) << " "
            << wavenumber << " "
            << pow(wavenumber,2) << " "
            << sExp_tauFit_interpolateFvalue(c,tcalc) << " "
            << rep.r2 << " "
            << sExp_tau_lo << " "
            << sExp_cut_hi << " "
            << "\n";
            
            tau_isfs.push_back
            ({wavenumber,log10(sExp_tauFit_interpolateFvalue(c,tcalc))});
            
            //---------------------------
            qfound.push_back(wavenumber);
            cfound.push_back(c);
            repfound.push_back(rep);
            //---------------------------
            
            double q=wavenumber;
            double q2=pow(wavenumber,2);
            betaKWW_sortincreasing.push_back({log10(q2),c[2]});
            
            write_sExp
            << Temp_d*pow(corrFacforT,-1) << " "
            << q  << " "
            << q2 << " ";
            for (int i=0; i<c.length(); ++i) write_sExp << c[i] << " ";
            write_sExp
            << rep.r2      << " "
            << sExp_tau_lo << " "
            << sExp_cut_hi << "\n";
            
            write_sExp_all
            << Temp_d*pow(corrFacforT,-1) << " "
            << q  << " "
            << q2 << " ";
            for (int i=0; i<c.length(); ++i) write_sExp_all << c[i] << " ";
            write_sExp_all
            << rep.r2      << " "
            << sExp_tau_lo << " "
            << sExp_cut_hi << "\n";
            
            double T=Temp_d*pow(corrFacforT,-1);
            bKWWvsq.push_back({log10(q2),c[2],T,q2});//(log(q2),beta,T,q2)
            
        } else {
            cout
            << "trial_00"<<n_trl<<"_T"+to_string((long long int)Temp_d)<<" "
            << "q="<<wavenumber<<", cut_hi="<<sExp_cut_hi<<", "
            << "r2="<<rep.r2<<"(<"<<r2_threshold<<")\n";
        }
    } write_sExp_all << "\n";//add a spece between different Ts
    
    /** NOTE:resume value **/
    sExp_tau_lo=sExp_tau_lo_org;
    sExp_cut_hi=sExp_cut_hi_org;
    
    /** smooth bKWW vs q2 curve by penalized spline **/
    //--------------------------------------------------------------------------
    alglib::spline1dinterpolant inp0,inp1,inph,inp2;
    build_penalizedspline_interpolant(bKWWvsq,0,1,50,0.0);
    inp0=interpolant;
    build_penalizedspline_interpolant(bKWWvsq,0,1,50,0.5);
    inph=interpolant;
    build_penalizedspline_interpolant(bKWWvsq,0,1,50,1.0);
    inp1=interpolant;
    build_penalizedspline_interpolant(bKWWvsq,0,1,50,2.0);
    inp2=interpolant;
    for (int i=0;i<(int)bKWWvsq.size();++i)
    {
        double T=bKWWvsq.at(i).at(2);
        double q2=bKWWvsq.at(i).at(3);
        double logq2=bKWWvsq.at(i).at(0);
        double b=bKWWvsq.at(i).at(1);
        write_bKWW_all
        << T                        << " " //T
        << q2                       << " " //q2
        << logq2                    << " " //log(q2)
        << b                        << " " //beta
        << spline1dcalc(inp0,logq2) << " " //beta_rho0
        << spline1dcalc(inph,logq2) << " " //beta_rho0.5
        << spline1dcalc(inp1,logq2) << " " //beta_rho1.0
        << spline1dcalc(inp2,logq2) << " " //beta_rho2.0
        << "\n";
    } write_bKWW_all << "\n";
    //--------------------------------------------------------------------------
    
    /** find matching q value of tau_isfs for corresponding tau_baf **/
    //--------------------------------------------------------------------------
    if (analysispart!="mode12")
    {
        if (tau_baf>0&&tau_isfs.size()>0)
        {
            build_penalizedspline_interpolant(tau_isfs,0,1);
            
            double dnow=0,dmin=0;
            double q_match=0,tau_match=0;
            double qnow=0;
            double qmin=tau_isfs.at(0).at(0);
            double qmax=tau_isfs.at(tau_isfs.size()-1).at(0);
            double dq=fabs(qmax-qmin)/10000.0;
            
            for (int i=0;i<10000;++i) {
                qnow=qmin+dq*i;
                if (i==0) {
                    dmin=fabs(spline1dcalc(interpolant,qnow)-tau_baf);
                    //dmin=fabs(tau_isfs.at(i).at(1)-tau_baf);
                } else {
                    dnow=fabs(spline1dcalc(interpolant,qnow)-tau_baf);
                    //dnow=fabs(tau_isfs.at(i).at(1)-tau_baf);
                }
                /** update dmin if a smaller one than current dmin is found **/
                if (dnow<dmin) {
                    if (i>0) dmin=dnow;
                    q_match=qnow;
                    tau_match=spline1dcalc(interpolant,qnow);
                    //q_match  =tau_isfs.at(i).at(0);
                    //tau_match=tau_isfs.at(i).at(1);
                }
            }
            
            /** find closest q to q_match **/
            int index=0;
            bool is_lt=(qfound.at(0)<q_match);
            for (int i=0;i<(int)qfound.size();++i) {
                if (is_lt) {
                    if (qfound.at(i)>=q_match) {
                        index=i; break;
                    }
                } else {
                    if (qfound.at(i)<q_match) {
                        index=i; break;
                    }
                }
            }
            double q_isfs=0;
            real_1d_array c_isfs;
            lsfitreport rep_isfs;
            if (index>0)
            {
                double pre=fabs(qfound.at(index-1)-q_match);
                double pos=fabs(qfound.at(index)-q_match);
                if (pre<pos) {
                    q_isfs=qfound.at(index-1);
                    c_isfs=cfound.at(index-1);
                    rep_isfs=repfound.at(index-1);
                } else {
                    q_isfs=qfound.at(index);
                    c_isfs=cfound.at(index);
                    rep_isfs=repfound.at(index);
                }
            } else {
                cout
                << "in FitData::fit_q_dependences():\n"
                << "q index = 0 \n"; exit(EXIT_FAILURE);
            }
            
            write_q_match
            << Temp_d*pow(corrFacforT,-1) << " "
            << q_match   << " "
            << tau_match << " "
            << tau_baf   << " "
            << pow(10,tau_match)/pow(10,tau_baf) << "\n";
            
            //T(0) q qfs Afs taufs bfs(5) r2fs Abaf taubaf bbaf(9) r2baf
            
            write_q_KWW << Temp_d*pow(corrFacforT,-1) << " " << q_match << " ";
            write_q_KWW << q_isfs << " ";
            for (int i=0;i<c_isfs.length();++i) {
                write_q_KWW << c_isfs[i] << " ";
            } write_q_KWW << rep_isfs.r2 << " ";
            for (int i=0;i<c_baf.length();++i) {
                write_q_KWW << c_baf[i] << " ";
            } write_q_KWW << rep_baf.r2 << "\n";
        }
    }
    //--------------------------------------------------------------------------
    
#ifdef NEVER
    if (betaKWW_sortincreasing.size()>0)
    {
        build_penalizedspline_interpolant(betaKWW_sortincreasing,0,1);
        
        double logq2=0,f=0,df=0,d2f=0;
        size_t s=betaKWW_sortincreasing.size();
        double logq2_min=betaKWW_sortincreasing.at(0).at(0);
        double logq2_max=betaKWW_sortincreasing.at(s-1).at(0);
        double ds=fabs(logq2_max-logq2_min)/10000;
        double dfpre=0,dfnow=0;
        double min_logq2=0,min_beta=0;
        double max_logq2=0,max_beta=0;
        int count_min=0,count_max=0;
        //consider interpolated data in [3rd from top, 3rd from bottom]
        for (int index=2;index<=9998;++index) {
            logq2=logq2_min+ds*index;
            spline1ddiff(interpolant,logq2,f,df,d2f);
            dfnow=df;
            if (dfpre*dfnow<0) {
                double turnx=((logq2-ds)+logq2)/2;
                double turny=spline1dcalc(interpolant,turnx);
                double pre1=spline1dcalc(interpolant,logq2-ds);
                double pre2=spline1dcalc(interpolant,logq2-ds*2);
                double pos1=spline1dcalc(interpolant,logq2+ds);
                double pos2=spline1dcalc(interpolant,logq2+ds*2);
                /* local minimum */
                if (pre1<pre2&&pos1<pos2) {
                    if (turny<pre1&&turny<pos1) {
                        if (count_min==0) {
                            min_logq2=logq2;
                            min_beta =spline1dcalc(interpolant,logq2);
                            ++count_min;
                        }
                    }
                }
                /* local maximum */
                if (pre1>pre2&&pos1>pos2) {
                    if (turny>pre1&&turny>pos1) {
                        if (count_max==0) {
                            max_logq2=logq2;
                            max_beta =spline1dcalc(interpolant,logq2);
                            ++count_max;
                        }
                    }
                }
            } dfpre=dfnow;
        }
        write_beta_all
        << Temp_d*pow(corrFacforT,-1)     << " "
        << min_logq2                      << " "
        << sqrt(pow(10,min_logq2))        << " "
        << 2*pi()/sqrt(pow(10,min_logq2)) << " "
        << min_beta                       << "   "
        << max_logq2                      << " "
        << sqrt(pow(10,max_logq2))        << " "
        << 2*pi()/sqrt(pow(10,max_logq2)) << " "
        << max_beta                       << "\n";
    }
#endif
    
}





bool FitData::check_is_retry(StructureClass& sysVar)
{
    const int n_trial        = sysVar.get_n_trial();
    const int n_system       = sysVar.get_n_system();
    const int n_sys_beg      = sysVar.get_n_sys_beg();
    const int n_sys_end      = sysVar.get_n_sys_end();
    const int current_regime = sysVar.get_current_regime();
    const int n_regimeTs     = sysVar.get_n_regime_temps().at(current_regime);
    
    double Tie_min=0;              // lowest in-equil.(ie) T
    double tie_min=0;              // tau_alpha at Tie_min
    double Toe_max=0;              // highest out-of-equil.(oe) T
    double toe_max=0;              // tau_alpha at Toe_max
    
    double Tmin=0;                 // T at current tau_alpha cutoff
    double tmin=log10(tauFit_cut); // log(tau_alpha) of Tmin
    double Tmax=0;                 // lowest Tie of previous regime
    double tmax=0;                 // log(tau_alpha) of Tmax
    
    double tmin_pre=0;             // tau_alpha cutoff of previous regime
    if (current_regime>0) {
        tmin_pre=
        log10(sysVar.get_equilibration_times().at(current_regime-1)/n_equ_blocks);
    }
    
    bool is_empty=false;
    int  count_neq=0;
    int  n_threshold=(int)(0.5*n_regimeTs);
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            index=n_trl*n_system+(n_sys-n_sys_beg);
            
            if (is_imposeStrictTeq) {
                //--------------------------------------------------------------
                /** check_is_retry() is accessed after write_equTs(), so the
                 ** "equilibratedTs" & "tauEqu" containers had been properly
                 ** adjusted if is_imposeStrictTeq was turned true **/
                //--------------------------------------------------------------
                Toe_max=0;
                for (size_t i=0; i<sysVar.get_tauFit().at(index).size(); ++i) {
                    if (sysVar.get_tauFit().at(index).at(i).at(1)>log10(tauFit_cut)) {
                        Toe_max=sysVar.get_tauFit().at(index).at(i).at(0); //expanded
                        break;
                    }
                }
                if (Toe_max==0) {
                    count_neq += sysVar.get_tauFit().at(index).size();
                } else {
                    for (size_t i=0; i<sysVar.get_equilibratedTs().at(index).size(); ++i) {
                        if (sysVar.get_equilibratedTs().at(index).at(i)>Toe_max) {
                            ++count_neq;
                        }
                    }
                }
            } else {
                count_neq += sysVar.get_tauEqu().at(index).size();
            }
            /** check if any equil. data container is empty **/
            if (sysVar.get_equilibratedTs().at(index).size()==0) {
                is_empty=true;
                cout <<"\n"
                <<sysVar.get_usic()<<"_00"<<n_trl<<"_"<<sysVar.get_nameString(n_sys)<<":\n"
                <<"No equilibrium data found for this species at current run\n";
            }
        }
    } sysVar.add_n_cum_taueq(count_neq);
    
    cout << "\n"
    << "n_cum_taueq    " << sysVar.get_n_cum_taueq() << "\n"
    << "n_threshold    " << n_threshold << "\n";
    
    /** NOTE:
     ** The retry algorithm makes sure the TOTAL number of equilibrium data of
     ** all simulated species combined is at least n_threshold at current regime
     ** .AND. each equil. data container is NOT empty if automation proceeds **/
    
    bool is_retry=false;
    if (sysVar.get_n_cum_taueq()<n_threshold || is_empty) {
        if (current_regime==0) {
            /** NOTE: retry is only attempted after first regime **/
            cout << "\n"
            << "in FitData::check_is_retry():\n"
            << "# of in-equil. data in first regime = "<<count_neq
            << ", which is too few. Start with higher T or longer equilibration.\n";
            exit(EXIT_FAILURE);
        } is_retry=true;
    } else {
        is_retry=false;
    }
    
    /** retry algorithm kernel **/
    //--------------------------------------------------------------------------
    if (is_retry && (sysVar.get_times_retry()<sysVar.get_max_times_retry()))
    {
        int n_retryTs=0; // counting retry temperatures of all trials & sys
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                
                index=n_trl*n_system+(n_sys-n_sys_beg);
                
                /** Find Tie_min **/
                //--------------------------------------------------------------
                read_individual_taueq_data(sysVar,n_trl,n_sys);
                Tie_min=std::round(taueq_sortincreasing[0][0]*corrFacforT);//expanded
                tie_min=taueq_sortincreasing[0][1];
                
                /** Find Tmax **/
                //--------------------------------------------------------------
                Tmax=sysVar.get_extrpTmax().at(index)[0];//expanded
                tmax=sysVar.get_extrpTmax().at(index)[1];
                
                /** Find first Toe **/
                //--------------------------------------------------------------
                for (size_t i=0; i<sysVar.get_tauFit().at(index).size(); ++i) {
                    if (sysVar.get_tauFit().at(index).at(i)[1]>log10(tauFit_cut)) {
                        Toe_max=sysVar.get_tauFit().at(index).at(i)[0];//expanded
                        toe_max=sysVar.get_tauFit().at(index).at(i)[1];
                        break;
                    }
                }
                
                /** All Tie should be g.t. the first Toe **/
                //--------------------------------------------------------------
                if (Toe_max>Tie_min) {
                    double Ti=0;
                    for (size_t i=0; i<taueq_sortdecreasing.size(); ++i) {
                        Ti=std::round(taueq_sortdecreasing.at(i)[0]*corrFacforT);//expanded
                        if (Ti>Toe_max) {
                            Tie_min=Ti;
                            tie_min=taueq_sortdecreasing.at(i)[1];
                        }
                    }
                }
                
                /** write simulated temperatures to file **/
                //--------------------------------------------------------------
                string o;
                o.append(return_AnalysisFolderPath(sysVar));
                o.append("/fit_data");
                o.append("/tempsInfo_");
                o.append(sysVar.get_usic());
                o.append("_00"+to_string((long long int)n_trl));
                o.append("_"+sysVar.get_nameString(n_sys));
                o.append(".dat");
                ofstream write_simtempsinfo(o.c_str(),ofstream::app);//"append"
                //--------------------------------------------------------------
                // NOTE:
                // temperatures in this file has been normalized to "expanded" format
                // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
                // so "setprecision" needs to be set to zero
                //--------------------------------------------------------------
                write_simtempsinfo << fixed << setprecision(0);
                for (size_t i=0; i<sysVar.get_temperaturesInfo().at(index).size(); ++i) {
                    write_simtempsinfo << sysVar.get_temperaturesInfo().at(index).at(i) << "\n";
                } write_simtempsinfo.close();
                
                // NOTE:
                // use actual value instead of expanded for calculating slope
                // on an "Arrhenius plot" (NOTE: inverse T scale)
                //--------------------------------------------------------------
                double Tie_actual  = Tie_min*pow(corrFacforT,-1);//actual
                double Toe_actual  = Toe_max*pow(corrFacforT,-1);//actual
                double Tmax_actual = Tmax*pow(corrFacforT,-1);   //actual
                double Tmin_actual = 0;
                double TminArr=0,TminVFT=0;
                
                /** NOTE: 1/T or 1000/T **/
                double A=convert; /** A(real)=1000, A(LJ)=1 **/
                double slope=get_2pt_slope({A/Tie_actual,tie_min},{A/Toe_actual,toe_max});
                if (slope<0) {
                    cout << "\n"
                    << "in FitData::check_is_retry():\n"
                    << "WARNINIG: slope of interpolation < 0\n"
                    << "Tie_actual = " << Tie_actual << "\n"
                    << "tie_min    = " << tie_min    << "\n"
                    << "Toe_actual = " << Toe_actual << "\n"
                    << "toe_max    = " << toe_max    << "\n";
                    exit(EXIT_FAILURE);
                }
                
                /** crossover T by Arrhenius **/
                //--------------------------------------------------------------
                Tmin_actual=A*pow((A/Tie_actual)+((tmin-tie_min)/slope),-1);
                TminArr=std::round(Tmin_actual*corrFacforT);//expanded
                
                /** crossover T by VFT **/
                //--------------------------------------------------------------
                set_fitParams("VFT");
                read_individual_taueq_data(sysVar,n_trl,n_sys);
                fitdata_processing_lsfit(taueq_sortdecreasing);
                alglib_lsfit("VFT");
                TminVFT=calc_x_given_y(c,compu_time,"VFT")[0];//actual
                TminVFT=std::round(TminVFT*corrFacforT);//expanded
                
                /** crossover T is the larger Tmin  **/
                //--------------------------------------------------------------
                bool is_byVFT=false;
                if (TminVFT>TminArr) {
                    is_byVFT=true;
                    Tmin=TminVFT;
                } else {
                    Tmin=TminArr;
                } Tmin_actual=Tmin*pow(corrFacforT,-1);//actual
                
                /** number of retry temperatures **/
                //--------------------------------------------------------------
                int n_retry_tmp=0;
                if (n_regimeTs<2) n_retry_tmp=2;
                else              n_retry_tmp=n_regimeTs;
                
                /** get interpolated temperatures between Tmax & Tmin **/
                //--------------------------------------------------------------
                vector<double> hightolow=
                get_interpolated_Ts(sysVar,n_trl,n_sys,n_retry_tmp,Tmax,Tmin);
                
                if (sysVar.get_is_updatetinfo())
                {
                    sysVar.get_temperaturesInfo().at(index)=hightolow;
                    
                    /** clear container for new quench temperatures **/
                    //**********************************************************
                    sysVar.get_quenchingTs().at(index).clear();
                    
                    sysVar.get_equilibratedTs().at(index).clear();//NOTE
                    sysVar.get_equilibratedTs().at(index).push_back({Tmax});
                    //**********************************************************
                    
                    /** quenchingT's should be < Tmax **/
                    //----------------------------------------------------------
                    int count=0;
                    do {
                        double Ti=sysVar.get_temperaturesInfo().at(index).at(count);
                        if (Ti<Tmax) { //NOTE: Tmax
                            sysVar.get_quenchingTs().at(index).push_back(Ti);
                            ++n_retryTs;
                        } ++count;
                    } while (count<sysVar.get_temperaturesInfo().at(index).size());
                    //----------------------------------------------------------
                }
                cout << "\n"
                <<sysVar.get_usic()<<"_00"<<n_trl<<"_"<<sysVar.get_nameString(n_sys)<<":\n"
                << "  (Tie_min,tie) ("<<Tie_actual<<","<<tie_min<<")\n"
                << "  (Toe_max,toe) ("<<Toe_actual<<","<<toe_max<<")\n"
                << "(Tmin_Arr,tmin) ("<<TminArr*pow(corrFacforT,-1)<<","<<tmin<<")\n"
                << "(Tmin_VFT,tmin) ("<<TminVFT*pow(corrFacforT,-1)<<","<<tmin<<")\n"
                << "(Tmax_pre,tmax) ("<<Tmax_actual<<","<<tmax<<")\n";
                if (is_byVFT) cout << "use VFT interpolation:" << "\n";
                else cout << "use Arrhenius interpolation:"    << "\n";
                cout << "retry Ts: ";
                for (size_t i=0;i<sysVar.get_quenchingTs().at(index).size();++i) {
                    cout << sysVar.get_quenchingTs().at(index).at(i) << " ";
                } cout << "\n";
            }
        } cout << "\nn_totalretryTs " << n_retryTs << "\n";
    } return is_retry;
}





void FitData::save_fit_coeffs(StructureClass& sysVar)
{
    // save well-fit KWW parameters
    // NOTE:
    // if any of the c parameters is less than 1e-3, don't use the parameters
    // because it would possibly become extremly small and lead to ALGLIB error
    
    vector<double> tmpCoeffs;
    for (size_t i=0; i<c.length(); ++i) {
        if (c[i]<1e-3) return;
        tmpCoeffs.push_back(c[i]);
    }
    sysVar.set_c_vD(tmpCoeffs);
    sysVar.set_rep(rep);
    sysVar.set_coeffs_bndl_vD(coeffs_bndl_vD);
    sysVar.set_coeffs_bndu_vD(coeffs_bndu_vD);
    sysVar.set_coeffs_bndl(coeffs_bndl);
    sysVar.set_coeffs_bndu(coeffs_bndu);
    sysVar.set_epsf(epsf);
    sysVar.set_epsx(epsx);
    sysVar.set_maxits(maxits);
    sysVar.set_diffstep(diffstep);
}





bool FitData::load_fit_coeffs(const StructureClass& sysVar)
{
    // load saved well-fit KWW parameters
    if (sysVar.get_c_vD().size()==0) {
        cout
        << "in FitData::load_fit_coeffs: "
        << "sysVar.c_vD.size()==0, "
        << "use default KWW coeffs.\n";
        return false;
    }
    set_coeffs(sysVar.get_c_vD());
    set_coeffs_scale(sysVar.get_c_vD());
    coeffs_bndl_vD = sysVar.get_coeffs_bndl_vD();
    coeffs_bndu_vD = sysVar.get_coeffs_bndu_vD();
    coeffs_bndl    = sysVar.get_coeffs_bndl();
    coeffs_bndu    = sysVar.get_coeffs_bndu();
    epsf           = sysVar.get_epsf();
    epsx           = sysVar.get_epsx();
    maxits         = sysVar.get_maxits();
    diffstep       = sysVar.get_diffstep();
    return true;
}





void FitData::fit_Arrhenius(StructureClass& sysVar,
                            const int n_trl,
                            const int n_sys)
{
    string model="Arrhenius";
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    if (index==0) {
        ArrheniusCoeffs.clear();
        r2_Arrhenius.clear();
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_Arrhenius_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fit_Arrhenius(o.c_str());
    
    //find_cutTforArrhenius(sysVar,n_sys);
    set_Theq(get_Theq(sysVar,n_sys)[0]);
    
    set_fitParams(model);
    read_individual_taueq_data(sysVar,n_trl,n_sys,model);
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    
    real_1d_array c_org=c;
    r2_Arrhenius.push_back(rep.r2);
    vector<double> tmpCoeffs;
    for (size_t i=0; i<c.length(); ++i) tmpCoeffs.push_back(c[i]);
    ArrheniusCoeffs.push_back(tmpCoeffs);
    
    write_fitModel(write_fit_Arrhenius,model);
    write_fitarrays(write_fit_Arrhenius);
    write_fitcutoff(write_fit_Arrhenius,sysVar,model);
    write_stopcond(write_fit_Arrhenius);
    write_fitinfo(write_fit_Arrhenius);
    write_errcurve(write_fit_Arrhenius,model);
    write_tauFitfit_vs_T(write_fit_Arrhenius,c_org,Theq,model);
    write_tauFit_vs_T(write_fit_Arrhenius);
    write_tauFitfit_vs_invT(write_fit_Arrhenius,c_org,Theq,model);
    write_tauFit_vs_invT(write_fit_Arrhenius);
    
    avgfitParams(sysVar,n_trl,n_sys,model);
    
    if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg))) {
        // last index
        write_fitavg(sysVar,n_sys,model);
        write_TA(sysVar,n_sys);
    } write_fit_Arrhenius.close();
}





void FitData::fit_tauFit(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys,
                         const string& extrp_model,
                         const string& presq_model)
{
    //string model=relmodel[0];
    //string extrp=relmodel[1];
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    if (index==0)
    {
        ModelCoeffs.clear();
        ExtrpCoeffs.clear();
        r2_Model.clear();
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_taueq_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    ofstream write_fit_tauFit(o.c_str());
    
    set_Theq(get_Theq(sysVar,n_sys)[0]);
    set_Tleq(get_Tleq(sysVar,n_trl,n_sys)[0]);
    
    set_fitParams(extrp_model);
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(extrp_model);
    
    real_1d_array c_org=c;
    r2_Model.push_back(rep.r2);
    vector<double> tmpCoeffs;
    for (size_t i=0; i<c.length(); ++i) tmpCoeffs.push_back(c[i]);
    ModelCoeffs.push_back(tmpCoeffs);
    compute_Tg_fragility(sysVar,c_org,extrp_model);
    
    write_fitModel(write_fit_tauFit,extrp_model);
    write_fitarrays(write_fit_tauFit);
    write_fitcutoff(write_fit_tauFit,sysVar,extrp_model);
    write_stopcond(write_fit_tauFit);
    write_Tg_fragility(write_fit_tauFit,c_org,extrp_model);
    write_fitinfo(write_fit_tauFit);
    write_errcurve(write_fit_tauFit,extrp_model);
    write_tauFitfit_vs_T(write_fit_tauFit,c_org,Theq,extrp_model);
    write_tauFit_vs_T(write_fit_tauFit);
    write_tauFitfit_vs_invT(write_fit_tauFit,c_org,Theq,extrp_model);
    write_tauFit_vs_invT(write_fit_tauFit);
    
    avgfitParams(sysVar,n_trl,n_sys,extrp_model);
    
    if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
    {
        write_fitavg(sysVar,n_sys,extrp_model);
        write_Tc(sysVar,n_sys);
        //write_tau0(sysVar,n_sys,extrp_model);
    }
    
    ////////////////////////////////////////////////////////////////////////
    //======================================================================
    if (get_is_shootForNext())
    {
        set_fitParams(presq_model);
        read_individual_taueq_data(sysVar,n_trl,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit(presq_model);
        
        real_1d_array c_extrp=c;
        vector<double> tmpCoeffs;
        for (size_t i=0; i<c.length(); ++i) tmpCoeffs.push_back(c[i]);
        ExtrpCoeffs.push_back(tmpCoeffs);
        
        find_largest_Tg_fragility(sysVar,n_trl,n_sys,c_extrp,presq_model);
        avgfitParams_extrp(sysVar,n_trl,n_sys,presq_model);
        
        shootForNewTemperatures(sysVar,n_trl,n_sys,c_extrp,presq_model);
        
        /** write fit info to file (VFT) **/
        //------------------------------------------------------------------
        set_fitParams("VFT");
        read_individual_taueq_data(sysVar,n_trl,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit("VFT");
        compute_VFT_Tg_fragility(sysVar,c);
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_taueq_VFT_");
        o.append(sysVar.get_usic());
        o.append("_00"+to_string((long long int)n_trl));
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append("_Regime"+to_string((long long int)current_regime));
        o.append(".dat");
        ofstream write_VFT(o.c_str());
        write_fitModel(write_VFT,"VFT");
        write_fitarrays(write_VFT);
        write_fitcutoff(write_VFT,sysVar,"VFT");
        write_stopcond(write_VFT);
        write_Tg_fragility(write_VFT,c,"VFT");
        write_fitinfo(write_VFT);
        write_errcurve(write_VFT,"VFT");
        write_tauFitfit_vs_T(write_VFT,c,Theq,"VFT");
        write_tauFit_vs_T(write_VFT);
        write_tauFitfit_vs_invT(write_VFT,c,Theq,"VFT");
        write_tauFit_vs_invT(write_VFT);
        write_VFT.close();
        //------------------------------------------------------------------
        
        /** write fit info to file (COOP) **/
        //------------------------------------------------------------------
        set_fitParams("COOP");
        read_individual_taueq_data(sysVar,n_trl,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit("COOP");
        compute_COOP_Tg_fragility(sysVar,c);
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_taueq_COOP_");
        o.append(sysVar.get_usic());
        o.append("_00"+to_string((long long int)n_trl));
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append("_Regime"+to_string((long long int)current_regime));
        o.append(".dat");
        ofstream write_COOP(o.c_str());
        write_fitModel(write_COOP,"COOP");
        write_fitarrays(write_COOP);
        write_fitcutoff(write_COOP,sysVar,"COOP");
        write_stopcond(write_COOP);
        write_Tg_fragility(write_COOP,c,"COOP");
        write_fitinfo(write_COOP);
        write_errcurve(write_COOP,"COOP");
        write_tauFitfit_vs_T(write_COOP,c,Theq,"COOP");
        write_tauFit_vs_T(write_COOP);
        write_tauFitfit_vs_invT(write_COOP,c,Theq,"COOP");
        write_tauFit_vs_invT(write_COOP);
        write_COOP.close();
        //------------------------------------------------------------------
    }
    //======================================================================
    ////////////////////////////////////////////////////////////////////////
    
    /** final regime **/
    if (current_regime==(n_regime-1))
    {
        /** Stickel derivative analysis **/
        //stickelanalysis(sysVar);
        
        if (!sysVar.get_is_aging()) {
            if (get_is_fitByEveryPoint()) {
                if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg))) {
                    fit_continuousRange(sysVar,n_sys,extrp_model);
                    fit_neighboringNpoints(sysVar,n_sys,extrp_model);
                    if (get_is_fit_by_Tc()||get_is_fit_by_TA()) {
                        fit_partial(sysVar,n_sys,extrp_model);
                    }
                }
            }
        }
    } write_fit_tauFit.close();
}





void FitData::alglibfit_routine(const StructureClass& sysVar)
{
    string in,out;
    /* input data path */
    in.append(return_AnalysisFolderPath(sysVar));
    in.append("/fit_data");
    in.append("/ain.dat");
    /* output data path */
    out.append(return_AnalysisFolderPath(sysVar));
    out.append("/fit_data");
    out.append("/aout.dat");
    
    bool is_sortX_ascending=true;
    bool is_show_pspline   =true;
    bool is_show_fitfunc   =true;
    
    //double cutoff_X=0;
    //double cutoff_Y=0;
    
    vector<vector<double>> datavvd;//read-in data ({X,Y1,Y2,...,Yn})
    
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        double xdata=0,ydata=0;
        while (getline(readFile,lineContent))
        {
            vector<double> datavd;
            istringstream iss(lineContent);
            /** avoid blank lines in file **/
            if (iss>>xdata) {
                datavd.push_back(xdata);
            } else {
                continue;
            }
            while (iss>>ydata) {
                datavd.push_back(ydata);
            } datavvd.push_back(datavd);
        } readFile.close();
    } else {
        cout
        << "in FitData::alglibfit_routine():\n"
        << "input file name should be 'ain.dat', please check.\n";
        exit(EXIT_FAILURE);
    }
    if (is_sortX_ascending) {
        sort(datavvd.begin(),datavvd.end(),sortIncreasing0);
    }
    int n_dataset=(int)datavvd.at(0).size()-1;//datavvd={X,Y1,Y2,...}
    
    if (false)
    {
        for (int i=0;i<(int)datavvd.size();++i) {
            cout << datavvd.at(i).at(0) << " ";
            for (int ii=1;ii<=n_dataset;++ii) {
                cout << datavvd.at(i).at(ii) << " ";
            } cout << "\n";
        }
    }
    
    ofstream writeFile(out.c_str());
    writeFile<<"X  {Y}\n";
    for (int ii=0;ii<datavvd.size();++ii) {
        for (int jj=0;jj<(int)datavvd.at(ii).size();++jj) {
            writeFile << datavvd[ii][jj] << " ";
        } writeFile << "\n";
    } writeFile << "\n";
    
    if (is_show_pspline)
    {
        /** build penalized spline **/
        vector<alglib::spline1dinterpolant> interpolantii;
        for (int ii=1;ii<=n_dataset;++ii) {
            build_penalizedspline_interpolant(datavvd,0,ii,50,2);
            interpolantii.push_back(interpolant);
        }
        
        writeFile<<"pSpline fit results:\n\n";
        if (n_dataset>1) {
            writeFile << "X  {Y}  \n";
        } else {
            writeFile << "X  Y  \n";
        }
        for (int i=0;i<(int)datavvd.size();++i) {
            writeFile << datavvd.at(i).at(0) << " ";
            for (int ii=1;ii<=n_dataset;++ii) {
                writeFile
                <<spline1dcalc(interpolantii.at(ii-1),datavvd.at(i).at(0))
                << " ";
            } writeFile << "\n";
        } writeFile << "\n\n";
    }
    if (is_show_fitfunc)
    {
        string func="HRM1";
        //string func="linear3";
        
        cout << "Fit to "<<n_dataset<<" dataset(s)\n";
        
        vector<real_1d_array> cii;
        vector<lsfitreport> repii;
        
        writeFile<<"Fit to <"<<func<<"> form\n\n";
        for (int ii=1;ii<=n_dataset;++ii) {//NOTE: ii=[1,dataset]
            vector<vector<double>> dataii;
            for (int jj=0;jj<(int)datavvd.size();++jj) {
                dataii.push_back
                ({
                    datavvd.at(jj).at(0),
                    datavvd.at(jj).at(ii)
                });
            }
            set_fitParams(func);
            fitdata_processing_lsfit(dataii);//NOTE
            alglib_lsfit(func);
            cii.push_back(c);
            repii.push_back(rep);
        }
        
        writeFile << "fit coeffs \n";
        for (int ii=0;ii<n_dataset;++ii) {//NOTE: ii=[0,dataset)
            for (int jj=0;jj<cii.at(ii).length();++jj) {
                writeFile << cii[ii][jj] << " ";
            } writeFile << "\n";
        } writeFile << "\n";
        
        writeFile << "fit r2 \n";
        for (int ii=0;ii<n_dataset;++ii) {//NOTE: ii=[0,dataset)
            writeFile << repii[ii].r2 << "\n";
        } writeFile << "\n";
        
    } writeFile << "\n";
}





void FitData::alglib_solverobj(const alglib::real_2d_array& x,
                               const alglib::real_1d_array& y,
                               const alglib::real_1d_array& c,
                               alglib::lsfitstate& state)
{
    if (is_use_FG) {
        lsfitcreatefg(x, y, c, true, state);
    } else {
        lsfitcreatef(x, y, c, diffstep, state);
    }
}





void FitData::alglib_lsfit_form(const std::string& model,
                                alglib::lsfitstate& state)
{
    if (model=="KWW")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, KWW_func, KWW_grad);
        } else {
            alglib::lsfitfit(state, KWW_func);
        }
    }
    else if (model=="KWW_lnFs")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, KWW_lnFs_func, KWW_lnFs_grad);
            alglib::lsfitfit(state, KWW_lnFs_func);
        } else {
            alglib::lsfitfit(state, KWW_lnFs_func);
        }
    }
    else if (model=="KWW_full")
    {
        if (is_use_FG) {
            //alglib::lsfitfit(state, KWW_full_func, KWW_full_grad);
            alglib::lsfitfit(state, KWW_full_func);
        } else {
            alglib::lsfitfit(state, KWW_full_func);
        }
    }
    else if (model=="KWW_full_lnFs")
    {
        if (is_use_FG) {
            //alglib::lsfitfit(state, KWW_full_lnFs_func, KWW_full_lnFs_grad);
            alglib::lsfitfit(state, KWW_full_lnFs_func);
        } else {
            alglib::lsfitfit(state, KWW_full_lnFs_func);
        }
    }
    else if (model=="KWW_pwr")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, KWW_pwr_func, KWW_pwr_grad);
        } else {
            alglib::lsfitfit(state, KWW_pwr_func);
        }
    }
    else if (model=="mKWW")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, mKWW_func, mKWW_grad);
        } else {
            alglib::lsfitfit(state, mKWW_func);
        }
    }
    else if (model=="mKWW_pwr")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, mKWW_pwr_func, mKWW_pwr_grad);
        } else {
            alglib::lsfitfit(state, mKWW_pwr_func);
        }
    }
    else if (model=="linear")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, linear_func, linear_grad);
        } else {
            alglib::lsfitfit(state, linear_func);
        }
    }
    else if (model=="diagonal")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, linear_func, linear_grad);
        } else {
            alglib::lsfitfit(state, linear_func);
        }
    }
    else if (model=="linear2")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, linear_func2, linear_grad2);
        } else {
            alglib::lsfitfit(state, linear_func2);
        }
    }
    else if (model=="linear3")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, linear_func3, linear_grad3);
        } else {
            alglib::lsfitfit(state, linear_func3);
        }
    }
    else if (model=="univu2T")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, univu2T_func, univu2T_grad);
        } else {
            alglib::lsfitfit(state, univu2T_func);
        }
    }
    else if (model=="secondpoly")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, secondpoly_func, secondpoly_grad);
        } else {
            alglib::lsfitfit(state, secondpoly_func);
        }
    }
    else if (model=="exp_string")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, exp_string_func, exp_string_grad);
        } else {
            alglib::lsfitfit(state, exp_string_func);
        }
    }
    else if (model=="pwr_exp_string")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, pwr_exp_string_func, pwr_exp_string_grad);
        } else {
            alglib::lsfitfit(state, pwr_exp_string_func);
        }
    }
    else if (model=="Arrhenius"||model=="Arrhenius_string")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, Arrhenius_func, Arrhenius_grad);
        } else {
            alglib::lsfitfit(state, Arrhenius_func);
        }
    }
    else if (model=="string2param")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, string2param_func, string2param_grad);
        } else {
            alglib::lsfitfit(state, string2param_func);
        }
    }
    else if (model=="string3param")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, string3param_func, string3param_grad);
        } else {
            alglib::lsfitfit(state, string3param_func);
        }
    }
    else if (model=="transitState")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, transitState_func, transitState_grad);
        } else {
            alglib::lsfitfit(state, transitState_func);
        }
    }
    else if (model=="hallwolynes")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, hallwolynes_func, hallwolynes_grad);
        } else {
            alglib::lsfitfit(state, hallwolynes_func);
        }
    }
    else if (model=="hallwolynes1param")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, hallwolynes1param_func, hallwolynes1param_grad);
        } else {
            alglib::lsfitfit(state, hallwolynes1param_func);
        }
    }
    else if (model=="leporiniref")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, leporiniref_func, leporiniref_grad);
        } else {
            alglib::lsfitfit(state, leporiniref_func);
        }
    }
    else if (model=="leporini_universal")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, leporini_universal_func, leporini_universal_grad);
        } else {
            alglib::lsfitfit(state, leporini_universal_func);
        }
    }
    else if (model=="GLMVFT1"||model=="GLMVFT3"||model=="GLMVFT4")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, GLMVFT_func, GLMVFT_grad);
        } else {
            alglib::lsfitfit(state, GLMVFT_func);
        }
    }
    else if (model=="VFT")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, VFT_func, VFT_grad);
        } else {
            alglib::lsfitfit(state, VFT_func);
        }
    }
    else if (model=="Mauro")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, Mauro_func, Mauro_grad);
        } else {
            alglib::lsfitfit(state, Mauro_func);
        }
    }
    else if (model=="AM")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, AM_func, AM_grad);
        } else {
            alglib::lsfitfit(state, AM_func);
        }
    }
    else if (model=="DG")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, DG_func, DG_grad);
        } else {
            alglib::lsfitfit(state, DG_func);
        }
    }
    else if (model=="MCT")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, MCT_func, MCT_grad);
        } else {
            alglib::lsfitfit(state, MCT_func);
        }
    }
    else if (model=="SOU")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, SOU_func, SOU_grad);
        } else {
            alglib::lsfitfit(state, SOU_func);
        }
    }
    else if (model=="ArrheniusII")
    {
        alglib::lsfitfit(state, ArrheniusII_func);
    }
    else if (model=="GLM")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, GLM_func, GLM_grad);
        } else {
            alglib::lsfitfit(state, GLM_func);
        }
    }
    else if (model=="GLM1"||model=="GLM1_hallwolynes")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, GLM1_func, GLM1_grad);
        } else {
            alglib::lsfitfit(state, GLM1_func);
        }
    }
    else if (model=="univGLM")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, univGLM_func, univGLM_grad);
        } else {
            alglib::lsfitfit(state, univGLM_func);
        }
    }
    else if (model=="univGT")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, univGT_func, univGT_grad);
        } else {
            alglib::lsfitfit(state, univGT_func);
        }
    }
    else if ((model=="COOP")||(model=="COOP_string")||(model=="COOP_DWF"))
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, COOP_func, COOP_grad);
        } else {
            alglib::lsfitfit(state, COOP_func);
        }
    }
    else if (model=="DEAG")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, DEAG_func, DEAG_grad);
        } else {
            alglib::lsfitfit(state, DEAG_func);
        }
    }
    else if (model=="CG")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, CG_func, CG_grad);
        } else {
            alglib::lsfitfit(state, CG_func);
        }
    }
    else if (model=="FourParamVFT")
    {
        if (is_use_FG) {
            alglib::lsfitfit(state, FourParamVFT_func, FourParamVFT_grad);
        } else {
            alglib::lsfitfit(state, FourParamVFT_func);
        }
    }
    else if (model=="ArrheniusIII")
    {
        alglib::lsfitfit(state, ArrheniusIII_func);
    }
    else if (model=="HRM1")
    {
        alglib::lsfitfit(state, HRM1_func);
    }
    else
    {
        cout
        << "in FitData::alglib_lsfit_form():\n"
        << "model specified ("<<model<<") not found!\n";
        exit(EXIT_FAILURE);
    }
}





void FitData::set_fitParams(const std::string& model)
{
    epsf     = 0;
    epsx     = 1e-10;
    maxits   = 1e+4;
    diffstep = 1e-5;
    
    if (model=="KWW"||model=="KWW_lnFs") // [A,tau,beta]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1e+5, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1e+2, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf});
        }
    }
    else if (model=="KWW_full"||model=="KWW_full_lnFs") // [A,tau,beta,tauf,p]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1e+3, 0.5, 1e+2, 0.5});
            set_coeffs_scale({ 1.0, 1e+5, 0.5, 1e+2, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1e+0, 0.5, 0.1, 0.5});
            set_coeffs_scale({ 1.0, 1e+2, 0.5, 0.1, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0, 0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf, inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0,  1.0, 0.5, 0.1, 0.5});
            set_coeffs_scale({ 1.0,  1.0, 0.5, 0.1, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0, 0.0, 0.0});
            set_coeffs_bndu({  1.0,  inf, inf, inf, inf});
        }
    }
    else if (model=="KWW_pwr") // [A,a,tau,beta]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1.0, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+5, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1.0, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0, 1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0,  inf, inf});
        }
    }
    else if (model=="mKWW") // [A,taulib,tau,beta]
    {
        // NOTE:
        // here restrict taulib = 1 tau or ps
        
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1e+3, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1e+3, 1e+5, 0.5});
            set_coeffs_bndl({  0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 1e+3,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1e+3, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1e+3, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 1e+3,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0,  1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0,  1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0,  1.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0,  1.0,  inf, inf});
        }
    }
    else if (model=="mKWW_pwr") // [A,a,taulib,tau,beta]
    {
        // NOTE:
        // here restrict taulib = 1 tau or ps
        
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1.0, 1e+3, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+3, 1e+5, 0.5});
            set_coeffs_bndl({  0.0, 0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0, 1e+3,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1.0, 1e+3, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+3, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0, 1e+3,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 1.0,  1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0, 1.0,  1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  1.0,  0.0, 0.0});
            set_coeffs_bndu({  1.0, 5.0,  1.0,  inf, inf});
        }
    }
    else if(model=="linear") // [slope,intercept]
    {
        set_coeffs({       1.0, 0.0});
        set_coeffs_scale({ 1.0, 1.0});
        set_coeffs_bndl({ -inf,-inf});
        set_coeffs_bndu({  inf, inf});
    }
    else if(model=="diagonal") // [slope,intercept]
    {
        set_coeffs({       1.0, 0.0});
        set_coeffs_scale({ 1.0, 1.0});
        set_coeffs_bndl({  1.0,-inf});
        set_coeffs_bndu({  1.0, inf});
    }
    else if(model=="linear2") // [slope,intercept]
    {
        set_coeffs({       1.5,  0.5});
        set_coeffs_scale({ 1.5,  0.5});
        set_coeffs_bndl({ -inf, -inf});
        set_coeffs_bndu({  inf,  inf});
    }
    else if(model=="linear3") // [intercept]
    {
        set_coeffs({       1.5});
        set_coeffs_scale({ 1.5});
        set_coeffs_bndl({ -inf});
        set_coeffs_bndu({  inf});
    }
    else if(model=="univu2T") // [u2A,slope,intercept]
    {
        set_coeffs({       u2A,  1.5, 0.5});
        set_coeffs_scale({ u2A,  1.5, 0.5});
        set_coeffs_bndl({  0.0,  1.5, 0.5});
        set_coeffs_bndu({  inf,  1.5, 0.5});
    }
    else if(model=="secondpoly")
    {
        set_coeffs({       1.0,  1.0,  1.0});
        set_coeffs_scale({ 1.0,  1.0,  1.0});
        set_coeffs_bndl({ -inf, -inf, -inf});
        set_coeffs_bndu({  inf,  inf,  inf});
    }
    else if(model=="exp_string") // [A,B]
    {
        set_coeffs({       1.0, 1.0});
        set_coeffs_scale({ 1.0, 1.0});
        set_coeffs_bndl({  0.0, 0.0});
        set_coeffs_bndu({  inf, inf});
    }
    else if (model=="pwr_exp_string") //[A,B,C]
    {
        set_coeffs({       1.0, 1.0, 1.0});
        set_coeffs_scale({ 1.0, 1.0, 1.0});
        set_coeffs_bndl({  0.0, 0.0, 0.0});
        set_coeffs_bndu({  inf, inf, inf});
    }
    else if (model=="Arrhenius") // [tau0_Arrhenius,Ea/k]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 2e+3});
            set_coeffs_scale({ 1e+3, 1e+3});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 2e+3});
            set_coeffs_scale({  1.0, 1e+3});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  2.0});
            set_coeffs_scale({  1.0,  1.0});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
    }
    else if (model=="Arrhenius_string") // [tau0_Arrhenius,Ea/k]
    {
        set_coeffs({       5.0, 1.0});
        set_coeffs_scale({ 1.0, 1.0});
        set_coeffs_bndl({  0.0, 0.0});
        set_coeffs_bndu({  inf, inf});
    }
    else if (model=="string2param") // [delH,delS]
    {
        set_coeffs({       1.0, -1.0});
        set_coeffs_scale({ 1.0,  1.0});
        set_coeffs_bndl({ -inf, -inf});
        set_coeffs_bndu({  inf,  inf});
    }
    else if (model=="string3param") // [tau0,delH,delS]
    {
        set_coeffs({        0.1,  1.0, -1.0});
        set_coeffs_scale({  0.1,  1.0,  1.0});
        set_coeffs_bndl({   0.0, -inf, -inf});
        set_coeffs_bndu({  1e+2,  inf,  inf});
    }
    else if (model=="transitState") // [delH,delS]
    {
        set_coeffs({       1.0, -1.0});
        set_coeffs_scale({ 1.0,  1.0});
        set_coeffs_bndl({ -inf, -inf});
        set_coeffs_bndu({  inf,  inf});
    }
    else if (model=="hallwolynes") // [tau0,u20]
    {
        if (systemUnit=="real")
        {
            set_coeffs({        1.0, 10.0});
            set_coeffs_scale({  1.0, 10.0});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1e-4, 10.0});
            set_coeffs_scale({ 1e-4, 10.0});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1e-4,  0.5});
            set_coeffs_scale({ 1e-4,  0.1});
            set_coeffs_bndl({   0.0,  0.0});
            set_coeffs_bndu({   inf,  inf});
        }
    }
    else if (model=="hallwolynes1param") // [u20]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       10.0});
            set_coeffs_scale({ 10.0});
            set_coeffs_bndl({   0.0});
            set_coeffs_bndu({   inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       10.0});
            set_coeffs_scale({ 10.0});
            set_coeffs_bndl({   0.0});
            set_coeffs_bndu({   inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.5});
            set_coeffs_scale({ 0.1});
            set_coeffs_bndl({  0.0});
            set_coeffs_bndu({  inf});
        }
    }
    else if (model=="leporiniref") // [A,B,C]
    {
        if (is_1paramLeporini) {
            set_coeffs({        1.0, 1.6285, 12.325});
            set_coeffs_scale({  1.0, 1.6285, 12.325});
            set_coeffs_bndl({  -inf, 1.6285, 12.325});
            set_coeffs_bndu({   inf, 1.6285, 12.325});
        } else {
            set_coeffs({        1.0,  1.0, 10.0});
            set_coeffs_scale({  1.0,  1.0, 10.0});
            set_coeffs_bndl({  -inf, -inf, -inf});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
    }
    else if (model=="leporini_universal") // [A,B,C]
    {
        if (systemUnit=="real")
        {
            set_coeffs({        2.576, 1.6285, 12.325});
            set_coeffs_scale({  2.576, 1.6285, 12.325});
            set_coeffs_bndl({   2.576, 1.6285, 12.325});
            set_coeffs_bndu({   2.576, 1.6285, 12.325});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       -0.424, 1.6285, 12.325});
            set_coeffs_scale({ -0.424, 1.6285, 12.325});
            set_coeffs_bndl({  -0.424, 1.6285, 12.325});
            set_coeffs_bndu({  -0.424, 1.6285, 12.325});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       -0.424, 1.6285, 12.325});
            set_coeffs_scale({ -0.424, 1.6285, 12.325});
            set_coeffs_bndl({  -0.424, 1.6285, 12.325});
            set_coeffs_bndu({  -0.424, 1.6285, 12.325});
        }
    }
    else if (model=="GLMVFT1") // [a]
    {
        set_coeffs({       tauA/exp(1), TA/3.0, 3.0, 1.0});
        set_coeffs_scale({ tauA/exp(1), TA/3.0, 3.0, 1.0});
        set_coeffs_bndl({  tauA/exp(1), TA/3.0, 0.0, 1.0});
        set_coeffs_bndu({  tauA/exp(1), TA/3.0, inf, 1.0});
    }
    else if (model=="GLMVFT3") // [tau0,T0,a]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 1e+2, 3.0, 1.0});
            set_coeffs_scale({ 1e+3, 1e+2, 3.0, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0, 1.0});
            set_coeffs_bndu({   inf,  inf, inf, 1.0});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 1e+2, 3.0, 1.0});
            set_coeffs_scale({  1.0, 1e+2, 3.0, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0, 1.0});
            set_coeffs_bndu({   inf,  inf, inf, 1.0});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.3, 3.0, 1.0});
            set_coeffs_scale({  1.0, 0.1, 3.0, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0, 1.0});
            set_coeffs_bndu({   inf, inf, inf, 1.0});
        }
    }
    else if (model=="GLMVFT4") // [tau0,T0,a,A]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 1e+2, 3.0, 1.0});
            set_coeffs_scale({ 1e+3, 1e+2, 3.0, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 1e+2, 3.0, 1.0});
            set_coeffs_scale({  1.0, 1e+2, 3.0, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.3, 3.0, 1.0});
            set_coeffs_scale({  1.0, 0.1, 3.0, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf, inf});
        }
    }
    else if (model=="VFT") // [tau0,D,T0]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2,  7.0, 1e+2});
            set_coeffs_scale({ 1e+3,  5.0, 1e+2});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1,  7.0, 1e+2});
            set_coeffs_scale({  1.0,  5.0, 1e+2});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  3.0,  0.3});
            set_coeffs_scale({  1.0,  1.0,  0.1});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
    }
    else if (model=="Mauro") // [tau0,K,C]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 8e+2, 5e+2});
            set_coeffs_scale({ 1e+3, 5e+2, 1e+2});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 8e+2, 5e+2});
            set_coeffs_scale({  1.0, 5e+2, 1e+2});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  0.3,  1.0});
            set_coeffs_scale({  1.0,  0.5,  1.0});
            set_coeffs_bndl({   0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
    }
    else if (model=="AM") // [tau0,p,alpha]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 8e+2, 2.5});
            set_coeffs_scale({ 1e+3, 1e+3, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 8e+2, 2.5});
            set_coeffs_scale({  1.0, 1e+3, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.5, 4.5});
            set_coeffs_scale({  1.0, 0.5, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
    }
    else if (model=="DG") // [tau0,A,B]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 2e+4, 5e-3});
            set_coeffs_scale({ 1e+3, 1e+4, 1e-3});
            set_coeffs_bndl({   0.0,  0.0, -inf});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 2e+4, 5e-3});
            set_coeffs_scale({  1.0, 1e+4, 1e-3});
            set_coeffs_bndl({   0.0,  0.0, -inf});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 1.0,   1.0});
            set_coeffs_scale({  1.0, 1.0,   1.0});
            set_coeffs_bndl({   0.0,  0.0, -inf});
            set_coeffs_bndu({   inf,  inf,  inf});
        }
    }
    else if (model=="MCT") // [A,Tc,r]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 2e+2, 2.5});
            set_coeffs_scale({ 1e+3, 1e+2, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 2e+2, 2.5});
            set_coeffs_scale({  1.0, 1e+2, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.5, 1.0});
            set_coeffs_scale({  1.0, 0.1, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
    }
    else if (model=="SOU") // [A,Tc,r]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 2e+2, 2.5});
            set_coeffs_scale({ 1e+3, 1e+2, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 2e+2, 2.5});
            set_coeffs_scale({  1.0, 1e+2, 1.0});
            set_coeffs_bndl({   0.0,  0.0, 0.0});
            set_coeffs_bndu({   inf,  inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.3, 5.0});
            set_coeffs_scale({  1.0, 0.1, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
    }
    else if (model=="ArrheniusII") // [tau0,A,B]
    {
        is_use_FG=false;
        
        if (systemUnit=="real")
        {
            //
        }
        else if (systemUnit=="metal")
        {
            //
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 1.0, 3.0});
            set_coeffs_scale({  1.0, 1.0, 1.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
    }
    else if (model=="GLM") // [tau0,u20,a]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 3.0, 3.0});
            set_coeffs_scale({ 1e+3, 1.0, 3.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 3.0, 3.0});
            set_coeffs_scale({  1.0, 1.0, 3.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 0.1, 3.0});
            set_coeffs_scale({  1.0, 0.1, 3.0});
            set_coeffs_bndl({   0.0, 0.0, 0.0});
            set_coeffs_bndu({   inf, inf, inf});
        }
    }
    else if (model=="GLM1")
    {
        if (GLM1mode=="mode1"||GLM1mode=="mode3") { // [alpha]
            set_coeffs(      { 3.0});
            set_coeffs_scale({ 3.0});
            set_coeffs_bndl( { 0.0});
            set_coeffs_bndu( { inf});
        } else if (GLM1mode=="mode2") { // [alpha,A]
            set_coeffs(      { 3.0, 1.0});
            set_coeffs_scale({ 3.0, 1.0});
            set_coeffs_bndl( { 0.0, 0.0});
            set_coeffs_bndu( { inf, inf});
        }
    }
    else if (model=="GLM1_hallwolynes") // alpha=2, A=1
    {
        set_coeffs(      { 2.0, 1.0});
        set_coeffs_scale({ 2.0, 1.0});
        set_coeffs_bndl( { 2.0, 1.0});
        set_coeffs_bndu( { 2.0, 1.0});
    }
    else if (model=="univGLM") // [alpha,b,pre]
    {
        if (true) //1-param, b=0.5, pre=1.5
        {
            set_coeffs(      { 3.0,0.5,1.5});
            set_coeffs_scale({ 3.0,0.5,1.5});
            set_coeffs_bndl( {-inf,0.5,1.5});
            set_coeffs_bndu( { inf,0.5,1.5});
        }
        else //3-param, b is shift-param, pre is preFac
        {
            set_coeffs(      { 3.0,0.5,1.5});
            set_coeffs_scale({ 3.0,0.5,1.5});
            set_coeffs_bndl( {-inf,  0,  0});
            set_coeffs_bndu( { inf,inf,inf});
        }
    }
    else if (model=="univGT") // [d]
    {
        set_coeffs(      { 3.0});
        set_coeffs_scale({ 3.0});
        set_coeffs_bndl( {-inf});
        set_coeffs_bndu( { inf});
    }
    else if (model=="COOP") // [tau0,E_inf,u,b]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 1e+3,  1.0,  0.1});
            set_coeffs_scale({ 1e+3, 1e+3,  1.0,  0.1});
            set_coeffs_bndl({   0.0,  0.0, -inf, -inf});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 1e+3,  1.0,  0.1});
            set_coeffs_scale({  1.0, 1e+3,  1.0,  0.1});
            set_coeffs_bndl({   0.0,  0.0, -inf, -inf});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  0.5,  1.0,  0.5});
            set_coeffs_scale({  1.0,  0.5,  1.0,  0.1});
            set_coeffs_bndl({   0.0,  0.0, -inf, -inf});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
    }
    else if (model=="COOP_string") // [tau0,E_inf,u,b]
    {
        set_coeffs({       5.0, 1.0, 10.0,  0.5});
        set_coeffs_scale({ 1.0, 1.0, 10.0,  0.1});
        set_coeffs_bndl({  0.0, 0.0, -inf, -inf});
        set_coeffs_bndu({  inf, inf,  inf,  inf});
    }
    else if (model=="COOP_DWF") // [tau0,E_inf,u,b]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 1.0, 10.0,  0.5});
            set_coeffs_scale({ 1e+3, 1.0, 10.0,  0.1});
            set_coeffs_bndl({   0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({   inf, inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 1.0, 10.0,  0.5});
            set_coeffs_scale({  1.0, 1.0, 10.0,  0.1});
            set_coeffs_bndl({   0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({   inf, inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 1.0, 10.0,  0.5});
            set_coeffs_scale({  1.0, 1.0, 10.0,  0.1});
            set_coeffs_bndl({   0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({   inf, inf,  inf,  inf});
        }
    }
    else if (model=="DEAG") // [tau0,A,B,C]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 1e+2, -1.0, 5e+2});
            set_coeffs_scale({ 1e+3, 1e+2,  1.0, 1e+2});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 1e+2, -1.0, 5e+2});
            set_coeffs_scale({  1.0, 1e+2,  1.0, 1e+2});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1, 1e-2, -1.0,  1.5});
            set_coeffs_scale({  1.0, 1e-2,  1.0,  1.0});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
    }
    else if (model=="CG") // [tau0,B,C,T0]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2, 3e+3,  0.5, 2e+2});
            set_coeffs_scale({ 1e+3, 1e+3,  0.1, 1e+2});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1, 3e+3,  0.5, 2e+2});
            set_coeffs_scale({  1.0, 1e+3,  0.1, 1e+2});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  2.0,  0.0,  0.3});
            set_coeffs_scale({  1.0,  1.0, 1e-2,  0.1});
            set_coeffs_bndl({   0.0,  0.0, -inf,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
    }
    else if (model=="FourParamVFT") // [tau0,D,T0,alpha]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1e+2,  5.0, 1e+2,  2.0});
            set_coeffs_scale({ 1e+3,  1.0, 1e+2,  1.0});
            set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({        0.1,  5.0, 1e+2,  2.0});
            set_coeffs_scale({  1.0,  1.0, 1e+2,  1.0});
            set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({        0.1,  3.0,  0.3,  1.0});
            set_coeffs_scale({  1.0,  5.0,  0.1,  1.0});
            set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
            set_coeffs_bndu({   inf,  inf,  inf,  inf});
        }
    }
    else if (model=="ArrheniusIII") // [tau0,A,B,C]
    {
        is_use_FG=false;
        
        if (systemUnit=="real")
        {
            //
        }
        else if (systemUnit=="metal")
        {
            //
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.1, 1.0, 3.0, 1.0});
            set_coeffs_scale({ 1.0, 1.0, 1.0, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 0.0, 0.0});
            set_coeffs_bndu({  inf, inf, inf, inf});
        }
    }
    else if (model=="HRM1") // [rse]
    {
        is_use_FG=false;
        set_coeffs({      1.0});
        set_coeffs_scale({1.0});
        set_coeffs_bndl({ 0.0});
        set_coeffs_bndu({ inf});
    }
    else
    {
        cout
        << "in FitData::set_fitParams(): \n"
        << "model("<<model<<") not found!\n";
        exit(EXIT_FAILURE);
    }
}





void FitData::fitParams_correction(const std::string& model)
{
    if (model=="KWW"||model=="KWW_lnFs")
    {
        // self-correct tau:
        // 1.15^100=1.17e+6
        
        if (cyclecount>0) {
            coeffs_vD[1] *= 1.15;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2]});
            c = get_coeffs().c_str();
        }
    }
    else if (model=="KWW_full"||model=="KWW_full_lnFs")
    {
        if (cyclecount>0) {
            coeffs_vD[1] *= 1.15;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3],coeffs_vD[4]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3],coeffs_vD[4]});
            c = get_coeffs().c_str();
        }
    }
    else if (model=="KWW_pwr")
    {
        if (cyclecount>0) {
            coeffs_vD[2] *= 1.15;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            c = get_coeffs().c_str();
        }
    }
    else if (model=="mKWW")
    {
        if (cyclecount>0) {
            coeffs_vD[2] *= 1.15;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            c = get_coeffs().c_str();
        }
    }
    else if (model=="mKWW_pwr")
    {
        if (cyclecount>0) {
            coeffs_vD[3] *= 1.15;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3],coeffs_vD[4]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3],coeffs_vD[4]});
            c = get_coeffs().c_str();
        }
    }
    else if (model=="Arrhenius")
    {
        // [tauA,Ea/k]
        
        // self-correct tauA:
        // 1.02^100=7.24
        
        if (cyclecount>0) {
            coeffs_vD[0] *= 1.02;
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2]});
            c = get_coeffs().c_str();
        }
    }
}





void FitData::alglib_lsfit(const std::string& model)
{
    Tc_correction(model);
    
    ////////////////////////////////////////////////////////////////////////
    //======================================================================
    // alglib nonlinear least-square fit routine
    //======================================================================
    x    = fit_xData.c_str();
    y    = fit_yData.c_str();
    c    = get_coeffs().c_str();
    s    = get_coeffs_scale().c_str();
    bndl = get_coeffs_bndl().c_str();
    bndu = get_coeffs_bndu().c_str();
    
    cyclecount=0;
    usedt.push_back(time(NULL));
    do
    {
        fitParams_correction(model);
        
        alglib_solverobj(x, y, c, state);
        lsfitsetcond(state, epsf, epsx, maxits);
        lsfitsetbc(state, bndl, bndu);
        lsfitsetscale(state, s);
        alglib_lsfit_form(model, state);
        lsfitresults(state, info, c, rep);
        
        ++cyclecount;
        
    } while((rep.r2<0.99)&&(cyclecount<100));
    usedt.push_back(time(NULL));
    //======================================================================
    ////////////////////////////////////////////////////////////////////////
}





void FitData::build_function_interpolant(const string& function,
                                         const vector<vector<double>>& raw,
                                         const int x_index,
                                         const int y_index)
{
    vector<vector<double>> data;
    for (int i=0; i<(int)raw.size(); ++i) {
        data.push_back
        ({
            raw.at(i).at(x_index),
            raw.at(i).at(y_index)
        });
    }
    set_fitParams(function);
    fitdata_processing_lsfit(data);
    alglib_lsfit(function);
}





void FitData::KWW_func(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       void *ptr)
{
    /*==========================================================================
     F = A*exp(-(t/tau)^beta)
     c[0]=A, c[1]=tau, c[2]=beta
     =========================================================================*/
    double A=c[0],tau=c[1],beta=c[2],time=x[0];
    func = A*exp(-pow(time/tau,beta));
}
void FitData::KWW_grad(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       real_1d_array &grad,
                       void *ptr)
{
    double A=c[0],tau=c[1],beta=c[2],time=x[0];
    func = A*exp(-pow(time/tau,beta));
    
    grad[0] = func/A;
    grad[1] = func*(beta/tau)*pow(time/tau,beta);
    grad[2] = -func*log(time/tau)*pow(time/tau,beta);
}

void FitData::KWW_lnFs_func(const real_1d_array &c,
                            const real_1d_array &x,
                            double &func,
                            void *ptr)
{
    /*==========================================================================
     lnF = lnA-(t/tau)^beta
     c[0]=A, c[1]=tau, c[2]=beta
     =========================================================================*/
    double A=c[0],tau=c[1],beta=c[2],time=x[0];
    func = log(A)-pow(time/tau,beta);
}
void FitData::KWW_lnFs_grad(const real_1d_array &c,
                            const real_1d_array &x,
                            double &func,
                            real_1d_array &grad,
                            void *ptr)
{
    double A=c[0],tau=c[1],beta=c[2],time=x[0];
    func = log(A)-pow(time/tau,beta);
    
    grad[0] = pow(A,-1);
    grad[1] = (beta/tau)*pow(time/tau,beta);
    grad[2] = -pow(time/tau,beta)*log(time/tau);
}





void FitData::KWW_full_func(const real_1d_array &c,
                            const real_1d_array &x,
                            double &func,
                            void *ptr)
{
    /*==========================================================================
     lnF = lnA*(1-exp(-(t/tauf)^(2p)))^(1/p)-(t/tau)^beta
     F   = exp(lnA*(1-exp(-(t/tauf)^(2p)))^(1/p)-(t/tau)^beta)
     
     c[0]=A, c[1]=tau, c[2]=beta, c[3]=tauf, c[4]=p
     =========================================================================*/
    double A=c[0],tau=c[1],beta=c[2],tauf=c[3],p=c[4],time=x[0];
    double fast=pow(1-exp(-pow(time/tauf,2*p)),pow(p,-1));
    double slow=-pow(time/tau,beta);
    double lnFs=log(A)*fast+slow;
    func = exp(lnFs);
}
void FitData::KWW_full_grad(const real_1d_array &c,
                            const real_1d_array &x,
                            double &func,
                            real_1d_array &grad,
                            void *ptr)
{
    //
}

void FitData::KWW_full_lnFs_func(const real_1d_array &c,
                                 const real_1d_array &x,
                                 double &func,
                                 void *ptr)
{
    /*==========================================================================
     lnF = lnA*(1-exp(-(t/tauf)^(2p)))^(1/p)-(t/tau)^beta
     F   = exp(lnA*(1-exp(-(t/tauf)^(2p)))^(1/p)-(t/tau)^beta)
     
     c[0]=A, c[1]=tau, c[2]=beta, c[3]=tauf, c[4]=p
     =========================================================================*/
    double A=c[0],tau=c[1],beta=c[2],tauf=c[3],p=c[4],time=x[0];
    double fast=pow(1-exp(-pow(time/tauf,2*p)),pow(p,-1));
    double slow=-pow(time/tau,beta);
    double lnFs=log(A)*fast+slow;
    func = lnFs;
}
void FitData::KWW_full_lnFs_grad(const real_1d_array &c,
                                 const real_1d_array &x,
                                 double &func,
                                 real_1d_array &grad,
                                 void *ptr)
{
    //
}





void FitData::KWW_pwr_func(const real_1d_array &c,
                           const real_1d_array &x,
                           double &func,
                           void *ptr)
{
    /*==========================================================================
     F = A*(t^-a)*exp(-(t/tau)^beta)
     c[0]=A, c[1]=a, c[2]=tau, c[3]=beta
     =========================================================================*/
    double A=c[0],a=c[1],tau=c[2],beta=c[3],time=x[0];
    func = A*pow(time,-a)*exp(-pow(time/tau,beta));
}
void FitData::KWW_pwr_grad(const real_1d_array &c,
                           const real_1d_array &x,
                           double &func,
                           real_1d_array &grad,
                           void *ptr)
{
    double A=c[0],a=c[1],tau=c[2],beta=c[3],time=x[0];
    func = A*pow(time,-a)*exp(-pow(time/tau,beta));
    
    // gradient (wrt c):
    // verified by Wolfram-Alpha
    grad[0] = func/A;
    grad[1] = -func*log(time);
    grad[2] = func*(beta/tau)*pow(time/tau,beta);
    grad[3] = -func*log(time/tau)*pow(time/tau,beta);
}





void FitData::mKWW_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    /*==========================================================================
     F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
     c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
     =========================================================================*/
    double A=c[0],taulib=c[1],tau=c[2],beta=c[3],time=x[0];
    func = (1-A)*exp(-time/taulib)+A*exp(-pow(time/tau,beta));
}
void FitData::mKWW_grad(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        alglib::real_1d_array &grad,
                        void *ptr)
{
    double A=c[0],taulib=c[1],tau=c[2],beta=c[3],time=x[0];
    func = (1-A)*exp(-time/taulib)+A*exp(-pow(time/tau,beta));
    
    // gradient (wrt c):
    // verified by Wolfram-Alpha
    grad[0] = exp(-pow(time/tau,beta))-exp(-time/taulib);
    grad[1] = (1-A)*exp(-time/taulib)*(time/pow(taulib,2));
    grad[2] = A*exp(-pow(time/tau,beta))*pow(time/tau,beta)*(beta/tau);
    grad[3] = -A*exp(-pow(time/tau,beta))*pow(time/tau,beta)*log(time/tau);
}





void FitData::mKWW_pwr_func(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr)
{
    /*==========================================================================
     F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
     c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
     =========================================================================*/
    double A=c[0],a=c[1],taulib=c[2],tau=c[3],beta=c[4],time=x[0];
    func = (1-A)*exp(-time/taulib)+A*pow(time,-a)*exp(-pow(time/tau,beta));
}
void FitData::mKWW_pwr_grad(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            alglib::real_1d_array &grad,
                            void *ptr)
{
    double A=c[0],a=c[1],taulib=c[2],tau=c[3],beta=c[4],time=x[0];
    func = (1-A)*exp(-time/taulib)+A*pow(time,-a)*exp(-pow(time/tau,beta));
    
    // gradient (wrt c):
    // verified by Wolfram-Alpha
    grad[0] =
    pow(time,-a)*exp(-pow(time/tau,beta))-exp(-time/taulib);
    
    grad[1] =
    -A*pow(time,-a)*log(time)*exp(-pow(time/tau,beta));
    
    grad[2] =
    (1-A)*exp(-time/taulib)*(time/pow(taulib,2));
    
    grad[3] =
    A*(beta/tau)*pow(time,-a)*exp(-pow(time/tau,beta))*pow(time/tau,beta);
    
    grad[4] =
    -A*pow(time,-a)*exp(-pow(time/tau,beta))*pow(time/tau,beta)*log(time/tau);
}





double FitData::sExp_tauFit_interpolateFvalue(const real_1d_array& c,
                                              const string& tcalc)
{
    double tauFit=0;
    
    if (tcalc=="KWW"||tcalc=="KWW_lnFs")
    {
        /*======================================================================
         F = A*exp(-(t/tau)^b)
         t = tau*(ln(A/F))^(1/b) -- analytic solution for t
         
         c[0]=A, c[1]=tau, c[2]=b
         =====================================================================*/
        double A=c[0],tau=c[1],b=c[2];
        tauFit=tau*pow(log(A/get_sExp_tauFit()),pow(b,-1));
    }
    else if (tcalc=="KWW_full"||tcalc=="KWW_full_lnFs")
    {
        cout << "No analytic solution for relxation time from "<<tcalc<<"\n";
        exit(EXIT_FAILURE);
    }
    else if (tcalc=="KWW_pwr")
    {
        /*======================================================================
         F = A*(t^-a)*exp(-(t/tau)^b)
         t = no analytic solution; use Newton's method for root finding
         
         c[0]=A, c[1]=a, c[2]=tau, c[3]=b
         =====================================================================*/
        double A=c[0],a=c[1],tau=c[2],b=c[3];
        /*------------------------------------------
         use Newton's method for root finding
         f  = F - A*(t^-a)*exp(-(t/tau)^b)
         df = A*(t^(-a-1))*exp(-(t/tau)^b)*(a+b*(t/tau)^b)
         ------------------------------------------*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double t_pre=tau;  // initial guess of time
        double t_now=tau;
        double abserr=0;
        int counter=0;
        do {
            double f =
            get_sExp_tauFit()-A*pow(t_pre,-a)*exp(-pow(t_pre/tau,b));
            
            double df =
            A*pow(t_pre,-a-1)*exp(-pow(t_pre/tau,b))*(a+b*pow(t_pre/tau,b));
            
            t_now=t_pre-(f/df);
            abserr=fabs(t_now-t_pre);
            t_pre=t_now;
            ++counter;
        } while( ((abserr>TOL)||(t_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        tauFit=t_now;
    }
    else if (tcalc=="mKWW")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
         t = no analytic solution; use Newton's method for root finding
         
         c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
         =====================================================================*/
        double A=c[0],taulib=c[1],tau=c[2],b=c[3];
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double t_pre=tau;  // initial guess of time
        double t_now=tau;
        double abserr=0;
        int counter=0;
        do {
            double f =
            get_sExp_tauFit() - ((1-A)*exp(-t_pre/taulib)+A*exp(-pow(t_pre/tau,b)));
            
            double df =
            A*exp(-pow(t_pre/tau,b))*(b/tau)*pow(t_pre/tau,b-1) +
            ((1-A)/taulib)*exp(-t_pre/taulib);
            
            t_now=t_pre-(f/df);
            abserr=fabs(t_now-t_pre);
            t_pre=t_now;
            ++counter;
        } while( ((abserr>TOL)||(t_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        tauFit=t_now;
    }
    else if (tcalc=="mKWW_pwr")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
         t = no analytic solution; use Newton's method for root finding
         
         c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
         =====================================================================*/
        double A=c[0],a=c[1],taulib=c[2],tau=c[3],b=c[4];
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double t_pre=tau;  // initial guess of time
        double t_now=tau;
        double abserr=0;
        int counter=0;
        do {
            double f =
            get_sExp_tauFit() -
            ((1-A)*exp(-t_pre/taulib)+A*pow(t_pre,-a)*exp(-pow(t_pre/tau,b)));
            
            double df =
            a*A*pow(t_pre,-a-1)*exp(-pow(t_pre/tau,b)) +
            A*(b/tau)*pow(t_pre,-a)*exp(-pow(t_pre/tau,b))*pow(t_pre/tau,b-1) +
            ((1-A)/taulib)*exp(-t_pre/taulib);
            
            t_now=t_pre-(f/df);
            abserr=fabs(t_now-t_pre);
            t_pre=t_now;
            ++counter;
        } while( ((abserr>TOL)||(t_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        tauFit=t_now;
    }
    else
    {
        cout
        << "in FitData::sExp_tauFit_interpolateFvalue():\n"
        << tcalc << " model not found! Please check.\n"; exit(EXIT_FAILURE);
    }
    return tauFit;
}





double FitData::sExp_tauFit_gammafunction(const real_1d_array& c,
                                          const string& tcalc)
{
    double tauFit=0;
    
    if (tcalc=="KWW"||tcalc=="KWW_lnFs")
    {
        /*======================================================================
         F = A*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=tau_kww, c[2]=beta_kww
         =====================================================================*/
        double tau=c[1],beta=c[2];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="KWW_pwr")
    {
        /*======================================================================
         F = A*(t^-a)*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=a, c[2]=tau_kww, c[3]=beta_kww
         =====================================================================*/
        double tau=c[2],beta=c[3];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="mKWW")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
         =====================================================================*/
        double tau=c[2],beta=c[3];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="mKWW_pwr")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
         =====================================================================*/
        double tau=c[3],beta=c[4];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else
    {
        cout
        << "in FitData::sExp_tauFit_gammafunction():\n"
        << tcalc << " model not found.\n"; exit(EXIT_FAILURE);
    }
    return tauFit;
}





double FitData::sExp_func_value(const real_1d_array& c,
                                const double time,
                                const string& tcalc)
{
    double Fvalue=0;
    
    if (tcalc=="KWW"||tcalc=="KWW_lnFs")
    {
        /*======================================================================
         F = A*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=tau, c[2]=beta
         =====================================================================*/
        double A=c[0],tau=c[1],beta=c[2];
        Fvalue=A*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="KWW_full"||tcalc=="KWW_full_lnFs")
    {
        double A=c[0],tau=c[1],beta=c[2],tauf=c[3],p=c[4];
        double fast=pow(1-exp(-pow(time/tauf,2*p)),pow(p,-1));
        double slow=-pow(time/tau,beta);
        double lnFs=log(A)*fast+slow;
        Fvalue=exp(lnFs);
    }
    else if (tcalc=="KWW_pwr")
    {
        /*======================================================================
         F = A*(t^-a)*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=a, c[2]=tau, c[3]=beta
         =====================================================================*/
        double A=c[0],a=c[1],tau=c[2],beta=c[3];
        Fvalue=A*pow(time,-a)*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="mKWW")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
         =====================================================================*/
        double A=c[0],taulib=c[1],tau=c[2],beta=c[3];
        Fvalue=(1-A)*exp(-time/taulib)+A*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="mKWW_pwr")
    {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
         =====================================================================*/
        double A=c[0],a=c[1],taulib=c[2],tau=c[3],beta=c[4];
        Fvalue=(1-A)*exp(-time/taulib)+A*pow(time,-a)*exp(-pow(time/tau,beta));
    }
    else
    {
        cout
        << "in FitData::in sExp_func_value():\n"
        << tcalc << " model not found.\n"; exit(EXIT_FAILURE);
    }
    return Fvalue;
}






void FitData::linear_func(const alglib::real_1d_array &c,
                          const alglib::real_1d_array &x,
                          double &func,
                          void *ptr)
{
    /*==========================================================================
     linear form:
     y(x) = c[0]*x + c[1]
     
     param: {c[0]=slope, c[1]=intercept}
     =========================================================================*/
    double slope=c[0],intercept=c[1],x_axis=x[0];
    func = slope*x_axis+intercept;
}
void FitData::linear_grad(const alglib::real_1d_array &c,
                          const alglib::real_1d_array &x,
                          double &func,
                          alglib::real_1d_array &grad,
                          void *ptr)
{
    double slope=c[0],intercept=c[1],x_axis=x[0];
    func = slope*x_axis+intercept;
    grad[0]=x_axis;
    grad[1]=1.0;
}





void FitData::linear_func2(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           void *ptr)
{
    /*==========================================================================
     linear form:
     y(x) = c[0]*x - c[1]
     
     param: {c[0]=slope, c[1]=intercept}
     =========================================================================*/
    double slope=c[0],intercept=c[1],x_axis=x[0];
    func = slope*x_axis-intercept;
}
void FitData::linear_grad2(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           alglib::real_1d_array &grad,
                           void *ptr)
{
    double slope=c[0],intercept=c[1],x_axis=x[0];
    func = slope*x_axis-intercept;
    grad[0]=x_axis;
    grad[1]=-1.0;
}





void FitData::linear_func3(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           void *ptr)
{
    /*==========================================================================
     linear form with constraint (x & y intercepts are same):
     y(x) = c[0] - x
     
     param: {c[0]=intercept}
     =========================================================================*/
    double intercept=c[0],x_axis=x[0];
    func = intercept-x_axis;
}
void FitData::linear_grad3(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           alglib::real_1d_array &grad,
                           void *ptr)
{
    double intercept=c[0],x_axis=x[0];
    func = intercept-x_axis;
    grad[0]=1.0;
}





void FitData::univu2T_func(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           void *ptr)
{
    /*==========================================================================
     <u2> = (1.5*(T/TA)-0.5)*<u2>A
     param: {c[0]=d}
     input: {T/TA}
     
     func: (with base 10)
     func = log10((1.5*(T/TA)-0.5)*u2A)
     =========================================================================*/
    double u2A=c[0],slope=c[1],intercept=c[2];
    double Tr=x[0];
    func = (slope*Tr-intercept)*u2A;
}
void FitData::univu2T_grad(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           alglib::real_1d_array &grad,
                           void *ptr)
{
    double u2A=c[0],slope=c[1],intercept=c[2];
    double Tr=x[0];
    func = (slope*Tr-intercept)*u2A;
    grad[0]= slope*Tr-intercept;
}





void FitData::secondpoly_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr)
{
    /*==========================================================================
     linear form:
     y(x) = c[0]+c[1]*x+c[2]*x2
     =========================================================================*/
    func = c[0]+c[1]*x[0]+c[2]*pow(x[0],2);
}
void FitData::secondpoly_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr)
{
    func = c[0]+c[1]*x[0]+c[2]*pow(x[0],2);
    grad[0]=1.0;
    grad[1]=x[0];
    grad[2]=pow(x[0],2);
}





void FitData::exp_string_func(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              void *ptr)
{
    /*==========================================================================
     f(n) = A*exp(-n/B)
     
     param: {c[0]=A, c[1]=B}
     =========================================================================*/
    double A=c[0],B=c[1],n=x[0];
    func = A*exp(-n/B);
}
void FitData::exp_string_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr)
{
    double A=c[0],B=c[1],n=x[0];
    func = A*exp(-n/B);
    grad[0]=func/A;
    grad[1]=(A*n*exp(-n/B))/pow(B,2);
}





void FitData::pwr_exp_string_func(const alglib::real_1d_array &c,
                                  const alglib::real_1d_array &x,
                                  double &func,
                                  void *ptr)
{
    /*==========================================================================
     f(n) = A*(n^-B)*exp(-n/C)
     
     param: {c[0]=A, c[1]=B, c[2]=C}
     =========================================================================*/
    double A=c[0],B=c[1],C=c[2],n=x[0];
    func = A*pow(n,-B)*exp(-n/C);
}
void FitData::pwr_exp_string_grad(const alglib::real_1d_array &c,
                                  const alglib::real_1d_array &x,
                                  double &func,
                                  alglib::real_1d_array &grad,
                                  void *ptr)
{
    double A=c[0],B=c[1],C=c[2],n=x[0];
    func = A*pow(n,-B)*exp(-n/C);
    grad[0]=func/A;
    grad[1]=-A*pow(n,-B)*exp(-n/C)*log(n);
    grad[2]=(A*pow(n,1-B)*exp(-n/C))/pow(C,2);
}





void FitData::Arrhenius_func(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             void *ptr)
{
    /*==========================================================================
     Arrhenius form:
     tau = tauA*exp(Ea/kT)
     param: {c[0]=tauA, c[1]=Ea/k}
     
     func: (with base 10)
     func = log10(tau) = (ln(tauA)+(Ea/kT))*log10(e)
     =========================================================================*/
    double tauA=c[0],Ea=c[1],T=x[0];
    func = (log(tauA)+(Ea/T))*log10(exp(1));
}
void FitData::Arrhenius_grad(const alglib::real_1d_array &c,
                             const alglib::real_1d_array &x,
                             double &func,
                             real_1d_array &grad,
                             void *ptr)
{
    /*==========================================================================
     Arrhenius form:
     tau = tauA*exp(Ea/kT)
     param: {c[0]=tauA, c[1]=Ea/k}
     
     func: (with base 10)
     func = log10(tau) = (ln(tauA)+(Ea/kT))*log10(e)
     
     gradient (wrt c):
     grad|tauA = (1/tauA)*log10(e)
     grad|Ea/K = (1/T)*log10(e)
     =========================================================================*/
    double tauA=c[0],Ea=c[1],T=x[0];
    func = (log(tauA)+(Ea/T))*log10(exp(1));
    grad[0]= (pow(tauA,-1))*log10(exp(1));
    grad[1]= (pow(T,-1))*log10(exp(1));
}





void FitData::string2param_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr)
{
    /** NOTE:
     ** it has 5 input parameters: [T,L(T),LA,tauA,TA] **/
    /*==========================================================================
     functional form:
     tau = tauA*exp((L/LA)*(delH-T*delS)/T-(delH-TA*delS)/TA)
     param: {c[0]=delH, c[1]=delS}
     input: {T,L(T),LA,tauA,TA,taueq}
     
     func: (with base 10)
     func = log10(tau) = (ln(tauA)+((L/LA)*(delH-T*delS)/T-(delH-TA*delS)/TA))*log10(e)
     =========================================================================*/
    double delH=c[0],delS=c[1];
    double T=x[0],L=x[1],LA=x[2],tauA=x[3],TA=x[4];
    double z=L/LA; //z=1.0;
    func = (log(tauA)+(z*(delH-T*delS)*pow(T,-1))-((delH-TA*delS)*pow(TA,-1)))*log10(exp(1));
}
void FitData::string2param_grad(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                real_1d_array &grad,
                                void *ptr)
{
    double delH=c[0],delS=c[1];
    double T=x[0],L=x[1],LA=x[2],tauA=x[3],TA=x[4];
    double z=L/LA; //z=1.0;
    func = (log(tauA)+(z*(delH-T*delS)*pow(T,-1))-((delH-TA*delS)*pow(TA,-1)))*log10(exp(1));
    grad[0]= (z*pow(T,-1)-pow(TA,-1))*log10(exp(1));
    grad[1]= (1.0-z)*log10(exp(1));
}





void FitData::string3param_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr)
{
    /*==========================================================================
     functional form:
     tau = tau0*exp((L/LA)*(delH-T*delS)/T)
     param: {c[0]=tau0, c[1]=delH; c[2]=delS}
     input: {T,L(T),LA}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+((L/LA)*(delH-T*delS)/T))*log10(e)
     =========================================================================*/
    double tau0=c[0],delH=c[1],delS=c[2];
    double T=x[0],L=x[1],LA=x[2];
    double z=L/LA;
    func = (log(tau0)+(z*(delH-T*delS)*pow(T,-1)))*log10(exp(1));
}
void FitData::string3param_grad(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                real_1d_array &grad,
                                void *ptr)
{
    double tau0=c[0],delH=c[1],delS=c[2];
    double T=x[0],L=x[1],LA=x[2];
    double z=L/LA;
    func = (log(tau0)+(z*(delH-T*delS)*pow(T,-1)))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (z*pow(T,-1))*log10(exp(1));
    grad[2]= (-z)*log10(exp(1));
}





void FitData::transitState_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr)
{
    /*==========================================================================
     functional form:
     tau = tau0*exp((delH-T*delS)/T)
     param: {c[0]=delH; c[1]=delS}
     input: {T,tau_vib}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+((delH-T*delS)/T))*log10(e)
     =========================================================================*/
    double delH=c[0],delS=c[1];
    double T=x[0],tau0=x[1];
    func = (log(tau0)+((delH-T*delS)*pow(T,-1)))*log10(exp(1));
}
void FitData::transitState_grad(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                real_1d_array &grad,
                                void *ptr)
{
    double delH=c[0],delS=c[1];
    double T=x[0],tau0=x[1];
    func = (log(tau0)+((delH-T*delS)*pow(T,-1)))*log10(exp(1));
    grad[0]= (pow(T,-1))*log10(exp(1));
    grad[1]= (-1.0)*log10(exp(1));
}





void FitData::hallwolynes_func(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               void *ptr)
{
    /*==========================================================================
     functional form:
     tau = tau0*exp(u20/u2)
     param: {c[0]=tau0, c[1]=u20}
     input: {u2,T,tau}
     =========================================================================*/
    double tau0=c[0],u20=c[1];
    double u2=x[0];
    func = (log(tau0)+(u20/u2))*log10(exp(1));
}
void FitData::hallwolynes_grad(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               real_1d_array &grad,
                               void *ptr)
{
    double tau0=c[0],u20=c[1];
    double u2=x[0];
    func = (log(tau0)+(u20/u2))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (pow(u2,-1))*log10(exp(1));
}





void FitData::hallwolynes1param_func(const alglib::real_1d_array &c,
                                     const alglib::real_1d_array &x,
                                     double &func,
                                     void *ptr)
{
    /*==========================================================================
     functional form:
     tau = tauA*exp[(u20/u2A)*(u2A/u2-1)]
     param: {c[0]=u20}
     input: {u2,u2A,tauA,T,tau}
     =========================================================================*/
    double u20=c[0];
    double u2=x[0],u2A=x[1],tauA=x[2];
    double u20r=u20/u2A;
    double u2r =u2/u2A;
    func = (log(tauA)+u20r*(pow(u2r,-1)-1))*log10(exp(1));
}
void FitData::hallwolynes1param_grad(const alglib::real_1d_array &c,
                                     const alglib::real_1d_array &x,
                                     double &func,
                                     real_1d_array &grad,
                                     void *ptr)
{
    double u20=c[0];
    double u2=x[0],u2A=x[1],tauA=x[2];
    double u20r=u20/u2A;
    double u2r =u2/u2A;
    func = (log(tauA)+u20r*(pow(u2r,-1)-1))*log10(exp(1));
    grad[0]= (pow(u2A,-1)*(pow(u2r,-1)-1))*log10(exp(1));
}





double FitData::get_Xr_leporini(const alglib::real_1d_array& c,
                                const double logtaur)
{
    double A=c[0],B=c[1],C=c[2];
    return (-B+sqrt(pow(B,2)-4*A*C+4*C*logtaur))*pow(2*C,-1);
}
void FitData::leporiniref_func(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               void *ptr)
{
    /*==========================================================================
     functional form:
     tau = 10^(A+B*Xr*Y+C*Xr2*Y2)
     param: {c[0]=A, c[1]=B=Const.; c[2]=C=Const.}
     input: {logtauref,u2ref,u2,T,tau}
     
     func: (with base 10)
     func = log10(tau) = A+B*Xr*Y+C*Xr2*Y2
     =========================================================================*/
    double A=c[0],B=c[1],C=c[2];
    double logtaur=x[0],u2r=x[1],u2=x[2];
    double Y=u2r/u2;
    double Xr=(-B+sqrt(pow(B,2)-4.0*A*C+4.0*C*logtaur))/(2.0*C);
    double X=Xr*Y;
    func = A+B*X+C*pow(X,2);
}
void FitData::leporiniref_grad(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               real_1d_array &grad,
                               void *ptr)
{
    double A=c[0],B=c[1],C=c[2];
    double logtaur=x[0],u2r=x[1],u2=x[2];
    double Y=u2r/u2;
    double Xr=(-B+sqrt(pow(B,2)-4.0*A*C+4.0*C*logtaur))/(2.0*C);
    double X=Xr*Y;
    func = A+B*X+C*pow(X,2);
    grad[0]= 1.0;
    grad[1]= X;
    grad[2]= pow(X,2);
}





void FitData::leporini_universal_func(const alglib::real_1d_array &c,
                                      const alglib::real_1d_array &x,
                                      double &func,
                                      void *ptr)
{
    /*==========================================================================
     functional form:
     tau = 10^(A+B*Xr*Y+C*Xr2*Y2)
     param: {c[0]=A=-0.424, c[1]=1.6285; c[2]=12.325}
     input: {logtaur,u2r,u2}
     
     func: (with base 10)
     func = log10(tau) = A+B*Xr*Y+C*Xr2*Y2
     =========================================================================*/
    double A=c[0],B=c[1],C=c[2];
    double logtaur=x[0],u2r=x[1],u2=x[2];
    double Y=u2r*pow(u2,-1);
    double Xr=(-B+sqrt(pow(B,2)-4*A*C+4*C*logtaur))*pow(2*C,-1);
    func = A+B*Xr*Y+C*pow(Xr,2)*pow(Y,2);
}
void FitData::leporini_universal_grad(const alglib::real_1d_array &c,
                                      const alglib::real_1d_array &x,
                                      double &func,
                                      real_1d_array &grad,
                                      void *ptr)
{
    double A=c[0],B=c[1],C=c[2];
    double logtaur=x[0],u2r=x[1],u2=x[2];
    double Y=u2r*pow(u2,-1);
    double Xr=(-B+sqrt(pow(B,2)-4*A*C+4*C*logtaur))*pow(2*C,-1);
    func = A+B*Xr*Y+C*pow(Xr,2)*pow(Y,2);
    grad[0]= 1.0;
    grad[1]= Xr*Y;
    grad[2]= pow(Xr,2)*pow(Y,2);
}





void FitData::MCT_func(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       void *ptr)
{
    /*==========================================================================
     MCT form:
     tau = A*((T-Tc)/Tc)^(-r)
     param: {c[0]=A, c[1]=Tc, c[2]=r}
     
     func: (with base 10)
     func = log10(tau) = log10(A)-r*log10((T/Tc)-1)
     =========================================================================*/
    double A=c[0],Tc=c[1],r=c[2],T=x[0];
    func = log10(A)-r*log10((T/Tc)-1.0);
}
void FitData::MCT_grad(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       real_1d_array &grad,
                       void *ptr)
{
    double A=c[0],Tc=c[1],r=c[2],T=x[0];
    func = log10(A)-r*log10((T/Tc)-1.0);
    grad[0]= pow(A*log(10),-1);
    grad[1]= (r*T)/(T*Tc*log(10)-log(10)*pow(Tc,2));
    grad[2]= -log((T/Tc)-1.0)/log(10);
}





void FitData::SOU_func(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       void *ptr)
{
    /*==========================================================================
     SOU form:
     tau = A*((T-Tc)/T)^(-r)
     param: {c[0]=A, c[1]=Tc, c[2]=r}
     
     func: (with base 10)
     func = log10(tau) = log10(A)-r*log10(1-(Tc/T))
     =========================================================================*/
    double A=c[0],Tc=c[1],r=c[2],T=x[0];
    func = log10(A)-r*log10(1.0-(Tc/T));
}
void FitData::SOU_grad(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       real_1d_array &grad,
                       void *ptr)
{
    double A=c[0],Tc=c[1],r=c[2],T=x[0];
    func = log10(A)-r*log10(1.0-(Tc/T));
    grad[0]= pow(A*log(10),-1);
    grad[1]= r/(log(10)*(T-Tc));
    grad[2]= -log(1.0-(Tc/T))/log(10);
}





void FitData::GLMVFT_func(const real_1d_array &c,
                          const real_1d_array &x,
                          double &func,
                          void *ptr)
{
    /*==========================================================================
     tau = tau0*exp(A*[(2*T0)/(T-T0)]^(a/2))
     param: {c[0]=tau0, c[1]=T0, c[2]=a, c[3]=A}
     input: {T}
     =========================================================================*/
    double tau0=c[0],T0=c[1],a=c[2],A=c[3],T=x[0];
    func = (log(tau0)+A*pow((2*T0)/(T-T0),a/2))*log10(exp(1));
}
void FitData::GLMVFT_grad(const real_1d_array &c,
                          const real_1d_array &x,
                          double &func,
                          real_1d_array &grad,
                          void *ptr)
{
    double tau0=c[0],T0=c[1],a=c[2],A=c[3],T=x[0];
    func = (log(tau0)+A*pow((2*T0)/(T-T0),a/2))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (a*pow(2,(a/2)-1)*pow(T0/(T-T0),(a/2)-1)*((T0*pow(T-T0,-2))+pow(T-T0,-1)))*log10(exp(1));
    grad[2]= ((pow(2,(a/2)-1)*pow(T0/(T-T0),a/2))*(log(T0/(T-T0))+log(2)))*log10(exp(1));
    grad[3]= (pow((2*T0)/(T-T0),a/2))*log10(exp(1));
}





void FitData::VFT_func(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       void *ptr)
{
    /*==========================================================================
     VFT form:
     tau = tau0*exp((D*T0)/(T-T0))
     param: {c[0]=tau0, c[1]=D, c[2]=T0}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+((D*T0)/(T-T0)))*log10(e)
     =========================================================================*/
    double tau0=c[0],D=c[1],T0=c[2],T=x[0];
    func = (log(tau0)+((D*T0)/(T-T0)))*log10(exp(1));
}
void FitData::VFT_grad(const real_1d_array &c,
                       const real_1d_array &x,
                       double &func,
                       real_1d_array &grad,
                       void *ptr)
{
    /*==========================================================================
     VFT form:
     tau = tau0*exp((D*T0)/(T-T0))
     param: {c[0]=tau0, c[1]=D, c[2]=T0}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+((D*T0)/(T-T0)))*log10(e)
     
     gradient (wrt c):
     grad|tau0 = (1/tau0)*log10(e)
     grad|D    = (T0/(T-T0))*log10(e)
     grad|T0   = (D*T/(T-T0)^2)*log10(e)
     =========================================================================*/
    double tau0=c[0],D=c[1],T0=c[2],T=x[0];
    func = (log(tau0)+((D*T0)/(T-T0)))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (T0/(T-T0))*log10(exp(1));
    grad[2]= ((D*T)/pow(T-T0,2))*log10(exp(1));
}





void FitData::Mauro_func(const alglib::real_1d_array &c,
                         const alglib::real_1d_array &x,
                         double &func,
                         void *ptr)
{
    /*==========================================================================
     Double exponential form:
     tau = tau0*exp((K/T)exp(C/T))
     param: {c[0]=tau0, c[1]=K, c[2]=C}
     
     ***************************************************************************
     Note:
     This functional form is derived from the Adam-Gibbs equation.
     (Mauro et al. 19780-19784 PNAS vol.106 no.47)
     
     --> S = f*Nkln(omega)
     
     S     == configurational entropy
     f     == topological degrees of freedom per atom
     omega == # of degenerate configurations per floppy mode
     
     --> f=3*exp(-H/kT)
     
     H == energy difference of intact and broken of a two-state constrained
     network
     
     Define:
     K=B/3Nkln(omega)
     C=H/K
     
     --> tau = tau0*exp((K/T)exp(C/T)) -- <double exponential form>
     ***************************************************************************
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(K/T)exp(C/T))*log10(e)
     =========================================================================*/
    double tau0=c[0],K=c[1],C=c[2],T=x[0];
    func = (log(tau0)+(K/T)*exp(C/T))*log10(exp(1));
}
void FitData::Mauro_grad(const alglib::real_1d_array &c,
                         const alglib::real_1d_array &x,
                         double &func,
                         real_1d_array &grad,
                         void *ptr)
{
    /*==========================================================================
     Double exponential form:
     tau = tau0*exp((K/T)exp(C/T))
     param: {c[0]=tau0, c[1]=K, c[2]=C}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(K/T)exp(C/T))*log10(e)
     
     gradient (wrt c):
     
     grad|tau0 = (1/tau0)*log(e)
     grad|K    = ((1/T)exp(C/T))*log(e)
     grad|C    = ((K/T^2)exp(C/T))*log(e)
     =========================================================================*/
    double tau0=c[0],K=c[1],C=c[2],T=x[0];
    func = (log(tau0)+(K/T)*exp(C/T))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (pow(T,-1)*exp(C/T))*log10(exp(1));
    grad[2]= (K*pow(T,-2)*exp(C/T))*log10(exp(1));
}





void FitData::AM_func(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      void *ptr)
{
    /*==========================================================================
     AM form:
     tau = tau0*exp((p/T)^alpha)
     param:{c[0]=tau0, c[1]=p, c[2]=alpha}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(p/T)^alpha))*log10(e)
     =========================================================================*/
    double tau0=c[0],p=c[1],alpha=c[2],T=x[0];
    func = (log(tau0)+pow(p/T,alpha))*log10(exp(1));
}
void FitData::AM_grad(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      real_1d_array &grad,
                      void *ptr)
{
    double tau0=c[0],p=c[1],alpha=c[2],T=x[0];
    func = (log(tau0)+pow(p/T,alpha))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= ((alpha/p)*pow(p/T,alpha))*log10(exp(1));
    grad[2]= (pow(p/T,alpha)*log(p/T))*log10(exp(1));
}





void FitData::DG_func(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      void *ptr)
{
    /*==========================================================================
     DG form:
     tau = tau0*exp((A/T)*exp(-B*T))
     param:{c[0]=tau0, c[1]=A, c[2]=B}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(A/T)*exp(-B*T))*log10(e)
     =========================================================================*/
    double tau0=c[0],A=c[1],B=c[2],T=x[0];
    func = (log(tau0)+(A/T)*exp(-B*T))*log10(exp(1));
}
void FitData::DG_grad(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      alglib::real_1d_array &grad,
                      void *ptr)
{
    double tau0=c[0],A=c[1],B=c[2],T=x[0];
    func = (log(tau0)+(A/T)*exp(-B*T))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (pow(T,-1)*exp(-B*T))*log10(exp(1));
    grad[2]= (-A*exp(-B*T))*log10(exp(1));
}





void FitData::ArrheniusII_func(const alglib::real_1d_array &c,
                               const alglib::real_1d_array &x,
                               double &func,
                               void *ptr)
{
    /*==========================================================================
     ArrheniusII form:
     tau = tau0*exp((A/T)+(B/T^2))
     param: {c[0]=tau0, c[1]=A, c[2]=B}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(A/T)+(B/T^2))*log10(e)
     =========================================================================*/
    func = (log(c[0])+(c[1]/x[0])+(c[2]/pow(x[0],2)))*log10(exp(1));
}





void FitData::GLM_func(const alglib::real_1d_array &c,
                       const alglib::real_1d_array &x,
                       double &func,
                       void *ptr)
{
    /*==========================================================================
     Generalized Localization Model:
     tau = tau0*exp((u20/u2)^(a/2))
     param: {c[0]=tau0, c[1]=u20, c[2]=a}
     input: {tauA,u2A,u2,T,tau}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(u20/u2)^(a/2))*log10(e)
     =========================================================================*/
    double tau0=c[0],u20=c[1],a=c[2];
    double u2=x[2];
    func = (log(tau0)+pow(u20/u2,a/2.0))*log10(exp(1));
}
void FitData::GLM_grad(const alglib::real_1d_array &c,
                       const alglib::real_1d_array &x,
                       double &func,
                       alglib::real_1d_array &grad,
                       void *ptr)
{
    double tau0=c[0],u20=c[1],a=c[2];
    double u2=x[2];
    func = (log(tau0)+pow(u20/u2,a/2.0))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= ((a/(2*u20))*pow(u20/u2,a/2.0))*log10(exp(1));
    grad[2]= (0.5*pow(u20/u2,a/2.0)*log(u20/u2))*log10(exp(1));
}





void FitData::GLM1_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    string GLM1mode="mode1";//check Fitdata initializer "GLM1mode" correctly synced
    /*==========================================================================
     mode1: assume u20=u2A
     tau =  tauA*exp((u2A/u2)^(a/2)-1)
     
     mode2: treat A as free param
     tau =  tauA*exp(A^(a/2)*[(u2A/u2)^(a/2)-1])
     A   =  u20/u2A
     
     mode3: introduce highT activation
     tau =  tauA*exp((EA/akT)*[(u2A/u2)^(a/2)-1])
     
     param: {c[0]=a,c[1]=A}
     input: {tauA,u2A,u2,Ea,TA,T,tau}
     =========================================================================*/
    double a=0,A=0;
    if (GLM1mode=="mode1"){a=c[0];A=1;}
    if (GLM1mode=="mode2"){a=c[0];A=pow(c[1],c[0]/2);}
    if (GLM1mode=="mode3"){a=c[0];A=4*x[3]*pow(c[0]*x[4],-1);}
    double tauA=x[0],u2A=x[1],u2=x[2];
    func = (log(tauA)+A*(pow(u2A/u2,a/2)-1))*log10(exp(1));
}
void FitData::GLM1_grad(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        alglib::real_1d_array &grad,
                        void *ptr)
{
    string GLM1mode="mode1";//check Fitdata initializer "GLM1mode" correctly synced
    
    double a=0,A=0;
    if (GLM1mode=="mode1"){a=c[0];A=1;}
    if (GLM1mode=="mode2"){a=c[0];A=pow(c[1],c[0]/2);}
    if (GLM1mode=="mode3"){a=c[0];A=4*x[3]*pow(c[0]*x[4],-1);}
    double tauA=x[0],u2A=x[1],u2=x[2];
    func = (log(tauA)+A*(pow(u2A/u2,a/2)-1))*log10(exp(1));
    
    if (GLM1mode=="mode1"||GLM1mode=="mode3") {
        grad[0]= (0.5*A*pow(u2A/u2,a/2)*log(u2A/u2))*log10(exp(1));
    } else if (GLM1mode=="mode2") {
        grad[0]= (0.5*A*pow(u2A/u2,a/2)*log(u2A/u2))*log10(exp(1));
        grad[1]= (0.5*(pow(u2A/u2,a/2)-1)*pow(c[1],a/2)*log(c[1]))*log10(exp(1));
    }
}





void FitData::univGLM_func(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           void *ptr)
{
    /*==========================================================================
     tau = tauA*exp((1.5*(T/TA)-0.5)^(-a/2)-1)
     param: {c[0]=a}
     input: {tauA,TA,T}
     
     func: (with base 10)
     func = log10(tau) = (ln(tauA)+(1.5*(T/TA)-0.5)^(-a/2)-1)*log10(e)
     =========================================================================*/
    double a=c[0],b=c[1],pre=c[2];
    double tauA=x[0],TA=x[1],T=x[2];
    double Tr=T/TA;
    func = (log(tauA)+pow(pre*Tr-b,-a/2.0)-1.0)*log10(exp(1));
}
void FitData::univGLM_grad(const alglib::real_1d_array &c,
                           const alglib::real_1d_array &x,
                           double &func,
                           alglib::real_1d_array &grad,
                           void *ptr)
{
    double a=c[0],b=c[1],pre=c[2];
    double tauA=x[0],TA=x[1],T=x[2];
    double Tr=T/TA;
    func = (log(tauA)+pow(pre*Tr-b,-a/2.0)-1.0)*log10(exp(1));
    grad[0]= (-0.5*pow(1.5*Tr-b,-a/2.0)*log(1.5*Tr-b))*log10(exp(1));
}





void FitData::univGT_func(const alglib::real_1d_array &c,
                          const alglib::real_1d_array &x,
                          double &func,
                          void *ptr)
{
    /*==========================================================================
     ln(taur)=(u2r)^(-d/2)+ln(tauA/tau0)*(1/Tr-1)-1
     param: {c[0]=<u2>A,c[1]=slope,c[2]=intercept}
     input: {TA,tauA,u2r,tau0}
     =========================================================================*/
    double d=c[0];
    double Tr=x[0],u2r=x[1],tauA=x[2],tau0=x[3];
    func = (log(tauA)+pow(u2r,-d/2)+log(tauA/tau0)*(pow(Tr,-1)-1)-1)*log(exp(1));
}
void FitData::univGT_grad(const alglib::real_1d_array &c,
                          const alglib::real_1d_array &x,
                          double &func,
                          alglib::real_1d_array &grad,
                          void *ptr)
{
    double d=c[0];
    double Tr=x[0],u2r=x[1],tauA=x[2],tau0=x[3];
    func = (log(tauA)+pow(u2r,-d/2)+log(tauA/tau0)*(pow(Tr,-1)-1)-1)*log(exp(1));
    grad[0]= -0.5*pow(u2r,-d/2)*log(u2r);
}





void FitData::COOP_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    /*==========================================================================
     COOP (Decomposition of activation E(T) into E_inf+E_coop(T)) form:
     tau = tau0*exp(E_inf*(1+exp(-u*((T/E_inf)-b)))/T)
     
     tau = tau0*exp(E(T)/T)
     E(T)= E_inf+E_coop(T)
     E_coop(T)=E_inf*exp(-u*((T/E_inf)-b))
     param:{c[0]=tau0, c[1]=E_inf, c[2]=u, c[3]=b}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+E(T)/T)*log10(e)
     =========================================================================*/
    double tau0=c[0],E_inf=c[1],u=c[2],b=c[3],T=x[0];
    double exponent=exp(-u*((T/E_inf)-b));
    double E_coop=E_inf*exponent;
    double E=E_inf+E_coop;
    func = (log(tau0)+E/T)*log10(exp(1));
}
void FitData::COOP_grad(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        alglib::real_1d_array &grad,
                        void *ptr)
{
    double tau0=c[0],E_inf=c[1],u=c[2],b=c[3],T=x[0];
    double exponent=exp(-u*((T/E_inf)-b));
    double E_coop=E_inf*exponent;
    double E=E_inf+E_coop;
    func = (log(tau0)+E/T)*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= ((((T*u*exponent)/E_inf)+exponent+1)/T)*log10(exp(1));
    grad[2]= ((E_inf*(b-T/E_inf)*exponent)/T)*log10(exp(1));
    grad[3]= ((E_inf*u*exponent)/T)*log10(exp(1));
}





void FitData::DEAG_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    /*==========================================================================
     DEAG (Double Exponential from Adam-Gibbs) form:
     tau = tau0*exp(((A-B*T)/T)*exp(C/T))
     param:{c[0]=tau0, c[1]=A, c[2]=B, c[3]=C}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+((A-B*T)/T)*exp(C/T))*log10(e)
     =========================================================================*/
    double tau0=c[0],A=c[1],B=c[2],C=c[3],T=x[0];
    func = (log(tau0)+((A-B*T)/T)*exp(C/T))*log10(exp(1));
}
void FitData::DEAG_grad(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        real_1d_array &grad,
                        void *ptr)
{
    double tau0=c[0],A=c[1],B=c[2],C=c[3],T=x[0];
    func = (log(tau0)+((A-B*T)/T)*exp(C/T))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (pow(T,-1)*exp(C/T))*log10(exp(1));
    grad[2]= (-exp(C/T))*log10(exp(1));
    grad[3]= ((A-B*T)*pow(T,-2)*exp(C/T))*log10(exp(1));
}





void FitData::CG_func(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      void *ptr)
{
    /*==========================================================================
     CG form:
     tau = tau0*exp(B/{(T-T0)+[(T-T0)^2+C*T]^0.5})
     param:{c[0]=tau0, c[1]=B, c[2]=C, c[3]=T0}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+B/{(T-T0)+[(T-T0)^2+C*T]}^0.5)*log10(e)
     =========================================================================*/
    double tau0=c[0],B=c[1],C=c[2],T0=c[3],T=x[0];
    func = (log(tau0)+B/((T-T0)+sqrt(pow(T-T0,2)+C*T)))*log10(exp(1));
}
void FitData::CG_grad(const alglib::real_1d_array &c,
                      const alglib::real_1d_array &x,
                      double &func,
                      real_1d_array &grad,
                      void *ptr)
{
    double tau0=c[0],B=c[1],C=c[2],T0=c[3],T=x[0];
    func = (log(tau0)+B/((T-T0)+sqrt(pow(T-T0,2)+C*T)))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (pow((T-T0)+sqrt(pow(T-T0,2)+C*T),-1))*log10(exp(1));
    grad[2]= (-B*T/(2*sqrt(pow(T-T0,2)+C*T)*pow(sqrt(pow(T-T0,2)+C*T)+T-T0,2)))*log10(exp(1));
    grad[3]= (-B*(-((T-T0)/sqrt(pow(T-T0,2)+C*T))-1)/pow(sqrt(pow(T-T0,2)+C*T)+T-T0,2))*log10(exp(1));
}





void FitData::FourParamVFT_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr)
{
    /*==========================================================================
     Four-Param VFT form:
     tau = tau0*exp(((D*T0)/(T-T0))^alpha)
     param: {c[0)=tau0, c[1]=D, c[2]=T0, c[3]=alpha}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+[(D*T0)/(T-T0)]^alpha)*log10(e)
     =========================================================================*/
    double tau0=c[0],D=c[1],T0=c[2],alpha=c[3],T=x[0];
    func = (log(tau0)+pow((D*T0)/(T-T0),alpha))*log10(exp(1));
}
void FitData::FourParamVFT_grad(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                real_1d_array &grad,
                                void *ptr)
{
    double tau0=c[0],D=c[1],T0=c[2],alpha=c[3],T=x[0];
    func = (log(tau0)+pow((D*T0)/(T-T0),alpha))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= ((alpha/D)*pow((D*T0)/(T-T0),alpha))*log10(exp(1));
    grad[2]= (((alpha*T)/(D*pow(T0,2)))*pow((D*T0)/(T-T0),alpha+1))*log10(exp(1));
    grad[3]= (pow((D*T0)/(T-T0),alpha)*log((D*T0)/(T-T0)))*log10(exp(1));
}





void FitData::ArrheniusIII_func(const alglib::real_1d_array &c,
                                const alglib::real_1d_array &x,
                                double &func,
                                void *ptr)
{
    /*==========================================================================
     ArrheniusIII form:
     tau = tau0*exp((A/T)+(B/T^2)+(C/T^3))
     param: {c[0]=tau0, c[1]=A, c[2]=B, c[3]=C}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(A/T)+(B/T^2)+(C/T^3))*log10(e)
     =========================================================================*/
    func=
    (log(c[0])+(c[1]/x[0])+(c[2]/pow(x[0],2.0))+(c[3]/pow(x[0],3.0)))*log10(exp(1));
}





void FitData::HRM1_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    double rse=c[0],n=x[0];
    func = log10(n*(1-(1-pow(n,-1))*(1-pow(rse,-1))));
}





double FitData::calc_y_given_x(const real_1d_array& c,
                               const double x,
                               const string& model)
{
    if (model=="linear"||model=="diagonal") {
        // y(x) = slope*x + intercept
        // param: {c[0]=slope, c[1]=intercept}
        return c[0]*x+c[1];
    }
    else if (model=="linear2") {
        // y(x) = slope*x - intercept
        // param: {c[0]=slope, c[1]=intercept}
        return c[0]*x-c[1];
    }
    else if (model=="linear3") {
        // y(x) = c[0] - x
        // param: {c[0]=intercept}
        return c[0]-x;
    }
    else if (model=="univu2T") {
        // y(x) = (1.5*x-0.5)*u2A
        // param: {c[0]=u2A, c[1]=slope, c[2]=intercept}
        return (c[1]*x-c[2])*c[0];
    }
    else if (model=="secondpoly") {
        // y(x) = c[0]+c[1]*x+c[2]*x2
        return c[0]+c[1]*x+c[2]*pow(x,2);
    }
    else if (model=="exp_string") {
        // f(n) = A*exp(-n/B)
        // param: {c[0]=A, c[1]=B}
        return c[0]*exp(-x/c[1]);
    }
    else if (model=="pwr_exp_string") {
        // f(n) = A*pow(n,-B)*exp(-n/C);
        // param: {c[0]=A,c[1]=B,c[2]=C}
        return c[0]*pow(x,-c[1])*exp(-x/c[2]);
    }
    else if (model=="Arrhenius"||model=="Arrhenius_string") {
        // tau = tauA*exp(Ea/kT)
        // param: {c[0]=tauA, c[1]=Ea/k}
        return c[0]*exp(c[1]/x);
    }
    else if (model=="MCT") {
        // tau = A*((T-Tc)/Tc)^(-r)
        // param: {c[0]=A, c[1]=Tc, c[2]=r}
        return c[0]*pow((x/c[1])-1.0,-c[2]);
    }
    else if (model=="SOU") {
        // tau = A*((T-Tc)/T)^(-r)
        // param: {c[0]=A, c[1]=Tc, c[2]=r}
        return c[0]*pow(1.0-(c[1]/x),-c[2]);
    }
    else if (model=="GLMVFT1"||model=="GLMVFT3"||model=="GLMVFT4") {
        // tau = tau0*exp([(2*T0)/(T-T0)]^(a/2))
        // param: {c[0]=tau0,c[1]=T0,c[2]=a}
        return c[0]*exp(c[3]*pow((2*c[1])/(x-c[1]),c[2]/2));
    }
    else if (model=="VFT") {
        // tau = tau0*exp((D*T0)/(T-T0))
        // param: {c[0]=tau0, c[1]=D, c[2]=T0}
        return c[0]*exp((c[1]*c[2])/(x-c[2]));
    }
    else if (model=="Mauro") {
        // tau = tau0*exp((K/T)exp(C/T))
        // param: {c[0]=tau0, c[1]=K, c[2]=C}
        return c[0]*exp((c[1]/x)*exp(c[2]/x));
    }
    else if (model=="AM") {
        // tau = tau0*exp((p/T)^alpha)
        // param:{c[0]=tau0, c[1]=p, c[2]=alpha}
        return c[0]*exp(pow((c[1]/x),c[2]));
    }
    else if (model=="DG") {
        // tau = tau0*exp((A/T)*exp(-B*T))
        // param:{c[0]=tau0, c[1]=A, c[2]=B}
        return c[0]*exp((c[1]/x)*exp(-c[2]*x));
    }
    else if (model=="ArrheniusII") {
        // tau = tau0*exp((A/T)+(B/T^2))
        // param: {c[0]=tau0, c[1]=A, c[2]=B}
        return c[0]*exp((c[1]/x)+(c[2]/pow(x,2.0)));
    }
    else if ((model=="COOP")||(model=="COOP_string")||(model=="COOP_DWF")) {
        // tau = tau0*exp(E_inf*(1+exp(-u*((T/E_inf)-b)))/T)
        // param:{c[0]=tau0, c[1]=E_inf, c[2]=u, c[3]=b}
        return c[0]*exp((c[1]*(1+exp(-c[2]*((x/c[1])-c[3]))))/x);
    }
    else if (model=="DEAG") {
        // tau = tau0*exp(((A-B*T)/T)*exp(C/T))
        // param:{c[0]=tau0, c[1]=A, c[2]=B, c[3]=C}
        return c[0]*exp(((c[1]-c[2]*x)/x)*exp(c[3]/x));
    }
    else if (model=="CG") {
        // tau = tau0*exp(D*T0/{(T-T0)+[(T-T0)^2+C*T]^0.5})
        // param:{c[0]=tau0, c[1]=B, c[2]=C, c[3]=T0}
        return c[0]*exp(c[1]/((x-c[3])+pow(pow(x-c[3],2)+c[2]*x,0.5)));
    }
    else if (model=="FourParamVFT") {
        // tau = tau0*exp([(D*T0)/(T-T0)]^alpha)
        // param: {c[0]=tau0, c[1]=D, c[2]=T0, c[3]=alpha}
        return c[0]*exp(pow((c[1]*c[2])/(x-c[2]),c[3]));
    }
    else if (model=="ArrheniusIII") {
        // tau = tau0*exp((A/T)+(B/T^2)+(C/T^3))
        // param: {c[0]=tau0, c[1]=A, c[2]=B, c[3]=C}
        return c[0]*exp((c[1]/x)+(c[2]/pow(x,2.0))+(c[3]/pow(x,3.0)));
    }
    else
    {
        cout
        << "in FitData::calc_y_given_x():\n"
        << "model provided ("<<model<<") not found.\n";
        return 0; exit(EXIT_FAILURE);
    }
}





double FitData::calc_time_given_vecinput(const alglib::real_1d_array& c,
                                         const std::vector<double>& v,
                                         const std::string& model)
{
    if (model=="string2param") {
        // tau = tauA*exp((L/LA)*(delH-T*delS)/T-(delH-TA*delS)/TA)
        // param: {c[0]=delH, c[1]=delS}
        // v={Teqtau,LT,LA,tauA,TA,taueq}
        return v[3]*exp((v[1]/v[2])*(c[0]-v[0]*c[1])*pow(v[0],-1)-(c[0]-v[4]*c[1])*pow(v[4],-1));
    }
    else if (model=="string3param") {
        // tau = tau0*exp((L/LA)*(delH-T*delS)/T)
        // param: {c[0]=tau0, c[1]=delH; c[2]=delS}
        // v={Teqtau,LT,LA,taueq}
        return c[0]*exp((v[1]/v[2])*(c[1]-v[0]*c[2])*pow(v[0],-1));
    }
    else if (model=="transitState") {
        // tau = tau0*exp((delH-T*delS)/T)
        // param: {c[0]=delH; c[1]=delS}
        // v={Teqtau,tau0_vib}
        return v[1]*exp((c[0]-v[0]*c[1])*pow(v[0],-1));
    }
    else if (model=="GLM") {
        // tau = tau0*exp((u20/u2)^(a/2))
        // param: {c[0]=tau0, c[1]=u20, c[2]=a}
        // v={logtauref,u2ref,u2,T,tau,tauA,u2A}
        return c[0]*exp(pow(c[1]/v[2],c[2]/2.0));
    }
    else if (model=="GLM1"||model=="GLM1_hallwolynes") {
        // tau = tauA*exp(A*[(u2A/u2)^(a/2)-1])
        // param: {c[0]=a,[c[1]=A]}
        // v={logtauref,u2ref,u2,T,tau,tauA,u2A,EA,TA}
        if (GLM1mode=="mode1") {
            return v[5]*exp(pow(v[6]/v[2],c[0]/2)-1.0);
        } else if (GLM1mode=="mode2") {
            return v[5]*exp(pow(c[1],c[0]/2)*(pow(v[6]/v[2],c[0]/2)-1.0));
        } else if (GLM1mode=="mode3") {
            return v[5]*exp((4*v[7]*pow(c[0]*v[8],-1))*(pow(v[6]/v[2],c[0]/2)-1.0));
        } else {
            cout
            << "in FitData::calc_time_given_vecinput(): mode used for GLM1 ("
            << GLM1mode << ") is not found.\n"; exit(EXIT_FAILURE);
        }
    }
    else if (model=="univGLM") {
        // tau = tauA*exp((pre*(T/TA)-b)^(-a/2)-1)
        // param: {c[0]=a,c[1]=b,c[2]=pre}
        // v=({tauA,TA,T,tau}
        return v[0]*exp(pow(c[2]*(v[2]/v[1])-c[1],-c[0]/2.0)-1.0);
    }
    else if (model=="univGT") {
        // tau = tauA*exp((pre*(T/TA)-b)^(-a/2)-1)
        // param: {c[0]=a,c[1]=b,c[2]=pre}
        // v={T/TA,u2/u2A,tauA,tau0_Arrhenius}
        return v[2]*exp(pow(v[1],-c[0]/2.0)+log(v[2]/v[3])*(pow(v[0],-1)-1)-1);
    }
    else if (model=="hallwolynes") {
        // tau = tau0*exp(u20/u2)
        // param: {c[0]=tau0,c[1]=u20}
        // v={u2,T,tau}
        return c[0]*exp(c[1]/v[0]);
    }
    else if (model=="hallwolynes1param") {
        // tau = tauA*exp((u20/u2A)*(u2A/u2-1))
        // param: {c[0]=u20}
        // v={u2,u2A,tauA,T,tau}
        return v[2]*exp((c[0]/v[1])*((v[1]/v[0])-1.0));
    }
    else if (model=="leporiniref"||model=="leporini_universal") {
        // tau = 10^(A+B*Xr*Y+C*Xr2*Y2)
        // param: {c[0]=A, c[1]=B; c[2]=C}
        // v={logtauref,u2ref,u2,T,tau}
        double Y=v[1]*pow(v[2],-1);
        double Xr=(-c[1]+sqrt(pow(c[1],2)-4*c[0]*c[2]+4*c[2]*v[0]))*pow(2*c[2],-1);
        return pow(10,c[0]+c[1]*Xr*Y+c[2]*pow(Xr,2)*pow(Y,2));
    }
    else
    {
        cout
        << "in FitData::calc_time_given_vecinput():\n"
        << "model provided ("<<model<<") not found.\n";
        return 0; exit(EXIT_FAILURE);
    }
}





vector<double> FitData::calc_x_given_y(const real_1d_array& c,
                                       const double time,
                                       const string& model)
{
    double lnTime=log(time/c[0]);//F=ln(time/tau0)
    
    vector<double> glass_data;
    
    if (model=="MCT")
    {
        double Tg =
        c[1]*pow(time/c[0],-1/c[2])*(1.0+pow(time/c[0],1/c[2]));
        glass_data.push_back(Tg);
        
        double m =
        (c[2]*Tg)/(log(10)*(Tg-c[1]));
        glass_data.push_back(m);
    }
    else if (model=="SOU")
    {
        double Tg =
        (c[1]*pow(time/c[0],1/c[2]))/(pow(time/c[0],1/c[2])-1.0);
        glass_data.push_back(Tg);
        
        double m =
        (c[2]*c[1])/(log(10)*(Tg-c[1]));
        glass_data.push_back(m);
    }
    else if (model=="GLMVFT1"||model=="GLMVFT3"||model=="GLMVFT4")
    {
        /*===========================================
         F=ln(time/tau0)
         Tg = T0*(1+((F*2^(-a/2))/A)^(-2/a))
         c[0]=tau0, c[1]=T0, c[2]=a, c[3]=A
         ===========================================*/
        double Tg=c[1]*(1+pow((lnTime*pow(2,-c[2]/2))/c[3],-2/c[2]));
        glass_data.push_back(Tg);
        
        /* ==========================================
         m = (a*2^((a/2)-1)*((T0/(Tg-T0))^(a/2)*(1+(T0/(Tg-T0))))*log10(e)
         c[0]=tau0, c[1]=T0, c[2]=a, c[3]=A
         ===========================================*/
        double m=(c[3]*c[2]*pow(2,(c[2]/2)-1)*pow(c[1]*pow(Tg-c[1],-1),c[2]/2)*
                  (1+c[1]*pow(Tg-c[1],-1)))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="VFT")
    {
        /*===========================================
         F=ln(time/tau0)
         Tg = (T0/F)*(D+F) and F!=0 and D*T0!=0
         c[0]=tau0, c[1]=D, c[2]=T0
         ===========================================*/
        double Tg=(c[2]/lnTime)*(c[1]+lnTime);
        glass_data.push_back(Tg);
        
        /* ==========================================
         m = ((D*T0*Tg)/(Tg-T0)^2)*log10(e)
         c[0]=tau0, c[1]=D, c[2]=T0
         ===========================================*/
        double m=((c[1]*c[2]*Tg)/pow(Tg-c[2],2))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="Mauro")
    {
        /*===========================================
         ln(time/tau0)=(K/T)exp(C/T)
         c[0]=tau0, c[1]=K, c[2]=C
         ------------------------------------------
         use Newton's method for root finding
         f  = (K/T)exp(C/T)-ln(time/tau0)
         df = (-K/T^2)exp(C/T)-(KC/T^3)exp(C/T)
         ===========================================*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double Tg_pre=NewtonGuessTemp/corrFacforT; // initial guess of Tg
        double Tg_now=NewtonGuessTemp/corrFacforT;
        double abserr=0;
        int counter=0;
        do {
            double f =
            (c[1]/Tg_pre)*exp(c[2]/Tg_pre)-lnTime;
            
            double df =
            -(c[1]/pow(Tg_pre,2))*exp(c[2]/Tg_pre)
            -(c[1]*c[2]/pow(Tg_pre,3))*exp(c[2]/Tg_pre);
            
            Tg_now=Tg_pre-(f/df);
            abserr=fabs(Tg_now-Tg_pre);
            Tg_pre=Tg_now;
            ++counter;
        } while( ((abserr>TOL)||(Tg_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        double Tg=Tg_now;
        glass_data.push_back(Tg);
        
        /*===========================================
         m = (1+(C/Tg))*(K/Tg)exp(C/Tg)*log10(e)
         c[0]=tau0, c[1]=K, c[2]=C
         ===========================================*/
        double m =
        ((1+(c[2]/Tg))*(c[1]/Tg)*exp(c[2]/Tg))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="AM")
    {
        /*===========================================
         Tg = p/[ln(time/tau0)]^(1/alpha)
         c[0]=tau0, c[1]=p, c[2]=alpha
         ===========================================*/
        double Tg=c[1]/(pow(lnTime,1/c[2]));
        glass_data.push_back(Tg);
        
        /*===========================================
         m = (alpha*(p/Tg)^alpha)*log10(e)
         c[0]=tau0, c[1]=p, c[2]=alpha
         ===========================================*/
        double m=(c[2]*pow(c[1]/Tg,c[2]))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="DG")
    {
        /*===========================================
         ln(time/tau0)=(A/T)exp(-B*T)
         c[0]=tau0, c[1]=A, c[2]=B
         ------------------------------------------
         use Newton's method for root finding
         (verified by Wolfram-Alpha)
         f  = (A/T)exp(-B*T)-ln(time/tau0)
         df = -(A*e^(-B*T)*(B*(T+1)))/T^2
         ===========================================*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double Tg_pre=NewtonGuessTemp/corrFacforT; // initial guess of Tg
        double Tg_now=NewtonGuessTemp/corrFacforT;
        double abserr=0;
        int counter=0;
        do {
            double f =
            ((c[1]/Tg_pre)*exp(-c[2]*Tg_pre))-lnTime;
            
            double df =
            -(c[1]/pow(Tg_pre,2))*(c[2]*Tg_pre+1.0)*exp(-c[2]*Tg_pre);
            
            Tg_now=Tg_pre-(f/df);
            abserr=fabs(Tg_now-Tg_pre);
            Tg_pre=Tg_now;
            ++counter;
        } while( ((abserr>TOL)||(Tg_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        double Tg=Tg_now;
        glass_data.push_back(Tg);
        
        /* Analytic form for fragility */
        // (verified by Wolfram-Alpha)
        double m =
        ((c[1]/Tg)*(c[2]*Tg+1.0)*exp(-c[2]*Tg))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if(model=="ArrheniusII")
    {
        /*===========================================
         Tg = (A+square(A^2+4*B*ln(time/tau0)))/2*ln(time/tau0)
         c[0]=tau0, c[1]=A, c[2]=B
         ===========================================*/
        double squareRoot=sqrt(pow(c[1],2)+4*c[2]*lnTime);
        
        double Tg=(c[1]+squareRoot)/(2*lnTime);
        glass_data.push_back(Tg);
        
        /* ==========================================
         m = ( A/Tg + 2B/Tg^2 )*log10(e)
         c[0]=tau0, c[1]=A, c[2]=B
         ===========================================*/
        double m=((c[1]/Tg)+(2*c[2]/pow(Tg,2)))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if ((model=="COOP")||(model=="COOP_string")||(model=="COOP_DWF"))
    {
        /*===========================================
         ln(time/tau0)=E(T)/T
         c[0]=tau0, c[1]=E_inf, c[2]=u, c[3]=b
         ------------------------------------------
         use Newton's method for root finding
         (verified by Wolfram-Alpha)
         f  = (E(T)/T)-ln(time/tau0)
         df = (E'(T)*T-E(T))/T^2
         ===========================================*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double Tg_pre=NewtonGuessTemp/corrFacforT; // initial guess of Tg
        double Tg_now=NewtonGuessTemp/corrFacforT;
        double abserr=0;
        
        double E_inf=c[1],u=c[2],b=c[3];
        double exponent=0;
        double E_coop=0;
        double E=0;
        
        int counter=0;
        do {
            exponent=exp(-u*((Tg_pre/E_inf)-b));
            E_coop=E_inf*exponent;
            E=E_inf+E_coop;
            
            double f  = (E/Tg_pre)-lnTime;
            double dE = -u*exponent;
            double df = (dE*Tg_pre-E)/pow(Tg_pre,2);
            
            Tg_now=Tg_pre-(f/df);
            abserr=fabs(Tg_now-Tg_pre);
            Tg_pre=Tg_now;
            ++counter;
        } while( ((abserr>TOL)||(Tg_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        double Tg=Tg_now;
        glass_data.push_back(Tg);
        
        exponent=exp(-u*((Tg/E_inf)-b));
        E_coop=E_inf*exponent;
        E=E_inf+E_coop;
        
        /* Analytic form for fragility */
        // (verified by Wolfram-Alpha)
        double m=(((E_inf*(1+exponent))/Tg)+(u*exponent))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="DEAG")
    {
        /*===========================================
         ln(time/tau0)=((A-B*T)/T)*exp(C/T)
         c[0]=tau0, c[1]=A, c[2]=B, c[3]=C
         ------------------------------------------
         use Newton's method for root finding
         (verified by Wolfram-Alpha)
         f  = ((A-BT)/T)*exp(C/T)-ln(time/tau0)
         df = (e^(C/x) (B C x-A (C+x)))/x^3
         ===========================================*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double Tg_pre=NewtonGuessTemp/corrFacforT; // initial guess of Tg
        double Tg_now=NewtonGuessTemp/corrFacforT;
        double abserr=0;
        int counter=0;
        do {
            double f =
            (((c[1]-c[2]*Tg_pre)/Tg_pre)*exp(c[3]/Tg_pre))-lnTime;
            
            double df =
            (exp(c[3]/Tg_pre)*(c[2]*c[3]*Tg_pre-c[1]*(c[3]+Tg_pre)))/pow(Tg_pre,3);
            
            Tg_now=Tg_pre-(f/df);
            abserr=fabs(Tg_now-Tg_pre);
            Tg_pre=Tg_now;
            ++counter;
        } while( ((abserr>TOL)||(Tg_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        double Tg=Tg_now;
        glass_data.push_back(Tg);
        
        /* Analytic form for fragility */
        // (verified by Wolfram-Alpha)
        double m =
        ((exp(c[3]/Tg)*(c[1]*(c[3]+Tg)-c[2]*c[3]*Tg))/pow(Tg,2))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if (model=="CG")
    {
        /*===========================================
         ln(time/tau0)=B/{(T-T0)+[(T-T0)^2+C*T]^0.5}
         c[0]=tau0, c[1]=B, c[2]=C, c[3]=T0
         ------------------------------------------
         Analytic form for Tg:
         (verified by Wolfram-Alpha)
         F=ln(time/tau0)
         Tg = (B*(B+2*F*T0))/(F*(2*B+C*F))
         ===========================================*/
        double Tg=(c[1]/lnTime)*(c[1]+2*lnTime*c[3])/(2*c[1]+c[2]*lnTime);
        glass_data.push_back(Tg);
        
        /* Analytic form for fragility */
        // (verified by Wolfram-Alpha)
        // m = (-(B*((-C*Tg-2*Tg*(Tg-T0))/(2*sqrt(C*Tg+(Tg-T0)^2))-Tg))/(sqrt(C*Tg+(Tg-T0)^2)+Tg-T0)^2)*log10(e)
        double
        delTg=Tg-c[3]; // Tg-T0
        double
        divSqrtTg = pow(pow(delTg,2)+c[2]*Tg,0.5);
        double
        m =
        (-c[1]*(((-c[2]*Tg-2*Tg*delTg)/(2*divSqrtTg))-Tg)/pow(divSqrtTg+delTg,2))
        *log10(exp(1));
        glass_data.push_back(m);
    }
    else if(model=="FourParamVFT")
    {
        /*===========================================
         Tg = T0*lnTime^(-1/alpha)*(D+lnTime^(1/alpha))
         c[0]=tau0, c[1]=D, c[2]=T0, c[3]=alpha
         ===========================================*/
        double Tg =
        c[2]*pow(lnTime,-1/c[3])*(c[1]+pow(lnTime,1/c[3]));
        glass_data.push_back(Tg);
        
        /*===========================================
         m = ((alpha*Tg*((D*T0)/(Tg-T0))^alpha)/(Tg-T0))*log10(e)
         c[0]=tau0, c[1]=D, c[2]=T0, c[3]=alpha
         ===========================================*/
        double m =
        ((c[3]*Tg*pow((c[1]*c[2])/(Tg-c[2]),c[3]))/(Tg-c[2]))*log10(exp(1));
        glass_data.push_back(m);
    }
    else if(model=="ArrheniusIII")
    {
        /*===========================================
         ln(time/tau0)=(A/T)+(B/T^2)+(C/T^3)
         c[0]=tau0, c[1]=A, c[2]=B, c[3]=C
         ------------------------------------------
         use Newton's method for root finding
         f  = ln(time/tau0)-(A/T)-(B/T^2)-(C/T^3)
         df = (A/T^2)+(2*B/T^3)+(3*C/T^4)
         ===========================================*/
        double TOL=1.0e-6; // tolerance for Newton's algorithm
        double Tg_pre=NewtonGuessTemp/corrFacforT; // initial guess of Tg
        double Tg_now=NewtonGuessTemp/corrFacforT;
        double abserr=0;
        int counter=0;
        do {
            double
            f =
            lnTime-
            (c[1]/pow(Tg_pre,1))-
            (c[2]/pow(Tg_pre,2))-
            (c[3]/pow(Tg_pre,3));
            
            double
            df=
            (1*c[1]/pow(Tg_pre,2))+
            (2*c[2]/pow(Tg_pre,3))+
            (3*c[3]/pow(Tg_pre,4));
            
            Tg_now=Tg_pre-(f/df);
            abserr=fabs(Tg_now-Tg_pre);
            Tg_pre=Tg_now;
            ++counter;
        } while( ((abserr>TOL)||(Tg_now<0)) && (counter<1e+4) );
        //cout << counter << "\n";
        //cout << abserr << "\n";
        double Tg=Tg_now;
        glass_data.push_back(Tg);
        
        /* ==========================================
         m = ( A/Tg + 2B/Tg^2 + 3C/Tg^3 )*log10(e)
         
         c[0]=tau0, c[1]=A, c[2]=B, c[3]=C
         ===========================================*/
        double m=
        ((c[1]/Tg)+(2*c[2]/pow(Tg,2))+(3*c[3]/pow(Tg,3)))*log10(exp(1));
        glass_data.push_back(m);
    }
    else
    {
        cout
        << "in FitData::calc_x_given_y():\n"
        << "model provided ("<<model<<") not found.\n";
        return glass_data; exit(EXIT_FAILURE);
    }
    return glass_data;
    //==========================================================================
    // 1st element: Tg (@extrapolated time) || T (@time=t)
    // 2nd element: m  (@Tg)
    //==========================================================================
}





vector<double> FitData::calc_Y_given_time(const real_1d_array& c,
                                          const double time,
                                          const vector<double>& ref,
                                          const string& model)
{
    double logtime=log10(time);
    double logtaur=ref.at(0);
    double u2r=ref.at(1);
    //double Tr=ref.at(2);
    double Yt=0,mt=0;
    
    if (model=="GLM")
    {
        Yt=pow(log(time/c[0]),2.0/c[2])*pow(c[1]/u2r,-1);
        mt=(c[2]*pow(2*Yt,-1)*pow((c[1]/u2r)*Yt,c[2]/2.0))*log10(exp(1));
    }
    else if (model=="GLM1"||model=="GLM1_hallwolynes")
    {
        Yt=pow(log(time/tauA)+1.0,2.0/c[0])*pow(u2A/u2r,-1);
        mt=(c[0]*pow(2*Yt,-1)*pow((u2A/u2r)*Yt,c[0]/2.0))*log10(exp(1));
    }
    else if (model=="univGLM")
    {
        //
    }
    else if (model=="univGT")
    {
        //
    }
    else if (model=="hallwolynes")
    {
        //
    }
    else if (model=="hallwolynes1param")
    {
        //
    }
    else if (model=="leporiniref" || model=="leporini_universal")
    {
        // tau = 10^(A+B*Xr*Y+C*Xr^2+Y^2)
        // param: {c[0]=A, c[1]=B; c[2]=C}
        double Xr  = (-c[1]+sqrt(pow(c[1],2)-4*c[0]*c[2]+4*c[2]*logtaur))*pow(2*c[2],-1);
        double A   = c[0];
        double Br  = c[1]*Xr;
        double Cr  = c[2]*pow(Xr,2);
        Yt = (-Br+sqrt(pow(Br,2)-4*Cr*(A-logtime)))*pow(2*Cr,-1);
        mt = Br+2*Cr*Yt;
    }
    else
    {
        cout
        << "in FitData::calc_Y_given_time():\n"
        << "model provided ("<<model<<") not found.\n";
        return {Yt,mt}; exit(EXIT_FAILURE);
    }
    return {Yt,mt};
    //==========================================================================
    // 1st element: Yt (@time=t)
    // 2nd element: mt (@Yt)
    //==========================================================================
}





double FitData::error_Tg(const string& model)
{
    if (model=="VFT")
    {
        //
        // F = Tg = T0+(D*T0/ln(time/tau0))
        //
        // elements: {tau0, D, T0}
        //
        
        vector<double> errori={rep.errpar[0],rep.errpar[1],rep.errpar[2]};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0;
        double Tg_error=0;
        double tau0=c[0],D=c[1],T0=c[2];
        double time=extrp_time; // time @ Tg
        double lnTime=log(time/tau0);
        
        // (verified by Wolfram-Alpha)
        dF0 = (D*T0)/(tau0*pow(lnTime,2));
        dF1 = T0/lnTime;
        dF2 = 1+(D/lnTime);
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        
        Tg_error=error_propagation(dFi,errori);
        
        return Tg_error;
    }
    else if (model=="AM")
    {
        //
        // F = Tg = p/lnTime^(1/alpha)
        //
        // elements: {tau0, p, alpha}
        //
        
        vector<double> errori={rep.errpar[0],rep.errpar[1],rep.errpar[2]};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0;
        double Tg_error=0;
        double time=extrp_time; // time @ Tg
        double lnTime=log(time/c[0]);
        
        // (verified by Wolfram-Alpha)
        dF0 = (c[1]/(c[2]*c[0]))*pow(lnTime,-(c[2]+1)/c[2]);
        dF1 = 1/pow(lnTime,1/c[2]);
        dF2 = (c[1]/pow(c[2],2))*pow(lnTime,-1/c[2])*log(lnTime);
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        
        Tg_error=error_propagation(dFi,errori);
        
        return Tg_error;
    }
    else if (model=="CG")
    {
        //
        // F = Tg = Tg = (B*(B+2*lnTime*T0))/(lnTime*(2*B+C*lnTime))
        //
        // elements: {tau0, B, C, T0}
        //
        
        vector<double> errori={rep.errpar[0],rep.errpar[1],rep.errpar[2],rep.errpar[3]};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0,dF3=0;
        double Tg_error=0;
        double time=extrp_time; // time @ Tg
        double lnTime=log(time/c[0]);
        
        // (verified by Wolfram-Alpha)
        dF0 = ((2*c[1])/(c[0]*pow(lnTime,2)))*(pow(c[1],2)+c[1]*c[2]*lnTime+c[2]*c[3]*lnTime)/pow(2*c[1]+c[2]*lnTime,2);
        dF1 = (2/lnTime)*(c[2]*lnTime*(lnTime*c[3]+c[1])+pow(c[1],2))/pow(c[2]*lnTime+2*c[1],2);
        dF2 = -c[1]*(c[1]+2*lnTime*c[3])/pow(2*c[1]+c[2]*lnTime,2);
        dF3 = (2*c[1])/(2*c[1]+c[2]*lnTime);
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        dFi.push_back(dF3);
        
        Tg_error=error_propagation(dFi,errori);
        
        return Tg_error;
    }
    else if (model=="FourParamVFT")
    {
        //
        // F = Tg = T0*lnTime^(-1/alpha)*(D+lnTime^(1/alpha))
        //
        // elements: {tau0, D, T0, alpha}
        //
        
        vector<double> errori={rep.errpar[0],rep.errpar[1],rep.errpar[2],rep.errpar[3]};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0,dF3=0;
        double Tg_error=0;
        double time=extrp_time; // time @ Tg
        double lnTime=log(time/c[0]);
        
        // (verified by Wolfram-Alpha)
        dF0 = ((c[1]*c[2])/(c[0]*c[3]))*pow(lnTime,-((c[3]+1)/c[3]));
        dF1 = c[2]*pow(lnTime,-1/c[3]);
        dF2 = pow(lnTime,-1/c[3])*(c[1]+pow(lnTime,1/c[3]));
        dF3 = ((c[1]*c[2])/pow(c[3],2))*pow(lnTime,-1/c[3])*log(lnTime);
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        dFi.push_back(dF3);
        
        Tg_error=error_propagation(dFi,errori);
        
        return Tg_error;
    }
    else
    {
        return 0;
    }
    
}





double FitData::error_fragility(const string& model,
                                const double Tg_value,
                                const double Tg_error)
{
    if (model=="VFT")
    {
        //
        // F = m = ((D*T0*Tg)/(Tg-T0)^2)*log10(e)
        //
        // elements: {D, T0, Tg}
        //
        
        vector<double> errori={rep.errpar[1],rep.errpar[2],Tg_error};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0;
        double m_error=0;
        double D=c[1],T0=c[2],Tg=Tg_value;
        
        // (verified by Wolfram-Alpha)
        dF0 = ((T0*Tg)/pow(Tg-T0,2))*log10(exp(1));
        dF1 = ((D*Tg*(T0+Tg))/pow(Tg-T0,3))*log10(exp(1));
        dF2 = ((D*T0*(T0+Tg))/pow(T0-Tg,3))*log10(exp(1));
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        
        m_error=error_propagation(dFi,errori);
        
        return m_error;
    }
    else if (model=="AM")
    {
        //
        // F = m = (alpha*(p/Tg)^alpha)*log10(e)
        //
        // elements: {p, alpha, Tg}
        //
        
        vector<double> errori={rep.errpar[1],rep.errpar[2],Tg_error};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0;
        double m_error=0;
        
        // (verified by Wolfram-Alpha)
        dF0 = ((pow(c[2],2)/c[1])*pow(c[1]/Tg_value,c[2]))*log10(exp(1));
        dF1 = ((pow(c[1]/Tg_value,c[2]))*(1+c[2]*log(c[1]/Tg_value)))*log10(exp(1));
        dF2 = (-(pow(c[2],2)/Tg_value)*pow(c[1]/Tg_value,c[2]))*log10(exp(1));
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        
        m_error=error_propagation(dFi,errori);
        
        return m_error;
    }
    else if (model=="CG")
    {
        //
        // F = m = (-(B*((-C*Tg-2*Tg*(Tg-T0))/(2*sqrt(C*Tg+(Tg-T0)^2))-Tg))/(sqrt(C*Tg+(Tg-T0)^2)+Tg-T0)^2)*log10(e)
        //
        // elements: {B, C, T0, Tg}
        //
        
        vector<double> errori={rep.errpar[1],rep.errpar[2],rep.errpar[3],Tg_error};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0,dF3=0;
        double m_error=0;
        
        double B=c[1],C=c[2],T0=c[3],T=Tg_value;
        
        // (verified by Wolfram-Alpha)
        dF0 = (-((-T+(-(C*T)-2*T*(T-T0))/(2*sqrt(C*T+pow((T-T0),2))))/pow((T+sqrt(C*T+pow((T-T0),2))-T0),2)))*log10(exp(1));
        dF1 = ((B*T*(-T+(-2*T*(T-T0)-T*C)/(2*sqrt(pow((T-T0),2)+T*C))))/(sqrt(pow((T-T0),2)+T*C)*pow((T-T0+sqrt(pow((T-T0),2)+T*C)),3))-(B*(-(T*(-2*T*(T-T0)-T*C))/(4*pow((pow((T-T0),2)+T*C),1.5))-T/(2*sqrt(pow((T-T0),2)+T*C))))/pow((T-T0+sqrt(pow((T-T0),2)+T*C)),2))*log10(exp(1));
        dF2 = ((2*B*(-T+(-(C*T)-2*T*(T-T0))/(2*sqrt(C*T+pow((T-T0),2))))*(-1-(T-T0)/sqrt(C*T+pow((T-T0),2))))/pow((T+sqrt(C*T+pow((T-T0),2))-T0),3)-(B*(T/sqrt(C*T+pow((T-T0),2))+((-(C*T)-2*T*(T-T0))*(T-T0))/(2*pow((C*T+pow((T-T0),2)),1.5))))/pow((T+sqrt(C*T+pow((T-T0),2))-T0),2))*log10(exp(1));
        dF3 = ((2*B*(1+(C+2*(-T0+T))/(2*sqrt(C*T+pow((-T0+T),2))))*(-T+(-(C*T)-2*T*(-T0+T))/(2*sqrt(C*T+pow((-T0+T),2)))))/pow((-T0+T+sqrt(C*T+pow((-T0+T),2))),3)-(B*(-1-((C+2*(-T0+T))*(-(C*T) -2*T*(-T0+T)))/(4*pow((C*T+pow((-T0+T),2)),1.5))+(-C-2*T-2*(-T0+T))/(2*sqrt(C*T+pow((-T0+T),2)))))/pow((-T0+T+sqrt(C*T+pow((-T0+T),2))),2))*log10(exp(1));
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        dFi.push_back(dF3);
        
        m_error=error_propagation(dFi,errori);
        
        return m_error;
    }
    else if (model=="FourParamVFT")
    {
        //
        // F = m = ((alpha*Tg*((D*T0)/(Tg-T0))^alpha)/(Tg-T0))*log10(e)
        //
        // elements: {D, T0, alpha, Tg}
        //
        
        vector<double> errori={rep.errpar[1],rep.errpar[2],rep.errpar[3],Tg_error};
        vector<double> dFi;
        
        double dF0=0,dF1=0,dF2=0,dF3=0;
        double m_error=0;
        
        // (verified by Wolfram-Alpha)
        dF0 = (((Tg_value*pow(c[3],2))/(c[2]*pow(c[1],2)))*pow((c[1]*c[2])/(Tg_value-c[2]),c[3]+1))*log10(exp(1));
        dF1 = ((c[3]*Tg_value*(c[3]*Tg_value+c[2]))*pow((c[1]*c[2])/(Tg_value-c[2]),c[3])/(c[2]*pow(Tg_value-c[2],2)))*log10(exp(1));
        dF2 = ((Tg_value/(Tg_value-c[2]))*(1+c[3]*log((c[1]*c[2])/(Tg_value-c[2])))*pow((c[1]*c[2])/(Tg_value-c[2]),c[3]))*log10(exp(1));
        dF3 = (-c[3]*(c[2]+c[3]*Tg_value)*pow((c[1]*c[2])/(Tg_value-c[2]),c[3])/pow(Tg_value-c[2],2))*log10(exp(1));
        
        dFi.push_back(dF0);
        dFi.push_back(dF1);
        dFi.push_back(dF2);
        dFi.push_back(dF3);
        
        m_error=error_propagation(dFi,errori);
        
        return m_error;
    }
    else
    {
        return 0;
    }
}





vector<double> FitData::get_interpolated_Ts(StructureClass& sysVar,
                                            const int n_trl,
                                            const int n_sys,
                                            const int n_temps,
                                            const double Tl,
                                            const double Th)
{
    vector<double> originalTs=sysVar.get_temperaturesInfo().at(index);
    vector<double> simtemps; // simulated temperatures
    vector<double> newtemperatures;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/tempsInfo_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    ifstream readfile(i.c_str());
    if (readfile.is_open()) {
        string lineContent;
        double Tdata=0;
        while (getline(readfile,lineContent)) {
            istringstream iss(lineContent);
            iss >> Tdata;
            simtemps.push_back(Tdata);
        }
    } else {
        cout
        << "in FitData::get_interpolated_Ts():\n"
        << "file ("+i+") cannot open.\n";
    } readfile.close();
    
    double dT=fabs(Tl-Th);
    double Ti=min(Tl,Th);
    
    for (int i=0; i<n_temps; ++i)
    {
        ////////////////////////////////////////////////////////////////////////
        /** Check repeatedness of temperatures **/
        //----------------------------------------------------------------------
        bool is_foundrepeat=false;
        do {
            vector<double>::iterator itr0,itr1,itr2;
            itr0=find(originalTs.begin(),originalTs.end(),Ti);
            itr1=find(simtemps.begin(),simtemps.end(),Ti);
            itr2=find(newtemperatures.begin(),newtemperatures.end(),Ti);
            //------------------------------------------------------------------
            // "find" returns iterators in the range: [first,last); it would
            // return a pointer to the "last" element if nothing is found.
            //------------------------------------------------------------------
            if ((itr0!=originalTs.end())||
                (itr1!=simtemps.end())||
                (itr2!=newtemperatures.end()))
            {
                Ti -= precision;
                is_foundrepeat=true;
            } else {
                is_foundrepeat=false;
            }
        } while (is_foundrepeat);
        //----------------------------------------------------------------------
        newtemperatures.push_back(Ti); // sotre Ti if it's unique (not repeated)
        ////////////////////////////////////////////////////////////////////////
        
        /** Temperature spacing **/
        if (get_shootratio()>0) {
            /** Ti's are spaced between Th and Tl based on shootratio **/
            dT *= get_shootratio();
            // NOTE:
            // this algorithm makes extrp. T's more crowded
            // in higher temperature end that the temperature spacing (dT)
            // is progressively shrinking in each iteration.
        } else {
            /** Ti's are evenly spaced between Th and Tl **/
            dT = fabs(Tl-Th)/(double)(n_temps);
            // NOTE: n_temps >= 2
        }
        /** New T's are increased from low to high **/
        Ti = std::round(Ti+dT);
    }
    //--------------------------------------------------------------------------
    // Now reverse the order of new temperatures;
    // by convention, T's in temperature containers are from high to low
    //--------------------------------------------------------------------------
    vector<double> hightolow;
    for (int i=(int)(newtemperatures.size()-1); i>=0; --i) {
        hightolow.push_back(newtemperatures.at(i));
    } return hightolow;
}





void FitData::read_raw_msd_data(const StructureClass& sysVar,
                                const int n_trl,
                                const int n_sys,
                                const double Temp_d,
                                const int frame)
{
    msd_sortincreasing.clear();       //(time,msd)
    loglogmsd_sortincreasing.clear(); //(log.time,log.msd)
    semilogmsd_sortincreasing.clear();//(log.time,msd)
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        i.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        i.append("/statistics");
    }
    i.append("/msd/msd_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    if (frame!=0) i.append("_frame"+to_string((long long int)frame));
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        string strData;
        double time=0,msd=0;
        /* read 1st line */
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in FitData::read_raw_msd_data():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        /* read 2nd to final line */
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> time; // time
            iss >> msd;  // msd (length^2)
            if (time==0) continue;//not taking time=0
            /** msd data already sorted by time (increasing) in amdat **/
            msd_sortincreasing.push_back({time,msd});
            loglogmsd_sortincreasing.push_back({log10(time),log10(msd)});
            semilogmsd_sortincreasing.push_back({log10(time),msd});
        } readFile.close();
    } else {
        cout
        << "in FitData::read_raw_msd_data():\n"
        << "file ("<<i<<") cannot open.\n";
    }
}





void FitData::write_raw_msd_data(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d)
{
    read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/MSD_raw_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream writeFile(o.c_str());
    
    writeFile << "t  <r2(t)> \n";
    for (int i=0;i<(int)msd_sortincreasing.size();++i) {
        writeFile
        << msd_sortincreasing.at(i).at(0) << " "
        << msd_sortincreasing.at(i).at(1) << "\n";
    } writeFile.close();
}





double FitData::get_DiffCoeff(const StructureClass& sysVar,
                              const int n_trl,
                              const int n_sys,
                              const double Temp_d)
{
    read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
    size_t n_data=msd_sortincreasing.size();
    double time=msd_sortincreasing.at(n_data-1).at(0);
    double r2fn=msd_sortincreasing.at(n_data-1).at(1);
    double dc=r2fn/(6*time);
    return dc;
}





void FitData::write_diffusion_coeff(const StructureClass& sysVar,
                                    const int n_trl,
                                    const int n_sys,
                                    const double Temp_d)
{
    read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/MSD_DiffCoeff_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream writeFile(o.c_str());
    
    double f=0,df=0,d2f=0;
    double dcf=0,dcdf=0,dcd2f=0;
    double r2f=0,r2df=0,r2d2f=0;
    vector<vector<double>> r2dc2d,r22dloglog,dc2dloglog;
    
    for (int i=0;i<(int)msd_sortincreasing.size();++i)
    {
        double time=msd_sortincreasing.at(i).at(0);
        double r2  =msd_sortincreasing.at(i).at(1);
        if (time>0) {
            double dc=r2/(6*time);//diffusion coefficient
            r2dc2d.push_back({time,r2,dc});
            r22dloglog.push_back({log10(time),log10(r2)});
            dc2dloglog.push_back({log10(time),log10(dc)});
        }
    }
    build_penalizedspline_interpolant(r22dloglog,0,1);
    auto intrp_msdt=interpolant;
    build_penalizedspline_interpolant(dc2dloglog,0,1);
    auto intrp_dc=interpolant;
    
    size_t n_data=r2dc2d.size();
    writeFile
    << "D     "<< r2dc2d.at(n_data-1).at(2)        << "\n"
    << "Log.D "<< log10(r2dc2d.at(n_data-1).at(2)) << "\n\n";;
    
    writeFile << "t r2 D L.t L.r2 L.r2(sp) dLr2 d2r2 L.D L.D(sp) dLD d2LD \n";
    for (int i=0;i<(int)r2dc2d.size();++i)
    {
        double logtime=log10(r2dc2d.at(i).at(0));
        spline1ddiff(intrp_msdt,logtime,f,df,d2f);
        r2f=f; r2df=df; r2d2f=d2f;
        spline1ddiff(intrp_dc,logtime,f,df,d2f);
        dcf=f; dcdf=df; dcd2f=d2f;
        
        writeFile
        << r2dc2d.at(i).at(0) << " "        // time
        << r2dc2d.at(i).at(1) << " "        // r2
        << r2dc2d.at(i).at(2) << " "        // dc
        << logtime << " "                   // log(time)
        << log10(r2dc2d.at(i).at(1)) << " " // log(r2)
        << r2f   << " "                     // log(r2)_pspline
        << r2df  << " "                     // D.log(r2)
        << r2d2f << " "                     // D2.log(r2)
        << log10(r2dc2d.at(i).at(2)) << " " // log(dc)
        << dcf   << " "                     // log(dc)_pspline
        << dcdf  << " "                     // D.log(dc)
        << dcd2f << " "                     // D2.log(dc)
        << "\n";
    } writeFile.close();
}





void FitData::write_loglogMSD_derivatives(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d)
{
    read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/msd_loglogdf_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_dmsd(o.c_str(),ofstream::app);//'append'
    
    /** build interpolant DS **/
    //--------------------------------------------------------------------------
    double logtime=0,f=0,df=0,d2f=0;
    build_cubicspline_interpolant(loglogmsd_sortincreasing,0,1);
    spline1dinterpolant loglogMSD_interpolant=interpolant;
    vector<vector<double>> logtimedlogMSD;
    for (indexi=0; indexi<(int)loglogmsd_sortincreasing.size(); ++indexi) {
        logtime=loglogmsd_sortincreasing.at(indexi).at(0);
        spline1ddiff(loglogMSD_interpolant,logtime,f,df,d2f);
        logtimedlogMSD.push_back({logtime,df});
    }
    /** noise reduction for dLog(MSD)/dt **/
    build_penalizedspline_interpolant(logtimedlogMSD,0,1);
    spline1dinterpolant logtimedlogMSD_interpolant=interpolant;
    
    write_dmsd << "Log(t) Log(MSD) dLog(MSD) dLog(MSD)_sfit d2Log(MSD)_sfit\n";
    for (indexi=0; indexi<(int)loglogmsd_sortincreasing.size(); ++indexi) {
        logtime=loglogmsd_sortincreasing.at(indexi).at(0);
        spline1ddiff(loglogMSD_interpolant,logtime,f,df,d2f);
        write_dmsd << logtime << " " << f << " " << df << " ";
        spline1ddiff(logtimedlogMSD_interpolant,logtime,f,df,d2f);
        write_dmsd << f << " " << df << "\n";
    } write_dmsd.close();
    
    /** find inflection point **/
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/msd_loglogdf_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //ofstream write_inf(o.c_str(),ofstream::app); // 'append'
}





void FitData::write_DWF(const StructureClass& sysVar,
                        const int n_trl,
                        const int n_sys,
                        const double Temp_d,
                        const int frame)
{
    read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d,frame);
    //find_balltocage_time(sysVar,n_trl,n_sys,Temp_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/DWF_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_DWF(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/DWF_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_DWFinv(o.c_str(),ofstream::app);//'append'
    
    bool flag=false;
    double u2=
    pow(10,tfunction_at_given_time(loglogmsd_sortincreasing,log10(DWF_time),flag));
    if (flag) return;
    
    double T_actual=Temp_d*pow(corrFacforT,-1);//NOTE: convert!
    
    // T DWF
    streamsize ss=write_DWF.precision();
    
    write_DWF
    << fixed << setprecision(10) //NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_DWF.precision(ss);
    write_DWF << resetiosflags(ios::fixed|ios::showpoint);
    
    write_DWF
    << u2
    << "\n";
    write_DWFinv
    << convert/T_actual
    << " "
    << u2
    << "\n";
    write_DWF.close();
    write_DWFinv.close();
}





void FitData::write_NGP(const StructureClass& sysVar,
                        const int n_trl,
                        const int n_sys,
                        const double Temp_d)
{
    write_ngp_peak_frame(sysVar,n_trl,n_sys,Temp_d);
}





void FitData::write_avg_NGP(const StructureClass& sysVar)
{
    write_avg_ngp_peak_frame(sysVar);
}





void FitData::read_individual_DWF(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys)
{
    dwf_sortdecreasing.clear();
    vector<vector<double>> DWF;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/DWF_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,DWFData=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent); // (T,DWF)
            iss >> Tdata;
            iss >> DWFData;
            DWF.push_back({Tdata,DWFData});
        } readFile.close();
    } else {
        cout
        << "in FitData::read_individual_DWF_data():\n"
        << i << " cannot open.\n";
    }
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_individual_DWF_data():" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
}





void FitData::read_individual_equ_DWF(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys)
{
    dwf_sortdecreasing.clear();
    dwf_tau_data.clear();
    vector<vector<double>> DWF;
    vector<vector<double>> dwf_tau;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/DWF_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,DWFData=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent); // (T,DWF)
            iss >> Tdata;
            iss >> DWFData;
            DWF.push_back({Tdata,DWFData});
        } readFile.close();
    } else {
        cout
        << "in FitData::read_individual_equ_DWF():\n"
        << "file ("+i+") cannot open.\n";
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_individual_equ_DWF():\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: individual */
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int dwfSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<dwfSize; ++ii) {
            if (taueq_sortdecreasing.at(i)[0]==dwf_sortdecreasing.at(ii)[0]) {
                taue=taueq_sortdecreasing.at(i)[1];
                temp=dwf_sortdecreasing.at(ii)[0];
                dwfr=dwf_sortdecreasing.at(ii)[1];
                DWF.push_back({temp,dwfr});
                dwf_tau.push_back({dwfr,taue});
                break;
            }
        }
    }
    dwf_sortdecreasing=DWF;
    dwf_tau_data=dwf_tau;
}





void FitData::read_all_DWF(const StructureClass& sysVar,
                           const int n_sys)
{
    dwf_sortdecreasing.clear();
    vector<vector<double>> DWF;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/DWF_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent); // (T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                DWF.push_back({Tdata,DWFData});
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_DWF():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_DWF():\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
}





void FitData::read_all_equ_DWF(const StructureClass& sysVar,
                               const int n_sys)
{
    dwf_sortdecreasing.clear();
    dwf_invT_sortdecreasing.clear();
    dwf_tau_data.clear();
    vector<vector<double>> DWF;
    vector<vector<double>> DWF_invT;
    vector<vector<double>> dwf_tau;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        //i.append("/DWF_all_");
        //i.append("/DWF_equ_avg_");
        i.append("/DWF_equ_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);//(T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                DWF.push_back({Tdata,DWFData});
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_equ_DWF():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_equ_DWF():\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: all */
    read_all_taueq_data(sysVar,n_sys);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int dwfSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<dwfSize; ++ii) {
            if (taueq_sortdecreasing.at(i)[0]==dwf_sortdecreasing.at(ii)[0]) {
                taue=taueq_sortdecreasing.at(i)[1];
                temp=dwf_sortdecreasing.at(ii)[0];
                dwfr=dwf_sortdecreasing.at(ii)[1];
                DWF.push_back({temp,dwfr});
                DWF_invT.push_back({convert/temp,dwfr});//NOTE
                dwf_tau.push_back({dwfr,taue});
                break;
            }
        }
    }
    dwf_sortdecreasing=DWF;
    dwf_invT_sortdecreasing=DWF_invT;
    dwf_tau_data=dwf_tau;
}





void FitData::read_all_equ_DWF(const StructureClass& sysVar,
                               const int n_sys,
                               const double TA_d)
{
    dwf_sortdecreasing.clear();
    dwf_invT_sortdecreasing.clear();
    dwf_tau_data.clear();
    vector<vector<double>> DWF;
    vector<vector<double>> DWF_invT;
    vector<vector<double>> dwf_tau;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        //i.append("/DWF_all_");
        //i.append("/DWF_equ_avg_");
        i.append("/DWF_equ_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent); // (T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    DWF.push_back({Tdata,DWFData});}
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_equ_DWF():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_equ_DWF():\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing0);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: all */
    read_all_taueq_data(sysVar,n_sys,TA_d);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int dwfSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<dwfSize; ++ii) {
            if (taueq_sortdecreasing.at(i)[0]==dwf_sortdecreasing.at(ii)[0]) {
                taue=taueq_sortdecreasing.at(i)[1];
                temp=dwf_sortdecreasing.at(ii)[0];
                dwfr=dwf_sortdecreasing.at(ii)[1];
                DWF.push_back({temp,dwfr});
                DWF_invT.push_back({convert/temp,dwfr});//NOTE
                dwf_tau.push_back({dwfr,taue});
                break;
            }
        }
    }
    dwf_sortdecreasing=DWF;
    dwf_invT_sortdecreasing=DWF_invT;
    dwf_tau_data=dwf_tau;
}





void FitData::find_ngp_peak_frame(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d)
{
    ngp_data_raw.clear();
    ngp_data_raw_smoothed.clear();
    ngp_data_sortincreasing.clear();
    vector<vector<double>> ngp_data;
    vector<vector<double>> ngp_data_logtime;
    vector<vector<double>> ngp_raw_logtime;
    vector<vector<double>> ngp_data_smoothed;
    vector<vector<double>> ngp_raw_smoothed;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        i.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        i.append("/statistics");
    }
    i.append("/ngp/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    ifstream read_ngp(i.c_str());
    if (read_ngp.is_open())
    {
        string lineContent;
        string strData;
        int    frame_index=0;//zero-based index
        double time=0,ngp=0;
        
        getline(read_ngp,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in FitData::find_ngp_peak_frame():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        while (getline(read_ngp,lineContent))
        {
            istringstream iss(lineContent);
            iss >> time;
            iss >>  ngp;
            if (time>0) {
                ngp_data_raw.push_back
                ({time,ngp,(double)frame_index});//raw
                ngp_raw_logtime.push_back
                ({log10(time),ngp,(double)frame_index});//raw in logtime
            }
            /** collect "in-block" data: 0 < time <= block time **/
            if ((time>0)&&(time<=compu_time)) {
                ngp_data.push_back
                ({time,ngp,(double)frame_index});//in-block
                ngp_data_logtime.push_back
                ({log10(time),ngp,(double)frame_index});//in-block in logtime
            } ++frame_index;
        } read_ngp.close();
        
        /** smoothed all data **/
        //----------------------------------------------------------------------
        double logtime=0,f=0;
        build_penalizedspline_interpolant(ngp_raw_logtime,0,1,50,rho_ngp);//use logtime
        for (int i=0;i<(int)ngp_raw_logtime.size();++i) {
            logtime=ngp_raw_logtime.at(i).at(0);//use logtime
            f=spline1dcalc(interpolant,logtime);
            ngp_data_raw_smoothed.push_back
            ({ngp_data_raw.at(i).at(0),f,ngp_data_raw.at(i).at(2)});
        }
        //----------------------------------------------------------------------
        
        /** smoothed in-block data **/
        //----------------------------------------------------------------------
        logtime=0,f=0;
        build_penalizedspline_interpolant(ngp_data_logtime,0,1,50,rho_ngp);//use logtime
        for (int i=0;i<(int)ngp_data_logtime.size();++i) {
            logtime=ngp_data_logtime.at(i).at(0);//use logtime
            f=spline1dcalc(interpolant,logtime);
            ngp_data_smoothed.push_back
            ({ngp_data.at(i).at(0),f,ngp_data.at(i).at(2)});
        }
        //----------------------------------------------------------------------
    } else {
        cout
        << "in FitData::find_ngp_peak_frame():\n"
        << i <<" cannot open.\n";
        exit(EXIT_FAILURE);
    }
    try {
        if (ngp_data.size()==0) throw 0;
        if (is_use_ngp_smoothed) {
            std::sort(ngp_data_smoothed.begin(),ngp_data_smoothed.end(),sortIncreasing0);
            std::sort(ngp_data.begin(),ngp_data.end(),sortIncreasing0);
            ngp_data_sortincreasing=ngp_data_smoothed;
        } else {
            std::sort(ngp_data.begin(),ngp_data.end(),sortIncreasing0);
            ngp_data_sortincreasing=ngp_data;
        }
    } catch (int i) {
        cout
        << "in FitData::find_ngp_peak_frame():" << "\n"
        << "ngp size = " << i << "\n";
        ngp_data.push_back({0,0});
    }
}





void FitData::find_avg_ngp_peak_frame(const StructureClass& sysVar)
{
    ngp_avg_data_sortDecreasing.clear();
    vector<vector<double>> data_tmp;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/NGP_peak_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append("."+analysispart+".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    dubVec.clear();
                    istringstream iss(lineContent);
                    while (iss>>dubVar) {
                        dubVec.push_back(dubVar);
                    } data_tmp.push_back(dubVec);
                }
            } else {
                cout
                << "in FitData::find_avg_ngp_peak_frame():\n"
                << i <<" cannot open.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    std::sort(data_tmp.begin(),data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(data_tmp,ngp_avg_data_sortDecreasing);
}





void FitData::write_ngp_peak_frame(const StructureClass& sysVar,
                                   const int n_trl,
                                   const int n_sys,
                                   const double Temp_d)
{
    find_ngp_peak_frame(sysVar,n_trl,n_sys,Temp_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/NGP_raw_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));//T
    o.append("."+analysispart+".dat");
    ofstream write_rawNGP(o.c_str());
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/NGP_peak_vs_T_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("."+analysispart+".dat");
    ofstream write_peakNGPvsT(o.c_str(),ofstream::app);
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        o.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        o.append("/statistics");
    }
    o.append("/ngp/peakframe_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));//T
    o.append("."+analysispart+".dat");
    ofstream write_rawtostats(o.c_str());
    
    write_rawtostats<<"time log(time) NGP_raw frame\n";
    write_rawNGP    <<"time log(time) NGP_raw frame\n";
    
    for (int i=0;i<(int)ngp_data_raw.size();++i)
    {
        write_rawtostats
        << ngp_data_raw.at(i).at(0)        << " "
        << log10(ngp_data_raw.at(i).at(0)) << " "
        << ngp_data_raw.at(i).at(1)        << " "
        << ngp_data_raw.at(i).at(2)        << "\n";
        
        write_rawNGP
        << ngp_data_raw.at(i).at(0)        << " "
        << log10(ngp_data_raw.at(i).at(0)) << " "
        << ngp_data_raw.at(i).at(1)        << " "
        << ngp_data_raw.at(i).at(2)        << "\n";
    }
    if (is_use_ngp_smoothed)
    {
        write_rawtostats<<"\ntime log(time) NGP_smoothed frame\n";
        write_rawNGP    <<"\ntime log(time) NGP_smoothed frame\n";
        for (int i=0;i<(int)ngp_data_raw_smoothed.size();++i)
        {
            write_rawtostats
            << ngp_data_raw_smoothed.at(i).at(0)        << " "
            << log10(ngp_data_raw_smoothed.at(i).at(0)) << " "
            << ngp_data_raw_smoothed.at(i).at(1)        << " "
            << ngp_data_raw_smoothed.at(i).at(2)        << "\n";
            write_rawNGP
            << ngp_data_raw_smoothed.at(i).at(0)        << " "
            << log10(ngp_data_raw_smoothed.at(i).at(0)) << " "
            << ngp_data_raw_smoothed.at(i).at(1)        << " "
            << ngp_data_raw_smoothed.at(i).at(2)        << "\n";
        }
    }
    
    vector<vector<double>> ngp_raw_sort=ngp_data_raw;
    vector<vector<double>> ngp_data_sort=ngp_data_sortincreasing;
    std::sort(ngp_raw_sort.begin(),ngp_raw_sort.end(),sortDecreasing0);//sort time
    std::sort(ngp_data_sort.begin(),ngp_data_sort.end(),sortDecreasing0);//sort time
    
    write_rawtostats << "\n"
    << "final in-block time = " << ngp_data_sort.at(0).at(0)<< "\n"
    << "final  overall time = " << ngp_raw_sort.at(0).at(0) << "\n";
    write_rawNGP << "\n"
    << "final in-block time = " << ngp_data_sort.at(0).at(0)<< "\n"
    << "final  overall time = " << ngp_raw_sort.at(0).at(0) << "\n";
    
    /** find peak of NGP from in-block data **/
    std::sort(ngp_data_sort.begin(),ngp_data_sort.end(),sortDecreasing1);//sort NGP
    
    double logtaueq=get_individual_taueq_at_T(sysVar,n_trl,n_sys,Temp_d);
    
    write_peakNGPvsT
    << Temp_d*pow(corrFacforT,-1)       << " " //T
    << logtaueq                         << " " //log(tau)
    << ngp_data_sort.at(0).at(0)        << " " //t*
    << log10(ngp_data_sort.at(0).at(0)) << " " //logt*
    << ngp_data_sort.at(0).at(1)        << " " //NGP_smoothed
    << ngp_data_sort.at(0).at(2)        << " " //frame
    << "\n";
    
    write_rawtostats.close();
    write_rawNGP.close();
}





void FitData::write_avg_ngp_peak_frame(const StructureClass& sysVar)
{
    find_avg_ngp_peak_frame(sysVar);
    
    vector<vector<double>> ngp=ngp_avg_data_sortDecreasing;
    for (int i=0;i<(int)ngp.size();++i) {
        ngp.at(i).at(0)=convert/ngp.at(i).at(0);//replace T by 1/T
    }
    alglib::spline1dinterpolant inp0,inp1,inph,inp2;
    /** interpolate log.t* vs T **/
    build_penalizedspline_interpolant(ngp,0,3,50,0.0);
    inp0=interpolant;
    build_penalizedspline_interpolant(ngp,0,3,50,0.5);
    inph=interpolant;
    build_penalizedspline_interpolant(ngp,0,3,50,1.0);
    inp1=interpolant;
    build_penalizedspline_interpolant(ngp,0,3,50,2.0);
    inp2=interpolant;
    
    vector<vector<double>> ngpcoop;
    double T=0,logngp=0;
    for (int i=0;i<(int)ngp_avg_data_sortDecreasing.size();++i) {
        T=ngp_avg_data_sortDecreasing.at(i).at(0);
        logngp=ngp_avg_data_sortDecreasing.at(i).at(3);
        ngpcoop.push_back({T,logngp});
    }
    set_fitParams("COOP");
    fitdata_processing_lsfit(ngpcoop);
    alglib_lsfit("COOP");
    alglib::lsfitreport rep_ngpcoop=rep;
    alglib::real_1d_array c_ngpcoop=c;
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/NGP_avg_peak_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile(o.c_str());
    
    for (int i=0;i<(int)ngp_avg_data_sortDecreasing.size();++i)
    {
        double T=ngpcoop.at(i).at(0);
        double invT=ngp.at(i).at(0);
        for (int ii=0;ii<(int)ngp_avg_data_sortDecreasing.at(i).size();++ii) {
            writeFile << ngp_avg_data_sortDecreasing.at(i).at(ii) << " ";
        }
        writeFile
        << spline1dcalc(inp0,invT) << " "
        << spline1dcalc(inph,invT) << " "
        << spline1dcalc(inp1,invT) << " "
        << spline1dcalc(inp2,invT) << " "
        << log10(calc_y_given_x(c_ngpcoop,T,"COOP")) << " "
        << "\n";
        
    } writeFile.close();
}





void FitData::write_avg_normal_modes_decoupling(const StructureClass& sysVar)
{
    vector<vector<double>> decoupling_data_tmp;
    vector<vector<double>> decoupling_data;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/decoupling_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append("."+analysispart+".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double T=0,invT=0,dp1=0,dp=0,logtisfs=0,logtbaf=0;
                string lineContent;
                while(getline(readFile,lineContent))
                {
                    istringstream iss(lineContent);
                    iss >> T;
                    iss >> invT;
                    iss >> dp1;
                    iss >> dp;
                    iss >> logtisfs;
                    iss >> logtbaf;
                    decoupling_data_tmp.push_back
                    ({T,invT,dp1,dp,logtisfs,logtbaf});
                }
            } else {
                cout
                << "in FitData::write_avg_normal_modes_decoupling():\n"
                << i <<" cannot open.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    
    std::sort(decoupling_data_tmp.begin(),decoupling_data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(decoupling_data_tmp,decoupling_data);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/dcoupling_avg_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile(o.c_str());
    
    writeFile<<"T 1/T D(p=1) D(p) Log(t.isfs) Log(t.baf)\n";
    
    for (int i=0;i<(int)decoupling_data.size();++i) {
        writeFile
        << decoupling_data.at(i).at(0) << " "
        << decoupling_data.at(i).at(1) << " "
        << decoupling_data.at(i).at(2) << " "
        << decoupling_data.at(i).at(3) << " "
        << decoupling_data.at(i).at(4) << " "
        << decoupling_data.at(i).at(5) << "\n";
    } writeFile.close();
}





void FitData::write_avg_qmatch(const StructureClass& sysVar)
{
    vector<vector<double>> data_tmp;
    vector<vector<double>> data;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_q_match_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    dubVec.clear();
                    istringstream iss(lineContent);
                    while (iss>>dubVar) {
                        dubVec.push_back(dubVar);
                    } data_tmp.push_back(dubVec);
                }
            } else {
                cout
                << "in FitData::write_avg_qmatch():\n"
                << i <<" cannot open.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    
    std::sort(data_tmp.begin(),data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(data_tmp,data);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_q_match_avg_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile(o.c_str());
    
    for (int i=0;i<(int)data.size();++i) {
        for (int ii=0;ii<(int)data.at(i).size();++ii) {
            writeFile << data.at(i).at(ii) << " ";
        } writeFile << "\n";
    } writeFile.close();
}





void FitData::write_avg_qmatch_KWW(const StructureClass& sysVar)
{
    vector<vector<double>> data_tmp;
    vector<vector<double>> data;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_q_match_KWW_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    dubVec.clear();
                    istringstream iss(lineContent);
                    while (iss>>dubVar) {
                        dubVec.push_back(dubVar);
                    } data_tmp.push_back(dubVec);
                }
            } else {
                cout
                << "in FitData::write_avg_qmatch_KWW():\n"
                << i <<" cannot open.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    
    std::sort(data_tmp.begin(),data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(data_tmp,data);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_q_match_avg_KWW_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile(o.c_str());
    
    for (int i=0;i<(int)data.size();++i) {
        for (int ii=0;ii<(int)data.at(i).size();++ii) {
            writeFile << data.at(i).at(ii) << " ";
        } writeFile << "\n";
    } writeFile.close();
}





void FitData::write_smo_qmatch_KWW(const StructureClass& sysVar)
{
    vector<vector<double>> data_tmp;
    vector<vector<double>> data;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_q_match_KWW_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    dubVec.clear();
                    istringstream iss(lineContent);
                    while (iss>>dubVar) {
                        dubVec.push_back(dubVar);
                    } data_tmp.push_back(dubVec);
                }
            } else {
                cout
                << "in FitData::write_avg_qmatch_KWW():\n"
                << i <<" cannot open.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    
    std::sort(data_tmp.begin(),data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(data_tmp,data);
    
    vector<vector<double>> invT=data;
    for (int i=0;i<(int)invT.size();++i) {
        invT.at(i).at(0)=pow(invT.at(i).at(0),-1);
    }
    
    /** smooth beta_isfs vs 1/T **/
    alglib::spline1dinterpolant inp0,inp1,inph,inp2;
    build_penalizedspline_interpolant(invT,0,5,50,0.0);
    inp0=interpolant;
    build_penalizedspline_interpolant(invT,0,5,50,0.5);
    inph=interpolant;
    build_penalizedspline_interpolant(invT,0,5,50,1.0);
    inp1=interpolant;
    build_penalizedspline_interpolant(invT,0,5,50,2.0);
    inp2=interpolant;
    
    /** smooth beta_baf vs 1/T **/
    alglib::spline1dinterpolant inp0b,inp1b,inphb,inp2b;
    build_penalizedspline_interpolant(invT,0,9,50,0.0);
    inp0b=interpolant;
    build_penalizedspline_interpolant(invT,0,9,50,0.5);
    inphb=interpolant;
    build_penalizedspline_interpolant(invT,0,9,50,1.0);
    inp1b=interpolant;
    build_penalizedspline_interpolant(invT,0,9,50,2.0);
    inp2b=interpolant;
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_q_match_smo_KWW_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile(o.c_str());
    
    for (int i=0;i<(int)data.size();++i) {
        double iT=invT.at(i).at(0);
        writeFile
        << iT                     << " " //1/T
        << spline1dcalc(inp0,iT)  << " " //beta_rho0
        << spline1dcalc(inph,iT)  << " " //beta_rho0.5
        << spline1dcalc(inp1,iT)  << " " //beta_rho1.0
        << spline1dcalc(inp2,iT)  << " " //beta_rho2.0
        << "   "
        << spline1dcalc(inp0b,iT) << " " //beta_rho0
        << spline1dcalc(inphb,iT) << " " //beta_rho0.5
        << spline1dcalc(inp1b,iT) << " " //beta_rho1.0
        << spline1dcalc(inp2b,iT) << " " //beta_rho2.0
        << "\n";
    } writeFile.close();
}





void FitData::write_bKWW_vs_T(const StructureClass& sysVar)
{
    vector<vector<double>> data_tmp;
    vector<vector<double>> data_allT;
    vector<vector<double>> data_bmin;
    vector<vector<double>> data_bmax;
    vector<vector<double>> data_min;
    vector<vector<double>> data_max;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_bKWW_vs_q2_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    if (lineContent=="") //T transition @ bottom of each T
                    {
                        if (data_tmp.size()>0)
                        {
                            double f=0,df=0,d2f=0;
                            double bnow=0,qnow=0;
                            double bmin=0,bmax=0;
                            double qMin=0,bMinVal=0;
                            double qMax=0,bMaxVal=0;
                            double qmin=data_tmp.at(0).at(2);//log.q2
                            double qmax=data_tmp.at(data_tmp.size()-1).at(2);//log.q2
                            double dq=fabs(qmax-qmin)/10000.0;
                            double T=data_tmp.at(0).at(0);
                            build_penalizedspline_interpolant(data_tmp,2,7,50,0.0);//(log.q2,bKWW)
                            for (int i=0;i<10000;++i) {
                                qnow=qmin+dq*i;
                                spline1ddiff(interpolant,qnow,f,df,d2f);
                                if (d2f>0) { //convave upward
                                    
                                    if (i==0) {
                                        bmin=fabs(df);
                                    } else {
                                        bnow=fabs(df);
                                    }
                                    if (bnow<bmin) {
                                        if (i>0) bmin=bnow;
                                        qMin=qnow;
                                        bMinVal=f;
                                    }
                                } else { //convave downward
                                    if (i==0) {
                                        bmax=fabs(df);
                                    } else {
                                        bnow=fabs(df);
                                    }
                                    if (bnow<bmax) {
                                        if (i>0) bmax=bnow;
                                        qMax=qnow;
                                        bMaxVal=f;
                                    }
                                }
                            } data_allT.push_back({T,qMin,bMinVal,qMax,bMaxVal});
                            
                        } data_tmp.clear();
                        
                    } else {
                        dubVec.clear();
                        istringstream iss(lineContent);
                        while (iss>>dubVar) {
                            dubVec.push_back(dubVar);
                        } data_tmp.push_back(dubVec);
                    }
                }
            } else {
                cout
                << "in FitData::write_avg_qmatch_KWW():\n"
                << i <<" cannot open.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    
    for (int i=0;i<(int)data_allT.size();++i)
    {
        double T   =data_allT.at(i).at(0);
        double qmin=data_allT.at(i).at(1);
        double bmin=data_allT.at(i).at(2);
        double qmax=data_allT.at(i).at(3);
        double bmax=data_allT.at(i).at(4);
        if (qmin!=0&&bmin!=0) {
            data_bmin.push_back({T,qmin,bmin});
        }
        if (qmax!=0&&bmax!=0&&qmax>0) {
            data_bmax.push_back({T,qmax,bmax});
        }
    }
    
    std::sort(data_bmin.begin(),data_bmin.end(),sortDecreasing0);
    averaging_sorted_data(data_bmin,data_min);
    
    std::sort(data_bmax.begin(),data_bmax.end(),sortDecreasing0);
    averaging_sorted_data(data_bmax,data_max);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_bKWW_min_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile_min(o.c_str());
    
    for (int i=0;i<(int)data_min.size();++i) {
        for (int ii=0;ii<(int)data_min.at(i).size();++ii) {
            writeFile_min << data_min.at(i).at(ii) << " ";
        } writeFile_min << "\n";
    } writeFile_min.close();
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/qtest_bKWW_max_vs_T_");
    o.append(sysVar.get_usic());
    o.append("."+analysispart+".dat");
    ofstream writeFile_max(o.c_str());
    
    for (int i=0;i<(int)data_max.size();++i) {
        for (int ii=0;ii<(int)data_max.at(i).size();++ii) {
            writeFile_max << data_max.at(i).at(ii) << " ";
        } writeFile_max << "\n";
    } writeFile_max.close();
}





void FitData::find_taualpha_frame(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d)
{
    double frame_time_pre=0;
    double frame_time_now=0;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/fit_sExp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        string lineContent,strVar;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> strVar;
            if (strVar=="tauFit") {
                iss >> strVar;
                if (strVar=="(@KWW=0.2)") {
                    iss >> strVar;
                    taualpha_time=atof(strVar.c_str());
                }
            }
        } readFile.close();
    } else {
        cout
        << "in FitData::find_taualpha_frame():\n"
        << i<<" cannot open.\n";
        exit(EXIT_FAILURE);
    }
    i.clear();
    i.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        i.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        i.append("/statistics");
    }
    i.append("/ngp/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    readFile.clear();
    readFile.open(i.c_str());
    if (readFile.is_open())
    {
        string lineContent,strVar;
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in FitData::find_taualpha_frame():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        int frame_index=0;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strVar;
            frame_time_now=atof(strVar.c_str());
            if (taualpha_time<frame_time_now) {
                if (taualpha_time>frame_time_pre) {
                    double dtnow=fabs(taualpha_time-frame_time_now);
                    double dtpre=fabs(taualpha_time-frame_time_pre);
                    if (frame_index>1) {
                        if (dtnow<dtpre) {
                            taualpha_frame=frame_index;
                            taualpha_frame_time=frame_time_now;
                        } else {
                            taualpha_frame=frame_index-1;
                            taualpha_frame_time=frame_time_pre;
                        }
                    } else {
                        cout
                        << "in FitData::find_taualpha_frame():\n"
                        << "frame_index of tau_alpha ("<<frame_index<<") is not correct. "
                        << "Please check.\n\n";
                        exit(EXIT_FAILURE);
                    }
                }
            } frame_time_pre=frame_time_now;
            ++frame_index;
        } readFile.close();
    } else {
        cout
        << "in FitData::find_taualpha_frame():\n"
        << i<<" cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void FitData::find_ngp_block_frame(const StructureClass& sysVar,
                                   const int n_trl,
                                   const int n_sys,
                                   const double Temp_d)
{
    find_ngp_peak_frame(sysVar,n_trl,n_sys,Temp_d);
    vector<vector<double>> ngp_data_tmp=ngp_data_sortincreasing;
    std::sort(ngp_data_tmp.begin(),ngp_data_tmp.end(),sortDecreasing0);//sort time
    ngp_block_frame=ngp_data_tmp.at(0).at(2);//(time,ngp,frame_index)
}





void FitData::write_MSDt_equ(const StructureClass& sysVar)
{
    /** This function tracks <r2> at a specified time by interpolating MSD
     ** raw data by penalized spline **/
    string o;
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/MSDt_equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream writeFile(o.c_str(),ofstream::app); // 'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_data=(int)taueq_sortdecreasing.size();
            
            streamsize ss=writeFile.precision();
            
            for (int i=0; i<n_data; ++i) {
                
                double Temp_d=std::round(taueq_sortdecreasing.at(i).at(0)*corrFacforT);
                read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
                double log10time=taueq_sortdecreasing.at(i)[1];//log10 value
                bool flag=false;
                
                /* <r2> at a specified time */
                double r2attime=
                pow(10,tfunction_at_given_time(loglogmsd_sortincreasing,log10time,flag));
                
                /* time at a specified <r2> */
                double r2val=1.0;
                if (false)
                {
                    if (sysVar.get_systemUnit()=="real") {
                        r2val=pow((2*pi())/1.19,2);//Angstrom^2
                    } else if (sysVar.get_systemUnit()=="lj") {
                        r2val=pow((2*pi())/7.07,2);//LJ^2 ~ nm^2
                    } else if (sysVar.get_systemUnit()=="metal") {
                        r2val=pow((2*pi())/1.19,2);//Angstrom^2
                    }
                }
                double timeatr2=
                pow(10,tfunction_at_given_valu(loglogmsd_sortincreasing,log10(r2val),flag));
                
                if (false) {
                    cout
                    << sysVar.get_usic()
                    << "_00"+to_string((long long int)n_trl)
                    << "_"+sysVar.get_nameString(n_sys)
                    << "_T"+to_string((long long int)Temp_d)
                    << "\n"
                    << "MSD(@"<<log10time<<") "<< r2attime
                    << "\n";
                }
                
                writeFile
                << fixed << setprecision(10)//NOTE:to be consistent in T format
                << taueq_sortdecreasing.at(i)[0] << "  ";
                
                writeFile.precision(ss);
                writeFile << resetiosflags(ios::fixed|ios::showpoint);
                
                writeFile << log10(timeatr2)               << "  ";//translation
                writeFile << taueq_sortdecreasing.at(i)[1] << "  ";//relaxation
                writeFile << r2attime                      << "\n";
                
            } writeFile.close();
        }
    }
}





void FitData::write_MSDt_equ_avg(const StructureClass& sysVar)
{
    vector<vector<double>> data_tmp;
    vector<vector<double>> data;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            string i;
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/MSDt_equ_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readFile(i.c_str());
            if (readFile.is_open())
            {
                double dubVar=0;
                string lineContent;
                vector<double> dubVec;
                while(getline(readFile,lineContent))
                {
                    dubVec.clear();
                    istringstream iss(lineContent);
                    while (iss>>dubVar) {
                        dubVec.push_back(dubVar);
                    } data_tmp.push_back(dubVec);
                }
            } else {
                cout
                << "in FitData::write_MSDt_equ_avg():\n"
                << i <<" cannot open.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    
    std::sort(data_tmp.begin(),data_tmp.end(),sortDecreasing0);
    averaging_sorted_data(data_tmp,data);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/MSDt_equ_avg_");
    o.append(sysVar.get_usic());
    o.append(".dat");
    ofstream writeFile(o.c_str());
    
    for (int i=0;i<(int)data.size();++i) {
        for (int ii=0;ii<(int)data.at(i).size();++ii) {
            writeFile << data.at(i).at(ii) << " ";
        } writeFile << "\n";
    } writeFile.close();
}





void FitData::write_DWF_equ(const StructureClass& sysVar)
{
    string o;
    
    /* write out "equilibrium" DWF to file (similar to taueq file) */
    //--------------------------------------------------------------------------
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/DWF_equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_DWF(o.c_str(),ofstream::app);//'append'
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/DWF_equ_inverseT_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_DWFinv(o.c_str(),ofstream::app);//'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            read_individual_DWF(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_data =(int)taueq_sortdecreasing.size();
            int dwfSize=(int)dwf_sortdecreasing.size();
            
            streamsize ss=write_DWF.precision();
            
            for (int i=0; i<n_data; ++i) {
                for (int ii=0; ii<dwfSize; ++ii) {
                    if (taueq_sortdecreasing.at(i)[0]==dwf_sortdecreasing.at(ii)[0])
                    {
                        double Temp_d=std::round(taueq_sortdecreasing.at(i)[0]*corrFacforT);
                        
                        /** find transition from ballistic to caging regimes **/
                        find_balltocage_time(sysVar,n_trl,n_sys,Temp_d);
                        
                        write_DWF
                        << fixed << setprecision(10) // NOTE: be consistent in temperature format
                        << dwf_sortdecreasing.at(ii)[0] << " ";
                        
                        write_DWF.precision(ss);
                        write_DWF << resetiosflags(ios::fixed|ios::showpoint);
                        
                        write_DWF
                        << dwf_sortdecreasing.at(ii)[1] << "\n";
                        write_DWFinv
                        << convert/dwf_sortdecreasing.at(ii)[0] << " "
                        << dwf_sortdecreasing.at(ii)[1] << "\n";
                        break;
                    }
                }
            } write_DWF.close(); write_DWFinv.close();
        }
    }
}





void FitData::write_DWF_equ_avg(const StructureClass& sysVar)
{
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/DWF_equ_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_DWF(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/DWF_equ_avg_inverseT_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_DWFinv(o.c_str(),ofstream::app); // 'append'
        
        /* NOTE: all */
        //-----------------------------------
        read_all_taueq_data(sysVar,n_sys);
        read_all_DWF(sysVar,n_sys);
        //-----------------------------------
        
        int n_data =(int)taueq_sortdecreasing.size();
        int dwfSize=(int)dwf_sortdecreasing.size();
        
        streamsize ss=write_DWF.precision();
        
        for (int i=0; i<n_data; ++i) {
            for (int ii=0; ii<dwfSize; ++ii) {
                if (taueq_sortdecreasing.at(i)[0]==dwf_sortdecreasing.at(ii)[0])
                {
                    write_DWF
                    << fixed << setprecision(10) // NOTE: be consistent in temperature format
                    << dwf_sortdecreasing.at(ii)[0] << " ";
                    
                    write_DWF.precision(ss);
                    write_DWF << resetiosflags(ios::fixed|ios::showpoint);
                    
                    write_DWF
                    << dwf_sortdecreasing.at(ii)[1] << "\n";
                    write_DWFinv
                    << convert/dwf_sortdecreasing.at(ii)[0] << " "
                    << dwf_sortdecreasing.at(ii)[1] << "\n";
                    break;
                }
            }
        } write_DWF.close(); write_DWFinv.close();
    }
}





bool FitData::skim_individual_correlation_data(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d,
                                               const int frame,
                                               const int waveind,
                                               const string& target)
{
    /** save a copy of original relaxation target **/
    string relaxation_target_org=relaxation_target;
    if (target!="") {
        relaxation_target=target;
    }
    
    /** read target relaxtion profile **/
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        i.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        i.append("/statistics");
    }
    i.append("/"+relaxation_target);
    i.append("/"+relaxation_target+"_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    if (frame!=0) i.append("_frame"+to_string((long long int)frame));
    if (is_qspectrum) {
        if (waveind>=0) {
            i.append("_waveindex"+to_string((long long int)waveind));
        } else {
            i.append("_waveindex"+to_string((long long int)waveindices.at(2)));//default q
        }
    } i.append("."+analysispart+".dat");
    
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        string strData;
        double time=0,corr=0;
        /* read 1st line */
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in FitData::skim_individual_correlation_data():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        /* read 2nd line */
        bool is_fullyRelaxed=false;
        int  firstpassage=0;
        getline(readFile,lineContent);
        
        //assign original values
        is_use_sExp_cut_hi=is_use_sExp_cut_hi_org;
        is_use_plateautime=is_use_plateautime_org;
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> time;//time
            iss >> corr;//Fs
            if (time>sExp_tau_lo&&firstpassage==0) {
                if (corr>sExp_cut_hi) {
                    is_use_sExp_cut_hi=false;
                    is_use_plateautime=true;
                } else {
                    is_use_sExp_cut_hi=true;
                    is_use_plateautime=false;
                } ++firstpassage;
            }
        } readFile.close();
        
        /** Regard as fully relaxed when final point is below a threshold **/
        if (corr>sExp_tauFit) {
            is_fullyRelaxed=false;
            cout
            << relaxation_target<<"_"
            << sysVar.get_usic()<<"_"
            << "00"<<n_trl<<"_"<<n_sys<<"_T"<<Temp_d<<" ";
            if (is_qspectrum)if(waveind>=0)cout<<"(waveind "<<waveind<<") ";
            cout << "is NOT fully relaxed. Data not counted.\n";
        } else {
            is_fullyRelaxed=true;
        }
        
        /** resume original value of relaxation target **/
        relaxation_target=relaxation_target_org;
        
        return is_fullyRelaxed;
    }
    else
    {
        cout
        << "in FitData::skim_individual_correlation_data():\n"
        << i << " cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void FitData::read_individual_correlation_data(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d,
                                               const int frame,
                                               const int waveind,
                                               const string& target)
{
    correlation_original.clear();        //(time,corr)
    correlation_semilog.clear();         //(log.time,corr)
    correlation_loglog.clear();          //(log.time,log.corr)
    correlation_sortincreasing.clear();  //(time,corr)
    lncorrelation_sortincreasing.clear();//(time,ln.corr)
    
    /** save a copy of original relaxation target **/
    string relaxation_target_org=relaxation_target;
    if (target!="") {
        relaxation_target=target;
    }
    
    /** read the target relaxtion profile **/
    string in;
    in.append(return_AnalysisFolderPath(sysVar));
    if (sysVar.get_is_aging()) {
        in.append("/statistics_Regime_"+to_string(sysVar.get_current_regime()));
    } else {
        in.append("/statistics");
    }
    in.append("/"+relaxation_target);
    in.append("/"+relaxation_target+"_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    if (frame!=0) in.append("_frame"+to_string((long long int)frame));
    if (is_qspectrum) {
        if (waveind>=0) {
            in.append("_waveindex"+to_string((long long int)waveind));
        } else {
            in.append("_waveindex"+to_string((long long int)waveindices.at(2)));//default q
        }
    } in.append("."+analysispart+".dat");
    
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        string strData;
        double time=0,corr=0;
        
        //auto start_position=readFile.tellg(); /* mark the beginning of file */
        //readFile.clear(); /* reset state flag of read file */
        //readFile.seekg(start_position); /* rewind to start of file */
        
        /* read 1st line */
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in FitData::read_individual_correlation_data():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        /* read 2nd line */
        getline(readFile,lineContent);
        istringstream iss1(lineContent);
        iss1 >> corr;
        wavenumber = corr;
        
        int cut_index_hi=0;
        int cut_index_lo=0;
        int cut_index_total=0;
        int first_passage=0;
        
        sExp_cut_index_lo=0;
        
        while (getline(readFile,lineContent))
        {
            ++cut_index_hi;
            ++cut_index_lo;
            ++cut_index_total;
            
            istringstream iss(lineContent);
            iss >> time;//time
            iss >> corr;//fs(q,time)
            
            /** these containers store raw relaxation data **/
            //------------------------------------------------------------------
            correlation_original.push_back({time,corr});
            if (time>0) {
                correlation_semilog.push_back({log10(time),corr});
                if (corr>0) {
                    correlation_loglog.push_back({log10(time),log10(corr)});
                }
            }
            //------------------------------------------------------------------
            
            
            /** cutoff data based on time or value **/
            //------------------------------------------------------------------
            if (is_use_sExp_cut_hi)
            {
                if (corr<=sExp_cut_hi)
                {
                    sExp_cut_index_hi = cut_index_hi;
                    --cut_index_hi;
                    // NOTE:
                    // cut_index_hi fixed at where the higher cutoff is
                }
            }
            if (is_use_plateautime)
            {
                if (time>sExp_tau_lo)
                {
                    sExp_cut_index_hi = cut_index_hi;
                    --cut_index_hi;
                    // NOTE:
                    // cut_index_hi fixed at where the higher cutoff is
                }
            }
            if (corr<sExp_cut_lo)
            {
                ++first_passage;
                if (first_passage==1)
                {
                    sExp_cut_index_lo = cut_index_lo-1;
                    // NOTE:
                    // cut_index_lo fixed at where the lower cutoff is
                }
            }
            //==================================================================
            // Correlation data for fitting:
            // 1st point below sExp_cut_hi to 1st point above sExp_cut_lo
            //==================================================================
            if (is_use_sExp_cut_hi)
            {
                if ((corr<=sExp_cut_hi)&&(corr>=sExp_cut_lo)&&(sExp_cut_index_lo==0))
                {
                    correlation_sortincreasing.push_back({time,corr});
                    lncorrelation_sortincreasing.push_back({time,log(corr)});
                }
            }
            //==================================================================
            // Correlation data for fitting:
            // 1st point past sExp_tau_lo to 1st point above sExp_cut_lo
            //==================================================================
            if (is_use_plateautime)
            {
                if ((time>sExp_tau_lo)&&(corr>=sExp_cut_lo)&&(sExp_cut_index_lo==0))
                {
                    correlation_sortincreasing.push_back({time,corr});
                    lncorrelation_sortincreasing.push_back({time,log(corr)});
                }
            }
        } readFile.close();
        
        /** resume original value of relaxation target **/
        relaxation_target=relaxation_target_org;
    }
    else
    {
        cout
        << "in FitData::read_individual_correlation_data():\n"
        << in << " cannot open.\n";
        exit(EXIT_FAILURE);
    }
    
}





void FitData::read_individual_taueq_data(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string in;
    in.append(return_AnalysisFolderPath(sysVar));
    in.append("/fit_data");
    in.append("/taueq_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append(".dat");
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);//(T,taueq)
            iss >> Tdata;
            iss >> tauFitdata;
            taueq.push_back({Tdata,tauFitdata});
        } readFile.close();
    } else {
        cout
        << "in FitData::read_individual_taueq_data():\n"
        << "file ("+in+") cannot open.\n";
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_individual_taueq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





double FitData::get_individual_taueq_at_T(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d)
{
    double T=0,taueq=0;
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
        T=std::round(taueq_sortdecreasing.at(i).at(0)*corrFacforT);
        if (T==Temp_d) {
            taueq=taueq_sortdecreasing.at(i).at(1);
            break;
        }
    } return taueq;
}





void FitData::read_individual_taueq_data(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys,
                                         const std::string& model)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string in;
    in.append(return_AnalysisFolderPath(sysVar));
    in.append("/fit_data");
    in.append("/taueq_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append(".dat");
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent); // (T,taueq)
            if (model=="Arrhenius")
            {
                iss >> Tdata;
                if (Tdata>=cutTforArrhenius) {
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            } else {
                iss >> Tdata;
                iss >> tauFitdata;
                taueq.push_back({Tdata,tauFitdata});
            }
        } readFile.close();
    } else {
        cout
        << "in FitData::read_individual_taueq_data():\n"
        << "file ("+in+") cannot open." << "\n";
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_individual_taueq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitData::read_individual_temps_data(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys)
{
    temps_sortdecreasing.clear();
    temps_sortincreasing.clear();
    vector<vector<double>> temps;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/tempsInfo_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    ifstream readtempsdata(i.c_str());
    if (readtempsdata.is_open()) {
        string lineContent;
        double Tdata=0;
        while (getline(readtempsdata,lineContent)) {
            istringstream iss(lineContent); // T
            iss >> Tdata;
            temps.push_back({Tdata});
        } readtempsdata.close();
    } else {
        cout
        << "in FitData::read_individual_temps_data():\n"
        << "file ("+i+") cannot open.\n";
    }
    try {
        if(temps.size()==0) throw 0;
        std::sort(temps.begin(),temps.end(),sortDecreasing0);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing0);
        averaging_sorted_data(temps,temps_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_individual_temps_data():\n"
        << "temps size = " << i << "\n";
        temps.push_back({0});
        std::sort(temps.begin(),temps.end(),sortDecreasing0);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing0);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
}





void FitData::read_all_temps_data(const StructureClass& sysVar,
                                  const int n_sys)
{
    temps_sortdecreasing.clear();
    temps_sortincreasing.clear();
    vector<vector<double>> temps;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/tempsInfo_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readtempsdata(i.c_str());
        if (readtempsdata.is_open()) {
            string lineContent;
            double Tdata=0;
            while (getline(readtempsdata,lineContent)) {
                istringstream iss(lineContent); // T
                iss >> Tdata;
                temps.push_back({Tdata});
            } readtempsdata.close();
        } else {
            cout
            << "in FitData::read_all_temps_data():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(temps.size()==0) throw 0;
        std::sort(temps.begin(),temps.end(),sortDecreasing0);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing0);
        averaging_sorted_data(temps,temps_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_temps_data():\n"
        << "temps size = " << i << "\n";
        temps.push_back({0});
        std::sort(temps.begin(),temps.end(),sortDecreasing0);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing0);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
}





void FitData::read_all_thermo_data(const StructureClass& sysVar,
                                   const int n_sys,
                                   const string& keyword)
{
    thermo_sortdecreasing.clear();
    vector<vector<double>> thermo;
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/thermo_"+keyword+"_");//NOTE: thermo style keyword
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,avg=0,stdev=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);//(T,avg.thermo,stdev.thermo)
                iss >> Tdata;
                iss >> avg;
                iss >> stdev;
                thermo.push_back({Tdata,avg,stdev});
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_thermo_data():\n"
            << "file ("+i+") cannot open.\n";
            exit(EXIT_FAILURE);
        }
    }
    try {
        if(thermo.size()==0) throw 0;
        std::sort(thermo.begin(),thermo.end(),sortDecreasing0);
        averaging_sorted_data(thermo,thermo_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_thermo_data():\n"
        << "thermo("+keyword+") data size = " << i << "\n";
        thermo.push_back({0,0});
        std::sort(thermo.begin(),thermo.end(),sortDecreasing0);
        averaging_sorted_data(thermo,thermo_sortdecreasing);
    }
    
}





void FitData::read_all_thermo_data(const StructureClass& sysVar,
                                   const int n_sys,
                                   const string& keyword,
                                   const double TA_d,
                                   bool is_ltTA)
{
    thermo_sortdecreasing.clear();
    vector<vector<double>> thermo;
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/thermo_"+keyword+"_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,avg=0,stdev=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);//(T,avg.thermo,stdev.thermo)
                iss >> Tdata;
                iss >> avg;
                iss >> stdev;
                if (is_ltTA) {
                    if (Tdata<=TA_d) {//only collect T<=TA
                        thermo.push_back({Tdata,avg,stdev});
                    }
                } else {
                    if (Tdata>TA_d) {//only collect T>TA
                        thermo.push_back({Tdata,avg,stdev});
                    }
                }
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_thermo_data():\n"
            << "file ("+i+") cannot open.\n";
            exit(EXIT_FAILURE);
        }
    }
    try {
        if(thermo.size()==0) throw 0;
        std::sort(thermo.begin(),thermo.end(),sortDecreasing0);
        averaging_sorted_data(thermo,thermo_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_thermo_data():\n"
        << "thermo("+keyword+") data size = " << i << "\n";
        thermo.push_back({0,0});
        std::sort(thermo.begin(),thermo.end(),sortDecreasing0);
        averaging_sorted_data(thermo,thermo_sortdecreasing);
    }
    
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys)
{
    taueq_sortdecreasing.clear();
    taueqinvT_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent); // (T,taueq)
                iss >> Tdata;
                iss >> tauFitdata;
                taueq.push_back({Tdata,tauFitdata});
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_taueq_data():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_taueq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
        double invT =convert/taueq_sortdecreasing.at(i).at(0);
        double taueq=taueq_sortdecreasing.at(i).at(1);
        taueqinvT_sortdecreasing.push_back({invT,taueq});
    }
}





void FitData::read_all_r2teq_data(const StructureClass& sysVar,
                                  const int n_sys)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string in;
    in.append(return_AnalysisFolderPath(sysVar));
    in.append("/fit_data");
    in.append("/MSDt_equ_avg_");
    in.append(sysVar.get_usic());
    in.append(".dat");
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,logr2teq=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);//(T,logr2teq)
            iss >> Tdata;
            iss >> logr2teq;
            taueq.push_back({Tdata,logr2teq});
        } readFile.close();
    } else {
        cout
        << "in FitData::read_all_r2teq_data():\n"
        << "file ("+in+") cannot open.\n";
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_r2teq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys,
                                  const double TA_d,
                                  bool is_ltTA)
{
    taueq_sortdecreasing.clear();
    taueqinvT_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);//(T,taueq)
                iss >> Tdata;
                iss >> tauFitdata;
                if (is_ltTA) {
                    if (Tdata<=TA_d) {//only collect T<=TA
                        taueq.push_back({Tdata,tauFitdata});
                    }
                } else {
                    if (Tdata>TA_d) {//only collect T>TA
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_taueq_data():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_taueq_data()\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
        double invT =convert/taueq_sortdecreasing.at(i).at(0);
        double taueq=taueq_sortdecreasing.at(i).at(1);
        taueqinvT_sortdecreasing.push_back({invT,taueq});
    }
}





void FitData::read_all_r2teq_data(const StructureClass& sysVar,
                                  const int n_sys,
                                  const double TA_d,
                                  bool is_ltTA)
{
    taueq_sortdecreasing.clear();
    taueqinvT_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/MSDt_equ_avg_");
    i.append(sysVar.get_usic());
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        string lineContent;
        double Tdata=0,logr2teq=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);//(T,r2teq)
            iss >> Tdata;
            iss >> logr2teq;
            if (is_ltTA) {
                if (Tdata<=TA_d) {//only collect T<=TA
                    taueq.push_back({Tdata,logr2teq});
                }
            } else {
                if (Tdata>TA_d) {//only collect T>TA
                    taueq.push_back({Tdata,logr2teq});
                }
            }
        } readFile.close();
    } else {
        cout
        << "in FitData::read_all_r2teq_data():\n"
        << "file ("+i+") cannot open.\n";
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_r2teq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
        double invT =convert/taueq_sortdecreasing.at(i).at(0);
        double taueq=taueq_sortdecreasing.at(i).at(1);
        taueqinvT_sortdecreasing.push_back({invT,taueq});
    }
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys,
                                  const std::string& cond)
{
    taueq_sortdecreasing.clear();
    taueqinvT_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent); // (T,taueq)
                iss >> Tdata;
                if (cond=="Arrhenius")
                {
                    if (Tdata>=cutTforArrhenius) {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                else if (cond=="continuous")
                {
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
                else if (cond=="neighboring")
                {
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
                else if (cond=="partial_lt") // T<Tx
                {
                    if (get_is_fit_by_Tc()) {
                        if (Tdata<get_Tc_MCT()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    } else if (get_is_fit_by_TA()) {
                        if (Tdata<get_TA_avg()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    } else {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                else if (cond=="partial_gt") // T>Tx
                {
                    if (get_is_fit_by_Tc()) {
                        if (Tdata>get_Tc_MCT()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    } else if (get_is_fit_by_TA()) {
                        if (Tdata>get_TA_avg()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    } else {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                else
                {
                    iss >> Tdata;
                    iss >> tauFitdata;
                    //if (tauFitdata<log10(tauFit_cut)) {
                    //    taueq.push_back({Tdata,tauFitdata});}
                    taueq.push_back({Tdata,tauFitdata});
                }
            } readFile.close();
        }
        else {
            cout
            << "in FitData::read_all_taueq_data():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_taueq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
        double invT =convert/taueq_sortdecreasing.at(i).at(0);
        double taueq=taueq_sortdecreasing.at(i).at(1);
        taueqinvT_sortdecreasing.push_back({invT,taueq});
    }
}





void FitData::read_all_sExp_params(const StructureClass& sysVar,
                                   const int n_sys)
{
    param_sortdecreasing.clear();
    vector<vector<double>> param;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/sExp_params_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open()) {
            string lineContent;
            double tmpdata;
            vector<double> rowdata;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);
                rowdata.clear();
                while (iss>>tmpdata) rowdata.push_back(tmpdata);
                param.push_back(rowdata);
            } readFile.close();
        } else {
            cout
            << "in FitData::read_all_sExp_params():\n"
            << "file ("+i+") cannot open.\n";
        }
    }
    try {
        if(param.size()==0) throw 0;
        std::sort(param.begin(),param.end(),sortDecreasing0);
        averaging_sorted_data(param,param_sortdecreasing);
    } catch (int i) {
        cout
        << "in FitData::read_all_sExp_params():\n"
        << "param size = " << i << "\n";
        param.push_back({0,0});
        std::sort(param.begin(),param.end(),sortDecreasing0);
        averaging_sorted_data(param,param_sortdecreasing);
    }
}





void FitData::Tc_correction(const std::string& model)
{
    if ((model=="MCT")||(model=="SOU")) {
        
        size_t vecSize=xdata.size();
        string tleq_str=xdata.at(vecSize-2);
        double tleq=stod(tleq_str);
        
        /* the original coefficients */
        double A  = get_coeffs_vD()[0];
        double Tc = get_coeffs_vD()[1]; Tc = tleq*Tc_percent;
        double r  = get_coeffs_vD()[2];
        set_coeffs({A,Tc,r});
        
        double A_bndu  = get_coeffs_bndu_vD()[0];
        double Tc_bndu = get_coeffs_bndu_vD()[1]; Tc_bndu = Tc;
        double r_bndu  = get_coeffs_bndu_vD()[2];
        set_coeffs_bndu({A_bndu,Tc_bndu,r_bndu});
    }
}





void FitData::write_fitcutoff(std::ofstream& outputFile,
                              const StructureClass& sysVar,
                              const std::string& model)
{
    if (model=="KWW"||model=="KWW_lnFs"||model=="KWW_pwr"||model=="mKWW"||model=="mKWW_pwr")
    {
        outputFile
        << amdat_version
        << "\n\n"
        << "Wavenumber " << wavenumber
        << "\n\n";
        
        if (is_use_sExp_cut_hi) {
            outputFile
            << "sExp_cut_hi "
            << sExp_cut_hi << " "
            << "\n";
        }
        if (is_use_plateautime) {
            outputFile
            << "sExp_tau_lo "
            << sExp_tau_lo << " "
            << "\n";
        }
        outputFile
        << "sExp_cut_lo "
        << sExp_cut_lo << " "
        << "\n\n";
        outputFile
        << "cut_index_hi(1-based) " << sExp_cut_index_hi
        << "\n"
        << "cut_index_lo(1-based) " << sExp_cut_index_lo
        << "\n\n";
        
    }
    else if (model=="KWW_full"||model=="KWW_full_lnFs")
    {
        outputFile
        << amdat_version
        << "\n\n"
        << "Wavenumber " << wavenumber
        << "\n\n";
    }
    else if (model=="Arrhenius")
    {
        //
    }
    else
    {
        outputFile << "tauFit_cut " << "\n";
        for (size_t i=0; i<sysVar.get_n_regime(); ++i) {
            double previous_teq,current_teq,teq;
            if (i>0) {
                previous_teq = sysVar.get_equilibration_times().at(i-1);
                current_teq  = sysVar.get_equilibration_times().at(i);
                teq = max(previous_teq,current_teq);
            } else {
                teq = sysVar.get_equilibration_times().at(i);
            }
            outputFile
            << "Regime_" << i << "\n"
            << teq/n_equ_blocks<< "\n";
        } outputFile << "\n";
    }
}





void FitData::write_spline_actual(ofstream& outputFile,
                                  const string& model,
                                  const vector<vector<double>>& vec2d_sortdecreasing)
{
    outputFile << "\n"
    << "x  y  y_"<<model<<"\n"
    << "=========================================== " << "\n";
    for (int i=0; i<vec2d_sortdecreasing.size(); ++i) {
        outputFile
        << vec2d_sortdecreasing.at(i)[0] << " "
        << vec2d_sortdecreasing.at(i)[1] << " "
        << spline1dcalc(interpolant,vec2d_sortdecreasing.at(i)[0])<<"\n";
    } outputFile << "\n\n";
}





void FitData::write_errcurve_actual(ofstream& outputFile,
                                    const string& model,
                                    const vector<vector<double>>& vec2d_sortdecreasing)
{
    outputFile << "\n"
    << "x  y  y_fit  ErrCurve  " << "\n"
    << "=========================================== " << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << vec2d_sortdecreasing.at(i)[0] << " "
        << vec2d_sortdecreasing.at(i)[1] << " "
        << calc_y_given_x(c,vec2d_sortdecreasing.at(i)[0],model) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n\n";
}





void FitData::write_errcurve_log(ofstream& outputFile,
                                 const string& model,
                                 const vector<vector<double>>& vec2d_sortdecreasing)
{
    outputFile << "\n"
    << "x  y  y_fit  ErrCurve  " << "\n"
    << "=========================================== " << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << vec2d_sortdecreasing.at(i)[0] << " "
        << vec2d_sortdecreasing.at(i)[1] << " "
        << log10(calc_y_given_x(c,vec2d_sortdecreasing.at(i)[0],model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n\n";
}





void FitData::write_errcurve2d_log(ofstream& outputFile,
                                   const string& model,
                                   const vector<vector<double>>& vec2d_sortdecreasing)
{
    outputFile << "\n"
    << "{x}  y  y_fit  ErrCurve  " << "\n"
    << "=========================================== " << "\n";
    vector<double> v;
    int n_column=(int)vec2d_sortdecreasing.at(0).size();
    for (int i=0; i<rep.errcurve.length(); ++i) {
        v.clear(); v=vec2d_sortdecreasing.at(i);
        for (int ii=0; ii<(n_column-1); ++ii) {
            outputFile << v.at(ii) << " ";
        }
        outputFile
        << v.at(n_column-1) << " "
        << log10(calc_time_given_vecinput(c,v,model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n\n";
}




void FitData::compute_Tg_fragility(StructureClass& sysVar,
                                   const alglib::real_1d_array& c,
                                   const std::string& model)
{
    /** extrapolated Tg and fragility **/
    //--------------------------------------------------------------------------
    double Tg_extrp = calc_x_given_y(c,extrp_time,model)[0];
    double m_extrp  = calc_x_given_y(c,extrp_time,model)[1];
    sysVar.get_Tg_extrp().at(current_regime).push_back(Tg_extrp);
    sysVar.get_m_extrp().at(current_regime).push_back(m_extrp);
    
    /** computational Tg and fragility **/
    //--------------------------------------------------------------------------
    double Tg_compu = calc_x_given_y(c,compu_time,model)[0];
    double m_compu  = calc_x_given_y(c,compu_time,model)[1];
    sysVar.get_Tg_compu().at(current_regime).push_back(Tg_compu);
    sysVar.get_m_compu().at(current_regime).push_back(m_compu);
}





void FitData::compute_VFT_Tg_fragility(StructureClass& sysVar,
                                       const alglib::real_1d_array& c)
{
    /** VFT Tg and fragility **/
    //--------------------------------------------------------------------------
    double Tg_VFT = calc_x_given_y(c,extrp_time,"VFT")[0];
    double m_VFT  = calc_x_given_y(c,extrp_time,"VFT")[1];
    sysVar.get_Tg_VFT().at(current_regime).push_back(Tg_VFT);
    sysVar.get_m_VFT().at(current_regime).push_back(m_VFT);
}




void FitData::compute_COOP_Tg_fragility(StructureClass& sysVar,
                                        const alglib::real_1d_array& c)
{
    /** COOP Tg and fragility **/
    //--------------------------------------------------------------------------
    double Tg_COOP = calc_x_given_y(c,extrp_time,"COOP")[0];
    double m_COOP  = calc_x_given_y(c,extrp_time,"COOP")[1];
    sysVar.get_Tg_COOP().at(current_regime).push_back(Tg_COOP);
    sysVar.get_m_COOP().at(current_regime).push_back(m_COOP);
}





void FitData::find_largest_Tg_fragility(StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys,
                                        const alglib::real_1d_array& c,
                                        const std::string& model)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    /* extrapolated Tg and fragility */
    //------------------------------------------------------------------
    double Tg_extrp = calc_x_given_y(c,extrp_time,model)[0];
    double m_extrp  = calc_x_given_y(c,extrp_time,model)[1];
    //------------------------------------------------------------------
    
    if (index==0) {
        set_largest_Tg(Tg_extrp);
        set_index_largest_Tg(index);
        set_largest_m(m_extrp);
        set_index_largest_m(index);
    } else {
        if (Tg_extrp>get_largest_Tg()) {
            set_largest_Tg(Tg_extrp);
            set_index_largest_Tg(index);
        }
        if (m_extrp>get_largest_m()) {
            set_largest_m(m_extrp);
            set_index_largest_m(index);
        }
    }
    /* calculate accumulative fragility avg and stdev */
    sysVar.set_calcvector(sysVar.get_m_extrp().at(current_regime));
    double m_accumAvg=sysVar.calc_mean();
    double m_accumStdev=sysVar.get_sample_stdev();
    relativeStdev=m_accumStdev/m_accumAvg;
}





vector<double> FitData::write_Tg_fragility(std::ofstream& outputFile,
                                           const alglib::real_1d_array& c,
                                           const std::string& model)
{
    vector<double> glass_extrp,glass_compu;
    
    /* extrapolated Tg and fragility */
    //------------------------------------------------------------------
    glass_extrp = calc_x_given_y(c,extrp_time,model);
    
    /* computational Tg and fragility */
    //------------------------------------------------------------------
    glass_compu = calc_x_given_y(c,compu_time,model);
    
    outputFile
    << "-------------------------------------\n"
    << "Tg_extrp [@("<<extrp_time<<")timeUnit] " << glass_extrp.at(0) << "\n"
    << "m_extrp  [@Tg_extrp] " << glass_extrp.at(1) << "\n\n"
    << "Tg_compu [@("<<compu_time<<")timeUnit] " << glass_compu.at(0) << "\n"
    << "m_compu  [@Tg_compu] " << glass_compu.at(1) << "\n"
    << "-------------------------------------\n\n";
    
    return
    { glass_extrp.at(0),glass_extrp.at(1),glass_compu.at(0),glass_compu.at(1)};
}





vector<double> FitData::write_Yg_fragility(std::ofstream& outputFile,
                                           const alglib::real_1d_array& c,
                                           const vector<double>& ref,
                                           const std::string& model)
{
    vector<double> glass_extrp;
    double u2g_extrp=0;
    //double logtaur=ref.at(0);
    double u2r=ref.at(1);
    //double Tr=ref.at(2);
    
    /* extrapolated Yg and fragility */
    //------------------------------------------------------------------
    glass_extrp = calc_Y_given_time(c,extrp_time,ref,model);
    u2g_extrp   = u2r*pow(glass_extrp.at(0),-1);
    //------------------------------------------------------------------
    
    outputFile
    << "-------------------------------------" << "\n"
    << "Yg_extrp [@("<<extrp_time<<")timeUnit] "
    << glass_extrp.at(0)
    << "\n"
    << "mg_extrp [@Tg_extrp] "
    << glass_extrp.at(1)
    << "\n"
    << "u2g_extrp[@("<<extrp_time<<")timeUnit] "
    << u2g_extrp
    << "\n"
    << "-------------------------------------" << "\n"
    << "\n";
    
    return {glass_extrp.at(0),glass_extrp.at(1),u2g_extrp};
}





void FitData::write_fit_correlation(std::ofstream& outputFile,
                                    const std::string& tcalc)
{
    if (is_use_gammafunc) {
        outputFile
        << "tauFit (integral by Gamma function) "
        << sExp_tauFit_gammafunction(c,tcalc)
        << "\n\n";
    } else {
        outputFile
        << "tauFit (@KWW=" << get_sExp_tauFit() << ") "
        << sExp_tauFit_interpolateFvalue(c,tcalc)
        << "\n\n";
    }
    
    outputFile
    << "time  correlation_fit            " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<correlation_sortincreasing.size(); ++i) {
        outputFile
        << correlation_sortincreasing.at(i)[0] << " "
        << sExp_func_value(c,correlation_sortincreasing.at(i)[0],tcalc)
        << "\n";
    } outputFile << "\n\n";
    
    outputFile
    << "time  correlation(time)          " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<correlation_original.size(); ++i) {
        outputFile
        << correlation_original.at(i)[0] << " "
        << correlation_original.at(i)[1]
        << "\n";
    }
}





void FitData::write_fit_correlation_spline(std::ofstream& outputFile)
{
    spline1dinterpolant org=interpolant;
    
    build_penalizedspline_interpolant(correlation_sortincreasing,1,0);//NOTE: time vs F
    outputFile
    << "tauFit (@Corr="<<get_sExp_tauFit()<<") "
    << spline1dcalc(interpolant,get_sExp_tauFit())
    << "\n\n";
    
    build_penalizedspline_interpolant(correlation_sortincreasing,0,1);//NOTE: F vs time
    outputFile
    << "time  correlation_fit            " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<correlation_sortincreasing.size(); ++i) {
        outputFile
        << correlation_sortincreasing.at(i).at(0) << " "
        << spline1dcalc(interpolant,correlation_sortincreasing.at(i).at(0))
        << "\n";
    } outputFile << "\n\n";
    
    outputFile
    << "time  correlation(time)          " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<correlation_original.size(); ++i) {
        outputFile
        << correlation_original.at(i).at(0) << " "
        << correlation_original.at(i).at(1) << "\n";
    } interpolant=org;
}





void FitData::write_sExp_params(std::ofstream& outputFile,
                                const double Temp_d)
{
    //if (outputFile.tellp()==0) {
    //    outputFile
    //    << "T A tau beta r2" << "\n";
    //}
    outputFile << Temp_d*pow(corrFacforT,-1) << " ";
    for (int i=0; i<c.length(); ++i) outputFile << c[i] << " ";
    outputFile << rep.r2 << "\n";
}





void FitData::write_sExp_params_avg(const StructureClass& sysVar,
                                    const int n_sys,
                                    const std::string& tcalc)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/sExp_params_avg_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_params_avg(o.c_str());
    
    read_all_sExp_params(sysVar,n_sys);
    
    if (tcalc=="KWW"||tcalc=="KWW_lnFs") {
        write_params_avg << "T A tau beta r2" << "\n";
    }
    else if (tcalc=="KWW_full"||tcalc=="KWW_full_lnFs") {
        write_params_avg << "T A tau beta tauf p r2" << "\n";
    }
    else if (tcalc=="KWW_pwr") {
        write_params_avg << "T A a tau beta r2" << "\n";
    }
    else if (tcalc=="mKWW") {
        write_params_avg << "T A taulib tau beta r2" << "\n";
    }
    else if (tcalc=="mKWW_pwr") {
        write_params_avg << "T A a taulib tau beta r2" << "\n";
    }
    
    /** smooth beta vs 1/T by penalized spline **/
    vector<vector<double>> invT=param_sortdecreasing;
    for (int i=0;i<(int)invT.size();++i) {
        invT.at(i).at(0)=pow(invT.at(i).at(0),-1);
    }
    alglib::spline1dinterpolant inp0,inp1,inph,inp2;
    build_penalizedspline_interpolant(invT,0,3,50,0.0);
    inp0=interpolant;
    build_penalizedspline_interpolant(invT,0,3,50,0.5);
    inph=interpolant;
    build_penalizedspline_interpolant(invT,0,3,50,1.0);
    inp1=interpolant;
    build_penalizedspline_interpolant(invT,0,3,50,2.0);
    inp2=interpolant;
    
    for (size_t i=0; i<param_sortdecreasing.size(); ++i) {
        for (size_t ii=0; ii<param_sortdecreasing.at(0).size(); ++ii) {
            write_params_avg << param_sortdecreasing.at(i).at(ii) << " ";
        }
        double iT=invT.at(i).at(0);
        write_params_avg
        << iT                    << " " //1/T
        << spline1dcalc(inp0,iT) << " " //beta_rho0
        << spline1dcalc(inph,iT) << " " //beta_rho0.5
        << spline1dcalc(inp1,iT) << " " //beta_rho1.0
        << spline1dcalc(inp2,iT) << " " //beta_rho2.0
        << "\n";
    }
}





void FitData::write_tauFitfit_vs_T(std::ofstream& outputFile,
                                   const alglib::real_1d_array& c,
                                   const double T_initial,
                                   const std::string& model)
{
    //double T_initial=Theq;
    double Ti=0;
    
    /* format: T log10(tauFit_fit) */
    //------------------------------------------------------------------
    if (model=="Arrhenius") {
        
        double res=res_Arrhenius/(corrFacforT/precision);//T resolution
        
        outputFile << "\n"
        << "T  log10(tauFit_fit)             " << "\n"
        << "=================================" << "\n";
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-i*res;
            
            outputFile
            << Ti << " "
            << log10(calc_y_given_x(c,Ti,model))
            << "\n";
        }
        
    } else {
        
        double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
        
        outputFile << "\n"
        << "T  log10(tauFit_fit)             " << "\n"
        << "=================================" << "\n";
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
            
            outputFile
            << Ti << " "
            << log10(calc_y_given_x(c,Ti,model))
            << "\n";
        }
    } outputFile << "\n\n";
}





void FitData::write_tauFitfit_vs_invT(std::ofstream& outputFile,
                                      const alglib::real_1d_array& c,
                                      const double T_initial,
                                      const std::string& model)
{
    //double T_initial=Theq;
    double Ti=0;
    
    /* format: 1/T log10(tauFit_fit) */
    //------------------------------------------------------------------
    if (model=="Arrhenius") {
        
        double res=res_Arrhenius/(corrFacforT/precision); // T resolution
        
        outputFile << "\n"
        << 1000.0/(corrFacforT/precision)
        << "/T  log10(tauFit_fit)            " << "\n"
        << "=================================" << "\n";
        
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-i*res;
            
            outputFile
            << convert/Ti << " "
            << log10(calc_y_given_x(c,Ti,model))
            << "\n";
        }
        
    } else {
        
        double Tg_extrp = calc_x_given_y(c,extrp_time,model)[0];
        
        outputFile << "\n"
        << 1000.0/(corrFacforT/precision)
        << "/T  log10(tauFit_fit)            " << "\n"
        << "=================================" << "\n";
        
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
            
            outputFile
            << convert/Ti << " "
            << log10(calc_y_given_x(c,Ti,model))
            << "\n";
        }
    } outputFile << "\n\n";
}





void FitData::write_tauFitfit_vs_dwf(std::ofstream& outputFile,
                                     const alglib::real_1d_array& c,
                                     const double dwf_initial,
                                     const std::string& model)
{
    //
}





void FitData::write_tauFitfit_vs_reduceddwf(std::ofstream& outputFile,
                                            const alglib::real_1d_array& c,
                                            const vector<vector<double>>& standard2d,
                                            const std::vector<double>& ref,
                                            const double Y_initial,
                                            const std::string& model)
{
    // standard2d.at(i): {logtaur,u2r,u2,T,tau}
    
    double Yi=0;
    double Yg_extrp=calc_Y_given_time(c,extrp_time,ref,model).at(0);
    double Xr=pow(Yg_extrp,-1);
    double logtaur=ref.at(0);
    double u2r=ref.at(1);
    double Tr=ref.at(2);
    double u2=0;
    
    vector<double> v;
    
    outputFile << "\n"
    << "Tr      " << Tr      << "\n"
    << "logtaur " << logtaur << "\n"
    << "u2r     " << u2r     << "\n\n";
    
    if (model=="leporiniref"||model=="leporini_universal") {
        outputFile
        << "Xr(u2g/u2r)       " << Xr  << "\n"
        << "Xr(A,B,C,logtaur) " << get_Xr_leporini(c,logtaur) << "\n\n";
        Xr=get_Xr_leporini(c,logtaur);
    } else {
        outputFile
        << "Xr(u2g/u2r) " << Xr << "\n\n";
    }
    
    outputFile
    << "T  u2  Xr  Y  log10(tauFit)" << "\n";
    for (size_t i=0; i<standard2d.size(); ++i) {
        outputFile
        << standard2d.at(i).at(3) << " "
        << standard2d.at(i).at(2) << " "
        << Xr << " "
        << standard2d.at(i).at(1)/standard2d.at(i).at(2) << " "
        << standard2d.at(i).at(4) << "\n";
    } outputFile << "\n\n";
    
    outputFile
    << "u2  Xr  Y  log10(tauFit_fit)" << "\n";
    for (int i=1; i<=equal_pieces; ++i) {
        
        Yi=Y_initial+(fabs(Yg_extrp-Y_initial)/(double)equal_pieces)*i;
        u2=u2r*pow(Yi,-1);
        v.clear(); v={logtaur,u2r,u2};
        
        outputFile
        << u2 << " "
        << Xr << " "
        << Yi << " "
        << log10(calc_time_given_vecinput(c,v,model))
        << "\n";
    }
}





void FitData::write_qchTs(StructureClass& sysVar,
                          const int n_trl,
                          const int n_sys)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tempsInfo_qchTs_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    ofstream write_qchTs(o.c_str(),ofstream::app); // use "append"
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to the "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero and with fixed format
    //--------------------------------------------------------------------------
    write_qchTs << fixed << setprecision(0);
    write_qchTs << "{";
    for (size_t i=0; i<sysVar.get_quenchingTs().at(index).size(); ++i)
    {
        if (i!=(sysVar.get_quenchingTs().at(index).size()-1)) {
            write_qchTs
            << sysVar.get_quenchingTs().at(index).at(i) << ",";
        } else {
            write_qchTs
            << sysVar.get_quenchingTs().at(index).at(i) << "}";
        }
    } write_qchTs << "\n";
    write_qchTs.close();
}





void FitData::write_equTs(StructureClass& sysVar,
                          const int n_trl,
                          const int n_sys)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tempsInfo_equTs_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    ofstream write_equTs(o.c_str(),ofstream::app); // use "append"
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to the "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero and with fixed format
    //--------------------------------------------------------------------------
    write_equTs << fixed << setprecision(0);
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_taueq(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_taueqinv(o.c_str(),ofstream::app); // 'append'
    
    if (is_imposeStrictTeq)
    {
        /** NOTE: impose strict monotonic increasing rule on taueq **/
        //------------------------------------------------------------------
        vector<double> inequTs;
        vector<vector<double>> tauEqus;
        /** Find first Toe **/
        double Toe_max=0;
        for (size_t i=0; i<sysVar.get_tauFit().at(index).size(); ++i) {
            if (sysVar.get_tauFit().at(index).at(i).at(1)>log10(tauFit_cut)) {
                Toe_max=sysVar.get_tauFit().at(index).at(i).at(0);
                break;
            }
        } /** NOTE: if all T's are in-equilibrium, Toe_max=0 **/
        
        int counter=0;
        double tau_big=0;
        for (size_t i=0; i<sysVar.get_tauEqu().at(index).size(); ++i)
        {
            double Ti  =sysVar.get_tauEqu().at(index).at(i).at(0);
            double taui=sysVar.get_tauEqu().at(index).at(i).at(1);
            if (counter==0) {
                if (Ti>Toe_max) {
                    // first T (of current regime) is g.t. Toe
                    inequTs.push_back(Ti);
                    tauEqus.push_back({Ti,taui});
                    tau_big=taui;
                    ++counter;
                }
            } else {
                if ((Ti>Toe_max)&&(taui>tau_big)) {
                    // tau_alpha shuold be monotonic increasing
                    inequTs.push_back(Ti);
                    tauEqus.push_back({Ti,taui});
                    tau_big=taui;
                }
            }
        }
        /** replace ie T of current regime by new ie T's **/
        sysVar.get_equilibratedTs().at(index).clear();
        sysVar.get_equilibratedTs().at(index)=inequTs;
        sysVar.get_tauEqu().at(index).clear();
        sysVar.get_tauEqu().at(index)=tauEqus;
    }
    
    /** store in-equilibrium T's into the contianer **/
    //----------------------------------------------------------------------
    write_equTs << "{";
    for (size_t i=0; i<sysVar.get_equilibratedTs().at(index).size(); ++i)
    {
        if (i!=(sysVar.get_equilibratedTs().at(index).size()-1)) {
            write_equTs
            << sysVar.get_equilibratedTs().at(index).at(i) << ",";
        } else {
            write_equTs
            << sysVar.get_equilibratedTs().at(index).at(i) << "}";
        }
    } write_equTs << "\n";
    
    /** write tau_eq data to file (NOTE: setprecision()) **/
    //----------------------------------------------------------------------
    for (size_t i=0; i<sysVar.get_tauEqu().at(index).size(); ++i)
    {
        double T_actual=sysVar.get_tauEqu().at(index).at(i).at(0)*pow(corrFacforT,-1);
        double log_tauFit=sysVar.get_tauEqu().at(index).at(i).at(1);
        
        write_taueq
        << setprecision(10)
        << T_actual << " " << log_tauFit << "\n";
        
        write_taueqinv
        << setprecision(10)
        << convert/T_actual << " "
        << log_tauFit
        << "\n";
    }
    write_equTs.close();
    write_taueq.close();
    write_taueqinv.close();
}





void FitData::write_fitavg(StructureClass& sysVar,
                           const int n_sys,
                           const std::string& model)
{
    real_1d_array c_avg=c;
    
    for (size_t i=0; i<c_avgVec.size(); ++i) {
        c_avg[i] = c_avgVec[i];
    }
    
    string o;
    
    if (model=="Arrhenius")
    {
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_Arrhenius_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_fit_Arrhenius_avg(o.c_str());
        
        write_fit_Arrhenius_avg
        << "Data fit to <" << model << "> functional form \n\n";
        
        write_fit_Arrhenius_avg << "avg fit coeffs" << "\n";
        for (size_t i=0; i<c_avg.length(); ++i) {
            write_fit_Arrhenius_avg << c_avg[i] << " ";
        } write_fit_Arrhenius_avg << "\n\n";
        
        sysVar.set_calcvector(r2_Arrhenius);
        
        write_fit_Arrhenius_avg
        << "Coefficient of Determination(R2)        " << "\n"
        << "----------------------------------------" << "\n"
        << "avg(R2)       " << sysVar.calc_mean() << "\n"
        << "stdev(R2)     " << sysVar.get_sample_stdev() << "\n"
        << "rel_stdev(R2) " << sysVar.get_sample_stdev()/sysVar.calc_mean()
        << "\n\n";
        
        for (size_t i=0; i<r2_Arrhenius.size(); ++i) {
            write_fit_Arrhenius_avg << r2_Arrhenius[i] << "\n";
        } write_fit_Arrhenius_avg << "\n\n";
        
        write_tauFitfit_vs_T(write_fit_Arrhenius_avg,c_avg,Theq,"Arrhenius");
        write_tauFitfit_vs_invT(write_fit_Arrhenius_avg,c_avg,Theq,"Arrhenius");
        
        write_fit_Arrhenius_avg.close();
        
    } else {
        
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_taueq_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append("_Regime"+to_string((long long int)current_regime));
        o.append(".dat");
        ofstream write_fit_tauFit_avg(o.c_str());
        
        write_fit_tauFit_avg
        << "Data fit to <" << model << "> functional form \n\n"
        << "(tauFit defined @KWW=" << get_sExp_tauFit() << ") \n\n";
        
        write_Tg_fragility(write_fit_tauFit_avg,c_avg,model);
        write_fitcutoff(write_fit_tauFit_avg,sysVar,model);
        
        write_fit_tauFit_avg << "avg fit coeffs" << "\n";
        for (size_t i=0; i<c_avg.length(); ++i) {
            write_fit_tauFit_avg << c_avg[i] << " ";
        } write_fit_tauFit_avg << "\n\n";
        
        sysVar.set_calcvector(r2_Model);
        
        write_fit_tauFit_avg
        << "Coefficient of Determination(R2)        " << "\n"
        << "----------------------------------------" << "\n"
        << "avg(R2)       " << sysVar.calc_mean() << "\n"
        << "stdev(R2)     " << sysVar.get_sample_stdev() << "\n"
        << "rel_stdev(R2) " << sysVar.get_sample_stdev()/sysVar.calc_mean()
        << "\n\n";
        
        for (size_t i=0; i<r2_Model.size(); ++i) {
            write_fit_tauFit_avg << r2_Model[i] << "\n";
        } write_fit_tauFit_avg << "\n\n";
        
        write_tauFitfit_vs_T(write_fit_tauFit_avg,c_avg,Theq,model);
        write_tauFitfit_vs_invT(write_fit_tauFit_avg,c_avg,Theq,model);
        
        write_fit_tauFit_avg.close();
    }
}





void FitData::avgfitParams(const StructureClass& sysVar,
                           const int n_trl,
                           const int n_sys,
                           const std::string& model)
{
    vector<vector<double>> coeffs;
    if (model=="Arrhenius") coeffs=ArrheniusCoeffs;
    else coeffs=ModelCoeffs;
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    vector<double> avgParams;
    for (size_t i=0; i<coeffs.at(index).size(); ++i) {
        avgParams.push_back(0);
    }
    
    int total_count=0; // total counts of all trials & systems
    for (int i=0; i<=n_trl; ++i) { // up to current trial
        for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) { // up to current sys
            
            ++total_count;
            
            int index_i=i*n_system+ii;
            
            for (size_t iii=0; iii<coeffs.at(index_i).size(); ++iii) {
                avgParams.at(iii) += coeffs.at(index_i).at(iii);
            }
        }
    }
    for (size_t i=0; i<coeffs.at(index).size(); ++i) {
        avgParams.at(i) /= (double)(total_count);
    }
    
    c_avgVec.clear(); // NOTE: vector renewed every time!
    
    for (size_t i=0; i<avgParams.size(); ++i) {
        c_avgVec.push_back(avgParams.at(i));
    }
}





void FitData::avgfitParams_extrp(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const std::string& extrp)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    vector<double> avgParams;
    for (size_t i=0; i<ExtrpCoeffs.at(index).size(); ++i) {
        avgParams.push_back(0);
    }
    //--------------------------------------------------------------------------
    // from second regime on, if stdev/avg < 0.1,
    // use avg fitParams to shoot new temperatures
    //--------------------------------------------------------------------------
    if ((current_regime>0)&&(relativeStdev<0.1)) {
        int total_count=0; // total counts of all trials & systems
        for (int i=0; i<=n_trl; ++i) { // up to current trial
            for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) { // up to current sys
                
                ++total_count;
                
                int index_i=i*n_system+ii;
                for (size_t iii=0; iii<ExtrpCoeffs.at(index_i).size(); ++iii) {
                    avgParams.at(iii) += ExtrpCoeffs.at(index_i).at(iii);
                }
            }
        }
        for (size_t i=0; i<ExtrpCoeffs.at(index).size(); ++i) {
            avgParams.at(i) /= (double)(total_count);
        }
    }
    //--------------------------------------------------------------------------
    // in first regime || relstdev >= 0.1 in higher regimes,
    // use the fitParams of the highest fragility system to shoot
    //--------------------------------------------------------------------------
    else {
        int ilm=get_index_largest_m();
        for (size_t i=0; i<ExtrpCoeffs.at(ilm).size(); ++i) {
            avgParams.at(i) = ExtrpCoeffs.at(ilm).at(i);
        }
    }
    c_avgVec.clear(); // NOTE: vector renewed every time!
    for (size_t i=0; i<avgParams.size(); ++i) {
        c_avgVec.push_back(avgParams.at(i));
    }
}





void FitData::find_cutTforArrhenius(const StructureClass& sysVar,
                                    const int n_sys)
{
    string model="Arrhenius";
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys,"continuous");
    taueq_sortincreasing=taueq_sortdecreasing;
    std::sort(taueq_sortincreasing.begin(),taueq_sortincreasing.end(),sortIncreasing0);
    fitdata_processing_lsfit(taueq_sortincreasing);
    int n_data=(int)taueq_sortincreasing.size();
    int n_iteration=n_data-n_fit_chopping;
    for (int i=0; i<n_iteration; i++) {
        alglib_lsfit(model);
        string tem = xdata.at(xdata.size()-2);
        //string tau = ydata.at(ydata.size()-2);
        double tem_double=stod(tem);
        //double tau_double=stod(tau);
        if (rep.r2>0.99) set_cutTforArrhenius(tem_double);
        fit_every_postprocessing();
    }
}





void FitData::find_TA(const StructureClass& sysVar,
                      const int n_sys,
                      const vector<vector<double>> tauVec)
{
    if (is_fixed_tauA)
    {
        cout
        << "in FitData::find_TA():\n"
        << "is_fixed_tauA=true, this functionality is disabled.\n";
        exit(EXIT_FAILURE);
        //read_all_taueq_data(sysVar,n_sys);
        /*
         vector<vector<double>> Tvstau;
         for (indexi=0;indexi<taueq_sortdecreasing.size();++indexi) {
         double T=taueq_sortdecreasing.at(indexi).at(0);
         double tau=taueq_sortdecreasing.at(indexi).at(1);
         Tvstau.push_back({tau,T});
         }
         */
        //build_penalizedspline_interpolant(taueq_sortdecreasing,0,1,50,2.0);//(time,T)
        //tauA=pow(10,tauA_fixed);
        //TA=spline1dcalc(interpolant,tauA_fixed);
        
    } else {
        
        vector<vector<double>> taueq_org=tauVec;
        vector<string> xdata_org=xdata;
        vector<string> ydata_org=ydata;
        
        /* fit highT tau(T) to Arrhenius */
        //----------------------------------------------------------------------
        set_fitParams("Arrhenius");
        read_all_taueq_data(sysVar,n_sys,"Arrhenius");
        fitdata_processing_lsfit(taueq_sortdecreasing);//NOTE: T >= cutoffT
        alglib_lsfit("Arrhenius");
        tau0_Arrhenius=c[0];
        Ea=c[1];
        //----------------------------------------------------------------------
        taueq_sortdecreasing=tauVec;
        /* use all taueq data to find TA
         * cf. read taueq with Arrhenius condition would only take relaxation times
         * whose T >= cutoff T */
        
        if (definitionofTA=="ArrFit") {
            use_ArrFit();
        } else if (definitionofTA=="ArrNor") {
            use_ArrNor();
        } else if (definitionofTA=="MsdCag") {
            use_MsdCag();
        } else if (definitionofTA=="FsqCag") {
            //use_FsqCag();
        }
        TA=TA_avg;
        tauA=tauA_avg;
        taueq_sortdecreasing=taueq_org;
        xdata=xdata_org;
        ydata=ydata_org;
    }
}





void FitData::find_u2A(const StructureClass& sysVar,
                       const int n_sys,
                       const double TA_d)
{
    /* interpolate <u2(T)> to get <u2(TA)> */
    double T_pre=0,T_pos=0;
    double DWF_pre=0,DWF_pos=0;
    double slope=0;
    //double Tr=0,u2r=0;
    
    string fit="pSpline";//COOP_DWF,pSpline
    
    read_all_equ_DWF(sysVar,n_sys);//NOTE: use "equ" data
    
    if (true)
    {
        if (fit=="COOP_DWF") {
            set_fitParams("COOP_DWF");
            fitdata_processing_lsfit(dwf_invT_sortdecreasing);
            alglib_lsfit("COOP_DWF");
            u2A=log10(calc_y_given_x(c,convert/TA_d,"COOP_DWF"));
            /* u2* modification */
            //Tr=xTA*TA_d;
            //u2r=log10(calc_y_given_x(c,convert/Tr,"COOP_DWF"));
            //u2A=u2r/(1.5*xTA-0.5);
        } else if (fit=="pSpline") {
            //cout
            //<< "in FitData::find_u2A():\n"
            //<< "currently not supporting using ("<<fit<<") fit for u2A.\n";
            //exit(EXIT_FAILURE);
            build_penalizedspline_interpolant(dwf_invT_sortdecreasing,0,1,50,rho_dwf);
            u2A=spline1dcalc(interpolant,convert/TA_d);
            /* u2* modification */
            //Tr=xTA*TA_d;
            //u2r=spline1dcalc(interpolant,convert/Tr);
            //u2A=u2r/(1.5*xTA-0.5);
        } else {
            cout
            << "in FitData::find_u2A(): "
            << "fitting model ("<<fit<<") not found.\n";
            exit(EXIT_FAILURE);
        }
    } else { //OLD: linear interpolation
        
        for (int i=0; i<(int)(dwf_sortdecreasing.size()-1); ++i) {
            
            T_pre = dwf_sortdecreasing[i][0];
            T_pos = dwf_sortdecreasing[i+1][0];
            
            if (((TA_d-T_pre)*(TA_d-T_pos))<0)
            {
                DWF_pre=dwf_sortdecreasing[i][1];
                DWF_pos=dwf_sortdecreasing[i+1][1];
                slope = (DWF_pos-DWF_pre)/(T_pos-T_pre);
                u2A   = DWF_pre+slope*(TA_d-T_pre);
                return;
            }
        }
    }
}





void FitData::find_rhoA(const StructureClass& sysVar,
                        const int n_sys,
                        const double TA_d)
{
    string model_d="pSpline";
    
    read_all_thermo_data(sysVar,n_sys,"Density");
    build_penalizedspline_interpolant(thermo_sortdecreasing,0,1);
    rhoA=spline1dcalc(interpolant,TA_d);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_rho_spline_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    fitdata_processing_lsfit(thermo_sortdecreasing);
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_spline_actual(writefile,model_d,thermo_sortdecreasing);
    writefile.close();
}





void FitData::find_balltocage_time(const StructureClass& sysVar,
                                   const int n_trl,
                                   const int n_sys,
                                   const double Temp_d)
{
    string which_method="isf";//msd,isf
    
    double T_actual=Temp_d*pow(corrFacforT,-1);
    double logtime=0,f=0,df=0,d2f=0;
    double logcutofftime=1;//cutoff time in ps
    if (sysVar.get_systemUnit()=="real") {
        logcutofftime+=3.0;
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/dMSD_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //ofstream writeFile_dMSD(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/dISF_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //ofstream writeFile_dISF(o.c_str(),ofstream::app);//'append'
    
    /** Derivative Analysis for MSD **/
    //--------------------------------------------------------------------------
    if (T_actual<TA)
    {
        read_raw_msd_data(sysVar,n_trl,n_sys,Temp_d);
        vector<vector<double>> r2cont=semilogmsd_sortincreasing;//loglogmsd_sortincreasing
        
        vector<vector<double>> cutmsd;//increase-time-sorted msd
        for (size_t i=0;i<r2cont.size();++i) {
            logtime=r2cont.at(i).at(0);
            f      =r2cont.at(i).at(1);
            if (logtime<logcutofftime) cutmsd.push_back({logtime,f});
        }
        build_penalizedspline_interpolant(cutmsd,0,1);
        alglib::spline1dinterpolant interpolant_r2=interpolant;
        vector<vector<double>> dmsd;
        for (size_t i=0;i<cutmsd.size();++i) {
            logtime=cutmsd.at(i).at(0);
            spline1ddiff(interpolant_r2,logtime,f,df,d2f);
            dmsd.push_back({d2f,logtime,f,df});//NOTE: 1d=d2f
        } std::sort(dmsd.begin(),dmsd.end(),sortIncreasing0);//NOTE: sort increasing
        
        //writeFile_dMSD
        //<< T_actual         << " "  //T
        //<< dmsd.at(0).at(1) << " "  //log.time
        //<< dmsd.at(0).at(2) << " "  //log.msd
        //<< dmsd.at(0).at(3) << " "  //dlog.msd
        //<< dmsd.at(0).at(0) << "\n";//d2log.msd
        
        //if (which_method=="msd") DWF_time=pow(10,dmsd.at(0).at(1));
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data/");
        o.append("/dMSDT_");
        o.append(sysVar.get_usic());
        //o.append("_00"+to_string((long long int)n_trl));
        //o.append("_"+sysVar.get_nameString(n_sys));
        //o.append("_T"+to_string((long long int)Temp_d));
        o.append(".dat");
        ofstream writemsdT(o.c_str(),ofstream::app);//'append'
        for (size_t i=0;i<cutmsd.size();++i) {
            logtime=cutmsd.at(i).at(0);
            spline1ddiff(interpolant_r2,logtime,f,df,d2f);
            if (sysVar.get_systemUnit()=="real") logtime-=3.0;
            writemsdT
            << T_actual<< " "
            << logtime << " "
            << f       << " "
            << df      << " "
            << d2f     << "\n";
        } writemsdT<<"\n";
    } //writemsdT.close();
    //--------------------------------------------------------------------------
    
    
    /** Derivative Analysis for ISF **/
    //--------------------------------------------------------------------------
    if (T_actual<TA)
    {
        read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d);
        vector<vector<double>> fscont=correlation_loglog;//correlation_semilog
        
        vector<vector<double>> cutisf;//increase-time-sorted isf
        for (size_t i=0;i<fscont.size();++i) {
            logtime=fscont.at(i).at(0);
            f      =fscont.at(i).at(1);
            if (logtime<logcutofftime) cutisf.push_back({logtime,f});
        }
        build_penalizedspline_interpolant(cutisf,0,1);
        alglib::spline1dinterpolant interpolant_fs=interpolant;
        vector<vector<double>> disf;
        for (size_t i=0;i<cutisf.size();++i) {
            logtime=cutisf.at(i).at(0);
            spline1ddiff(interpolant_fs,logtime,f,df,d2f);
            disf.push_back({d2f,logtime,f,df});//NOTE: 1d=d2f
        } std::sort(disf.begin(),disf.end(),sortDecreasing0);//NOTE: sort decreasing
        
        //writeFile_dISF
        //<< T_actual         << " "  //T
        //<< disf.at(0).at(1) << " "  //log.time
        //<< disf.at(0).at(2) << " "  //isf
        //<< disf.at(0).at(3) << " "  //disf
        //<< disf.at(0).at(0) << "\n";//d2isf
        
        //if (which_method=="isf") DWF_time=pow(10,disf.at(0).at(1));
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data/");
        o.append("/dISFT_");
        o.append(sysVar.get_usic());
        //o.append("_00"+to_string((long long int)n_trl));
        //o.append("_"+sysVar.get_nameString(n_sys));
        //o.append("_T"+to_string((long long int)Temp_d));
        o.append(".dat");
        ofstream writeisfT(o.c_str(),ofstream::app);//'append'
        for (size_t i=0;i<cutisf.size();++i) {
            logtime=cutisf.at(i).at(0);
            spline1ddiff(interpolant_fs,logtime,f,df,d2f);
            if (sysVar.get_systemUnit()=="real") logtime-=3.0;
            writeisfT
            << T_actual<< " "
            << logtime << " "
            << f       << " "
            << df      << " "
            << d2f     << "\n";
        } writeisfT<<"\n";
    } //writeisfT.close();
    //--------------------------------------------------------------------------
}





void FitData::stickelanalysis(const StructureClass& sysVar)
{
    /* penalized spline params */
    //const int    m=50;
    //const double rho=0;
    
    const int n_sys_beg=sysVar.get_n_sys_beg();
    //const int n_sys_end=sysVar.get_n_sys_end();
    //for (int n_sys=n_sys_beg;n_sys<=n_sys_end;++n_sys) {
    read_all_taueq_data(sysVar,n_sys_beg);
    //}
    vector<vector<double>> taueq=taueq_sortdecreasing;
    vector<vector<double>> dF,result;
    double T=0,tau=0,tauf=0,stickel=0,f=0,df=0,d2f=0;
    
    /** first Stickel derivative analysis **/
    //build_penalizedspline_interpolant(taueq,0,1,m,rho);
    build_cubicspline_interpolant(taueq,0,1);
    alglib::spline1dinterpolant interpdF=interpolant;
    for (size_t i=0;i<taueq.size();++i) {
        T=   taueq.at(i).at(0);
        spline1ddiff(interpdF,T,f,df,d2f);
        stickel=pow(-df,-0.5);//(-dLog.tau/dT)^(-1/2)
        dF.push_back({T,stickel});
    }
    
    /** second Stickel derivative analysis **/
    //build_penalizedspline_interpolant(dF,0,1,m,rho);
    build_cubicspline_interpolant(dF,0,1);
    alglib::spline1dinterpolant interpd2F=interpolant;
    for (size_t i=0;i<dF.size();++i) {
        T=   taueq.at(i).at(0);
        tau= taueq.at(i).at(1);
        tauf=spline1dcalc(interpdF,T);
        spline1ddiff(interpd2F,T,f,df,d2f);
        result.push_back({1/T,T,tau,tauf,dF[i][1],df});
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stickelAnalysis_");
    o.append(sysVar.get_usic());
    o.append(".dat");
    ofstream writeFile(o.c_str());
    writeFile << "1/T T tau tau_pslpine dF d2F \n";
    for (size_t i=0;i<result.size();++i) {
        writeFile
        << result.at(i).at(0) << " "
        << result.at(i).at(1) << " "
        << result.at(i).at(2) << " "
        << result.at(i).at(3) << " "
        << result.at(i).at(4) << " "
        << result.at(i).at(5) << "\n";
    } writeFile.close();
}





void FitData::use_ArrFit()
{
    deviationscale="raw";
    
    double tauA_calc=0;
    double Ti=0,T_hi=0,T_lo=0;
    double tau_hi=0,tau_lo=0;
    double slope=0,dT=0;
    int equal_parts=1000; // discretize by 1000 equal parts
    int whichT=0;
    
    if (deviationscale=="log") { // log10(tau_alpha) scale
        
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            whichT         = indexi;
            Ti             = taueq_sortdecreasing.at(indexi).at(0);
            tauA_calc      = taueq_sortdecreasing.at(indexi).at(1);        //log10
            tauA_Arrhenius = log10(calc_y_given_x(c,Ti,"Arrhenius"));//log10
            tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
            if (tauA_deviation>deviate_ArrFit) break;
        }
        /** linear interpolation **/
        if (whichT>0) {
            T_hi   = taueq_sortdecreasing.at(whichT-1).at(0);
            tau_hi = taueq_sortdecreasing.at(whichT-1).at(1);//log10
            T_lo   = taueq_sortdecreasing.at(whichT).at(0);
            tau_lo = taueq_sortdecreasing.at(whichT).at(1);  //log10
            slope  = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT     = fabs(T_hi-T_lo)/(double)equal_parts;
            /** Interpolation on log10 scale **/
            for (int i=0; i<=equal_parts; ++i) {
                Ti             = T_hi - i*dT;
                tauA_calc      = tau_hi+slope*(Ti-T_hi);
                tauA_Arrhenius = log10(calc_y_given_x(c,Ti,"Arrhenius"));//log10
                tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
                if (tauA_deviation>deviate_ArrFit) break;
            }
        }
        TA_avg         = Ti;                // raw temperature
        tauA_avg       = pow(10,tauA_calc); // raw tau
        tauA_Arrhenius = calc_y_given_x(c,Ti,"Arrhenius");
        
    } else if (deviationscale=="ln") { // ln(tau_alpha) scale
        
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            whichT         = indexi;
            Ti             = taueq_sortdecreasing.at(indexi).at(0);
            tauA_calc      = log(10)*taueq_sortdecreasing.at(indexi).at(1);//ln
            tauA_Arrhenius = log(calc_y_given_x(c,Ti,"Arrhenius"));  //ln
            tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
            if (tauA_deviation>deviate_ArrFit) break;
        }
        /** linear interpolation **/
        if (whichT>0) {
            T_hi   = taueq_sortdecreasing.at(whichT-1).at(0);
            tau_hi = log(10)*taueq_sortdecreasing.at(whichT-1).at(1);//ln
            T_lo   = taueq_sortdecreasing.at(whichT).at(0);
            tau_lo = log(10)*taueq_sortdecreasing.at(whichT).at(1);  //ln
            slope  = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT     = fabs(T_hi-T_lo)/(double)equal_parts;
            /** Interpolation on loge scale **/
            for (int i=0; i<=equal_parts; ++i) {
                Ti             = T_hi - i*dT;
                tauA_calc      = tau_hi + slope*(Ti-T_hi);
                tauA_Arrhenius = log(calc_y_given_x(c,Ti,"Arrhenius"));//ln
                tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
                if (tauA_deviation>deviate_ArrFit) break;
            }
        }
        TA_avg         = Ti;                    // raw temperature
        tauA_avg       = pow(exp(1),tauA_calc); // raw tau
        tauA_Arrhenius = calc_y_given_x(c,Ti,"Arrhenius");
        
    } else if (deviationscale=="raw") { // raw tau_alpha scale
        
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            whichT         = indexi;
            Ti             = taueq_sortdecreasing.at(indexi).at(0);
            tauA_calc      = pow(10,taueq_sortdecreasing.at(indexi).at(1));
            tauA_Arrhenius = calc_y_given_x(c,Ti,"Arrhenius");
            tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
            if (tauA_deviation>deviate_ArrFit) break;
        }
        /** linear interpolation **/
        if (whichT>0) {
            T_hi   = taueq_sortdecreasing.at(whichT-1).at(0);
            tau_hi = pow(10,taueq_sortdecreasing.at(whichT-1).at(1));//raw
            T_lo   = taueq_sortdecreasing.at(whichT).at(0);
            tau_lo = pow(10,taueq_sortdecreasing.at(whichT).at(1));  //raw
            slope  = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT     = fabs(T_hi-T_lo)/(double)equal_parts;
            /** Interpolation on raw tau scale **/
            for (int i=0; i<=equal_parts; ++i) {
                Ti             = T_hi - i*dT;
                tauA_calc      = tau_hi + slope*(Ti-T_hi);
                tauA_Arrhenius = calc_y_given_x(c,Ti,"Arrhenius");//raw
                tauA_deviation = (tauA_calc-tauA_Arrhenius)/tauA_Arrhenius;
                if (tauA_deviation>deviate_ArrFit) break;
            }
        }
        TA_avg   = Ti;        // raw temperature
        tauA_avg = tauA_calc; // raw tau
    }
}





void FitData::use_ArrNor()
{
    double norm=0,tauA_calc=0;
    double Ti=0,Tinv=0,T_hi=0,T_lo=0;
    double tau_hi=0,tau_lo=0;
    double norm_hi=0,norm_lo=0;
    double slope_norm=0,slope_tau=0,dT=0;
    int equal_parts=5000;//discretization resolution
    int whichT=0;
    
    if (true)
    {
        vector<vector<double>> arrNorm;
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            norm  = taueq_sortdecreasing.at(indexi).at(0);
            norm *= log(pow(10,taueq_sortdecreasing.at(indexi).at(1))/tau0_Arrhenius);
            norm /= Ea;
            Ti=taueq_sortdecreasing.at(indexi).at(0);
            Tinv=convert/Ti;//NOTE: use invT
            arrNorm.push_back({Tinv,norm});//(invT,norm)
        }
        /* use pspline to interpolate normalized Arrhenius relaxation data */
        build_penalizedspline_interpolant(arrNorm,0,1,50,rho_arrNorm);//(invT,norm)
        alglib::spline1dinterpolant arrNorm_spline=interpolant;
        T_hi=arrNorm.at(0).at(0);//NOTE:invT
        T_lo=arrNorm.at(arrNorm.size()-1).at(0);//NOTE:invT
        dT=fabs(T_hi-T_lo)/equal_parts;
        /* use discretized interval to get close value to TA above threshold */
        for (size_t i=0;i<=equal_parts;++i) {
            Tinv=T_hi+i*dT;
            norm=spline1dcalc(interpolant,Tinv);
            if ((norm-1.0)>deviate_ArrNor) break;
        } Ti=convert/Tinv;
        
        string fit="pSpline";//COOP,Arrhenius,pSpline
        
        if (fit=="COOP") {
            alglib::real_1d_array c_org=c;
            //------------------------------------------------
            set_fitParams("COOP");
            //read_all_taueq_data(sysVar,n_sys);//already read
            fitdata_processing_lsfit(taueq_sortdecreasing);
            alglib_lsfit("COOP");
            tauA_calc=calc_y_given_x(c,Ti,"COOP");
            //------------------------------------------------
            c=c_org;
        } else if (fit=="Arrhenius") {
            tauA_calc=calc_y_given_x(c,Ti,"Arrhenius");
        } else if (fit=="pSpline") {
            //NOTE:
            //taueqinvT_sortdecreasing not read before entering use_ArrNor()
            //so {invT,taueq} data container needs to be created here
            vector<vector<double>> invTtaueq;
            double invT=0,taueq=0;
            for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
                invT=convert/taueq_sortdecreasing.at(i).at(0);
                taueq=taueq_sortdecreasing.at(i).at(1);
                invTtaueq.push_back({invT,taueq});
            }
            //build_penalizedspline_interpolant(taueq_sortdecreasing,0,1,50,rho_tau);
            build_penalizedspline_interpolant(invTtaueq,0,1,50,rho_tau);
            tauA_calc=pow(10,spline1dcalc(interpolant,convert/Ti));
        } else {
            cout
            << "in FitData::use_ArrNor(): "
            << "fitting model ("<<fit<<") not found.\n";
            exit(EXIT_FAILURE);
        } interpolant=arrNorm_spline;//NOTE:use original
        
    } else { //OLD: linear interpolation between two adjecent points
        
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            whichT= indexi;
            norm  = taueq_sortdecreasing.at(indexi).at(0);
            norm *= log(pow(10,taueq_sortdecreasing.at(indexi).at(1))/tau0_Arrhenius);
            norm /= Ea;
            if ((norm-1.0)>deviate_ArrNor) break;
        } cout << "\n";
        
        /** linear interpolation for more precise X% deviation **/
        if (whichT>0) {
            T_hi   = taueq_sortdecreasing.at(whichT-1).at(0);
            tau_hi = pow(10,taueq_sortdecreasing.at(whichT-1).at(1));//raw
            norm_hi= T_hi*log(tau_hi/tau0_Arrhenius)/Ea;
            T_lo   = taueq_sortdecreasing.at(whichT).at(0);
            tau_lo = pow(10,taueq_sortdecreasing.at(whichT).at(1));  //raw
            norm_lo= T_lo*log(tau_lo/tau0_Arrhenius)/Ea;
            
            slope_norm = (norm_lo-norm_hi)/(T_lo-T_hi);
            slope_tau  = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT         = fabs(T_hi-T_lo)/(double)equal_parts;
            
            /** Interpolation on log10 scale **/
            for (int i=0; i<=equal_parts; ++i) {
                Ti        = T_hi - i*dT;
                norm      = norm_hi+slope_norm*(Ti-T_hi);
                tauA_calc = tau_hi +slope_tau*(Ti-T_hi);
                if ((norm-1.0)>deviate_ArrNor) break;
            }
        }
    }
    cout << "Deviation (Normalized Arrhenius) " << (norm-1.0)*100 << "%\n";
    cout << "TA       " << Ti << "\n";
    cout << "Log.tauA " << log10(tauA_calc) << "\n";
    TA_avg         = Ti;       //raw TA
    tauA_avg       = tauA_calc;//raw tauA
    tauA_deviation = norm-1.0; //deviation from normalized Arrhenius behavior
    tauA_Arrhenius = calc_y_given_x(c,Ti,"Arrhenius");//tau_alpha by Arrhenius
}





void FitData::use_MsdCag()
{
    //TODO
}





void FitData::use_FsqCag()
{
    //TODO
}





void FitData::write_TA(const StructureClass& sysVar,
                       const int n_sys)
{
    const int current_regime = sysVar.get_current_regime();
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tempsInfo_TA_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    ofstream outfile(o.c_str());
    
    read_all_taueq_data(sysVar,n_sys);
    find_TA(sysVar,n_sys,taueq_sortdecreasing);
    
    outfile
    << "TA defined by <"<<definitionofTA<<"> method"
    << "\n"
    << "TA raw.tauA raw.tauA(Arrhenius) %deviation"
    << "\n"
    << "---------------------------------------------------"
    << "\n";
    
    outfile
    << TA_avg << " "
    << tauA_avg << " "
    << tauA_Arrhenius << " "
    << tauA_deviation*100 << "\n\n";
    
    if (definitionofTA=="ArrFit") {
        outfile
        << 1000.0/(corrFacforT/precision)
        << "/T  T  log10.tau_alpha  log10.tau_alpha(Arrhenius)" << "\n";
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            double Ti=taueq_sortdecreasing.at(indexi).at(0);
            outfile
            << convert/Ti   << " "
            << taueq_sortdecreasing.at(indexi).at(0) << " "
            << taueq_sortdecreasing.at(indexi).at(1) << " "
            << log10(calc_y_given_x(c,Ti,"Arrhenius")) << "\n";
        }
    } else if (definitionofTA=="ArrNor") {
        outfile
        << 1000.0/(corrFacforT/precision)
        << "/T  T  log10.tau_alpha  ArrNorm  spline(rho="<<rho_arrNorm<<")\n";
        for (indexi=0; indexi<(int)taueq_sortdecreasing.size(); ++indexi) {
            double Ti=taueq_sortdecreasing.at(indexi).at(0);
            double tau=pow(10,taueq_sortdecreasing.at(indexi).at(1));
            double Tinv=convert/Ti;
            outfile
            << convert/Ti   << " "
            << taueq_sortdecreasing.at(indexi).at(0) << " "
            << taueq_sortdecreasing.at(indexi).at(1) << " "
            << Ti*log(tau/tau0_Arrhenius)*pow(Ea,-1) << " "
            << spline1dcalc(interpolant,Tinv)        << "\n";
        }
    } outfile.close();
}





void FitData::find_tau0(StructureClass& sysVar,
                        const int n_sys,
                        const std::string& model)
{
    set_fitParams("Arrhenius");
    read_all_taueq_data(sysVar,n_sys,"Arrhenius");
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit("Arrhenius");
    tau0_Arrhenius=c[0];
    
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys);
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    tau0_model=c[0];
}





void FitData::write_tau0(StructureClass& sysVar,
                         const int n_sys,
                         const std::string& model)
{
    const int current_regime = sysVar.get_current_regime();
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tempsInfo_tau0_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau0(o.c_str());
    
    find_tau0(sysVar,n_sys,model);
    double deviation=fabs(tau0_model-tau0_Arrhenius)/fabs(tau0_Arrhenius);
    
    write_tau0
    << "tau0("<<model<<")  tau0(Arrhenius)  %deviation     \n"
    << "---------------------------------------------------\n";
    write_tau0
    << tau0_model << " "
    << tau0_Arrhenius << " "
    << deviation*100
    << "\n";
    write_tau0.close();
}





double FitData::write_Tc(const StructureClass& sysVar,
                         const int n_sys)
{
    const int current_regime = sysVar.get_current_regime();
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_Tc_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_Regime"+to_string((long long int)current_regime));
    o.append(".dat");
    ofstream write_findTc(o.c_str());
    
    /** fit MCT **/
    set_fitParams("MCT");
    read_all_taueq_data(sysVar,n_sys);
    fitdata_processing_lsfit(taueq_sortdecreasing);
    write_fitModel(write_findTc,"MCT");
    write_fitarrays(write_findTc);
    alglib_lsfit("MCT");
    write_stopcond(write_findTc);
    
    if (rep.r2>0.9) Tc_MCT = c[1]; // NOTE!
    
    write_Tg_fragility(write_findTc,c,"MCT");
    write_fitinfo(write_findTc);
    write_errcurve(write_findTc,"MCT");
    
    if (false) {
        /** fit SOU **/
        set_fitParams("SOU");
        read_all_taueq_data(sysVar,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        write_fitModel(write_findTc,"SOU");
        write_fitarrays(write_findTc);
        alglib_lsfit("SOU");
        write_stopcond(write_findTc);
        
        if (rep.r2>0.9) Tc_SOU = c[1]; // NOTE!
        
        write_Tg_fragility(write_findTc,c,"SOU");
        write_fitinfo(write_findTc);
        write_errcurve(write_findTc,"SOU");
    }
    return Tc_MCT;
}





const vector<double>& FitData::get_Theq(StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys)
{
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    return taueq_sortdecreasing[0];//sortDecreasing0
}

const vector<double>& FitData::get_Theq(StructureClass& sysVar,
                                        const int n_sys)
{
    read_all_taueq_data(sysVar,n_sys);
    return taueq_sortdecreasing[0];//sortDecreasing0
}




const vector<double>& FitData::get_Tleq(StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys)
{
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    return taueq_sortincreasing[0];//sortIncreasing0
}

const vector<double>& FitData::get_Tleq(StructureClass& sysVar,
                                        const int n_sys)
{
    read_all_taueq_data(sysVar,n_sys);
    return taueq_sortincreasing[0];//sortIncreasing0
}





void FitData::fit_every_preprocessing(const StructureClass& sysVar,
                                      const int n_sys,
                                      const std::string& cond)
{
    read_all_taueq_data(sysVar,n_sys,cond);
    fitdata_processing_lsfit(taueq_sortdecreasing);
}





void FitData::fit_every_postprocessing()
{
    fit_xData.clear();
    fit_yData.clear();
    
    taueq_sortdecreasing.pop_back(); //
    
    for (int i=0; i<3; ++i) {
        // pop  "]]" & "value" & "],[" appending the real_2d_array
        xdata.pop_back();
        // pop "]" & "value" & "," appending the real_1d_array
        ydata.pop_back();
    }
    xdata.push_back("]]");
    ydata.push_back("]");
    
    for (size_t i=0; i<xdata.size(); ++i) {
        fit_xData += xdata.at(i);
        fit_yData += ydata.at(i);
    }
    
    //cout << "\n" << fit_xData << "\n";
    //cout << "\n" << fit_yData << "\n";
}





void FitData::fit_every_processing(const int iteration,
                                   const int n_threshold)
{
    if ((iteration+n_threshold)>taueq_sortdecreasing.size()) {
        cout << "fit_every_processing: out of vector boundary!" << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    xdata.clear();     ydata.clear();
    fit_xData.clear(); fit_yData.clear();
    
    xdata.push_back("[[");
    ydata.push_back("[");
    
    for (int i=iteration; i<iteration+n_threshold; ++i) {
        
        double x = taueq_sortdecreasing.at(i)[0];
        double y = taueq_sortdecreasing.at(i)[1];
        
        // x_data to fit
        xdata.push_back(to_string((long double)x));
        xdata.push_back("],[");
        
        // y_data to fit
        ydata.push_back(to_string((long double)y));
        ydata.push_back(",");
        
    }
    xdata.pop_back(); // rm appended "],["
    ydata.pop_back(); // rm appended ","
    
    xdata.push_back("]]");
    ydata.push_back("]");
    for (size_t i=0; i<xdata.size(); ++i) {
        fit_xData += xdata.at(i);
        fit_yData += ydata.at(i);
    }
    //cout << "\n" << fit_xData << "\n";
    //cout << "\n" << fit_yData << "\n";
}





void FitData::fit_continuousRange(const StructureClass& sysVar,
                                  const int n_sys,
                                  const std::string& model)
{
    set_fitParams(model);
    
    int n_data=0;
    string cond="continuous";
    
    fit_every_preprocessing(sysVar,n_sys,cond);
    n_data=(int)taueq_sortdecreasing.size();
    
    if (n_data<n_fit_chopping) {
        cout
        << "FitData::fit_continuousRange: n_data < " << n_fit_chopping
        << " points. returned." << "\n";
        return;
    }
    
    int n_iterations=n_data-n_fit_chopping+1;
    
    /* Define paths to output files */
    //--------------------------------------------------------------------------
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full_Tg_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_Tg_fit_every(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full_m_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_m_fit_every(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full_param_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fit_every_param(o.c_str(),ofstream::app); // 'append'
    //--------------------------------------------------------------------------
    static int count_access=0;
    ++count_access;
    
    if (n_data<=n_fit_chopping) {
        
        if (count_access==1) {
            
            write_Tg_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
            write_m_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
            write_fit_every_param
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
        }
        
    } else {
        
        if (count_access==1) {
            
            write_Tg_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  Tg_extrp  error"
            << "\n";
            
            write_m_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  m_extrp error"
            << "\n";
            
            write_fit_every_param
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  fitParams  R^2"
            << "\n";
        }
        
        for (int i=0; i<n_iterations; ++i) {
            
            alglib_lsfit(model);
            
            if(i==0) {
                string o;
                o.append(return_AnalysisFolderPath(sysVar));
                o.append("/fit_data");
                o.append("/fit_tau_full_");
                o.append(sysVar.get_usic());
                o.append("_"+sysVar.get_nameString(n_sys));
                o.append(".dat");
                //cout << o << "\n";
                ofstream write_fit_taueq_avg(o.c_str());
                write_fitModel(write_fit_taueq_avg,model);
                write_fitarrays(write_fit_taueq_avg);
                write_stopcond(write_fit_taueq_avg);
                write_Tg_fragility(write_fit_taueq_avg,c,model);
                write_fitinfo(write_fit_taueq_avg);
                write_errcurve(write_fit_taueq_avg,model);
                write_tauFitfit_vs_T(write_fit_taueq_avg,c,Theq,model);
                write_tauFitfit_vs_invT(write_fit_taueq_avg,c,Theq,model);
                write_fit_taueq_avg.close();
            }
            
            double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
            double m_extrp =calc_x_given_y(c,extrp_time,model)[1];
            
            string tem = xdata.at(xdata.size()-2);
            string tau = ydata.at(ydata.size()-2);
            
            double tem_double=stod(tem);
            double tau_double=stod(tau);
            double logtau=tau_double;
            if (sysVar.get_systemUnit()=="real") {
                logtau=tau_double-3.0;
            }
            
            /*------ write parameters to file ------*/
            write_fit_every_param
            << tem_double << " "
            << logtau     << " ";
            for (int i=0; i<get_coeffs_vD().size(); ++i) {
                write_fit_every_param << c[i] << " ";
            } write_fit_every_param << rep.r2 << "\n";
            
            /*------ write Tg_extrp, m_extrp to file ------*/
            write_Tg_fit_every
            << tem_double << " "
            << logtau     << " "
            << Tg_extrp   << " "
            << error_Tg(model)
            << "\n";
            
            write_m_fit_every
            << tem_double << " "
            << logtau     << " "
            << m_extrp    << " "
            << error_fragility(model,Tg_extrp,error_Tg(model))
            << "\n";
            
            /* data processing */
            fit_every_postprocessing();
        }
    }
    write_Tg_fit_every.close();
    write_m_fit_every.close();
    write_fit_every_param.close();
}





void FitData::fit_neighboringNpoints(const StructureClass& sysVar,
                                     const int n_sys,
                                     const std::string& model)
{
    set_fitParams(model);
    
    int n_data=0;
    string cond="neighboring";
    
    fit_every_preprocessing(sysVar,n_sys,cond);
    
    n_data=(int)taueq_sortdecreasing.size();
    
    if (n_data<n_fit_sliding) {
        cout
        << "FitData::fit_neighboringNpoints: n_data < " << n_fit_sliding
        << " points. returned." << "\n";
        return;
    }
    
    int n_iterations=n_data-n_fit_sliding+1;
    
    /* Define paths to output files */
    //--------------------------------------------------------------------------
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full2_Tg_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_Tg_fit_every(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full2_m_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_m_fit_every(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_tau_full2_param_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fit_every_param(o.c_str(),ofstream::app); // 'append'
    //--------------------------------------------------------------------------
    static int count_access=0;
    ++count_access;
    
    if (n_data<n_fit_sliding) {
        
        if (count_access==1) {
            
            write_Tg_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_sliding
            << "\n";
            write_m_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_sliding
            << "\n";
            write_fit_every_param
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_sliding
            << "\n";
        }
        
    } else {
        
        if (count_access==1) {
            
            write_Tg_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  Tg_extrp  error"
            << "\n";
            
            write_m_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  m_extrp error"
            << "\n";
            
            write_fit_every_param
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  fitParams  R^2"
            << "\n";
        }
        
        for (int i=0; i<n_iterations; ++i)
        {
            fit_every_processing(i,n_fit_sliding);
            
            alglib_lsfit(model);
            
            double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
            double m_extrp =calc_x_given_y(c,extrp_time,model)[1];
            
            string tem = xdata.at(xdata.size()-2);
            string tau = ydata.at(ydata.size()-2);
            
            double tem_double=stod(tem);
            double tau_double=stod(tau);
            double logtau=tau_double;
            if (sysVar.get_systemUnit()=="real") {
                logtau=tau_double-3.0;
            }
            
            /*------ write parameters to file ------*/
            write_fit_every_param
            << tem_double << " "
            << logtau     << " ";
            for (size_t i=0; i<get_coeffs_vD().size(); ++i) {
                write_fit_every_param << c[i] << " ";
            } write_fit_every_param << rep.r2 << "\n";
            
            /*------ write Tg_extrp, m_extrp to file ------*/
            write_Tg_fit_every
            << tem_double << " "
            << logtau     << " "
            << Tg_extrp   << " "
            << error_Tg(model)
            << "\n";
            
            write_m_fit_every
            << tem_double << " "
            << logtau     << " "
            << m_extrp    << " "
            << error_fragility(model,Tg_extrp,error_Tg(model))
            << "\n";
        }
    }
    write_Tg_fit_every.close();
    write_m_fit_every.close();
    write_fit_every_param.close();
}





void FitData::fit_partial(const StructureClass& sysVar,
                          const int n_sys,
                          const std::string& model)
{
    set_fitParams(model);
    
    int n_data=0;
    string cond;
    
    for (int partial=0; partial<=1; ++partial) {
        
        if (partial==0) cond="partial_gt";
        else if (partial==1) cond="partial_lt";
        
        fit_every_preprocessing(sysVar,n_sys,cond);
        
        n_data=(int)taueq_sortdecreasing.size();
        
        if (n_data<n_fit_chopping) {
            cout
            << "fit_partial("
            << cond << ") "
            << "n_data < " << n_fit_chopping
            << " points. returned." << "\n";
            continue;
        }
        
        int n_iterations=n_data-n_fit_chopping+1;
        
        /* Define paths to output files */
        //----------------------------------------------------------------------
        string o1,o2,o3,o4;
        
        if (partial==1) {
            
            o1.append(return_AnalysisFolderPath(sysVar));
            o1.append("/fit_data");
            o1.append("/fit_tau_TltTx_Tg_");
            o1.append(sysVar.get_usic());
            o1.append("_"+sysVar.get_nameString(n_sys));
            o1.append(".dat");
            
            o2.append(return_AnalysisFolderPath(sysVar));
            o2.append("/fit_data");
            o2.append("/fit_tau_TltTx_m_");
            o2.append(sysVar.get_usic());
            o2.append("_"+sysVar.get_nameString(n_sys));
            o2.append(".dat");
            
            o3.append(return_AnalysisFolderPath(sysVar));
            o3.append("/fit_data");
            o3.append("/fit_tau_TltTx_param_");
            o3.append(sysVar.get_usic());
            o3.append("_"+sysVar.get_nameString(n_sys));
            o3.append(".dat");
            
        } else if (partial==0) {
            
            o1.append(return_AnalysisFolderPath(sysVar));
            o1.append("/fit_data");
            o1.append("/fit_tau_TgtTx_Tg_");
            o1.append(sysVar.get_usic());
            o1.append("_"+sysVar.get_nameString(n_sys));
            o1.append(".dat");
            
            o2.append(return_AnalysisFolderPath(sysVar));
            o2.append("/fit_data");
            o2.append("/fit_tau_TgtTx_m_");
            o2.append(sysVar.get_usic());
            o2.append("_"+sysVar.get_nameString(n_sys));
            o2.append(".dat");
            
            o3.append(return_AnalysisFolderPath(sysVar));
            o3.append("/fit_data");
            o3.append("/fit_tau_TgtTx_param_");
            o3.append(sysVar.get_usic());
            o3.append("_"+sysVar.get_nameString(n_sys));
            o3.append(".dat");
            
            o4.append(return_AnalysisFolderPath(sysVar));
            o4.append("/fit_data");
            o4.append("/fit_tau_TgtTx_Tc_");
            o4.append(sysVar.get_usic());
            o4.append("_"+sysVar.get_nameString(n_sys));
            o4.append(".dat");
        }
        ofstream write_Tg_fit_every(o1.c_str(),ofstream::app);    // 'append'
        ofstream write_m_fit_every(o2.c_str(),ofstream::app);     // 'append'
        ofstream write_fit_every_param(o3.c_str(),ofstream::app); // 'append'
        ofstream write_fit_every_Tc(o4.c_str(),ofstream::app);    // 'append'
        //----------------------------------------------------------------------
        
        if (n_data<=n_fit_chopping) {
            
            write_Tg_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
            write_m_fit_every
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
            write_fit_every_param
            << "n_data (=" << n_data << ")"
            << " for fitting is less than " << n_fit_chopping
            << "\n";
            
        } else {
            
            if (get_is_fit_by_Tc()) {
                write_Tg_fit_every
                << "Tx=Tc=" << get_Tc_MCT() << "\n";
                write_m_fit_every
                << "Tx=Tc=" << get_Tc_MCT() << "\n";
                write_fit_every_param
                << "Tx=Tc=" << get_Tc_MCT() << "\n";
            }
            else if (get_is_fit_by_TA()) {
                write_Tg_fit_every    << "Tx=TA=" << get_TA_avg() << "\n";
                write_m_fit_every     << "Tx=TA=" << get_TA_avg() << "\n";
                write_fit_every_param << "Tx=TA=" << get_TA_avg() << "\n";
            }
            
            write_Tg_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  Tg_extrp  error"
            << "\n";
            
            write_m_fit_every
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  m_extrp error"
            << "\n";
            
            write_fit_every_param
            << "Data fit to <" << model << "> functional form"
            << "\n"
            << "Tmin  tau_alpha  fitParams  R^2"
            << "\n";
            
            if (partial==0) {
                write_fit_every_Tc
                << "Data fit to <" << model << "> functional form"
                << "\n"
                << "Tmin  tau_alpha  Tc  TA  tauA  R2(Arr)"
                << "\n";
            }
            
            for (int i=0; i<n_iterations; ++i)
            {
                set_fitParams(model);
                alglib_lsfit(model);
                
                if(i==0) {
                    
                    if (partial==1) {
                        
                        o1.clear();
                        o1.append(return_AnalysisFolderPath(sysVar));
                        o1.append("/fit_data");
                        o1.append("/fit_tau_TltTx_");
                        o1.append(sysVar.get_usic());
                        o1.append("_"+sysVar.get_nameString(n_sys));
                        o1.append(".dat");
                    }
                    else if (partial==0) {
                        
                        o1.clear();
                        o1.append(return_AnalysisFolderPath(sysVar));
                        o1.append("/fit_data");
                        o1.append("/fit_tau_TgtTx_");
                        o1.append(sysVar.get_usic());
                        o1.append("_"+sysVar.get_nameString(n_sys));
                        o1.append(".dat");
                    }
                    ofstream write_fit_taueq_avg(o1.c_str());
                    
                    if (write_fit_taueq_avg.is_open()) {
                        
                        if (get_is_fit_by_Tc()) {
                            write_fit_taueq_avg
                            << "Tx=Tc=" << get_Tc_MCT() << "\n"
                            << "log(tau|Tx)=" <<
                            log10(calc_y_given_x(c,get_Tc_MCT(),model))
                            <<"\n";
                        }
                        else if (get_is_fit_by_TA()) {
                            write_fit_taueq_avg
                            << "Tx=TA=" << get_TA_avg() << "\n"
                            << "log(tau|Tx)=" <<
                            log10(calc_y_given_x(c,get_TA_avg(),model))
                            <<"\n";
                        }
                        
                        write_fitModel(write_fit_taueq_avg,model);
                        write_fitarrays(write_fit_taueq_avg);
                        write_stopcond(write_fit_taueq_avg);
                        write_Tg_fragility(write_fit_taueq_avg,c,model);
                        write_fitinfo(write_fit_taueq_avg);
                        write_errcurve(write_fit_taueq_avg,model);
                        write_tauFitfit_vs_T(write_fit_taueq_avg,c,Theq,model);
                        write_tauFitfit_vs_invT(write_fit_taueq_avg,c,Theq,model);
                        
                        write_fit_taueq_avg.close();
                    }
                }
                double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
                double m_extrp =calc_x_given_y(c,extrp_time,model)[1];
                
                string tem = xdata.at(xdata.size()-2);
                string tau = ydata.at(ydata.size()-2);
                
                double tem_double=stod(tem);
                double tau_double=stod(tau);
                double logtau=tau_double;
                if (sysVar.get_systemUnit()=="real") {
                    logtau=tau_double-3.0;
                }
                
                /*------ write parameters to file ------*/
                write_fit_every_param
                << tem_double << " "
                << logtau     << " ";
                for (size_t i=0; i<get_coeffs_vD().size(); ++i) {
                    write_fit_every_param << c[i] << " ";
                } write_fit_every_param << rep.r2 << "\n";
                
                /*------ write Tg_extrp, m_extrp to file ------*/
                write_Tg_fit_every
                << tem_double << " "
                << logtau     << " "
                << Tg_extrp   << " "
                << error_Tg(model)
                << "\n";
                
                write_m_fit_every
                << tem_double << " "
                << logtau     << " "
                << m_extrp    << " "
                << error_fragility(model,Tg_extrp,error_Tg(model))
                << "\n";
                
                /*------ Tc ------*/
                if (false)
                {
                    if (partial==0) {
                        if ((model=="MCT")||(model=="SOU")) {
                            
                            double Tc=c[1];
                            find_TA(sysVar,n_sys,taueq_sortdecreasing);
                            
                            write_fit_every_Tc
                            << tem_double << " "
                            << logtau     << " "
                            << Tc         << " "
                            << TA_avg     << " "<< tauA_avg <<" "<< rep.r2
                            << "\n";
                        }
                    }
                }
                
                /* data postprocessing */
                fit_every_postprocessing();
            }
        }
        write_Tg_fit_every.close();
        write_m_fit_every.close();
        write_fit_every_param.close();
    }
}





void FitData::shootForNewTemperatures(StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const real_1d_array& c_extrp,
                                      const string& extrp)
{
    /** write fit info to file (extrapolation model) **/
    //--------------------------------------------------------------------------
    string o;
    if (false) {
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_taueq_extrp_");
        o.append(sysVar.get_usic());
        o.append("_00"+to_string((long long int)n_trl));
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append("_Regime"+to_string((long long int)current_regime));
        o.append(".dat");
        ofstream write_extrp(o.c_str());
        write_fitModel(write_extrp,extrp);
        write_fitarrays(write_extrp);
        write_fitcutoff(write_extrp,sysVar,extrp);
        write_stopcond(write_extrp);
        write_Tg_fragility(write_extrp,c_extrp,extrp);
        write_fitinfo(write_extrp);
        write_errcurve(write_extrp,extrp);
        write_tauFitfit_vs_T(write_extrp,c_extrp,Theq,extrp);
        write_tauFit_vs_T(write_extrp);
        write_tauFitfit_vs_invT(write_extrp,c_extrp,Theq,extrp);
        write_tauFit_vs_invT(write_extrp);
        write_extrp.close();
    }
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tempsInfo_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tempsinfo(o.c_str(),ofstream::app); //append
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero
    //--------------------------------------------------------------------------
    write_tempsinfo << fixed << setprecision(0);
    
    bool is_partial_init=sysVar.get_max_initialized_regime()!=n_regime;
    
    if ((current_regime==0)||(!sysVar.get_is_updatetinfo())||
        (is_partial_init&&current_regime==sysVar.get_max_initialized_regime()))
    {
        /** write simulated temperatures to file **/
        for (size_t i=0; i<sysVar.get_temperaturesInfo().at(index).size(); ++i) {
            write_tempsinfo << sysVar.get_temperaturesInfo().at(index).at(i) << "\n";
        }
    }
    
    /**=========================================================================
     ** Temperature Extrapolation: n_regime-1 times, not including final regime
     ** number of extrapolated temperatures is given by ${n_next_temps}
     ** Extrapolated temperatures will not repeat any existing temperature
     ========================================================================**/
    if (current_regime<(n_regime-1))
    {
        double teql=0,taueql=0;     // lower  cutoff time (i.e. higher T)
        double teqh=0,taueqh=0;     // higher cutoff time (i.e. lower  T)
        double Tl=0,Th=0;           // T calculated at "taueql" & "taueqh"
        double Tl_calc=0,Th_calc=0;
        
        teql   = sysVar.get_equilibration_times().at(current_regime);
        teqh   = sysVar.get_equilibration_times().at(current_regime+1);
        taueql = teql/n_equ_blocks;
        taueqh = teqh/n_equ_blocks;
        
        //----------------------------------------------------------------------
        /** Th is the T at hihger tau_alpha cutoff   (ie. lower  T) **/
        /** Tl is the lowest simulated equilibrium T (ie. higher T) **/
        //----------------------------------------------------------------------
        if (get_is_avgFitParams()) {
            real_1d_array c_avg=c_extrp;
            for (size_t i=0; i<c_avgVec.size(); ++i) {
                c_avg[i]=c_avgVec[i];
            }
            /** T's calculated by averaging over trials **/
            Tl_calc = calc_x_given_y(c_avg,taueql,extrp)[0];
            Th_calc = calc_x_given_y(c_avg,taueqh,extrp)[0];
            Th      = Th_calc;
            Tl      = Tleq;
        } else {
            /** T's calculated by each individul trial **/
            Tl_calc = calc_x_given_y(c_extrp,taueql,extrp)[0];
            Th_calc = calc_x_given_y(c_extrp,taueqh,extrp)[0];
            Th      = Th_calc;
            Tl      = Tleq;
        }
        
        /** process value to have correct precision **/
        //----------------------------------------------------------------------
        Tl = std::round(Tl*corrFacforT);//expanded
        Th = std::round(Th*corrFacforT);//expanded
        
        /** update extrpTmax **/
        //----------------------------------------------------------------------
        sysVar.get_extrpTmax().at(index)={Tl,get_Tleq(sysVar,n_trl,n_sys)[1]};//expanded
        
        /** "n_next_temps" (should >= 2) **/
        //----------------------------------------------------------------------
        int n_next_temps=sysVar.get_n_regime_temps().at(current_regime+1);
        if (n_next_temps<2) {n_next_temps=2;}
        
        /** get temperatures for next regime **/
        //----------------------------------------------------------------------
        vector<double> hightolow=
        get_interpolated_Ts(sysVar,n_trl,n_sys,n_next_temps,Tl,Th);
        
        /** replace current temperatures with temperatures of next regime **/
        //----------------------------------------------------------------------
        if (sysVar.get_is_updatetinfo())
        {
            sysVar.get_temperaturesInfo().at(index)=hightolow;
            
            /** clear container for new quench temperatures **/
            //**************************************************************
            sysVar.get_quenchingTs().at(index).clear();
            //**************************************************************
            
            size_t equSize = sysVar.get_equilibratedTs().at(index).size();
            double equTmin = sysVar.get_equilibratedTs().at(index).at(equSize-1);
            
            /** quenchingT's should be < lowest equilibrated T **/
            //--------------------------------------------------------------
            int count=0;
            do {
                double Ti=sysVar.get_temperaturesInfo().at(index).at(count);
                if (Ti<equTmin) { //NOTE: equTmin
                    sysVar.get_quenchingTs().at(index).push_back(Ti);
                } ++count;
            } while (count<sysVar.get_temperaturesInfo().at(index).size());
        }
        
        /** write new temperatures to file **/
        //------------------------------------------------------------------
        if (get_is_avgFitParams())
        {
            if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
            { // only at final index
                for (int i=0; i<n_trial; ++i) {
                    for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) {
                        o.clear();
                        o.append(return_AnalysisFolderPath(sysVar));
                        o.append("/fit_data");
                        o.append("/tempsInfo_");
                        o.append(sysVar.get_usic());
                        o.append("_00"+to_string((long long int)i));
                        o.append("_"+sysVar.get_nameString(ii));
                        o.append(".dat");
                        ofstream write_tempsavg(o.c_str(),ofstream::app); //append
                        //--------------------------------------------------
                        // NOTE:
                        // temperatures in this file has been normalized to the "expanded" format
                        // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
                        // so "setprecision" needs to be set to zero and with fixed format
                        //--------------------------------------------------
                        write_tempsavg << fixed << setprecision(0);
                        
                        int index_i=i*n_system+ii;
                        
                        if (sysVar.get_is_updatetinfo())
                        {
                            sysVar.get_temperaturesInfo().at(index_i)=hightolow;
                            
                            /** clear container for new quench temperatures **/
                            //**********************************************
                            sysVar.get_quenchingTs().at(index_i).clear();
                            //**********************************************
                            
                            size_t equSize=sysVar.get_equilibratedTs().at(index_i).size();
                            double equTmin=sysVar.get_equilibratedTs().at(index_i).at(equSize-1);
                            
                            /** quenchingT's should be < lowest equilibrated T **/
                            //----------------------------------------------
                            int count=0;
                            do {
                                double Ti=sysVar.get_temperaturesInfo().at(index_i).at(count);
                                if (Ti<equTmin) { //NOTE: equTmin
                                    sysVar.get_quenchingTs().at(index_i).push_back(Ti);
                                } ++count;
                            } while (count<sysVar.get_temperaturesInfo().at(index_i).size());
                            //----------------------------------------------
                            
                            for (size_t iii=0; iii<sysVar.get_temperaturesInfo().at(index_i).size(); ++iii) {
                                write_tempsavg << sysVar.get_temperaturesInfo().at(index_i).at(iii) << "\n";
                            }
                            
                        } write_tempsavg.close();
                    }
                }
            }
            
        } else {
            
            if (sysVar.get_is_updatetinfo())
            {
                for (size_t i=0; i<sysVar.get_temperaturesInfo().at(index).size(); ++i) {
                    write_tempsinfo << sysVar.get_temperaturesInfo().at(index).at(i) << "\n";
                }
            }
        }
        
    } /** End if (current_regime<(n_regime-1)) **/
    
    write_tempsinfo.close();
}





/** public setters **/
//------------------------------------------------------------------------------
/* bool */
void FitData::set_is_fit_sExp(const bool b){is_fit_sExp=b;}
void FitData::set_is_fit_lnFs(const bool b){is_fit_lnFs=b;}
void FitData::set_is_1stmoment_tau(const bool b){is_1stmoment_tau=b;}
void FitData::set_is_find_DWF(const bool b){is_find_DWF=b;}
void FitData::set_is_find_NGP(const bool b){is_find_NGP=b;}
void FitData::set_is_fit_Arrhenius(const bool b){is_fit_Arrhenius=b;}
void FitData::set_is_fit_tauFit(const bool b){is_fit_tauFit=b;}
void FitData::set_is_shootForNext(const bool b){is_shootForNext=b;}
void FitData::set_is_avgFitParams(const bool b){is_avgFitParams=b;}
void FitData::set_is_fit_by_TA(const bool b){is_fit_by_TA=b;}
void FitData::set_is_fitByEveryPoint(const bool b){is_fitByEveryPoint=b;}
void FitData::set_is_fit_by_Tc(const bool b){is_fit_by_Tc=b;}
void FitData::set_is_applytauFitcut(const bool b)
{
    is_applytauFitcut=b;
    if (is_applytauFitcut) {
        tauFit_cut=compu_time;
    } else {
        tauFit_cut=extrp_time;
    }
}
void FitData::set_is_fit_Fs_by_spline(const bool b){is_fit_Fs_by_spline=b;}
void FitData::set_is_fit_full_alpha(const bool b){is_fit_full_alpha=b;}
void FitData::set_is_use_gammafunc(const bool b){is_use_gammafunc=b;}
void FitData::set_is_imposeStrictTeq(const bool b){is_imposeStrictTeq=b;}
void FitData::set_is_calc_thermoData(const bool b){is_calc_thermoData=b;}
void FitData::set_is_use_KWWassist(const bool b){is_use_KWWassist=b;}
void FitData::set_is_normalModes(const bool b){is_normalModes=b;}
void FitData::set_is_qspectrum(const bool b){is_qspectrum=b;}
void FitData::set_is_use_ngp_peak_frame(const bool b){is_use_ngp_peak_frame=b;}
void FitData::set_is_use_ngp_smoothed(const bool b){is_use_ngp_smoothed=b;}
/* int */
void FitData::set_index_largest_Tg(const int i){index_largest_Tg=i;}
void FitData::set_index_largest_m(const int i){index_largest_m=i;}
/* double */
void FitData::set_sExp_tauFit(const double d){sExp_tauFit=d;}
void FitData::set_shootratio(const double d){shootratio=d;}
void FitData::set_cutTforArrhenius(const double d){cutTforArrhenius=d;}
void FitData::set_largest_Tg(const double d){largest_Tg=d;}
void FitData::set_largest_m(const double d){largest_m=d;}
void FitData::set_Theq(const double d){Theq=d;}
void FitData::set_Tleq(const double d){Tleq=d;}
void FitData::set_TA_avg(const double d){TA_avg=d;}
void FitData::set_Tc_MCT(const double d){Tc_MCT=d;}
void FitData::set_Tc_SOU(const double d){Tc_SOU=d;}
void FitData::set_tauFit_cut(const double d){tauFit_cut=d;}
void FitData::set_xtauA(const double d){xtauA=d;}
void FitData::set_xTA(const double d){xTA=d;}
/* string */
void FitData::set_relaxation_target(const string& str){relaxation_target=str;}
void FitData::set_definitionofTA(const string& str){definitionofTA=str;}
void FitData::set_fcorr_model(const string& str){fcorr_model=str;}
void FitData::set_extrp_model(const string& str){extrp_model=str;}
void FitData::set_presq_model(const string& str){presq_model=str;}
/* STL */
void FitData::set_waveindices(const vector<int>& vi){waveindices=vi;}
void FitData::set_logtselect(const vector<double>& vd){logtselect=vd;}




