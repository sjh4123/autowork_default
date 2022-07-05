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

#include "sysinfo.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

SysInfo::SysInfo():

/* bool */
is_makeFileSystem(false),
is_defaultLmpData(false),
is_moltempLmpData(false),
is_LammpsPhase(false),
is_use_prepScripts(false),
is_Simulations(false),
is_AMDAT(false),
is_fitData(false),
is_directSub(false),
is_watch_hold(false),
is_GPU(true),              // true: enables the GPU features
is_fullquench(false),      // true: quench from highest to lowest T or by regime
is_singleTempTest(false),  // true: single temperature test: gen-->prd-->amdat
is_use_viscosity(false),
is_set_CheckPoint(false),
is_updatetinfo(true),      // true: update temperaturesInfo and quenchingTs
is_theorytest(false),
is_fit_dielectric(false),
is_specify_GPU(false),
is_cutoff_traject(false),
is_alglibfit(false),
is_use_named_nameString(false),
is_makeNewAnanlysisFolder(false),
is_backupfitdata(false),
is_backupstatistics(false),
is_return(false),
is_amdatinp(false),
is_copy_from_to(false),
is_aging(false),
is_cancelRetryAlgo(false),

/* int */
indexi(0),
indexii(0),
indexiii(0),
n_digits(3),
n_trial(4),
n_system(1),
n_sys_beg(0),
n_sys_end(0),
n_startTemp(8),
n_Temp(8),
regime_beg(0),
regime_end(0),
current_regime(0),
n_regime(1),
times_retry(0),
max_times_retry(2),
n_cum_taueq(0),
which_GPU(0),
sim_restart(0),
gen_restart(0),
qch_restart(0),
equ_restart(0),
res_restart(0),
prd_restart(0),
ana_restart(0),
fd_res(0),
tt_res(0),
max_initialized_regime(0),
n_kuhnsegs(1),

/* double */
n_equ_blocks(10.0),
n_prd_blocks(10.0),
n_relaxation(n_equ_blocks),
startTemp(1000.0),
crossTemp(700.0),
finalTemp(200.0),
hTres(75.0),
lTres(50.0),
corrFacforT(1.0),
precision(1000.0),
cutTforArrhenius(1000.0),
dynamicRange_tau0(4.0),
dynamicRange_tauA(4.0),

/* string */
userName("jh148"),
masterFolder("AUTOWORK"),
Path(""),
simType(""),
year(""),
usicID(""),
systemUnit(""),
sys_targets(""),
computeCluster("Mycroft"),
copy_source(""),
copy_destin(""),

/** STL containers **/
//------------------------------------------------------------------------------
n_regime_temps(),
ts_regime(),
equilibration_times(),
equilibration_times_org(),
extrpTmax(),

/** simple stats **/
//------------------------------------------------------------------------------
max(0),
mean(0),
value(),

/** variables used in fitData for KWW fit **/
//------------------------------------------------------------------------------
epsf(0),
epsx(0),
diffstep(0),
coeffs_bndl_vD(),
coeffs_bndu_vD(),
coeffs_bndl(),
coeffs_bndu(),
c_vD(),
rep(),
maxits(),

/** variables used in the fianl reports **/
//------------------------------------------------------------------------------
is_node0_info(false),
is_node1_info(false),
is_node2_info(false),
extrp_time(0),
compu_time(0),
timestep_size(0),
quenchRate(0),
timediff_overall(0),
wavenumber(0),
startdatetime(),
message_dynRange(),
elapsedtime(),
extrp_model(),
presq_model(),
Tg_extrp(),
Tg_compu(),
m_extrp(),
m_compu(),
Tg_VFT(),
Tg_COOP(),
m_VFT(),
m_COOP(),
node0_info(),
node1_info(),
simulation_retry(),

/** temperature containers **/
//------------------------------------------------------------------------------
n_equ_blocks_plus(),
n_prd_blocks_plus(),
n_relaxation_plus(),
tautargets(),
initialtemps(),
temperaturesInfo(),
equilibratedTs(),
quenchingTs(),
qchTdiff(),
tauFit(),
tauEqu(),
inptinfo()
{
    /** Constructor Assignment of SysInfo **/
}


/** public setters **/
//------------------------------------------------------------------------------
/* bool */
void SysInfo::set_is_makeFileSystem(const bool b){is_makeFileSystem=b;}
void SysInfo::set_is_defaultLmpData(const bool b){is_defaultLmpData=b;}
void SysInfo::set_is_moltempLmpData(const bool b){is_moltempLmpData=b;}
void SysInfo::set_is_LammpsPhase(const bool b){is_LammpsPhase=b;}
void SysInfo::set_is_use_prepScripts(const bool b){is_use_prepScripts=b;}
void SysInfo::set_is_Simulations(const bool b){is_Simulations=b;}
void SysInfo::set_is_AMDAT(const bool b){is_AMDAT=b;}
void SysInfo::set_is_fitData(const bool b){is_fitData=b;}
void SysInfo::set_is_directSub(const bool b){is_directSub=b;}
void SysInfo::set_is_watch_hold(const bool b){is_watch_hold=b;}
void SysInfo::set_is_GPU(const bool b){is_GPU=b;}
void SysInfo::set_is_fullquench(const bool b){is_fullquench=b;}
void SysInfo::set_is_singleTempTest(const bool b){is_singleTempTest=b;}
void SysInfo::set_is_use_viscosity(const bool b){is_use_viscosity=b;}
void SysInfo::set_is_set_CheckPoint(const bool b){is_set_CheckPoint=b;}
void SysInfo::set_is_updatetinfo(const bool b){is_updatetinfo=b;}
void SysInfo::set_is_theorytest(const bool b){is_theorytest=b;}
void SysInfo::set_is_fit_dielectric(const bool b){is_fit_dielectric=b;}
void SysInfo::set_is_cutoff_traject(const bool b){is_cutoff_traject=b;}
void SysInfo::set_is_alglibfit(const bool b){is_alglibfit=b;}
void SysInfo::set_is_use_named_nameString(const bool b){is_use_named_nameString=b;}
void SysInfo::set_is_makeNewAnanlysisFolder(const bool b){is_makeNewAnanlysisFolder=b;}
void SysInfo::set_is_backupfitdata(const bool b){is_backupfitdata=b;}
void SysInfo::set_is_backupstatistics(const bool b){is_backupstatistics=b;}
void SysInfo::set_is_return(const bool b){is_return=b;}
void SysInfo::set_is_amdatinp(const bool b){is_amdatinp=b;}
void SysInfo::set_is_copy_from_to(const bool b){is_copy_from_to=b;}
void SysInfo::set_is_aging(const bool b){is_aging=b;}
void SysInfo::set_is_cancelRetryAlgo(const bool b){is_cancelRetryAlgo=b;}
/* int */
void SysInfo::set_n_digits(const int i){n_digits=i;}
void SysInfo::set_n_trial(const int i){n_trial=i;}
void SysInfo::set_n_system(const int i){n_system=i;}
void SysInfo::set_n_sys_beg(const int i){n_sys_beg=i;}
void SysInfo::set_n_sys_end(const int i){n_sys_end=i;}
void SysInfo::set_n_startTemp(const int i){n_startTemp=i;}
void SysInfo::set_regime_beg(const int i){regime_beg=i;}
void SysInfo::set_regime_end(const int i){regime_end=i;}
void SysInfo::set_n_Temp(const int i){n_Temp=i;}
void SysInfo::set_current_regime(const int i){current_regime=i;}
void SysInfo::set_n_regime(const int i){n_regime=i;}
void SysInfo::set_times_retry(const int i){times_retry=i;}
void SysInfo::set_max_times_retry(const int i){max_times_retry=i;}
void SysInfo::set_n_cum_taueq(const int i){n_cum_taueq=i;}
void SysInfo::add_n_cum_taueq(const int i){n_cum_taueq+=i;}//NOTE
void SysInfo::set_which_GPU(const int i)
{
    which_GPU=i;
    is_specify_GPU=true;
}
void SysInfo::set_sim_restart(const int i){sim_restart=i;}
void SysInfo::set_gen_restart(const int i){gen_restart=i;}
void SysInfo::set_qch_restart(const int i){qch_restart=i;}
void SysInfo::set_equ_restart(const int i){equ_restart=i;}
void SysInfo::set_res_restart(const int i){res_restart=i;}
void SysInfo::set_prd_restart(const int i){prd_restart=i;}
void SysInfo::set_ana_restart(const int i){ana_restart=i;}
void SysInfo::set_fd_res(const int i){fd_res=i;}
void SysInfo::set_tt_res(const int i){tt_res=i;}
void SysInfo::set_max_initialized_regime(const int i){max_initialized_regime=i;}
void SysInfo::set_n_kuhnsegs(const int i){n_kuhnsegs=i;}
/* double */
void SysInfo::set_n_equ_blocks(const double d){n_equ_blocks=d;}
void SysInfo::set_n_prd_blocks(const double d){n_prd_blocks=d;}
void SysInfo::set_n_relaxation(const double d){n_relaxation=d;}
void SysInfo::set_startTemp(const double f){startTemp=f;}
void SysInfo::set_crossTemp(const double f){crossTemp=f;}
void SysInfo::set_finalTemp(const double f){finalTemp=f;}
void SysInfo::set_hTres(const double f){hTres=f;}
void SysInfo::set_lTres(const double f){lTres=f;}
void SysInfo::set_corrFacforT(const double d){corrFacforT=d;}
void SysInfo::set_precision(const double d){precision=d;}
void SysInfo::set_cutTforArrhenius(const double d){cutTforArrhenius=d;}
void SysInfo::set_dynamicRange_tau0(const double d){dynamicRange_tau0=pow(10.0,d);}
void SysInfo::set_dynamicRange_tauA(const double d){dynamicRange_tauA=pow(10.0,d);}
/* string */
void SysInfo::set_userName(const std::string& s){userName=s;}
void SysInfo::set_masterFolder(const std::string& s){masterFolder=s;}
void SysInfo::set_Path(const string& s){Path=s;}
void SysInfo::set_simType(const string& s){simType=s;}
void SysInfo::set_year(const string& s){year=s;}
void SysInfo::set_usicID(const string& s){usicID=s;}
void SysInfo::set_systemUnit(const std::string& s){systemUnit=s;}
void SysInfo::set_sys_targets(const string& str){sys_targets=str;}
void SysInfo::set_computeCluster(const std::string& str){computeCluster=str;}
void SysInfo::set_copy_source(const std::string& str){copy_source=str;}
void SysInfo::set_copy_destin(const std::string& str){copy_destin=str;}
/* vector<int> */
void SysInfo::set_n_regime_temps(const vector<int>& vi){n_regime_temps=vi;}
/* vector<double> */
void SysInfo::set_n_equ_blocks_plus(const std::vector<double>& vd){n_equ_blocks_plus=vd;}
void SysInfo::set_n_prd_blocks_plus(const std::vector<double>& vd){n_prd_blocks_plus=vd;}
void SysInfo::set_n_relaxation_plus(const std::vector<double>& vd){n_relaxation_plus=vd;}
void SysInfo::set_tautargets(const std::vector<double>& vd){tautargets=vd;}
void SysInfo::set_ts_regime(const std::vector<double>& vd){ts_regime=vd;}
void SysInfo::set_equilibration_times(const vector<double>& vd)
{equilibration_times=vd;}
void SysInfo::set_equilibration_times_org(const std::vector<double>& vd)
{equilibration_times_org=vd;}
/* vector<vector<double>> */
void SysInfo::set_initialtemps(const std::vector<std::vector<double>>& vvd)
{initialtemps=vvd;}
/* vector<vector<vector<double>>> */
void SysInfo::set_inptinfo(const vector<vector<vector<double>>>& vvvd)
{inptinfo=vvvd;}
void SysInfo::set_regime_tinfo(const int regime,const vector<vector<double>>& tinfo)
{inptinfo.at(regime)=tinfo;}
void SysInfo::show_inptinfo()
{
    for (size_t i=0;i<inptinfo.size();++i) {//regime
        cout << "regime "<<i<<"\n";
        for (size_t ii=0;ii<inptinfo.at(i).size();++ii) {//trial
            for (size_t iii=0;iii<inptinfo.at(i).at(ii).size();++iii) {//T
                cout << (int)inptinfo.at(i).at(ii).at(iii) << " ";
            } cout << "\n";
        } //cout << "\n";
    }
}
void SysInfo::show_regime_tinfo(const int regime)
{
    for (size_t i=regime;i<regime+1;++i) {//regime
        cout << "regime "<<i<<"\n";
        for (size_t ii=0;ii<inptinfo.at(i).size();++ii) {//trial
            for (size_t iii=0;iii<inptinfo.at(i).at(ii).size();++iii) {//T
                cout << (int)inptinfo.at(i).at(ii).at(iii) << " ";
            } cout << "\n";
        } cout << "\n";
    }
}

/** variables used in fitData for KWW fit **/
//------------------------------------------------------------------------------
void SysInfo::set_coeffs_bndl_vD(const std::vector<double>& vd){coeffs_bndl_vD=vd;}
void SysInfo::set_coeffs_bndu_vD(const std::vector<double>& vd){coeffs_bndu_vD=vd;}
void SysInfo::set_coeffs_bndl(const std::string& str){coeffs_bndl=str;}
void SysInfo::set_coeffs_bndu(const std::string& str){coeffs_bndu=str;}
void SysInfo::set_c_vD(const std::vector<double>& vd){c_vD=vd;}
void SysInfo::set_rep(const alglib::lsfitreport& rep_d){rep=rep_d;}
void SysInfo::set_epsf(const double d){epsf=d;}
void SysInfo::set_epsx(const double d){epsx=d;}
void SysInfo::set_diffstep(const double d){diffstep=d;}
void SysInfo::set_maxits(const alglib::ae_int_t& maxits_d){maxits=maxits_d;}

/** variables used in fianl report **/
//------------------------------------------------------------------------------
/* string */
void SysInfo::set_startdatetime(const std::string& s){startdatetime=s;}
void SysInfo::set_message_dynRange(const std::string& s){message_dynRange=s;}
/* double */
void SysInfo::set_extrp_time(const double d){extrp_time=d;}
void SysInfo::set_compu_time(const double d){compu_time=d;}
void SysInfo::set_timestep_size(const double d){timestep_size=d;}
void SysInfo::set_quenchRate(const double d){quenchRate=d;}
void SysInfo::set_timediff_overall(const double d){timediff_overall=d;}
void SysInfo::set_wavenumber(const double d){wavenumber=d;}
/* bool */
void SysInfo::set_is_node0_info(const bool b){is_node0_info=b;}
void SysInfo::set_is_node1_info(const bool b){is_node1_info=b;}
void SysInfo::set_is_node2_info(const bool b){is_node2_info=b;}
/* vector<int> */
void SysInfo::set_node0_info(const vector<int>& vi){node0_info=vi;}
void SysInfo::set_node1_info(const vector<int>& vi){node1_info=vi;}
void SysInfo::set_node2_info(const vector<int>& vi){node2_info=vi;}


/** public getters **/
//------------------------------------------------------------------------------
const string SysInfo::get_year() const
{
    //return autoWork::return_year();
    return year;
}
const string SysInfo::get_usic() const
{
    //string usic=get_simType()+get_year()+get_usicID();
    string usic=get_simType()+get_usicID();
    return usic;
}
int SysInfo::get_n_system() const
{
    int first_sys=get_n_sys_beg();
    int final_sys=get_n_sys_end();
    return (final_sys-first_sys+1); // inclusive counting
}
void SysInfo::set_calcvector(const vector<double>& p)
{
    max = (int)p.size();
    value.clear();
    for(size_t i=0; i<p.size(); ++i) value.push_back(p[i]);
}

double SysInfo::calc_mean()
{
    double sum = 0;
    for(int i = 0; i < max; ++i) sum += value[i];
    return (sum/max);
}

double SysInfo::calc_variance()
{
    mean = calc_mean();
    double temp = 0;
    for(int i = 0; i < max; ++i)
    {
        temp += (value[i] - mean) * (value[i] - mean) ;
    }
    return (temp/max);
}
double SysInfo::calc_sampleVariance()
{
    mean = calc_mean();
    double temp = 0;
    for(int i = 0; i < max; ++i)
    {
        temp += (value[i] - mean) * (value[i] - mean) ;
    }
    return (temp/(max-1));
}




