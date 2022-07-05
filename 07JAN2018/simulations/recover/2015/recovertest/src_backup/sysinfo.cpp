//
//  Created by Sean Jh H. on 7/4/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "sysinfo.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

SysInfo::SysInfo():

/* string */
Path(),
simType(),
year(),
usicID(),
systemUnit(),
sys_targets(),

/* bool */
is_makeFileSystem(false),
is_defaultLmpData(false),
is_moltempLmpData(false),
is_LammpsPhase(false),
is_Simulations(false),
is_AMDAT(false),
is_fitData(false),
is_directSub(false),
is_watch_hold(false),
is_GPU(true),
is_fullquench(false),
is_singleTempTest(false),
is_useDynamicRange(false),
is_use_viscosity(false),
is_set_CheckPoint(false),
is_updatetinfo(true),

/* int */
indexi(0),
indexii(0),
indexiii(0),
n_digits(3),
n_trial(),
n_system(),
n_sys_beg(),
n_sys_end(),
n_startTemp(),
n_Temp(),
current_regime(),
n_regime(),
quenchMode(),
times_retry(0),
times_write_tempsinfo(0),
max_times_retry(0),

/* double */
n_equ_blocks(),
n_prd_blocks(),
startTemp(1000.0),
crossTemp(700.0),
finalTemp(200.0),
hTres(75.0),
lTres(50.0),
corrFacforT(1.0),
precision(1e+2),
dynamicRange_tau0(),
dynamicRange_tauA(),

/* vector<int> */
n_regime_temps(),

/* vector<double> */
ts_regime(),
equilibration_times(),
equilibration_times_org(),

/* temperature containers */
temperaturesInfo(),
equilibratedTs(),
quenchingTs(),
qchTdiff(),

/* variables used in the fianl reports */
startdatetime(),
relaxationModel(),
extrp_time(),
compu_time(),
timestep_size(),
quenchRate(),
timediff_overall(),
wavenumber(),
elapsedtime(),
Tg_extrp(),
Tg_compu(),
m_extrp(),
m_compu(),
is_node0_info(),
is_node1_info(),
node0_info(),
node1_info(),

/* simple stats */
max(),
mean(),
value()
{}


/*==( public setters )==*/

/* string */
void SysInfo::set_Path(const string& s){Path=s;}
void SysInfo::set_simType(const string& s){simType=s;}
void SysInfo::set_year(const string& s){year=s;}
void SysInfo::set_usicID(const string& s){usicID=s;}
void SysInfo::set_systemUnit(const std::string& s){systemUnit=s;}
void SysInfo::set_sys_targets(const string& str){sys_targets=str;}

/* bool */
void SysInfo::set_is_makeFileSystem(const bool b){is_makeFileSystem=b;}
void SysInfo::set_is_defaultLmpData(const bool b){is_defaultLmpData=b;}
void SysInfo::set_is_moltempLmpData(const bool b){is_moltempLmpData=b;}
void SysInfo::set_is_LammpsPhase(const bool b){is_LammpsPhase=b;}
void SysInfo::set_is_Simulations(const bool b){is_Simulations=b;}
void SysInfo::set_is_AMDAT(const bool b){is_AMDAT=b;}
void SysInfo::set_is_fitData(const bool b){is_fitData=b;}
void SysInfo::set_is_directSub(const bool b){is_directSub=b;}
void SysInfo::set_is_watch_hold(const bool b){is_watch_hold=b;}
void SysInfo::set_is_GPU(const bool b){is_GPU=b;}
void SysInfo::set_is_fullquench(const bool b){is_fullquench=b;}
void SysInfo::set_is_singleTempTest(const bool b){is_singleTempTest=b;}
void SysInfo::set_is_useDynamicRange(const bool b){is_useDynamicRange=b;}
void SysInfo::set_is_use_viscosity(const bool b){is_use_viscosity=b;}
void SysInfo::set_is_set_CheckPoint(const bool b){is_set_CheckPoint=b;}
void SysInfo::set_is_updatetinfo(const bool b){is_updatetinfo=b;}
void SysInfo::set_is_theorytest(const bool b){is_theorytest=b;}
void SysInfo::set_is_fit_dielectric(const bool b){is_fit_dielectric=b;}

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
void SysInfo::set_quenchMode(const int i){quenchMode=i;}

/* double */
void SysInfo::set_n_equ_blocks(const double d){n_equ_blocks=d;}
void SysInfo::set_n_prd_blocks(const double d){n_prd_blocks=d;}
void SysInfo::set_startTemp(const double f){startTemp=f;}
void SysInfo::set_crossTemp(const double f){crossTemp=f;}
void SysInfo::set_finalTemp(const double f){finalTemp=f;}
void SysInfo::set_hTres(const double f){hTres=f;}
void SysInfo::set_lTres(const double f){lTres=f;}
void SysInfo::set_corrFacforT(const double d){corrFacforT=d;}
void SysInfo::set_precision(const double d){precision=d;}
void SysInfo::set_cutTforArrhenius(const double d){cutTforArrhenius=d;}

void SysInfo::set_dynamicRange_tau0(const double d)
{
    dynamicRange_tau0=pow(10.0,d);
}

void SysInfo::set_dynamicRange_tauA(const double d)
{
    dynamicRange_tauA=pow(10.0,d);
}

/* vector<int> */
void SysInfo::set_n_regime_temps(const vector<int>& vi)
{n_regime_temps=vi;}

/* vector<double> */
void SysInfo::set_ts_regime(const std::vector<double>& vd)
{ts_regime=vd;}
void SysInfo::set_equilibration_times(const vector<double>& vd)
{equilibration_times=vd;}
void SysInfo::set_equilibration_times_org(const std::vector<double>& vd)
{equilibration_times_org=vd;}

/* vector<vector<double>> */
void SysInfo::set_initialtemps(const std::vector<std::vector<double>>& vvd)
{initialtemps=vvd;}

/* variables used in the fianl reports */
//------------------------------------------------------------------------------
/* string */
void SysInfo::set_startdatetime(const std::string& s){startdatetime=s;}
void SysInfo::set_relaxationModel(const std::string& s){relaxationModel=s;}
void SysInfo::set_extrpolateModel(const std::string& s){extrpolateModel=s;}
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
/* vector<int> */
void SysInfo::set_node0_info(const vector<int>& vi){node0_info=vi;}
void SysInfo::set_node1_info(const vector<int>& vi){node1_info=vi;}
//------------------------------------------------------------------------------



/*==( public getters )==*/

/* string */
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

/* int */
int SysInfo::get_n_system() const
{
    int first_sys=get_n_sys_beg();
    int final_sys=get_n_sys_end();
    return (final_sys-first_sys+1); // inclusive counting
}

/* simple stats */
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
