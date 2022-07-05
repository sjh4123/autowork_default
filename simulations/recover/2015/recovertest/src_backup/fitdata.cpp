//
//  Created by SJH on 8/11/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "fitdata.h"
#include "functions.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

/*
 
 TODO:
 make tau prefactor log scale
 
 */

FitData::FitData(const StructureClass& sysVar,
                 const WorkScripts& ws,
                 const AmdatAnalysis& aa):
/* base */
AlglibFittingKernel(),

/* bool */
is_fit_sExp(false),
is_fit_Arrhenius(false),
is_fit_tauFit(false),
is_shootForNext(true),
is_avgFitParams(false),
is_fit_by_TA(false),
is_fitByEveryPoint(false),
is_fit_by_Tc(false),
is_use_viscosity(false),
is_applytauFitcut(true),
is_fit_full_alpha(false),
is_use_sExp_cut_hi(false),
is_use_plateautime(false),
is_use_gammafunc(false),
is_imposeStrictTeq(false),

/* int */
n_fit_sExp(5),
n_fit_chopping(5),
n_fit_sliding(10),

/* double */
DWF(0),
deviateTA(0.1),
sExp_tau_lo(1.0),
sExp_cut_hi(0.6),
sExp_cut_lo(0.01),
sExp_tauFit(0.2),
shootratio(0),
cutTforArrhenius(1.0),
Tc_percent(0.99),
res_Arrhenius(1.1),
r2_threshold(0.95)
{
    
    /* cutoff T in Arrhenius fit */
    cutTforArrhenius=sysVar.get_cutTforArrhenius();
    
    if (sysVar.get_systemUnit()=="real")
    {
        sExp_tau_lo = 1e+3; // fs
    }
    else if (sysVar.get_systemUnit()=="lj")
    {
        sExp_tau_lo = 1.0;  // tau
    }
    else if (sysVar.get_systemUnit()=="metal")
    {
        sExp_tau_lo = 1.0;  // ps
    }
    
    sExp_cut_hi_org  = sExp_cut_hi;
    
    systemUnit       = sysVar.get_systemUnit();
    corrFacforT      = sysVar.get_corrFacforT();
    precision        = sysVar.get_precision();
    extrp_time       = sysVar.get_extrp_time();
    compu_time       = sysVar.get_compu_time();
    
    current_regime   = sysVar.get_current_regime();
    n_regime         = sysVar.get_n_regime();
    n_trial          = sysVar.get_n_trial();
    n_system         = sysVar.get_n_system();
    n_sys_beg        = sysVar.get_n_sys_beg();
    n_sys_end        = sysVar.get_n_sys_end();
    
    is_use_viscosity = sysVar.get_is_use_viscosity();
    
    NewtonGuessTemp  *= precision;
    
    n_equ_blocks     = ws.get_n_equ_blocks();
    n_prd_blocks     = ws.get_n_prd_blocks();
    
    analysispart     = aa.get_analysispart();
    DWF_time         = aa.get_DWF_time();
    
    if (get_is_applytauFitcut()) {
        tauFit_cut=ws.get_time_equ()/ws.get_n_equ_blocks();
    }
    else {
        tauFit_cut=sysVar.get_extrp_time();
    }
}





void FitData::fit_sExp(StructureClass& sysVar,
                       const int n_trl,
                       const int n_sys,
                       const double Temp_d,
                       const string& tcalc)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    //==========================================================================
    // tau_Fit file stores relaxation time of all fit temoerature
    //==========================================================================
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/tauFit_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tauFit(o.c_str(),ofstream::app); // 'append'
    
    //==========================================================================
    // taueq file stores the relaxation time of in-equilibrium temperatures
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_taueq(o.c_str(),ofstream::app); // 'append'
    
    //==========================================================================
    // taueq_inverseT stores the same content as taueq file except inverse T's
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_taueqinv(o.c_str(),ofstream::app); // 'append'
    
    //==========================================================================
    // sExp_params stores the fitting parameters at different fit temperatures
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/sExp_params_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_params(o.c_str(),ofstream::app); // 'append'
    
    //==========================================================================
    // fit_sExp file stores the overall fit information
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fit_sExp_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_sExp(o.c_str());
    
    
    /* setting up initial values of fitting parameters */
    //----------------------------------------------------------------------
    bool is_use_assist=true;
    if (is_use_assist)
    {
        // NOTE:
        // if previous fit r2 > r2_threshold, use params from previous fit
        // if < r2_threshold, use default params
        if (sysVar.rep.r2<r2_threshold) set_fitParams(tcalc);
        else load_fit_coeffs(sysVar);
    }
    else
    {
        set_fitParams(tcalc);
    }
    
    
    /* read in relaxation profile at each temperature */
    //----------------------------------------------------------------------
    if (is_fit_full_alpha)
    {
        /* fit full alpha relaxation */
        skim_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d);
        read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d);
    }
    else
    {
        /* data interpolated inbetween "sExp_cut_lo < corr < sExp_cut_hi" */
        is_use_sExp_cut_hi=true;
        is_use_plateautime=false;
        read_individual_correlation_data(sysVar,n_trl,n_sys,Temp_d);
    }
    sysVar.set_wavenumber(wavenumber);
    
    
    /* Data processing to ALGLIB readable form */
    //----------------------------------------------------------------------
    fitdata_processing_lsfit(correlation_sortincreasing);
    
    
    /* Curve-fitting Kernel */
    //----------------------------------------------------------------------
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
        
        if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
        {
            // last index
            size_t tmpsize=sysVar.get_quenchingTs()[index].size();
            if (Temp_d==sysVar.get_quenchingTs()[index][tmpsize-1])
            {
                // last T
                write_sExp_params_avg(sysVar,n_sys,tcalc);
            }
        }
    }
    if (rep.r2>r2_threshold) save_fit_coeffs(sysVar);
    
    
    if ((rep.r2<r2_threshold)||(correlation_sortincreasing.size()<n_fit_sExp))
    {
        double faketau=log10(tauFit_cut)+1; // NOTE!
        
        //------------------------------------------------------------------
        write_tauFit
        << Temp_d*pow(corrFacforT,-1) << " " << faketau << "\n";
        sysVar.get_tauFit()[index].push_back({Temp_d,faketau}); // use expanded form
        //------------------------------------------------------------------
        
        if (rep.r2<r2_threshold)
        {
            write_badR2(write_sExp,r2_threshold);
        }
        if (correlation_sortincreasing.size()<n_fit_sExp)
        {
            write_insufficientDara
            (write_sExp,correlation_sortincreasing.size(),n_fit_sExp);
        }
        return;
    }
    
    ////////////////////////////////////////////////////////////////////////
    //======================================================================
    // NOTE:
    // writing out in-equilibrium T's to file has been moved to
    // "write_equTs()" function
    //======================================================================
    bool is_write_to_file=false;
    
    
    double T_actual   = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    double log_tauFit = 0;
    
    if (is_use_gammafunc)
    {
        tauFit_calc=sExp_tauFit_gammafunction(c,tcalc);
    }
    else
    {
        tauFit_calc=sExp_tauFit_interpolateFvalue(c,tcalc);
    }
    log_tauFit = log10(tauFit_calc);
    
    
    /* write tauFit data to file */
    //----------------------------------------------------------------------
    write_tauFit << T_actual << " " << log_tauFit << "\n";
    sysVar.get_tauFit()[index].push_back({Temp_d,log_tauFit});
    //----------------------------------------------------------------------
    
    
    /* write tau_eq data to file */
    //----------------------------------------------------------------------
    if (sysVar.get_is_useDynamicRange())
    {
        if (current_regime==0)
        {
            if (is_write_to_file)
            {
                write_taueq
                << T_actual << " " << log_tauFit << "\n";
                write_taueqinv
                << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                << log_tauFit
                << "\n";
            }
            sysVar.get_equilibratedTs()[index].push_back(Temp_d);
            sysVar.get_tauEqu()[index].push_back({Temp_d,log_tauFit});
        }
        else
        {
            if (tauFit_calc<tauFit_cut)
            {
                // NOTE: apply cutoff
                if (is_write_to_file)
                {
                    write_taueq
                    << T_actual << " " << log_tauFit << "\n";
                    write_taueqinv
                    << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                    << log_tauFit
                    << "\n";
                }
                sysVar.get_equilibratedTs()[index].push_back(Temp_d);
                sysVar.get_tauEqu()[index].push_back({Temp_d,log_tauFit});
            }
        }
    }
    else
    {
        if (tauFit_calc<tauFit_cut)
        {
            // NOTE: apply cutoff
            if (is_write_to_file)
            {
                write_taueq
                << T_actual << " " << log_tauFit << "\n";
                write_taueqinv
                << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                << log_tauFit
                << "\n";
            }
            sysVar.get_equilibratedTs()[index].push_back(Temp_d);
            sysVar.get_tauEqu()[index].push_back({Temp_d,log_tauFit});
        }
    }
    //======================================================================
    ////////////////////////////////////////////////////////////////////////
    
    
    write_tauFit.close();
    write_taueq.close();
    write_taueqinv.close();
    write_params.close();
    write_sExp.close();
}





bool FitData::check_is_retry(StructureClass& sysVar)
{
    const int n_trial        = sysVar.get_n_trial();
    const int n_system       = sysVar.get_n_system();
    const int n_sys_beg      = sysVar.get_n_sys_beg();
    const int n_sys_end      = sysVar.get_n_sys_end();
    const int current_regime = sysVar.get_current_regime();
    const int n_regimeTs     = sysVar.get_n_regime_temps()[current_regime];
    
    double Tie_min=0; // lowest in-equilibrium(ie) T
    double tie_min=0; // tau_alpha at Tie_min
    double Toe_max=0; // highest out-of-equilibrium(oe) T
    double toe_max=0; // tau_alpha at Toe_max
    double slope=0;   // slope on the Arrhenius plot (linking Tie_min & Toe_max)
    
    double Tmin=0;    // interpolated T at higher tau_alpha cutoff
    double tmin=log10(tauFit_cut); // higher tau_alpha cutoff
    
    double fraction=0.5;
    int n_threshold=(int)(fraction*n_regimeTs);
    
    int count_neq=0;
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            int index=n_trl*n_system+(n_sys-n_sys_beg);
            
            /** Find first ooe T **/
            Toe_max=0;
            for (int i=0; i<sysVar.get_tauFit()[index].size(); ++i) {
                if (sysVar.get_tauFit()[index][i][1]>=log10(tauFit_cut)) {
                    Toe_max=sysVar.get_tauFit()[index][i][0]; // use expanded form
                    break;
                }
            }
            
            if (is_imposeStrictTeq) {
                
                /** all in-equilibrium T's should be greater than first
                 out-of-equilibrium T **/
                
                if (Toe_max==0) {
                    // all T's are in-equilibrium
                    count_neq += sysVar.get_tauFit()[index].size();
                } else {
                    for (int i=0; i<sysVar.get_equilibratedTs()[index].size(); ++i) {
                        if (sysVar.get_equilibratedTs()[index][i]>Toe_max) {
                            ++count_neq;
                        }
                    }
                }
            } else {
                count_neq += sysVar.get_tauEqu()[index].size();
            }
            
        }
    }
    cout
    << "\n"
    << "count_neq   " << count_neq   << "\n"
    << "n_threshold " << n_threshold << "\n";
    
    bool is_retry=false;
    if (count_neq<n_threshold) is_retry=true;
    else is_retry=false;
    
    //--------------------------------------------------------------------------
    // if number of ie temperatures < n_threshold:
    // Interpolate new temperatures between the loest ie temperature
    // and the first, also the highest, oe temperature.
    //--------------------------------------------------------------------------
    if (is_retry)
    {
        vector<double> hightolow;
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl)
        {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
            {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                /* Find lowest ie T */
                //--------------------------------------------------------------
                read_individual_taueq_data(sysVar,n_trl,n_sys);
                Tie_min=std::trunc(taueq_sortincreasing[0][0]*corrFacforT); // use expanded form
                tie_min=taueq_sortincreasing[0][1];
                
                /* Find first ooe T */
                //--------------------------------------------------------------
                for (int i=0; i<sysVar.get_tauFit()[index].size(); ++i)
                {
                    if (sysVar.get_tauFit()[index][i][1]>=log10(tauFit_cut))
                    {
                        Toe_max=sysVar.get_tauFit()[index][i][0]; // use expanded form
                        toe_max=sysVar.get_tauFit()[index][i][1];
                        break;
                    }
                }
                
                /* Every ie T should be higher than the first ooe T */
                //--------------------------------------------------------------
                if (Toe_max>Tie_min)
                {
                    double Ti=0;
                    for (int i=0; i<taueq_sortdecreasing.size(); ++i)
                    {
                        Ti=std::trunc(taueq_sortdecreasing[i][0]*corrFacforT); // use expanded form;
                        if (Ti>Toe_max)
                        {
                            Tie_min=Ti;
                            tie_min=taueq_sortdecreasing[i][1];
                        }
                    }
                }
                
                /* if ooe T is higher than all originally stored ie T of current
                 regime, use the lowest ie T simulated */
                //--------------------------------------------------------------
                if (sysVar.get_equilibratedTs()[index].size()==0)
                {
                    sysVar.get_equilibratedTs()[index].push_back(Tie_min);
                }
                //--------------------------------------------------------------
                
                
                
                //--------------------------------------------------------------
                // Record simulated Ts before retried Ts are interpolated
                //--------------------------------------------------------------
                string folder,file,output,input;
                folder.append(return_AnalysisFolderPath(sysVar));
                folder.append("/fit_data");
                file.append("/tempsInfo_");
                file.append(sysVar.get_usic());
                file.append("_00"+to_string((long long int)n_trl));
                file.append("_"+sysVar.get_nameString(n_sys));
                file.append(".dat");
                output=folder+file;
                ofstream write_tempsinfo(output.c_str(),ofstream::app); // use "append"
                //--------------------------------------------------------------
                // NOTE:
                // temperatures in this file has been normalized to "expanded" format
                // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
                // so "setprecision" needs to be set to zero
                //--------------------------------------------------------------
                write_tempsinfo << fixed << setprecision(0);
                for (size_t i=0; i<sysVar.get_temperaturesInfo()[index].size(); ++i)
                {
                    write_tempsinfo << sysVar.get_temperaturesInfo()[index][i] << "\n";
                }
                sysVar.times_write_tempsinfo++;
                //--------------------------------------------------------------
                
                
                
                //--------------------------------------------------------------
                // NOTE:
                // use actual value, instead of expanded, for calculating slope
                // on an Arrhenius plot.
                //--------------------------------------------------------------
                double Tie_actual=Tie_min*pow(corrFacforT,-1);
                double Toe_actual=Toe_max*pow(corrFacforT,-1);
                double Tmin_actual=0;
                slope=get_2pt_slope({pow(Tie_actual,-1),tie_min},{pow(Toe_actual,-1),toe_max});
                if (slope<0)
                {
                    cout
                    << "\n"
                    << "WARNINIG: slope of interpolation < 0; something is wrong."
                    << "\n";
                    cout
                    << "Tie_actual = " << Tie_actual << "\n"
                    << "tie_min    = " << tie_min    << "\n"
                    << "Toe_actual = " << Toe_actual << "\n"
                    << "toe_max    = " << toe_max    << "\n";
                }
                Tmin_actual=pow(pow(Tie_actual,-1)+((tmin-tie_min)/slope),-1);
                Tmin=std::trunc(Tmin_actual*corrFacforT); // use expanded form
                hightolow=get_interpolated_Ts(sysVar,n_trl,n_sys,n_regimeTs,Tie_min,Tmin);
                if (sysVar.get_is_updatetinfo()) {
                    sysVar.get_temperaturesInfo()[index]=hightolow;}
                sysVar.get_quenchingTs()[index].clear();
                
                cout
                << "\n"
                << "Tie_min " << Tie_actual  << "\n"
                << "tie_min " << tie_min     << "\n"
                << "Toe_max " << Toe_actual  << "\n"
                << "toe_max " << toe_max     << "\n"
                << "Tmin    " << Tmin_actual << "\n"
                << "tmin    " << tmin        << "\n"
                << "\n";
                
                if (sysVar.times_retry<sysVar.max_times_retry)
                {
                    // NOTE: "<"
                    int count=0;
                    do
                    {
                        double Ti=sysVar.get_temperaturesInfo()[index][count];
                        if (Ti<Tie_min)
                        {
                            // NOTE!
                            sysVar.get_quenchingTs()[index].push_back(Ti);
                        }
                        ++count;
                    } while (count<sysVar.get_temperaturesInfo()[index].size());
                }
            }
        }
    }
    
    return is_retry;
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
    
    sysVar.c_vD           = tmpCoeffs;
    sysVar.rep            = rep;
    
    sysVar.coeffs_bndl_vD = coeffs_bndl_vD;
    sysVar.coeffs_bndu_vD = coeffs_bndu_vD;
    sysVar.coeffs_bndl    = coeffs_bndl;
    sysVar.coeffs_bndu    = coeffs_bndu;
    sysVar.epsf           = epsf;
    sysVar.epsx           = epsx;
    sysVar.maxits         = maxits;
    sysVar.diffstep       = diffstep;
}





void FitData::load_fit_coeffs(const StructureClass& sysVar)
{
    // load saved well-fit KWW parameters
    
    set_coeffs(sysVar.c_vD);
    set_coeffs_scale(sysVar.c_vD);
    
    coeffs_bndl_vD = sysVar.coeffs_bndl_vD;
    coeffs_bndu_vD = sysVar.coeffs_bndu_vD;
    coeffs_bndl    = sysVar.coeffs_bndl;
    coeffs_bndu    = sysVar.coeffs_bndu;
    epsf           = sysVar.epsf;
    epsx           = sysVar.epsx;
    maxits         = sysVar.maxits;
    diffstep       = sysVar.diffstep;
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
    //cout << o << "\n";
    ofstream write_fit_Arrhenius(o.c_str());
    
    if (write_fit_Arrhenius.is_open())
    {
        //find_cutTforArrhenius(sysVar,n_sys);
        set_Theq(get_Theq(sysVar,n_sys));
        
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
        write_tauFitfit_vs_T(write_fit_Arrhenius,c_org,model);
        write_tauFit_vs_T(write_fit_Arrhenius);
        write_tauFitfit_vs_invT(write_fit_Arrhenius,c_org,model);
        write_tauFit_vs_invT(write_fit_Arrhenius);
        
        avgfitParams(sysVar,n_trl,n_sys,model);
        
        if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
        {
            // last index
            write_fitavg(sysVar,n_sys,model);
            write_TA(sysVar,n_sys);
        }
        
        write_fit_Arrhenius.close();
    }
    else {
        cout << "fit_Arrhenius: 'fit_Arrhenius.dat' cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
    }
}





void FitData::fit_tauFit(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys,
                         const vector<string>& relmodel)
{
    string model=relmodel[0];
    string extrp=relmodel[1];
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    if (index==0)
    {
        ModelCoeffs.clear();
        ExtrpCoeffs.clear();
        r2_Model.clear();
        
        /* change dynamics range based on tau0_Arrhenius */
        if (sysVar.get_is_useDynamicRange())
        {
            find_tau0(sysVar,n_sys,model);
            useDynamicsRange(sysVar,n_sys);
        }
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
    //cout << o << "\n";
    ofstream write_fit_tauFit(o.c_str());
    
    if (write_fit_tauFit.is_open())
    {
        set_Theq(get_Theq(sysVar,n_sys));
        set_Tleq(get_Tleq(sysVar,n_trl,n_sys));
        
        set_fitParams(model);
        read_individual_taueq_data(sysVar,n_trl,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit(model);
        
        real_1d_array c_org=c;
        r2_Model.push_back(rep.r2);
        vector<double> tmpCoeffs;
        for (size_t i=0; i<c.length(); ++i) tmpCoeffs.push_back(c[i]);
        ModelCoeffs.push_back(tmpCoeffs);
        compute_Tg_fragility(sysVar,c_org,model);
        
        write_fitModel(write_fit_tauFit,model);
        write_fitarrays(write_fit_tauFit);
        write_fitcutoff(write_fit_tauFit,sysVar,model);
        write_stopcond(write_fit_tauFit);
        write_Tg_fragility(write_fit_tauFit,c_org,model);
        write_fitinfo(write_fit_tauFit);
        write_errcurve(write_fit_tauFit,model);
        write_tauFitfit_vs_T(write_fit_tauFit,c_org,model);
        write_tauFit_vs_T(write_fit_tauFit);
        write_tauFitfit_vs_invT(write_fit_tauFit,c_org,model);
        write_tauFit_vs_invT(write_fit_tauFit);
        
        avgfitParams(sysVar,n_trl,n_sys,model);
        
        if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
        {
            // last index
            write_fitavg(sysVar,n_sys,model);
            write_Tc(sysVar,n_sys);
            write_tau0(sysVar,n_sys,model);
        }
        
        
        ////////////////////////////////////////////////////////////////////////
        //======================================================================
        if (get_is_shootForNext())
        {
            set_fitParams(extrp);
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            fitdata_processing_lsfit(taueq_sortdecreasing);
            alglib_lsfit(extrp);
            
            real_1d_array c_extrp=c;
            vector<double> tmpCoeffs;
            for (size_t i=0; i<c.length(); ++i) tmpCoeffs.push_back(c[i]);
            ExtrpCoeffs.push_back(tmpCoeffs);
            
            find_largest_Tg_fragility(sysVar,n_trl,n_sys,c_extrp,extrp);
            avgfitParams_extrp(sysVar,n_trl,n_sys,extrp);
            
            /* Shooting Algorithm Kernel */
            shootForNewTemperatures(sysVar,n_trl,n_sys,c_extrp,extrp);
        }
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        
        
        if (current_regime==(n_regime-1))
        {
            // final regime
            if (get_is_fitByEveryPoint())
            {
                if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg)))
                {
                    // last index
                    fit_continuousRange(sysVar,n_sys,model);
                    fit_neighboringNpoints(sysVar,n_sys,model);
                    
                    if (get_is_fit_by_Tc()||get_is_fit_by_TA())
                    {
                        fit_partial(sysVar,n_sys,model);
                    }
                }
            }
        }
        write_fit_tauFit.close();
    }
    else
    {
        cout << "FitData::fit_tauFit: 'fit_tauFit.dat' cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
    }
}





void FitData::alglib_solverobj(const alglib::real_2d_array& x,
                               const alglib::real_1d_array& y,
                               const alglib::real_1d_array& c,
                               alglib::lsfitstate& state)
{
    if (is_use_FG) {
        lsfitcreatefg(x, y, c, true, state);
    }
    else {
        lsfitcreatef(x, y, c, diffstep, state);
    }
}





void FitData::alglib_lsfit_form(const std::string& model,
                                alglib::lsfitstate& state)
{
    if (model=="KWW")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, KWW_func, KWW_grad);
        }
        else
        {
            alglib::lsfitfit(state, KWW_func);
        }
    }
    else if (model=="KWW_pwr")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, KWW_pwr_func, KWW_pwr_grad);
        }
        else
        {
            alglib::lsfitfit(state, KWW_pwr_func);
        }
    }
    else if (model=="mKWW")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, mKWW_func, mKWW_grad);
        }
        else
        {
            alglib::lsfitfit(state, mKWW_func);
        }
    }
    else if (model=="mKWW_pwr")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, mKWW_pwr_func, mKWW_pwr_grad);
        }
        else
        {
            alglib::lsfitfit(state, mKWW_pwr_func);
        }
    }
    else if (model=="linear")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, linear_func, linear_grad);
        }
        else
        {
            alglib::lsfitfit(state, linear_func);
        }
    }
    else if (model=="exp_string")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, exp_string_func, exp_string_grad);
        }
        else
        {
            alglib::lsfitfit(state, exp_string_func);
        }
    }
    else if (model=="pwr_exp_string")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, pwr_exp_string_func, pwr_exp_string_grad);
        }
        else
        {
            alglib::lsfitfit(state, pwr_exp_string_func);
        }
    }
    else if (model=="Arrhenius")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, Arrhenius_func, Arrhenius_grad);
        }
        else
        {
            alglib::lsfitfit(state, Arrhenius_func);
        }
    }
    else if (model=="VFT")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, VFT_func, VFT_grad);
        }
        else
        {
            alglib::lsfitfit(state, VFT_func);
        }
    }
    else if (model=="Mauro")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, Mauro_func, Mauro_grad);
        }
        else
        {
            alglib::lsfitfit(state, Mauro_func);
        }
    }
    else if (model=="AM")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, AM_func, AM_grad);
        }
        else
        {
            alglib::lsfitfit(state, AM_func);
        }
    }
    else if (model=="DG")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, DG_func, DG_grad);
        }
        else
        {
            alglib::lsfitfit(state, DG_func);
        }
    }
    else if (model=="MCT")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, MCT_func, MCT_grad);
        }
        else {
            alglib::lsfitfit(state, MCT_func);
        }
    }
    else if (model=="SOU")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, SOU_func, SOU_grad);
        }
        else
        {
            alglib::lsfitfit(state, SOU_func);
        }
    }
    else if (model=="ArrheniusII")
    {
        alglib::lsfitfit(state, ArrheniusII_func);
    }
    else if (model=="GLM")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, GLM_func, GLM_grad);
        }
        else
        {
            alglib::lsfitfit(state, GLM_func);
        }
    }
    else if (model=="GLM1")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, GLM1_func, GLM1_grad);
        }
        else
        {
            alglib::lsfitfit(state, GLM1_func);
        }
    }
    else if ((model=="COOP")||(model=="COOP_string")||(model=="COOP_DWF"))
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, COOP_func, COOP_grad);
        }
        else
        {
            alglib::lsfitfit(state, COOP_func);
        }
    }
    else if (model=="DEAG")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, DEAG_func, DEAG_grad);
        }
        else
        {
            alglib::lsfitfit(state, DEAG_func);
        }
    }
    else if (model=="CG")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, CG_func, CG_grad);
        }
        else
        {
            alglib::lsfitfit(state, CG_func);
        }
    }
    else if (model=="FourParamVFT")
    {
        if (is_use_FG)
        {
            alglib::lsfitfit(state, FourParamVFT_func, FourParamVFT_grad);
        }
        else
        {
            alglib::lsfitfit(state, FourParamVFT_func);
        }
    }
    else if (model=="ArrheniusIII")
    {
        alglib::lsfitfit(state, ArrheniusIII_func);
    }
}





void FitData::set_fitParams(const std::string& model)
{
    
    epsf     = 0;
    epsx     = 1e-10;
    maxits   = 1e+3;
    diffstep = 1e-5;
    
    
    if (model=="KWW") // [A,tau,beta]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1e+5, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0,  inf, 1.0});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1e+2, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0,  inf, 1.0});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0,  inf, 1.0});
        }
    }
    else if (model=="KWW_pwr") // [A,a,tau,beta]
    {
        if (systemUnit=="real")
        {
            set_coeffs({       1.0, 1.0, 1e+3, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+5, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 5.0,  inf, 1.0});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1.0, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 5.0,  inf, 1.0});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0, 1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 5.0,  inf, 1.0});
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
            set_coeffs_bndu({ 10.0, 1e+3,  inf, 1.0});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1e+3, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1e+3, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 1e+3,  inf, 1.0});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0,  1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0,  1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0,  1.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0,  1.0,  inf, 1.0});
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
            set_coeffs_bndu({ 10.0, 5.0, 1e+3,  inf, 1.0});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       1.0, 1.0, 1e+3, 1e+0, 0.5});
            set_coeffs_scale({ 1.0, 1.0, 1e+3, 1e+2, 0.5});
            set_coeffs_bndl({  0.0, 0.0, 1e+3,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 5.0, 1e+3,  inf, 1.0});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 1.0,  1.0,  1.0, 0.5});
            set_coeffs_scale({ 1.0, 1.0,  1.0,  1.0, 0.5});
            set_coeffs_bndl({  0.0, 0.0,  1.0,  0.0, 0.0});
            set_coeffs_bndu({ 10.0, 5.0,  1.0,  inf, 1.0});
        }
    }
    else if(model=="linear") // [slope,intercept]
    {
        set_coeffs({       1.0, 0.0});
        set_coeffs_scale({ 1.0, 1.0});
        set_coeffs_bndl({  1.0,-inf});
        set_coeffs_bndu({  1.0, inf});
        
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
            if (is_use_viscosity)
            {
                set_coeffs({       1e-6, 2.5e+3});
                set_coeffs_scale({ 1e-6,   1e+3});
                set_coeffs_bndl({   0.0,    0.0});
                set_coeffs_bndu({  1e+6,   1e+6});
            }
            else
            {
                set_coeffs({       2e+2, 2.5e+3});
                set_coeffs_scale({ 1e+2,   1e+3});
                set_coeffs_bndl({   0.0,    0.0});
                set_coeffs_bndu({  1e+6,   1e+6});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-6, 2.5e+3});
                set_coeffs_scale({ 1e-6,   1e+3});
                set_coeffs_bndl({   0.0,    0.0});
                set_coeffs_bndu({  1e+6,   1e+6});
            }
            else
            {
                set_coeffs({       2e+2, 2.5e+3});
                set_coeffs_scale({ 1e+2,   1e+3});
                set_coeffs_bndl({   0.0,    0.0});
                set_coeffs_bndu({  1e+6,   1e+6});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.1,  2.0});
            set_coeffs_scale({ 0.1,  1.0});
            set_coeffs_bndl({  0.0,  0.0});
            set_coeffs_bndu({ 1e+2, 1e+2});
        }
    }
    else if (model=="VFT") // [tau0,D,T0]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-6, 10.0, 2e+2});
                set_coeffs_scale({ 1e-6, 10.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3});
            }
            else
            {
                set_coeffs({       1e+2,  7.0, 1e+2});
                set_coeffs_scale({ 1e+2,  5.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-6, 10.0, 2e+2});
                set_coeffs_scale({ 1e-6, 10.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3});
            }
            else
            {
                set_coeffs({       1e+2,  7.0, 1e+2});
                set_coeffs_scale({ 1e+2,  5.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.2,  3.0,  0.3});
            set_coeffs_scale({ 0.1,  1.0,  0.1});
            set_coeffs_bndl({  0.0,  0.0,  0.0});
            set_coeffs_bndu({ 1e+2, 1e+3, 10.0});
        }
    }
    else if (model=="Mauro") // [tau0,K,C]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 8e+2, 5.5e+2});
                set_coeffs_scale({ 1e-5, 5e+2,   1e+2});
                set_coeffs_bndl({   0.0,  0.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,   1e+4});
            }
            else
            {
                set_coeffs({       6e+2, 8e+2, 5.5e+2});
                set_coeffs_scale({ 1e+2, 5e+2,   1e+2});
                set_coeffs_bndl({   0.0,  0.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,   1e+4});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 8e+2, 5.5e+2});
                set_coeffs_scale({ 1e-5, 5e+2,   1e+2});
                set_coeffs_bndl({   0.0,  0.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,   1e+4});
            }
            else {
                set_coeffs({       6e+2, 8e+2, 5.5e+2});
                set_coeffs_scale({ 1e+2, 5e+2,   1e+2});
                set_coeffs_bndl({   0.0,  0.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,   1e+4});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.3,  0.3,  1.0});
            set_coeffs_scale({ 0.1,  0.5,  1.0});
            set_coeffs_bndl({  0.0,  0.0,  0.0});
            set_coeffs_bndu({ 1e+2, 10.0, 10.0});
        }
    }
    else if (model=="AM") // [tau0,p,alpha]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 8.5e+2,  2.5});
                set_coeffs_scale({ 1e-5,   1e+3,  1.0});
                set_coeffs_bndl({   0.0,    0.0,  0.0});
                set_coeffs_bndu({  1e+6,   1e+4, 10.0});
            }
            else
            {
                set_coeffs({      1.3e+3, 8.5e+2,  2.5});
                set_coeffs_scale({  1e+3,   1e+3,  1.0});
                set_coeffs_bndl({    0.0,    0.0,  0.0});
                set_coeffs_bndu({   1e+6,   1e+4, 10.0});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 8.5e+2,  2.5});
                set_coeffs_scale({ 1e-5,   1e+3,  1.0});
                set_coeffs_bndl({   0.0,    0.0,  0.0});
                set_coeffs_bndu({  1e+6,   1e+4, 10.0});
            }
            else
            {
                set_coeffs({      1.3e+3, 8.5e+2,  2.5});
                set_coeffs_scale({  1e+3,   1e+3,  1.0});
                set_coeffs_bndl({    0.0,    0.0,  0.0});
                set_coeffs_bndu({   1e+6,   1e+4, 10.0});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.5,  0.5,  4.5});
            set_coeffs_scale({ 0.5,  0.5,  1.0});
            set_coeffs_bndl({  0.0,  0.0,  0.0});
            set_coeffs_bndu({ 1e+2, 10.0, 10.0});
        }
    }
    else if (model=="DG") // [tau0,A,B]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       2e+3, 2e+4,  5e-3});
                set_coeffs_scale({ 1e+3, 1e+4,  1e-3});
                set_coeffs_bndl({   0.0,  0.0, -1e+3});
                set_coeffs_bndu({  1e+6,  inf,  1e+3});
            }
            else
            {
                set_coeffs({       2e+3, 2e+4,  5e-3});
                set_coeffs_scale({ 1e+3, 1e+4,  1e-3});
                set_coeffs_bndl({   0.0,  0.0, -1e+3});
                set_coeffs_bndu({  1e+6,  inf,  1e+3});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       2e+3, 2e+4,  5e-3});
                set_coeffs_scale({ 1e+3, 1e+4,  1e-3});
                set_coeffs_bndl({   0.0,  0.0, -1e+3});
                set_coeffs_bndu({  1e+6,  inf,  1e+3});
            }
            else
            {
                set_coeffs({       2e+3, 2e+4,  5e-3});
                set_coeffs_scale({ 1e+3, 1e+4,  1e-3});
                set_coeffs_bndl({   0.0,  0.0, -1e+3});
                set_coeffs_bndu({  1e+6,  inf,  1e+3});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.5, 1.0,   1.0});
            set_coeffs_scale({ 1.0, 1.0,   1.0});
            set_coeffs_bndl({  0.0, 0.0, -1e+3});
            set_coeffs_bndu({ 1e+2, inf,  1e+3});
        }
    }
    else if (model=="MCT") // [A,Tc,r]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       1e+4, 2.5e+2, 2.5});
                set_coeffs_scale({ 1e+4,   1e+2, 1.0});
                set_coeffs_bndl({   0.0,    0.0, 0.0});
                set_coeffs_bndu({   inf,    inf, inf});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       1e+4, 2.5e+2, 2.5});
                set_coeffs_scale({ 1e+4,   1e+2, 1.0});
                set_coeffs_bndl({   0.0,    0.0, 0.0});
                set_coeffs_bndu({   inf,    inf, inf});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 0.5, 1.0});
            set_coeffs_scale({ 0.1, 0.1, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 0.0});
            set_coeffs_bndu({  inf, inf, inf});
        }
    }
    else if (model=="SOU") // [A,Tc,r]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       1e+4, 2.5e+2, 2.5});
                set_coeffs_scale({ 1e+4,   1e+2, 1.0});
                set_coeffs_bndl({   0.0,    0.0, 0.0});
                set_coeffs_bndu({   inf,    inf, inf});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       1e+4, 2.5e+2, 2.5});
                set_coeffs_scale({ 1e+4,   1e+2, 1.0});
                set_coeffs_bndl({   0.0,    0.0, 0.0});
                set_coeffs_bndu({   inf,    inf, inf});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.1, 0.3, 5.0});
            set_coeffs_scale({ 0.1, 0.1, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 0.0});
            set_coeffs_bndu({  inf, inf, inf});
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
            set_coeffs({       0.3, 1.0, 3.0});
            set_coeffs_scale({ 0.1, 1.0, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 0.0});
            set_coeffs_bndu({  inf, inf, inf});
        }
    }
    else if (model=="GLM") // [tau0,u0^2,a]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       2e+4, 3.0, 3.0});
                set_coeffs_scale({ 1e+4, 1.0, 1.0});
                set_coeffs_bndl({   0.0, 0.0, 3.0});
                set_coeffs_bndu({   inf, inf, inf});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs({       20.0, 3.0, 3.0});
                set_coeffs_scale({ 10.0, 1.0, 1.0});
                set_coeffs_bndl({   0.0, 0.0, 3.0});
                set_coeffs_bndu({   inf, inf, inf});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.5, 0.1, 3.0});
            set_coeffs_scale({ 1.0, 0.1, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 3.0});
            set_coeffs_bndu({  inf, inf, inf});
        }
    }
    else if (model=="GLM1") // [tauA,uA^2,a]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs(      {tauA, uA2,  3.0});
                set_coeffs_scale({tauA, uA2,  1.0});
                set_coeffs_bndl( {tauA, uA2, -inf});
                set_coeffs_bndu( {tauA, uA2,  inf});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                //
            }
            else
            {
                set_coeffs(      {tauA, uA2,  3.0});
                set_coeffs_scale({tauA, uA2,  1.0});
                set_coeffs_bndl( {tauA, uA2, -inf});
                set_coeffs_bndu( {tauA, uA2,  inf});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs(      {tauA, uA2,  3.0});
            set_coeffs_scale({tauA, uA2,  1.0});
            set_coeffs_bndl( {tauA, uA2, -inf});
            set_coeffs_bndu( {tauA, uA2,  inf});
        }
    }
    else if (model=="COOP") // [tau0,E_inf,u,b]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-4, 1e+3,  1.0,  0.1});
                set_coeffs_scale({ 1e-4, 1e+3,  1.0,  0.1});
                set_coeffs_bndl({  1e-7,  0.0,-5e+2,-1e+3});
                set_coeffs_bndu({   inf, 1e+6, 5e+2, 1e+3});
            }
            else
            {
                set_coeffs({       10.0, 1e+3,  1.0,  0.1});
                set_coeffs_scale({ 10.0, 1e+3,  1.0,  0.1});
                set_coeffs_bndl({   1.0,  0.0,-5e+2,-1e+3});
                set_coeffs_bndu({   inf, 1e+6, 5e+2, 1e+3});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-4, 1e+3,  1.0,  0.1});
                set_coeffs_scale({ 1e-4, 1e+3,  1.0,  0.1});
                set_coeffs_bndl({  1e-7,  0.0,-5e+2,-1e+3});
                set_coeffs_bndu({   inf, 1e+6, 5e+2, 1e+3});
            }
            else
            {
                set_coeffs({       1.0, 1e+3,  1.0,  0.1});
                set_coeffs_scale({ 1.0, 1e+3,  1.0,  0.1});
                set_coeffs_bndl({  0.0,  0.0,-5e+2,-1e+5});
                set_coeffs_bndu({  inf, 1e+6, 5e+2, 1e+5});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.5, 0.5,  1.0,  0.5});
            set_coeffs_scale({ 0.1, 0.5,  1.0,  0.1});
            set_coeffs_bndl({  0.0, 0.0,-5e+2,-1e+3});
            set_coeffs_bndu({  inf,1e+3, 5e+2, 1e+3});
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
            set_coeffs({       0.5, 1.0, 10.0,  0.5});
            set_coeffs_scale({ 1.0, 1.0, 10.0,  0.1});
            set_coeffs_bndl({  0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({  inf, inf,  inf,  inf});
        }
        else if (systemUnit=="metal")
        {
            set_coeffs({       0.5, 1.0, 10.0,  0.5});
            set_coeffs_scale({ 1.0, 1.0, 10.0,  0.1});
            set_coeffs_bndl({  0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({  inf, inf,  inf,  inf});
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.5, 1.0, 10.0,  0.5});
            set_coeffs_scale({ 1.0, 1.0, 10.0,  0.1});
            set_coeffs_bndl({  0.0, 0.0, -inf, -inf});
            set_coeffs_bndu({  inf, inf,  inf,  inf});
        }
    }
    else if (model=="DEAG") // [tau0,A,B,C]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 1e+2, -1.0, 5.5e+2});
                set_coeffs_scale({ 1e-5, 1e+2, 1.0, 1e+2});
                set_coeffs_bndl({   0.0, 0.0,  -1e+3, 0.0});
                set_coeffs_bndu({  1e+6, 1e+4, 1e+3, 1e+4});
            }
            else
            {
                set_coeffs({       2e+2, 1e+2,  -1.0, 5.5e+2});
                set_coeffs_scale({ 1e+2, 1e+2,   1.0,   1e+2});
                set_coeffs_bndl({   0.0,  0.0, -10.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  10.0,   1e+4});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 1e+2, -1.0, 5.5e+2});
                set_coeffs_scale({ 1e-5, 1e+2, 1.0, 1e+2});
                set_coeffs_bndl({   0.0, 0.0,  -1e+3, 0.0});
                set_coeffs_bndu({  1e+6, 1e+4, 1e+3, 1e+4});
            }
            else
            {
                set_coeffs({       2e+2, 1e+2,  -1.0, 5.5e+2});
                set_coeffs_scale({ 1e+2, 1e+2,   1.0,   1e+2});
                set_coeffs_bndl({   0.0,  0.0, -10.0,    0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  10.0,   1e+4});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       1.0, 1e-2,  -1.0,  1.5});
            set_coeffs_scale({ 1.0, 1e-2,   1.0,  1.0});
            set_coeffs_bndl({  0.0,  0.0, -10.0,  0.0});
            set_coeffs_bndu({  inf, 10.0,  10.0, 10.0});
        }
    }
    else if (model=="CG") // [tau0,B,C,T0]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 5e+2,   5.0, 2e+2});
                set_coeffs_scale({ 1e-5, 1e+3,   1.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0, -1e+3,  0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  1e+3, 1e+3});
            }
            else
            {
                set_coeffs({       3e+2, 3e+3,   0.5, 2e+2});
                set_coeffs_scale({ 1e+2, 1e+3,   0.1, 1e+2});
                set_coeffs_bndl({   0.0,  0.0, -1e+3,  0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  1e+3, 1e+3});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 5e+2,   5.0, 2e+2});
                set_coeffs_scale({ 1e-5, 1e+3,   1.0, 1e+2});
                set_coeffs_bndl({   0.0,  0.0, -1e+3,  0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  1e+3, 1e+3});
            }
            else
            {
                set_coeffs({       3e+2, 3e+3,   0.5, 2e+2});
                set_coeffs_scale({ 1e+2, 1e+3,   0.1, 1e+2});
                set_coeffs_bndl({   0.0,  0.0, -1e+3,  0.0});
                set_coeffs_bndu({  1e+6, 1e+4,  1e+3, 1e+3});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.2,  2.0,   0.0,  0.3});
            set_coeffs_scale({ 0.1,  1.0,  1e-2,  0.1});
            set_coeffs_bndl({  0.0,  0.0, -1e+3,  0.0});
            set_coeffs_bndu({ 1e+2, 10.0,  1e+3, 10.0});
        }
    }
    else if (model=="FourParamVFT") // [tau0,D,T0,alpha]
    {
        if (systemUnit=="real")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 1e+1, 1e+2, 10.0});
                set_coeffs_scale({ 1e-6, 1e+1, 1e+2,  1.0});
                set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3, 10.0});
            }
            else
            {
                set_coeffs({       1e+5,  5.0, 1e+2,  2.0});
                set_coeffs_scale({ 1e+5,  1.0, 1e+2,  1.0});
                set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3, 10.0});
            }
        }
        else if (systemUnit=="metal")
        {
            if (is_use_viscosity)
            {
                set_coeffs({       1e-5, 1e+1, 1e+2, 10.0});
                set_coeffs_scale({ 1e-6, 1e+1, 1e+2,  1.0});
                set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3, 10.0});
            }
            else
            {
                set_coeffs({       1e+5,  5.0, 1e+2,  2.0});
                set_coeffs_scale({ 1e+5,  1.0, 1e+2,  1.0});
                set_coeffs_bndl({   0.0,  0.0,  0.0,  0.0});
                set_coeffs_bndu({  1e+6, 1e+3, 1e+3, 10.0});
            }
        }
        else if (systemUnit=="lj")
        {
            set_coeffs({       0.2,  3.0,  0.3,  1.0});
            set_coeffs_scale({ 0.1,  5.0,  0.1,  1.0});
            set_coeffs_bndl({  0.0,  0.0,  0.0,  0.0});
            set_coeffs_bndu({ 1e+2, 1e+3, 10.0, 10.0});
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
            set_coeffs({       0.3, 1.0, 3.0, 1.0});
            set_coeffs_scale({ 0.1, 1.0, 1.0, 1.0});
            set_coeffs_bndl({  0.0, 0.0, 0.0, 0.0});
            set_coeffs_bndu({  inf, inf, inf, inf});
        }
    }
}





void FitData::fitParams_correction(const std::string& model)
{
    if (model=="KWW")
    {
        // [A,tau,beta]
        
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
    else if (model=="KWW_pwr")
    {
        // [A,a,tau,beta]
        
        // self-correct tau:
        // 1.15^100=1.17e+6
        
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
        // [A,taulib,tau,beta]
        
        // self-correct tau:
        // 1.15^100=1.17e+6
        
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
        // [A,a,taulib,tau,beta]
        
        // self-correct tau:
        // 1.15^100=1.17e+6
        
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
    
    // gradient (wrt c):
    // verified by Wolfram-Alpha
    grad[0] = func/A;
    grad[1] = func*(beta/tau)*pow(time/tau,beta);
    grad[2] = -func*log(time/tau)*pow(time/tau,beta);
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
    
    if (tcalc=="KWW") {
        /*======================================================================
         F = A*exp(-(t/tau)^b)
         t = tau*(ln(A/F))^(1/b) -- analytic solution for t
         
         c[0]=A, c[1]=tau, c[2]=b
         =====================================================================*/
        double A=c[0],tau=c[1],b=c[2];
        tauFit=tau*pow(log(A/get_sExp_tauFit()),pow(b,-1));
    }
    else if (tcalc=="KWW_pwr") {
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
    else if (tcalc=="mKWW") {
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
    else if (tcalc=="mKWW_pwr") {
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
    else {
        cout
        << "in FitData::sExp_tauFit_interpolateFvalue, " << tcalc << " model not found!" << "\n"
        << "Program aborted." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    return tauFit;
}





double FitData::sExp_tauFit_gammafunction(const real_1d_array& c,
                                          const string& tcalc)
{
    double tauFit=0;
    
    if (tcalc=="KWW") {
        /*======================================================================
         F = A*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=tau_kww, c[2]=beta_kww
         =====================================================================*/
        double tau=c[1],beta=c[2];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="KWW_pwr") {
        /*======================================================================
         F = A*(t^-a)*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=a, c[2]=tau_kww, c[3]=beta_kww
         =====================================================================*/
        double tau=c[2],beta=c[3];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="mKWW") {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
         =====================================================================*/
        double tau=c[2],beta=c[3];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else if (tcalc=="mKWW_pwr") {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
         t = (tau_kww/beta_kww)*Gammafunc(1/beta_kww)
         
         c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
         =====================================================================*/
        double tau=c[3],beta=c[4];
        double sgngam=0;
        tauFit=(tau/beta)*exp(alglib::lngamma(pow(beta,-1),sgngam));
    }
    else {
        cout
        << "in FitData::sExp_tauFit_gammafunction, " << tcalc << " model not found!" << "\n"
        << "Program aborted." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    return tauFit;
}





double FitData::sExp_func_value(const real_1d_array& c,
                                const double time,
                                const string& tcalc)
{
    double Fvalue=0;
    
    if (tcalc=="KWW") {
        /*======================================================================
         F = A*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=tau, c[2]=beta
         =====================================================================*/
        double A=c[0],tau=c[1],beta=c[2];
        Fvalue=A*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="KWW_pwr") {
        /*======================================================================
         F = A*(t^-a)*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=a, c[2]=tau, c[3]=beta
         =====================================================================*/
        double A=c[0],a=c[1],tau=c[2],beta=c[3];
        Fvalue=A*pow(time,-a)*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="mKWW") {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=taulib, c[2]=tau, c[3]=beta
         =====================================================================*/
        double A=c[0],taulib=c[1],tau=c[2],beta=c[3];
        Fvalue=(1-A)*exp(-time/taulib)+A*exp(-pow(time/tau,beta));
    }
    else if (tcalc=="mKWW_pwr") {
        /*======================================================================
         F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
         
         c[0]=A, c[1]=a, c[2]=taulib, c[3]=tau, c[4]=beta
         =====================================================================*/
        double A=c[0],a=c[1],taulib=c[2],tau=c[3],beta=c[4];
        Fvalue=(1-A)*exp(-time/taulib)+A*pow(time,-a)*exp(-pow(time/tau,beta));
    }
    else {
        cout
        << "FitData::in sExp_func_value, " << tcalc << " model not found!" << "\n"
        << "Program aborted." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
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
    double slope=c[0],intercept=c[1],dependet=x[0];
    func = slope*dependet+intercept;
}
void FitData::linear_grad(const alglib::real_1d_array &c,
                          const alglib::real_1d_array &x,
                          double &func,
                          alglib::real_1d_array &grad,
                          void *ptr)
{
    double slope=c[0],intercept=c[1],dependet=x[0];
    func = slope*dependet+intercept;
    grad[0]=dependet;
    grad[1]=1.0;
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
     tau = tau0*exp((u0^2/u^2)^(a/2))
     param: {c[0]=tau0, c[1]=u0^2, c[2]=a}
     
     func: (with base 10)
     func = log10(tau) = (ln(tau0)+(u0^2/u^2)^(a/2))*log10(e)
     =========================================================================*/
    double tau0=c[0],u02=c[1],a=c[2],u2=x[0];
    func = (log(tau0)+pow(u02/u2,a/2))*log10(exp(1));
}
void FitData::GLM_grad(const alglib::real_1d_array &c,
                       const alglib::real_1d_array &x,
                       double &func,
                       alglib::real_1d_array &grad,
                       void *ptr)
{
    double tau0=c[0],u02=c[1],a=c[2],u2=x[0];
    func = (log(tau0)+pow(u02/u2,a/2))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= ((a/(2*u02))*pow(u02/u2,a/2))*log10(exp(1));
    grad[2]= (0.5*pow(u02/u2,a/2)*log(u02/u2))*log10(exp(1));
}





void FitData::GLM1_func(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        void *ptr)
{
    /*==========================================================================
     Generalized Localization Model:
     tau = tauA*exp((uA^2/u^2)^(a/2)-1)
     param: {c[0]=tauA, c[1]=uA^2, c[2]=a}
     
     func: (with base 10)
     func = log10(tau) = (ln(tauA)+(uA^2/u^2)^(a/2)-1)*log10(e)
     =========================================================================*/
    double tauA=c[0],uA2=c[1],a=c[2],u2=x[0];
    func = (log(tauA)+pow(uA2/u2,a/2)-1)*log10(exp(1));
}
void FitData::GLM1_grad(const alglib::real_1d_array &c,
                        const alglib::real_1d_array &x,
                        double &func,
                        alglib::real_1d_array &grad,
                        void *ptr)
{
    double tauA=c[0],uA2=c[1],a=c[2],u2=x[0];
    func = (log(tauA)+pow(uA2/u2,a/2)-1)*log10(exp(1));
    grad[0]= (pow(tauA,-1))*log10(exp(1));
    grad[1]= ((a/(2*uA2))*pow(uA2/u2,a/2))*log10(exp(1));
    grad[2]= (0.5*pow(uA2/u2,a/2)*log(uA2/u2))*log10(exp(1));
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





double FitData::calc_time_given_temp(const real_1d_array& c,
                                     const double x,
                                     const string& model)
{
    if (model=="linear") {
        // y(x) = slope*x + intercept
        // param: {c[0]=slope, c[1]=intercept}
        return c[0]*x+c[1];
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
    else if (model=="Arrhenius") {
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
    else if (model=="GLM") {
        // tau = tau0*exp((u0^2/u^2)^(a/2))
        // param: {c[0]=tau0, c[1]=u0^2, c[2]=a}
        return c[0]*exp(pow(c[1]/x,c[2]/2.0));
    }
    else if (model=="GLM1") {
        // tau = tauA*exp((uA^2/u^2)^(a/2)-1)
        // param: {c[0]=tauA, c[1]=uA^2, c[2]=a}
        return c[0]*exp(pow(c[1]/x,c[2]/2)-1);
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
        return 0;
}





vector<double> FitData::calc_temp_given_time(const real_1d_array& c,
                                             const double time,
                                             const string& model)
{
    
    //==========================================================================
    // NOTE:
    //
    // 1st element: Tg (@ extrapolated time) || T (@ any time)
    // 2nd element: m  (@ Tg)
    //==========================================================================
    vector<double> glass_data;
    
    double lnTime=log(time/c[0]); // F=ln(time/tau0)
    
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
    
    return glass_data;
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
        return 0;
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
        return 0;
    return 0;
}





vector<double> FitData::get_interpolated_Ts(StructureClass& sysVar,
                                            const int n_trl,
                                            const int n_sys,
                                            const int n_temps,
                                            const double Tl,
                                            const double Th)
{
    vector<double> originalTs=sysVar.get_temperaturesInfo()[index];
    vector<double> simtemps; // simulated temperatures
    vector<double> newtemperatures,hightolow;
    
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
    }
    else {
        cout
        << "FitData::get_interpolated_Ts: readtempinfo cannot open."
        << "\n";
    }
    
    double dT=fabs(Tl-Th);
    double Ti=min(Tl,Th);
    
    for (int i=0; i<n_temps; ++i) {
        
        ////////////////////////////////////////////////////////////////////////
        /*----------------------------------------------------------------------
         check if Ti repeats with any existing T's;
         if found any that repeats, Ti is "decreased by one temperature unit".
         ---------------------------------------------------------------------*/
        bool is_foundrepeat=false;
        do {
            vector<double>::iterator itr0,itr1,itr2;
            
            itr0=find(originalTs.begin(),originalTs.end(),Ti);
            itr1=find(simtemps.begin(),simtemps.end(),Ti);
            itr2=find(newtemperatures.begin(),newtemperatures.end(),Ti);
            
            // NOTE:
            //------------------------------------------------------------------
            // "find" returns iterators in the range: [first,last); it would
            // return a pointer to the "last" element if nothing is found.
            //------------------------------------------------------------------
            if ((itr0!=originalTs.end())||
                (itr1!=simtemps.end())||
                (itr2!=newtemperatures.end())) {
                
                Ti -= precision;
                is_foundrepeat=true;
            }
            else {
                is_foundrepeat=false;
            }
            
        } while (is_foundrepeat);
        
        newtemperatures.push_back(Ti); // sotre new T's in this container
        //----------------------------------------------------------------------
        ////////////////////////////////////////////////////////////////////////
        
        /* Decide Temperature Spacing */
        //----------------------------------------------------------------------
        if (get_shootratio()>0) {
            
            /* Ti's are spaced between Th and Tl based on shootratio */
            
            dT *= get_shootratio();
            
            // NOTE:
            // this algorithm makes extrp. T's more crowded
            // in higher temperature end that the temperature spacing (dT)
            // is progressively shrinking in each iteration.
        }
        else {
            
            /* Ti's are evenly spaced between Th and Tl */
            
            dT = fabs(Tl-Th)/(double)(n_temps);
            
            // NOTE: n_temps >= 2
        }
        //----------------------------------------------------------------------
        
        /* New T's of next regime are added from low to high */
        //----------------------------------------------------------------------
        Ti = std::trunc(Ti+dT);
        //----------------------------------------------------------------------
    }
    
    // NOTE:
    //--------------------------------------------------------------------------
    // reverse the order of new temperatures;
    // by convention, T's in temperature containers are from high to low
    //--------------------------------------------------------------------------
    for (int i=(int)(newtemperatures.size()-1); i>=0; --i) {
        hightolow.push_back(newtemperatures[i]);
    }
    
    return hightolow;
}





void FitData::find_individual_msd_data(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d)
{
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/msd/");
    i.append(get_analysispart());
    i.append("/msd_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        
        string lineContent;
        string strData;
        double time=0,msd=0;
        vector<double> vec_time,vec_msd;
        
        /* read 1st line */
        getline(readFile,lineContent);
        amdat_version = lineContent;
        
        int current_index=0; // zero-based
        double t_hi=0,t_lo=0;
        double msd_hi=0,msd_lo=0;
        double slope=0;
        
        while (getline(readFile,lineContent)) {
            
            istringstream iss(lineContent);
            iss >> time; // time
            iss >> msd;  // msd (length^2)
            vec_time.push_back(time);
            vec_msd.push_back(msd);
            
            if(time>DWF_time) {
                
                if (current_index==0) current_index=1;
                
                t_lo   = vec_time[current_index-1];
                t_hi   = vec_time[current_index];
                msd_lo = vec_msd[current_index-1];
                msd_hi = vec_msd[current_index];
                slope  = (msd_hi-msd_lo)/(t_hi-t_lo);
                DWF    = slope*(DWF_time-t_lo)+msd_lo;
                break;
            }
            
            ++current_index;
        }
        readFile.close();
    }
    else {
        cout
        << "FitData::find_individual_msd_data: msd file cannot open."
        << "\n";
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/DWF_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_DWF(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/DWF_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_DWFinv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T DWF
    streamsize ss=write_DWF.precision();
    
    write_DWF
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_DWF.precision(ss);
    write_DWF << resetiosflags(ios::fixed|ios::showpoint);
    
    write_DWF
    << DWF
    << "\n";
    
    write_DWFinv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << DWF
    << "\n";
    
    write_DWF.close();
    write_DWFinv.close();
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
    //cout << i << "\n";
    ifstream readDWFdata(i.c_str());
    
    if (readDWFdata.is_open()) {
        string lineContent;
        double Tdata=0,DWFData=0;
        while (getline(readDWFdata,lineContent)) {
            istringstream iss(lineContent); // (T,DWF)
            iss >> Tdata;
            iss >> DWFData;
            DWF.push_back({Tdata,DWFData});
        }
        readDWFdata.close();
    }
    else {
        cout
        << "FitData::read_individual_DWF_data: DWF_all file cannot open."
        << "\n";
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_DWF_data:" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
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
    //cout << i << "\n";
    ifstream readDWFdata(i.c_str());
    
    if (readDWFdata.is_open()) {
        string lineContent;
        double Tdata=0,DWFData=0;
        while (getline(readDWFdata,lineContent)) {
            istringstream iss(lineContent); // (T,DWF)
            iss >> Tdata;
            iss >> DWFData;
            DWF.push_back({Tdata,DWFData});
        }
        readDWFdata.close();
    }
    else {
        cout
        << "FitData::read_individual_equ_DWF: DWF_all file cannot open."
        << "\n";
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_equ_DWF:" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: individual */
    //-----------------------------------
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int strSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==dwf_sortdecreasing[ii][0]) {
                taue=taueq_sortdecreasing[i][1];
                temp=dwf_sortdecreasing[ii][0];
                dwfr=dwf_sortdecreasing[ii][1];
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
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/DWF_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readDWFdata(i.c_str());
        
        if (readDWFdata.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readDWFdata,lineContent)) {
                istringstream iss(lineContent); // (T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                DWF.push_back({Tdata,DWFData});
            }
            readDWFdata.close();
        }
        else {
            cout
            << "read_all_DWF: DWF_all file cannot open."
            << "\n";
        }
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_DWF:" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
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
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/DWF_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readDWFdata(i.c_str());
        
        if (readDWFdata.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readDWFdata,lineContent)) {
                istringstream iss(lineContent); // (T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                DWF.push_back({Tdata,DWFData});
            }
            readDWFdata.close();
        }
        else {
            cout
            << "FitData::read_all_equ_DWF: DWF_all file cannot open."
            << "\n";
        }
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_string_equ_length:" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int strSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==dwf_sortdecreasing[ii][0]) {
                taue=taueq_sortdecreasing[i][1];
                temp=dwf_sortdecreasing[ii][0];
                dwfr=dwf_sortdecreasing[ii][1];
                DWF.push_back({temp,dwfr});
                DWF_invT.push_back({pow(temp,-1),dwfr});
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
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/DWF_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readDWFdata(i.c_str());
        
        if (readDWFdata.is_open()) {
            string lineContent;
            double Tdata=0,DWFData=0;
            while (getline(readDWFdata,lineContent)) {
                istringstream iss(lineContent); // (T,DWF)
                iss >> Tdata;
                iss >> DWFData;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    DWF.push_back({Tdata,DWFData});}
            }
            readDWFdata.close();
        }
        else {
            cout
            << "FitData::read_all_equ_DWF: DWF_all file cannot open."
            << "\n";
        }
    }
    
    try {
        if(DWF.size()==0) throw 0;
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_string_equ_length:" << "\n"
        << "DWF size = " << i << "\n";
        DWF.push_back({0,0});
        std::sort(DWF.begin(),DWF.end(),sortDecreasing);
        averaging_sorted_data(DWF,dwf_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys,TA_d);
    
    int n_data =(int)taueq_sortdecreasing.size();
    int strSize=(int)dwf_sortdecreasing.size();
    double temp=0;
    double taue=0;
    double dwfr=0;
    
    DWF.clear();
    
    for (int i=0; i<n_data; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==dwf_sortdecreasing[ii][0]) {
                taue=taueq_sortdecreasing[i][1];
                temp=dwf_sortdecreasing[ii][0];
                dwfr=dwf_sortdecreasing[ii][1];
                DWF.push_back({temp,dwfr});
                DWF_invT.push_back({pow(temp,-1),dwfr});
                dwf_tau.push_back({dwfr,taue});
                break;
            }
        }
    }
    dwf_sortdecreasing=DWF;
    dwf_invT_sortdecreasing=DWF_invT;
    dwf_tau_data=dwf_tau;
}





void FitData::write_DWF_equ(const StructureClass& sysVar)
{
    string o;
    
    /* write out "equilibrium" DWF to file (similar to taueq file) */
    //----------------------------------------------------------------------
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
            //cout << o << "\n";
            ofstream write_DWF(o.c_str(),ofstream::app); // 'append'
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/DWF_equ_inverseT_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_DWFinv(o.c_str(),ofstream::app); // 'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            read_individual_DWF(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_data =(int)taueq_sortdecreasing.size();
            int strSize=(int)dwf_sortdecreasing.size();
            
            for (int i=0; i<n_data; ++i) {
                for (int ii=0; ii<strSize; ++ii) {
                    if (taueq_sortdecreasing[i][0]==dwf_sortdecreasing[ii][0]) {
                        
                        write_DWF
                        << dwf_sortdecreasing[ii][0] << " "
                        << dwf_sortdecreasing[ii][1] << "\n";
                        write_DWFinv
                        << (1000.0/(corrFacforT/precision))/dwf_sortdecreasing[ii][0] << " "
                        << dwf_sortdecreasing[ii][1] << "\n";
                        break;
                    }
                }
            }
            write_DWF.close();
            write_DWFinv.close();
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
        //cout << o << "\n";
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
        int strSize=(int)dwf_sortdecreasing.size();
        
        for (int i=0; i<n_data; ++i) {
            for (int ii=0; ii<strSize; ++ii) {
                if (taueq_sortdecreasing[i][0]==dwf_sortdecreasing[ii][0]) {
                    
                    write_DWF
                    << dwf_sortdecreasing[ii][0] << " "
                    << dwf_sortdecreasing[ii][1] << "\n";
                    write_DWFinv
                    << (1000.0/(corrFacforT/precision))/dwf_sortdecreasing[ii][0] << " "
                    << dwf_sortdecreasing[ii][1] << "\n";
                    break;
                }
            }
        }
        write_DWF.close();
        write_DWFinv.close();
    }
}





void FitData::skim_individual_correlation_data(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d)
{
    /** read the target relaxtion profile **/
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/");
    i.append(relaxation_target+"/");
    i.append(get_analysispart());
    i.append("/"+relaxation_target+"_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    if(relaxation_target=="isfs") {
        i.append("."+get_analysispart()+".dat");
    } else if (relaxation_target=="baf") {
        i.append(".segmental.dat");
    }
    //cout << i << "\n";
    
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        
        string lineContent;
        string strData;
        double time=0,corr=0;
        
        /* read 1st line */
        getline(readFile,lineContent);
        
        /* read 2nd line */
        getline(readFile,lineContent);
        
        while (getline(readFile,lineContent)) {
            
            istringstream iss(lineContent);
            iss >> time; // time
            iss >> corr; // Fs(k,t)
            
            if (time>sExp_tau_lo) {
                
                if (corr>sExp_cut_hi_org) { // NOTE: sExp_cut_hi_org
                    is_use_sExp_cut_hi=false;
                    is_use_plateautime=true;
                    sExp_cut_hi=sExp_cut_hi_org;
                }
                else {
                    is_use_sExp_cut_hi=true;
                    is_use_plateautime=false;
                    //sExp_cut_hi=0.8; // NOTE!
                    sExp_cut_hi=sExp_cut_hi_org;
                }
                return;
            }
        }
        readFile.close();
    }
    else {
        cout
        << "FitData::skim_individual_correlation_data: correlation file cannot open."
        << "\n";
    }
}





void FitData::read_individual_correlation_data(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d)
{
    correlation_original.clear();
    correlation_sortincreasing.clear();
    vector<vector<double>> correlation;
    
    /** read the target relaxtion profile **/
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/");
    i.append(relaxation_target+"/");
    i.append(get_analysispart());
    i.append("/"+relaxation_target+"_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    if(relaxation_target=="isfs") {
        i.append("."+get_analysispart()+".dat");
    } else if (relaxation_target=="baf") {
        i.append(".segmental.dat");
    }
    //cout << i << "\n";
    
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        
        string lineContent;
        string strData;
        double time=0,corr=0;
        
        //auto start_position=readFile.tellg(); /* mark the beginning of file */
        //readFile.clear(); /* reset state flag of read file */
        //readFile.seekg(start_position); /* rewind to start of file */
        
        /* read 1st line */
        getline(readFile,lineContent);
        amdat_version = lineContent;
        
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
        
        while (getline(readFile,lineContent)) {
            
            ++cut_index_hi;
            ++cut_index_lo;
            ++cut_index_total;
            
            istringstream iss(lineContent);
            iss >> time; // time
            iss >> corr; // Fs(k,t)
            
            correlation_original.push_back({time,corr});
            
            if (is_use_sExp_cut_hi) {
                if (corr<=sExp_cut_hi) {
                    sExp_cut_index_hi = cut_index_hi;
                    --cut_index_hi;
                    // NOTE:
                    // cut_index_hi fixed at where the higher cutoff is
                }
            }
            
            if (is_use_plateautime) {
                if (time>sExp_tau_lo) {
                    sExp_cut_index_hi = cut_index_hi;
                    --cut_index_hi;
                    // NOTE:
                    // cut_index_hi fixed at where the higher cutoff is
                }
            }
            
            if (corr<sExp_cut_lo) {
                ++first_passage;
                if (first_passage==1) {
                    sExp_cut_index_lo = cut_index_lo-1;
                    // NOTE:
                    // cut_index_lo fixed at where the lower cutoff is
                }
            }
            
            //==================================================================
            // Correlation data for fitting:
            // 1st point below sExp_cut_hi to 1st point above sExp_cut_lo
            //==================================================================
            if (is_use_sExp_cut_hi) {
                if ((corr<=sExp_cut_hi) &&
                    (corr>=sExp_cut_lo) &&
                    (sExp_cut_index_lo==0)) {
                    correlation.push_back({time,corr});
                }
            }
            
            //==================================================================
            // Correlation data for fitting:
            // 1st point past sExp_tau_lo to 1st point above sExp_cut_lo
            //==================================================================
            if (is_use_plateautime) {
                if ((time>sExp_tau_lo) &&
                    (corr>=sExp_cut_lo) &&
                    (sExp_cut_index_lo==0)) {
                    correlation.push_back({time,corr});
                }
            }
        }
        readFile.close();
    }
    else {
        cout
        << "FitData::read_individual_correlation_data: correlation file cannot open."
        << "\n";
    }
    
    try {
        if (correlation.size()==0) throw 0;
        std::sort(correlation.begin(),correlation.end(),sortIncreasing);
        averaging_sorted_data(correlation,correlation_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_correlation_data:" << "\n"
        << "correlation size = " << i << "\n";
        correlation.push_back({0,0});
        std::sort(correlation.begin(),correlation.end(),sortIncreasing);
        averaging_sorted_data(correlation,correlation_sortincreasing);
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
    //cout << i << "\n";
    ifstream readtaueqdata(in.c_str());
    
    if (readtaueqdata.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readtaueqdata,lineContent)) {
            istringstream iss(lineContent); // (T,taueq)
            iss >> Tdata;
            iss >> tauFitdata;
            taueq.push_back({Tdata,tauFitdata});
        }
        readtaueqdata.close();
    }
    else {
        cout
        << "FitData::read_individual_taueq_data: taueq datafile cannot open." << "\n";
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
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
    //cout << i << "\n";
    ifstream readtaueqdata(in.c_str());
    
    if (readtaueqdata.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readtaueqdata,lineContent)) {
            istringstream iss(lineContent); // (T,taueq)
            
            if (model=="Arrhenius") {
                
                iss >> Tdata;
                
                if (Tdata>=cutTforArrhenius) {
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            }
            else {
                
                iss >> Tdata;
                iss >> tauFitdata;
                taueq.push_back({Tdata,tauFitdata});
            }
        }
        readtaueqdata.close();
    }
    else {
        cout
        << "FitData::read_individual_taueq_data: taueq datafile cannot open." << "\n";
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
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
    //cout << i << "\n";
    ifstream readtempsdata(i.c_str());
    
    if (readtempsdata.is_open()) {
        string lineContent;
        double Tdata=0;
        while (getline(readtempsdata,lineContent)) {
            istringstream iss(lineContent); // T
            iss >> Tdata;
            temps.push_back({Tdata});
        }
        readtempsdata.close();
    }
    else {
        cout
        << "FitData::read_individual_temps_data: temps datafile cannot open." << "\n";
    }
    
    try {
        if(temps.size()==0) throw 0;
        std::sort(temps.begin(),temps.end(),sortDecreasing);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_individual_temps_data:" << "\n"
        << "temps size = " << i << "\n";
        temps.push_back({0});
        std::sort(temps.begin(),temps.end(),sortDecreasing);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
}





void FitData::read_all_temps_data(const StructureClass& sysVar,
                                  const int n_sys)
{
    temps_sortdecreasing.clear();
    temps_sortincreasing.clear();
    vector<vector<double>> temps;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/tempsInfo_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readtempsdata(i.c_str());
        
        if (readtempsdata.is_open()) {
            string lineContent;
            double Tdata=0;
            while (getline(readtempsdata,lineContent)) {
                istringstream iss(lineContent); // T
                iss >> Tdata;
                temps.push_back({Tdata});
            }
            readtempsdata.close();
        }
        else {
            cout
            << "read_all_T_data: temps datafile cannot open." << "\n";
        }
    }
    
    try {
        if(temps.size()==0) throw 0;
        std::sort(temps.begin(),temps.end(),sortDecreasing);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_T_data:" << "\n"
        << "temps size = " << i << "\n";
        temps.push_back({0});
        std::sort(temps.begin(),temps.end(),sortDecreasing);
        averaging_sorted_data(temps,temps_sortdecreasing);
        std::sort(temps.begin(),temps.end(),sortIncreasing);
        averaging_sorted_data(temps,temps_sortincreasing);
    }
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readtaueqdata(i.c_str());
        
        if (readtaueqdata.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readtaueqdata,lineContent)) {
                istringstream iss(lineContent); // (T,taueq)
                iss >> Tdata;
                iss >> tauFitdata;
                taueq.push_back({Tdata,tauFitdata});
            }
            readtaueqdata.close();
        }
        else {
            cout
            << "read_all_taueq_data: taueq datafile cannot open." << "\n";
        }
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys,
                                  const double TA_d)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readtaueqdata(i.c_str());
        
        if (readtaueqdata.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readtaueqdata,lineContent)) {
                istringstream iss(lineContent); // (T,taueq)
                iss >> Tdata;
                iss >> tauFitdata;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    taueq.push_back({Tdata,tauFitdata});}
            }
            readtaueqdata.close();
        }
        else {
            cout
            << "read_all_taueq_data: taueq datafile cannot open." << "\n";
        }
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitData::read_all_taueq_data(const StructureClass& sysVar,
                                  const int n_sys,
                                  const std::string& cond)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/taueq_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readtaueqdata(i.c_str());
        
        if (readtaueqdata.is_open()) {
            string lineContent;
            double Tdata=0,tauFitdata=0;
            while (getline(readtaueqdata,lineContent)) {
                istringstream iss(lineContent); // (T,taueq)
                
                iss >> Tdata;
                
                if (cond=="Arrhenius") {
                    
                    if (Tdata>=cutTforArrhenius) {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                
                else if (cond=="continuous") {
                    
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
                
                else if (cond=="neighboring") {
                    
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
                
                else if (cond=="partial_lt") { // T<Tx
                    
                    if (get_is_fit_by_Tc()) {
                        if (Tdata<get_Tc_MCT()) {
                            
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    }
                    else if (get_is_fit_by_TA()) {
                        if (Tdata<get_TA_avg()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    }
                    else {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                
                else if (cond=="partial_gt") { // T>Tx
                    
                    if (get_is_fit_by_Tc()) {
                        if (Tdata>get_Tc_MCT()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    }
                    else if (get_is_fit_by_TA()) {
                        if (Tdata>get_TA_avg()) {
                            iss >> tauFitdata;
                            taueq.push_back({Tdata,tauFitdata});
                        }
                    }
                    else {
                        iss >> tauFitdata;
                        taueq.push_back({Tdata,tauFitdata});
                    }
                }
                
                else {
                    
                    iss >> Tdata;
                    iss >> tauFitdata;
                    //if (tauFitdata<log10(tauFit_cut)) {
                    //    taueq.push_back({Tdata,tauFitdata});}
                    taueq.push_back({Tdata,tauFitdata});
                }
            }
            readtaueqdata.close();
        }
        else {
            cout
            << "read_all_taueq_data: taueq datafile cannot open." << "\n";
        }
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitData::read_all_sExp_params(const StructureClass& sysVar,
                                   const int n_sys)
{
    param_sortdecreasing.clear();
    vector<vector<double>> param;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/sExp_params_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readparamdata(i.c_str());
        
        if (readparamdata.is_open()) {
            string lineContent;
            double tmpdata;
            vector<double> rowdata;
            while (getline(readparamdata,lineContent)) {
                istringstream iss(lineContent);
                rowdata.clear();
                
                while(iss>>tmpdata) rowdata.push_back(tmpdata); // iss>>tmpdata
                
                param.push_back(rowdata);
            }
            readparamdata.close();
        }
        else {
            cout
            << "read_all_sExp_params: sExp_Params file cannot open." << "\n";
        }
    }
    
    try {
        if(param.size()==0) throw 0;
        std::sort(param.begin(),param.end(),sortDecreasing);
        averaging_sorted_data(param,param_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitData::read_all_sExp_params:" << "\n"
        << "param size = " << i << "\n";
        param.push_back({0,0});
        std::sort(param.begin(),param.end(),sortDecreasing);
        averaging_sorted_data(param,param_sortdecreasing);
    }
}





void FitData::Tc_correction(const std::string& model)
{
    if ((model=="MCT")||(model=="SOU")) {
        
        size_t vecSize=xdata.size();
        string tleq_str=xdata[vecSize-2];
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
    if (outputFile.is_open()) {
        
        if ((model=="KWW") ||
            (model=="KWW_pwr") ||
            (model=="mKWW")    ||
            (model=="mKWW_pwr") ) {
            
            outputFile
            << amdat_version
            << "\n\n"
            
            << "Wavenumber "
            << wavenumber
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
            << "\n\n"
            
            << "cut_index_hi(1-based) " << sExp_cut_index_hi
            << "\n"
            << "cut_index_lo(1-based) " << sExp_cut_index_lo
            << "\n\n";
        }
        else if (model=="Arrhenius") {
            
            //
            
        }
        else {
            
            outputFile << "tauFit_cut " << "\n";
            
            for (size_t i=0; i<sysVar.get_n_regime(); ++i) {
                
                double previous_teq,current_teq,teq;
                if (i>0) {
                    previous_teq = sysVar.get_equilibration_times()[i-1];
                    current_teq  = sysVar.get_equilibration_times()[i];
                    teq = max(previous_teq,current_teq);
                }
                else {
                    teq = sysVar.get_equilibration_times()[i];
                }
                
                outputFile
                << "Regime_" << i << "\n"
                << teq/get_n_equ_blocks()<< "\n";
            }
            outputFile << "\n";
        }
    }
    else {
        cout
        << "FitData::write_taucutoff: file cannot open."
        << "\n";
    }
}





void FitData::write_errcurve_actual(ofstream& outputFile,
                                    const string& model,
                                    const vector<vector<double>>& vec2d_sortdecreasing)
{
    if (outputFile.is_open()) {
        
        outputFile
        << "\n"
        << "x  y  y_fit  ErrCurve  " << "\n"
        << "=========================================== " << "\n";
        for (int i=0; i<rep.errcurve.length(); ++i) {
            
            outputFile
            << vec2d_sortdecreasing[i][0] << " "
            << vec2d_sortdecreasing[i][1] << " "
            << calc_time_given_temp(c,vec2d_sortdecreasing[i][0],model) << " "
            << rep.errcurve[i] << "\n";
        }
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitData::write_errcurve_actual: file cannot open."
        << "\n";
    }
}





void FitData::write_errcurve_log(ofstream& outputFile,
                                 const string& model,
                                 const vector<vector<double>>& vec2d_sortdecreasing)
{
    if (outputFile.is_open()) {
        
        outputFile
        << "\n"
        << "x  y  y_fit  ErrCurve  " << "\n"
        << "=========================================== " << "\n";
        for (int i=0; i<rep.errcurve.length(); ++i) {
            
            outputFile
            << vec2d_sortdecreasing[i][0] << " "
            << vec2d_sortdecreasing[i][1] << " "
            << log10(calc_time_given_temp(c,vec2d_sortdecreasing[i][0],model)) << " "
            << rep.errcurve[i] << "\n";
        }
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitData::write_errcurve_log: file cannot open."
        << "\n";
    }
}





void FitData::compute_Tg_fragility(StructureClass& sysVar,
                                   const alglib::real_1d_array& c,
                                   const std::string& model)
{
    /* extrapolated Tg and fragility */
    //------------------------------------------------------------------
    double Tg_extrp = calc_temp_given_time(c,extrp_time,model)[0];
    double m_extrp  = calc_temp_given_time(c,extrp_time,model)[1];
    sysVar.get_Tg_extrp()[current_regime].push_back(Tg_extrp);
    sysVar.get_m_extrp()[current_regime].push_back(m_extrp);
    //------------------------------------------------------------------
    
    /* computational Tg and fragility */
    //------------------------------------------------------------------
    double Tg_compu = calc_temp_given_time(c,compu_time,model)[0];
    double m_compu  = calc_temp_given_time(c,compu_time,model)[1];
    sysVar.get_Tg_compu()[current_regime].push_back(Tg_compu);
    sysVar.get_m_compu()[current_regime].push_back(m_compu);
    //------------------------------------------------------------------
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
    double Tg_extrp = calc_temp_given_time(c,extrp_time,model)[0];
    double m_extrp  = calc_temp_given_time(c,extrp_time,model)[1];
    //------------------------------------------------------------------
    
    if (index==0) {
        set_largest_Tg(Tg_extrp);
        set_index_largest_Tg(index);
        set_largest_m(m_extrp);
        set_index_largest_m(index);
    }
    else {
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
    sysVar.set_calcvector(sysVar.get_m_extrp()[current_regime]);
    double m_accumAvg=sysVar.calc_mean();
    double m_accumStdev=sysVar.get_sample_stdev();
    relativeStdev=m_accumStdev/m_accumAvg;
}





void FitData::write_Tg_fragility(std::ofstream& outputFile,
                                 const alglib::real_1d_array& c,
                                 const std::string& model)
{
    if (outputFile.is_open()) {
        
        /* extrapolated Tg and fragility */
        //------------------------------------------------------------------
        double Tg_extrp = calc_temp_given_time(c,extrp_time,model)[0];
        double m_extrp  = calc_temp_given_time(c,extrp_time,model)[1];
        //------------------------------------------------------------------
        
        /* computational Tg and fragility */
        //------------------------------------------------------------------
        double Tg_compu = calc_temp_given_time(c,compu_time,model)[0];
        double m_compu  = calc_temp_given_time(c,compu_time,model)[1];
        //------------------------------------------------------------------
        
        outputFile
        << "-------------------------------------" << "\n"
        << "Tg_extrp [@("<<extrp_time<<")timeUnit] "
        << Tg_extrp
        << "\n"
        << "m_extrp [@Tg_extrp] "
        << m_extrp
        << "\n\n"
        
        << "Tg_compu [@("<<compu_time<<")timeUnit] "
        << Tg_compu
        << "\n"
        << "m_compu [@Tg_compu] "
        << m_compu
        << "\n"
        << "-------------------------------------" << "\n"
        << "\n";
        
    }
    else {
        cout
        << "FitData::write_Tg_fragility: file cannot open."
        << "\n";
    }
}





void FitData::write_fit_correlation(std::ofstream& outputFile,
                                    const std::string& tcalc)
{
    if (outputFile.is_open()) {
        
        if (is_use_gammafunc) {
            outputFile
            << "tauFit (integral by Gamma function) "
            << sExp_tauFit_gammafunction(c,tcalc)
            << "\n\n";
        }
        else {
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
            << correlation_sortincreasing[i][0] << " "
            << sExp_func_value(c,correlation_sortincreasing[i][0],tcalc)
            << "\n";
        }
        outputFile << "\n\n";
        
        outputFile
        << "time  correlation(time)          " << "\n"
        << "=================================" << "\n";
        for (size_t i=0; i<correlation_original.size(); ++i) {
            outputFile
            << correlation_original[i][0] << " "
            << correlation_original[i][1]
            << "\n";
        }
    }
    else {
        cout
        << "FitData::write_fit_correlation: file cannot open."
        << "\n";
    }
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
    //cout << o << "\n";
    ofstream write_params_avg(o.c_str());
    
    read_all_sExp_params(sysVar,n_sys);
    
    if (tcalc=="KWW") {
        write_params_avg << "T A tau beta r2" << "\n";
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
    
    for (size_t i=0; i<param_sortdecreasing.size(); ++i) {
        for (size_t ii=0; ii<param_sortdecreasing[0].size(); ++ii) {
            write_params_avg << param_sortdecreasing[i][ii] << " ";
        }
        write_params_avg << "\n";
    }
}





void FitData::write_tauFitfit_vs_T(std::ofstream& outputFile,
                                   const alglib::real_1d_array& c,
                                   const std::string& model)
{
    if (outputFile.is_open()) {
        
        double T_initial=Theq;
        
        int equal_pieces=1000;
        double Ti=0;
        
        /* format: T log10(tauFit_fit) */
        //------------------------------------------------------------------
        if (model=="Arrhenius") {
            
            double res=res_Arrhenius/(corrFacforT/precision); // T resolution
            
            outputFile << "\n"
            << "T  log10(tauFit_fit)             " << "\n"
            << "=================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                
                Ti=T_initial-i*res;
                
                outputFile
                << Ti << " "
                << log10(calc_time_given_temp(c,Ti,model))
                << "\n";
            }
        }
        
        else {
            
            double Tg_extrp=calc_temp_given_time(c,extrp_time,model)[0];
            
            outputFile << "\n"
            << "T  log10(tauFit_fit)             " << "\n"
            << "=================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                
                Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
                
                outputFile
                << Ti << " "
                << log10(calc_time_given_temp(c,Ti,model))
                << "\n";
            }
            
            if (get_is_fit_by_TA()) {
                
                if (current_regime==(n_regime-1)) {
                    
                    T_initial=TA_avg;
                    
                    outputFile << "\n"
                    << "fit T<TA" << "\n\n"
                    
                    << "T  log10(tauFit_fit)            " << "\n"
                    << "=================================" << "\n";
                    
                    for (int i=0; i<=equal_pieces; ++i) {
                        
                        Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
                        
                        outputFile
                        << Ti << " "
                        << log10(calc_time_given_temp(c,Ti,model))
                        << "\n";
                    }
                }
            }
        }
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitData::write_tauFitfit_vs_T: file cannot open."
        << "\n";
    }
    
}





void FitData::write_tauFitfit_vs_invT(std::ofstream& outputFile,
                                      const alglib::real_1d_array& c,
                                      const std::string& model)
{
    if (outputFile.is_open()) {
        
        double T_initial=Theq;
        
        int equal_pieces=1000;
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
                << (1000.0/(corrFacforT/precision))/Ti << " "
                << log10(calc_time_given_temp(c,Ti,model))
                << "\n";
            }
        }
        
        else {
            
            double Tg_extrp = calc_temp_given_time(c,extrp_time,model)[0];
            
            outputFile << "\n"
            << 1000.0/(corrFacforT/precision)
            << "/T  log10(tauFit_fit)            " << "\n"
            << "=================================" << "\n";
            
            for (int i=0; i<=equal_pieces; ++i) {
                
                Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
                
                outputFile
                << (1000.0/(corrFacforT/precision))/Ti << " "
                << log10(calc_time_given_temp(c,Ti,model))
                << "\n";
            }
            
            if (get_is_fit_by_TA()) {
                
                if (current_regime==(n_regime-1)) {
                    
                    T_initial=TA_avg;
                    
                    outputFile << "\n"
                    << "fit T<TA" << "\n\n"
                    
                    << 1000.0/(corrFacforT/precision)
                    << "/T  log10(tauFit_fit)            " << "\n"
                    << "=================================" << "\n";
                    
                    for (int i=0; i<=equal_pieces; ++i) {
                        
                        Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
                        
                        outputFile
                        << (1000.0/(corrFacforT/precision))/Ti << " "
                        << log10(calc_time_given_temp(c,Ti,model))
                        << "\n";
                    }
                }
            }
        }
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitData::write_tauFitfit_vs_invT: file cannot open."
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
    //cout << o << "\n";
    ofstream write_qchTs(o.c_str(),ofstream::app); // use "append"
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to the "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero and with fixed format
    //--------------------------------------------------------------------------
    write_qchTs << fixed << setprecision(0);
    
    if (write_qchTs.is_open()) {
        
        write_qchTs << "{";
        for (size_t i=0; i<sysVar.get_quenchingTs()[index].size(); ++i) {
            
            if (i!=(sysVar.get_quenchingTs()[index].size()-1)) {
                write_qchTs
                << sysVar.get_quenchingTs()[index][i] << ",";
            }
            else {
                write_qchTs
                << sysVar.get_quenchingTs()[index][i] << "}";
            }
        }
        write_qchTs << "\n";
        
        write_qchTs.close();
    }
    else {
        cout
        << "FitData::write_qchTs: 'tempsInfo_qchTs' cannot open."
        << "\n";
    }
    
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
    //cout << o << "\n";
    ofstream write_equTs(o.c_str(),ofstream::app); // use "append"
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to the "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero and with fixed format
    //--------------------------------------------------------------------------
    write_equTs << fixed << setprecision(0);
    
    //==========================================================================
    // WRITE: The path to the file storing {T,tau_eq} info
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_taueq(o.c_str(),ofstream::app); // 'append'
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has the actual values;
    // (ex. T=1234.567 for real unit; T=1.234567 for lj unit)
    // so "setprecision" is set to "n_digits" according to the unit adopted
    //--------------------------------------------------------------------------
    if (0) {
        if ((sysVar.get_systemUnit()=="real")||(sysVar.get_systemUnit()=="metal")) {
            write_taueq
            << fixed << setprecision(sysVar.get_n_digits());
        }
        else if (sysVar.get_systemUnit()=="lj") {
            write_taueq
            << fixed << setprecision(sysVar.get_n_digits()+3);
        }
    }
    
    
    //==========================================================================
    // WRITE: The path to the file storing {1/T,tau_eq} info
    //==========================================================================
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/taueq_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_taueqinv(o.c_str(),ofstream::app); // 'append'
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has the actual values;
    // (ex. T=1234.567 for real unit; T=1.234567 for lj unit)
    // so "setprecision" is set to "n_digits" according to the unit adopted
    //--------------------------------------------------------------------------
    if (0) {
        if ((sysVar.get_systemUnit()=="real")||(sysVar.get_systemUnit()=="metal")) {
            write_taueqinv
            << fixed << setprecision(sysVar.get_n_digits());
        }
        else if (sysVar.get_systemUnit()=="lj") {
            write_taueqinv
            << fixed << setprecision(sysVar.get_n_digits()+3);
        }
    }
    
    if (write_equTs.is_open()&&
        write_taueq.is_open()&&
        write_taueqinv.is_open()) {
        
        
        // if below block is enabled:
        // all in-equilibrium T's shuold be higher than the first
        // out-of-equilibrium T (Toe)
        //----------------------------------------------------------------------
        if (is_imposeStrictTeq) {
            
            vector<double> inequTs;
            vector<vector<double>> tauEqus;
            
            /** Find first ooe T **/
            double Toe_max=0;
            for (int i=0; i<sysVar.get_tauFit()[index].size(); ++i) {
                if (sysVar.get_tauFit()[index][i][1]>=log10(tauFit_cut)) {
                    Toe_max=sysVar.get_tauFit()[index][i][0];
                    break;
                }
            }
            /** Every ie T should be higher than the first ooe T **/
            for (int i=0; i<sysVar.get_tauEqu()[index].size(); ++i) {
                double Ti  =sysVar.get_tauEqu()[index][i][0];
                double taui=sysVar.get_tauEqu()[index][i][1];
                if (Ti>Toe_max) {
                    inequTs.push_back(Ti);
                    tauEqus.push_back({Ti,taui});
                }
            }
            /* replace ie T of current regime by new ie T's based on this feature */
            sysVar.get_equilibratedTs()[index].clear();
            sysVar.get_equilibratedTs()[index]=inequTs;
            sysVar.get_tauEqu()[index].clear();
            sysVar.get_tauEqu()[index]=tauEqus;
        }
        //----------------------------------------------------------------------
        
        
        /* store in-equilibrium T's into the contianer */
        //----------------------------------------------------------------------
        write_equTs << "{";
        for (int i=0; i<sysVar.get_equilibratedTs()[index].size(); ++i) {
            
            if (i!=(sysVar.get_equilibratedTs()[index].size()-1)) {
                write_equTs
                << sysVar.get_equilibratedTs()[index][i] << ",";
            }
            else {
                write_equTs
                << sysVar.get_equilibratedTs()[index][i] << "}";
            }
        }
        write_equTs << "\n";
        //----------------------------------------------------------------------
        
        
        /* write tau_eq data to file */
        //----------------------------------------------------------------------
        for (int i=0; i<sysVar.get_tauEqu()[index].size(); ++i) {
            
            double T_actual  =sysVar.get_tauEqu()[index][i][0]*pow(corrFacforT,-1); //NOTE!
            double log_tauFit=sysVar.get_tauEqu()[index][i][1];
            
            write_taueq
            << setprecision(10)
            << T_actual << " " << log_tauFit << "\n";
            
            write_taueqinv
            << setprecision(10)
            << (1000.0/(corrFacforT/precision))/(T_actual) << " "
            << log_tauFit
            << "\n";
        }
        //----------------------------------------------------------------------
        
        write_equTs.close();
        write_taueq.close();
        write_taueqinv.close();
    }
    else {
        cout
        << "FitData::write_equTs: 'tempsInfo_equTs' cannot open."
        << "\n";
    }
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
    
    if (model=="Arrhenius") {
        
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_Arrhenius_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        //cout << o << "\n";
        ofstream write_fit_Arrhenius_avg(o.c_str());
        
        if (write_fit_Arrhenius_avg.is_open()) {
            
            write_fit_Arrhenius_avg
            << "Data fit to <" << model << "> functional form"
            << "\n\n";
            
            write_fit_Arrhenius_avg << "avg fit coeffs" << "\n";
            for (size_t i=0; i<c_avg.length(); ++i) {
                write_fit_Arrhenius_avg << c_avg[i] << " ";
            }
            write_fit_Arrhenius_avg << "\n\n";
            
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
            }
            write_fit_Arrhenius_avg << "\n\n";
            
            write_tauFitfit_vs_T(write_fit_Arrhenius_avg,c_avg,"Arrhenius");
            write_tauFitfit_vs_invT(write_fit_Arrhenius_avg,c_avg,"Arrhenius");
            
            write_fit_Arrhenius_avg.close();
        }
        else {
            cout
            << "fit_Arrhenius_avg: 'fit_Arrhenius_avg.dat' cannot open."
            << "\n";
        }
    }
    
    else {
        
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fit_taueq_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append("_Regime"+to_string((long long int)current_regime));
        o.append(".dat");
        //cout << o << "\n";
        ofstream write_fit_tauFit_avg(o.c_str());
        
        if (write_fit_tauFit_avg.is_open()) {
            
            write_fit_tauFit_avg
            << "Data fit to <" << model << "> functional form"
            << "\n\n"
            
            << "(tauFit defined @KWW=" << get_sExp_tauFit() << ")"
            << "\n\n";
            
            write_Tg_fragility(write_fit_tauFit_avg,c_avg,model);
            write_fitcutoff(write_fit_tauFit_avg,sysVar,model);
            
            write_fit_tauFit_avg << "avg fit coeffs" << "\n";
            for (size_t i=0; i<c_avg.length(); ++i) {
                write_fit_tauFit_avg << c_avg[i] << " ";
            }
            write_fit_tauFit_avg << "\n\n";
            
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
            }
            write_fit_tauFit_avg << "\n\n";
            
            write_tauFitfit_vs_T(write_fit_tauFit_avg,c_avg,model);
            write_tauFitfit_vs_invT(write_fit_tauFit_avg,c_avg,model);
            
            write_fit_tauFit_avg.close();
        }
        else {
            cout
            << "fit_tauFit_avg: 'fit_taueq_avg.dat' cannot open."
            << "\n";
            //system("read -t5");exit(EXIT_FAILURE);
        }
        
    }
}





void FitData::useDynamicsRange(StructureClass& sysVar,
                               const int n_sys)
{
    /* change dynamics range based on tau0 */
    const int current_regime = sysVar.get_current_regime();
    
    if (current_regime==0) {  // NOTE!
        
        vector<double> teq_org=sysVar.get_equilibration_times_org();
        
        size_t size_org=teq_org.size();
        double equh_org=teq_org[size_org-1];
        double tauh_org=equh_org/get_n_equ_blocks();
        double range_tau0=sysVar.get_dynamicRange_tau0();
        
        double tauh_new=0; // new highest tau_alpha cut//
        double ratio=0;    // rescaling ratio
        
        tauh_new = range_tau0*tau0_Arrhenius; // use tau0 from Arrhenius fit
        ratio = tauh_new/tauh_org;
        
        /* rescale the original tau_alpha cutoffs */
        
        // NOTE:
        // rule of thumb for rescaling is that the first 3 regimes
        // (which usually only cover 1 decade dynamic range)
        // are rescaled relative to the high temperatures defined
        // by user; and rest of the regimes are rescaled by applying
        // the rescaling ratio calculated above
        
        read_all_taueq_data(sysVar,n_sys); // averaged taueq data
        
        double tauhiavg = taueq_sortdecreasing[0][1];
        
        teq_org[0]=get_n_equ_blocks()*tauhiavg;
        
        // NOTE::
        // if new dynamic range is smaller than 1 decade of hiT tau_alpha,
        // make new range based on highT alone
        // in case of insufficient space for extrapolation
        
        if (tauh_new<(tauhiavg*pow(10,1))) {
            for (size_t i=1; i<teq_org.size(); ++i) {
                teq_org[i] = teq_org[0]*pow(10,0.5*i);
            }
            sysVar.set_message_dynRange
            ("dynamic range: type1");
        }
        else {
            
            if (teq_org.size()>3) {
                
                for (size_t i=1; i<3; ++i) {
                    teq_org[i]=teq_org[0]*pow(10.0,0.5*i);
                }
                for (size_t i=3; i<teq_org.size(); ++i) {
                    teq_org[i] *= ratio;
                    if (teq_org[i]<(teq_org[0]*pow(10.0,1.0))) {
                        teq_org[i]=teq_org[0]*pow(10.0,0.5*i);
                        sysVar.set_message_dynRange
                        ("dynamic range: type2");
                    }
                }
            }
            else {
                
                for (size_t i=1; i<teq_org.size(); ++i) {
                    teq_org[i] *= ratio;
                    if (teq_org[i]<teq_org[0]) {
                        teq_org[i]=teq_org[0]*pow(10.0,0.5*i);
                        sysVar.set_message_dynRange
                        ("dynamic range: type3");
                    }
                }
            }
        }
        sysVar.set_equilibration_times(teq_org);
        
        //----------------------------------------------------------
        // Record original and rescaled equilibration times
        //----------------------------------------------------------
        string o;
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/new_teq_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        //cout << o << "\n";
        ofstream write_new_teq(o.c_str());
        
        if (write_new_teq.is_open()) {
            
            vector<double> teq=sysVar.get_equilibration_times_org();
            
            write_new_teq
            << "original teq" << "\n";
            for (size_t i=0; i<teq_org.size(); ++i) {
                write_new_teq << teq[i] << "\n";
            }
            
            write_new_teq
            << "\n"
            << "original tau_alpha" << "\n";
            for (size_t i=0; i<teq_org.size(); ++i) {
                write_new_teq << teq[i]/get_n_equ_blocks() << "\n";
            }
            write_new_teq << "\n\n";
            
            
            write_new_teq
            << "teq after rescaled to cover "
            << log10(sysVar.get_dynamicRange_tau0()) << " "
            << "decades of tau0" << "\n";
            for (size_t i=0; i<teq_org.size(); ++i) {
                write_new_teq
                << teq_org[i] << "\n";
            }
            
            write_new_teq
            << "\n"
            << "tau_alpha after rescaled" << "\n";
            for (size_t i=0; i<teq_org.size(); ++i) {
                write_new_teq
                << teq_org[i]/get_n_equ_blocks() << "\n";
            }
            
            write_new_teq
            << "\n"
            << "log(tau_alpha) after rescaled" << "\n";
            for (size_t i=0; i<teq_org.size(); ++i) {
                write_new_teq
                << log10(teq_org[i]/get_n_equ_blocks()) << "\n";
            }
            
            
            write_new_teq.close();
        }
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
    for (size_t i=0; i<coeffs[index].size(); ++i) {
        avgParams.push_back(0);
    }
    
    int total_count=0; // total counts of all trials & systems
    for (int i=0; i<=n_trl; ++i) { // up to current trial
        for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) { // up to current sys
            
            ++total_count;
            
            int index_i=i*n_system+ii;
            
            for (size_t iii=0; iii<coeffs[index_i].size(); ++iii) {
                avgParams[iii] += coeffs[index_i][iii];
            }
        }
    }
    for (size_t i=0; i<coeffs[index].size(); ++i) {
        avgParams[i] /= (double)(total_count);
    }
    
    c_avgVec.clear(); // NOTE: vector renewed every time!
    
    for (size_t i=0; i<avgParams.size(); ++i) {
        c_avgVec.push_back(avgParams[i]);
    }
}





void FitData::avgfitParams_extrp(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const std::string& extrp)
{
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    vector<double> avgParams;
    for (size_t i=0; i<ExtrpCoeffs[index].size(); ++i) {
        avgParams.push_back(0);
    }
    
    //--------------------------------------------------------------------------
    // from second regime on, if stdev/avg < 0.1,
    // use avg fitParams to shoot new temperatures
    //
    // NOTE:
    // thus desig because curve-fitting is usually bad in
    // predicting lowT behavior based on highT data;
    // so in first regime, use highest fragility case to ensure safe
    // extrapolation; and use avg fit coeffs in lowT regimes, where prediction
    // gets progressively improved as more lower T's are simulated;
    //--------------------------------------------------------------------------
    if ((current_regime>0)&&(relativeStdev<0.1)) {
        
        int total_count=0; // total counts of all trials & systems
        for (int i=0; i<=n_trl; ++i) { // up to current trial
            for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) { // up to current sys
                
                ++total_count;
                
                int index_i=i*n_system+ii;
                
                for (size_t iii=0; iii<ExtrpCoeffs[index_i].size(); ++iii) {
                    avgParams[iii] += ExtrpCoeffs[index_i][iii];
                }
            }
        }
        for (size_t i=0; i<ExtrpCoeffs[index].size(); ++i) {
            avgParams[i] /= (double)(total_count);
        }
    }
    //--------------------------------------------------------------------------
    // in first regime || relstdev >= 0.1 in higher regimes,
    // use the fitParams of the highest fragility system to shoot
    //--------------------------------------------------------------------------
    else {
        
        int ilm=get_index_largest_m();
        
        for (size_t i=0; i<ExtrpCoeffs[ilm].size(); ++i) {
            avgParams[i] = ExtrpCoeffs[ilm][i];
        }
    }
    
    c_avgVec.clear(); // NOTE: vector renewed every time!
    
    for (size_t i=0; i<avgParams.size(); ++i) {
        c_avgVec.push_back(avgParams[i]);
    }
}





void FitData::find_cutTforArrhenius(const StructureClass& sysVar,
                                    const int n_sys)
{
    string model="Arrhenius";
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys,"continuous");
    taueq_sortincreasing=taueq_sortdecreasing;
    std::sort(taueq_sortincreasing.begin(),taueq_sortincreasing.end(),sortIncreasing);
    fitdata_processing_lsfit(taueq_sortincreasing);
    int n_data=(int)taueq_sortincreasing.size();
    int n_iteration=n_data-n_fit_chopping;
    for (int i=0; i<n_iteration; i++) {
        alglib_lsfit(model);
        string tem = xdata[xdata.size()-2];
        //string tau = ydata[ydata.size()-2];
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
    vector<vector<double>> taueq_org=tauVec;
    vector<string> xdata_org=xdata;
    vector<string> ydata_org=ydata;
    
    double tauA_calc=0;
    double deviation=0;
    int which_T=0;
    
    /* use only highT's to fit to Arrhenius */
    set_fitParams("Arrhenius");
    read_all_taueq_data(sysVar,n_sys,"Arrhenius");
    fitdata_processing_lsfit(taueq_sortdecreasing); // NOTE: using taueq
    alglib_lsfit("Arrhenius");
    
    /* use all taueq data to find TA */
    taueq_sortdecreasing=tauVec;
    
    double Ti=0;
    
    if (false) { // use log-scale tau
        
        // use real data to find first deviating T from Arrhenius
        //--------------------------------------------------------------
        for (int i=0; i<taueq_sortdecreasing.size(); ++i) {
            
            which_T=i;
            
            Ti = taueq_sortdecreasing[i][0];
            
            tauA_calc = taueq_sortdecreasing[i][1];
            
            tauA_Arrhenius = log10(calc_time_given_temp(c,Ti,"Arrhenius"));
            
            deviation=(tauA_calc-tauA_Arrhenius)/(tauA_Arrhenius);
            
            if (deviation>deviateTA) break;
        }
        //--------------------------------------------------------------
        
        
        //--------------------------------------------------------------
        // do linear interpolation between two T's found below and
        // above critical deviation in real data and find the T's
        // that is more exact of devating magnitude
        //--------------------------------------------------------------
        if (which_T>0) {
            
            double T_hi=0,T_lo=0;
            double tau_hi=0,tau_lo=0;
            double slope=0;
            double dT=0;
            int equal_parts=1000; // discretize by 1000 equal parts
            
            T_hi   = taueq_sortdecreasing[which_T-1][0];
            tau_hi = taueq_sortdecreasing[which_T-1][1];
            T_lo   = taueq_sortdecreasing[which_T][0];
            tau_lo = taueq_sortdecreasing[which_T][1];
            
            slope = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT=fabs(T_hi-T_lo)/(double)equal_parts;
            
            // Interpolation
            for (int i=0; i<=equal_parts; ++i) {
                
                Ti = T_hi - i*dT;
                
                tauA_calc = tau_hi + slope*(Ti-T_hi);
                
                tauA_Arrhenius = log10(calc_time_given_temp(c,Ti,"Arrhenius"));
                
                deviation=(tauA_calc-tauA_Arrhenius)/(tauA_Arrhenius);
                
                if (deviation>deviateTA) break;
            }
        }
        //--------------------------------------------------------------
        
        TA_avg   = Ti;                // NOTE: actual temperature
        tauA_avg = pow(10,tauA_calc); // NOTE: actual tau
    }
    
    
    else { // use actual tau
        
        // use real data to find first deviating T from Arrhenius
        //--------------------------------------------------------------
        for (int i=0; i<taueq_sortdecreasing.size(); ++i) {
            
            which_T=i;
            
            Ti = taueq_sortdecreasing[i][0];
            
            tauA_calc = pow(10,taueq_sortdecreasing[i][1]); // NOTE: actual tau (not log.tau)
            
            tauA_Arrhenius = calc_time_given_temp(c,Ti,"Arrhenius");
            
            deviation=(tauA_calc-tauA_Arrhenius)/(tauA_Arrhenius);
            
            if (deviation>deviateTA) break;
        }
        //--------------------------------------------------------------
        
        
        //--------------------------------------------------------------
        // do linear interpolation between two T's found below and
        // above critical deviation in real data and find the T's
        // that is more exact of devating magnitude
        //--------------------------------------------------------------
        if (which_T>0) {
            
            double T_hi=0,T_lo=0;
            double tau_hi=0,tau_lo=0;
            double slope=0;
            double dT=0;
            int equal_parts=1000; // discretize by 1000 equal parts
            
            T_hi   = taueq_sortdecreasing[which_T-1][0];
            tau_hi = pow(10,taueq_sortdecreasing[which_T-1][1]); // NOTE: actual tau
            T_lo   = taueq_sortdecreasing[which_T][0];
            tau_lo = pow(10,taueq_sortdecreasing[which_T][1]); // NOTE: actual tau
            
            slope = (tau_lo-tau_hi)/(T_lo-T_hi);
            dT=fabs(T_hi-T_lo)/(double)equal_parts;
            
            // Interpolation
            for (int i=0; i<=equal_parts; ++i) {
                
                Ti = T_hi - i*dT;
                
                tauA_calc = tau_hi + slope*(Ti-T_hi);
                
                tauA_Arrhenius = calc_time_given_temp(c,Ti,"Arrhenius");
                
                deviation=(tauA_calc-tauA_Arrhenius)/(tauA_Arrhenius);
                
                if (deviation>deviateTA) break;
            }
        }
        //--------------------------------------------------------------
        
        TA_avg   = Ti;        // NOTE: actual temperature
        tauA_avg = tauA_calc; // NOTE: actual tau
    }
    
    taueq_sortdecreasing=taueq_org;
    xdata=xdata_org;
    ydata=ydata_org;
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
    //cout << o << "\n";
    ofstream write_TA(o.c_str());
    
    if (write_TA.is_open()) {
        
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);
        double deviation=(tauA_avg-tauA_Arrhenius)/(tauA_Arrhenius);
        
        write_TA
        << "TA tauA tauA(Arrhenius) %deviation"
        << "\n"
        << "---------------------------------------------------"
        << "\n";
        
        write_TA
        << TA_avg << " "
        << tauA_avg << " "
        << tauA_Arrhenius << " "
        << deviation*100 << "\n";
        
        write_TA.close();
    }
    else {
        cout
        << "FitData::findTA: 'tempsInfo_TA' cannot open."
        << "\n";
    }
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
    
    if (write_tau0.is_open()) {
        
        find_tau0(sysVar,n_sys,model);
        double deviation=fabs(tau0_model-tau0_Arrhenius)/fabs(tau0_Arrhenius);
        
        write_tau0
        << "tau0("<< model <<")  tau0(Arrhenius)  %deviation    "
        << "\n"
        << "---------------------------------------------------"
        << "\n";
        
        write_tau0
        << tau0_model << " "
        << tau0_Arrhenius << " "
        << deviation*100
        << "\n";
        
        write_tau0.close();
    }
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
    //cout << o << "\n";
    ofstream write_findTc(o.c_str());
    
    if (write_findTc.is_open()) {
        
        /* fit MCT */
        
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
        
        
        if (0) {
            
            /* fit SOU */
            
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
        
    }
    else {
        cout << "findTc: 'fit_Tc.dat' cannot open." << "\n";
    }
    
    write_findTc.close();
    
    return Tc_MCT;
}





double FitData::get_Theq(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys)
{
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    return taueq_sortdecreasing[0][0];
}


double FitData::get_Theq(StructureClass& sysVar,
                         const int n_sys)
{
    read_all_taueq_data(sysVar,n_sys);
    return taueq_sortdecreasing[0][0];
}




double FitData::get_Tleq(StructureClass& sysVar,
                         const int n_trl,
                         const int n_sys)
{
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    size_t length=taueq_sortdecreasing.size();
    return taueq_sortdecreasing[length-1][0];
}


double FitData::get_Tleq(StructureClass& sysVar,
                         const int n_sys)
{
    read_all_taueq_data(sysVar,n_sys);
    size_t length=taueq_sortdecreasing.size();
    return taueq_sortdecreasing[length-1][0];
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
        fit_xData += xdata[i];
        fit_yData += ydata[i];
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
        
        double x = taueq_sortdecreasing[i][0];
        double y = taueq_sortdecreasing[i][1];
        
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
        fit_xData += xdata[i];
        fit_yData += ydata[i];
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
    
    
    if (write_Tg_fit_every.is_open()&&
        write_m_fit_every.is_open() &&
        write_fit_every_param.is_open()) {
        
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
        }
        else {
            
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
                    
                    if (write_fit_taueq_avg.is_open()) {
                        
                        write_fitModel(write_fit_taueq_avg,model);
                        write_fitarrays(write_fit_taueq_avg);
                        write_stopcond(write_fit_taueq_avg);
                        write_Tg_fragility(write_fit_taueq_avg,c,model);
                        write_fitinfo(write_fit_taueq_avg);
                        write_errcurve(write_fit_taueq_avg,model);
                        write_tauFitfit_vs_T(write_fit_taueq_avg,c,model);
                        write_tauFitfit_vs_invT(write_fit_taueq_avg,c,model);
                        
                        write_fit_taueq_avg.close();
                    }
                }
                
                double Tg_extrp=calc_temp_given_time(c,extrp_time,model)[0];
                double m_extrp =calc_temp_given_time(c,extrp_time,model)[1];
                
                string tem = xdata[xdata.size()-2];
                string tau = ydata[ydata.size()-2];
                
                double tem_double=stod(tem);
                double tau_double=stod(tau);
                
                /*------ write parameters to file ------*/
                write_fit_every_param
                << tem_double << " "
                << pow(10,tau_double) << " ";
                for (int i=0; i<get_coeffs_vD().size(); ++i) {
                    write_fit_every_param << c[i] << " ";
                }
                write_fit_every_param << rep.r2 << "\n";
                
                /*------ write Tg_extrp, m_extrp to file ------*/
                write_Tg_fit_every
                << tem_double << " "
                << pow(10,tau_double) << " "
                << Tg_extrp << " "
                << error_Tg(model)
                << "\n";
                
                write_m_fit_every
                << tem_double << " "
                << pow(10,tau_double) << " "
                << m_extrp << " "
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
    else {
        cout << "fit_full: 'fit_Tg_vs_tau.dat' cannot open." << "\n";
    }
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
    
    
    if (write_Tg_fit_every.is_open()&&
        write_m_fit_every.is_open() &&
        write_fit_every_param.is_open()) {
        
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
        }
        else {
            
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
                
                fit_every_processing(i,n_fit_sliding);
                
                alglib_lsfit(model);
                
                double Tg_extrp=calc_temp_given_time(c,extrp_time,model)[0];
                double m_extrp =calc_temp_given_time(c,extrp_time,model)[1];
                
                string tem = xdata[xdata.size()-2];
                string tau = ydata[ydata.size()-2];
                
                double tem_double=stod(tem);
                double tau_double=stod(tau);
                
                /*------ write parameters to file ------*/
                write_fit_every_param
                << tem_double << " "
                << pow(10,tau_double) << " ";
                for (int i=0; i<get_coeffs_vD().size(); ++i) {
                    write_fit_every_param << c[i] << " ";
                }
                write_fit_every_param << rep.r2 << "\n";
                
                /*------ write Tg_extrp, m_extrp to file ------*/
                write_Tg_fit_every
                << tem_double << " "
                << pow(10,tau_double) << " "
                << Tg_extrp << " "
                << error_Tg(model)
                << "\n";
                
                write_m_fit_every
                << tem_double << " "
                << pow(10,tau_double) << " "
                << m_extrp << " "
                << error_fragility(model,Tg_extrp,error_Tg(model))
                << "\n";
            }
        }
        write_Tg_fit_every.close();
        write_m_fit_every.close();
        write_fit_every_param.close();
    }
    else {
        cout << "fit_full: 'fit_Tg_vs_tau.dat' cannot open." << "\n";
    }
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
        }
        else if (partial==0) {
            
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
        
        
        if (write_Tg_fit_every.is_open() &&
            write_m_fit_every.is_open()  &&
            write_fit_every_param.is_open()) {
            
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
            }
            else {
                
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
                
                
                for (int i=0; i<n_iterations; ++i) {
                    
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
                                log10(calc_time_given_temp(c,get_Tc_MCT(),model))
                                <<"\n";
                            }
                            else if (get_is_fit_by_TA()) {
                                write_fit_taueq_avg
                                << "Tx=TA=" << get_TA_avg() << "\n"
                                << "log(tau|Tx)=" <<
                                log10(calc_time_given_temp(c,get_TA_avg(),model))
                                <<"\n";
                            }
                            
                            write_fitModel(write_fit_taueq_avg,model);
                            write_fitarrays(write_fit_taueq_avg);
                            write_stopcond(write_fit_taueq_avg);
                            write_Tg_fragility(write_fit_taueq_avg,c,model);
                            write_fitinfo(write_fit_taueq_avg);
                            write_errcurve(write_fit_taueq_avg,model);
                            write_tauFitfit_vs_T(write_fit_taueq_avg,c,model);
                            write_tauFitfit_vs_invT(write_fit_taueq_avg,c,model);
                            
                            write_fit_taueq_avg.close();
                        }
                    }
                    
                    double Tg_extrp=calc_temp_given_time(c,extrp_time,model)[0];
                    double m_extrp =calc_temp_given_time(c,extrp_time,model)[1];
                    
                    string tem = xdata[xdata.size()-2];
                    string tau = ydata[ydata.size()-2];
                    
                    double tem_double=stod(tem);
                    double tau_double=stod(tau);
                    
                    /*------ write parameters to file ------*/
                    write_fit_every_param
                    << tem_double << " "
                    << pow(10,tau_double) << " ";
                    for (int i=0; i<get_coeffs_vD().size(); ++i) {
                        write_fit_every_param << c[i] << " ";
                    }
                    write_fit_every_param << rep.r2 << "\n";
                    
                    /*------ write Tg_extrp, m_extrp to file ------*/
                    write_Tg_fit_every
                    << tem_double << " "
                    << pow(10,tau_double) << " "
                    << Tg_extrp << " "
                    << error_Tg(model)
                    << "\n";
                    
                    write_m_fit_every
                    << tem_double << " "
                    << pow(10,tau_double) << " "
                    << m_extrp << " "
                    << error_fragility(model,Tg_extrp,error_Tg(model))
                    << "\n";
                    
                    /*------ Tc ------*/
                    if (partial==0) {
                        if ((model=="MCT")||(model=="SOU")) {
                            
                            double Tc=c[1];
                            find_TA(sysVar,n_sys,taueq_sortdecreasing);
                            
                            write_fit_every_Tc
                            << tem_double << " "
                            << pow(10,tau_double) << " "
                            << Tc << " "
                            << TA_avg << " " << tauA_avg << " " << rep.r2
                            << "\n";
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
        else {
            cout << "fit_partial: 'fit_Tg_vs_tau.dat' cannot open." << "\n";
        }
    }
}





void FitData::shootForNewTemperatures(StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const real_1d_array& c_extrp,
                                      const string& extrp)
{
    // record fit info by the extrapolation model to file
    //--------------------------------------------------------------------------
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/extrp_taueq_");
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
    write_tauFitfit_vs_T(write_extrp,c_extrp,extrp);
    write_tauFit_vs_T(write_extrp);
    write_tauFitfit_vs_invT(write_extrp,c_extrp,extrp);
    write_tauFit_vs_invT(write_extrp);
    //--------------------------------------------------------------------------
    
    index=n_trl*n_system+(n_sys-n_sys_beg);
    
    string folder,file,output,input;
    folder.append(return_AnalysisFolderPath(sysVar));
    folder.append("/fit_data");
    file.append("/tempsInfo_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append(".dat");
    output=folder+file;
    ofstream write_tempsinfo(output.c_str(),ofstream::app); // use "append"
    //--------------------------------------------------------------------------
    // NOTE:
    // temperatures in this file has been normalized to "expanded" format
    // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
    // so "setprecision" needs to be set to zero
    //--------------------------------------------------------------------------
    write_tempsinfo << fixed << setprecision(0);
    
    if (write_tempsinfo.is_open()) {
        
        
        if (current_regime==0) {
            if (sysVar.get_is_updatetinfo()) {
                if (sysVar.times_write_tempsinfo<(sysVar.max_times_retry+1)) {
                    for (size_t i=0; i<sysVar.get_temperaturesInfo()[index].size(); ++i) {
                        write_tempsinfo << sysVar.get_temperaturesInfo()[index][i] << "\n";
                    }
                }
            }
        }
        
        /*======================================================================
         Extrapolation of New Temperatures
         
         Invoke extrapolation (n_regime-1) times
         (not in the final regime)
         
         NOTE:
         1.
         number of extrapolated temperatures is given by
         ${n_next_temps} --> # of T's in the next regime;
         don't mix up with # of T's in the current regime.
         
         2.
         Any extrapolated temperature should not be repeated
         with any of the existing/simulated temperatures.
         =====================================================================*/
        if (current_regime<(n_regime-1)) {
            
            /* Collect all the simulated temperatures from file */
            vector<double> simtemps;
            string input=folder+file;
            ifstream readfile(input.c_str());
            if (readfile.is_open()) {
                string lineContent;
                double Tdata=0;
                while (getline(readfile,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    simtemps.push_back(Tdata);
                }
            }
            else {
                cout
                << "shootForNewTemperatures: readtempinfo cannot open."
                << "\n";
            }
            
            double teql,taueql;   // "lower  cutoff time" (i.e. higher T)
            double teqh,taueqh;   // "higher cutoff time" (i.e. lower  T)
            double Tl,Th;         // T calculated at "taueql" & "taueqh"
            double Tl_calc,Th_calc;
            
            teql   = sysVar.get_equilibration_times()[current_regime];
            teqh   = sysVar.get_equilibration_times()[current_regime+1];
            taueql = teql/get_n_equ_blocks();
            taueqh = teqh/get_n_equ_blocks();
            
            
            // NOTE: Determine Temperature Range in Next Simulation Regime
            //------------------------------------------------------------------
            // Th  : T @ hihger tau_alpha cut (lower T)
            // Tl  : lowest in-equilibrium T simulated
            //------------------------------------------------------------------
            if (get_is_avgFitParams()) {
                
                real_1d_array c_avg=c_extrp;
                
                for (size_t i=0; i<c_avgVec.size(); ++i) {
                    c_avg[i]=c_avgVec[i];
                }
                
                /* T's calculated by avg fitParams */
                Tl_calc = calc_temp_given_time(c_avg,taueql,extrp)[0];
                Th_calc = calc_temp_given_time(c_avg,taueqh,extrp)[0];
                
                Th = Th_calc;
                Tl = Tleq;
            }
            else {
                
                /* T's calculated by each individul trial */
                Tl_calc = calc_temp_given_time(c_extrp,taueql,extrp)[0];
                Th_calc = calc_temp_given_time(c_extrp,taueqh,extrp)[0];
                
                Th = Th_calc;
                Tl = Tleq;
            }
            
            /* process value to have correct precision */
            //------------------------------------------------------------------
            Tl = std::trunc(Tl*corrFacforT); // NOTE: convert!
            Th = std::trunc(Th*corrFacforT);
            //------------------------------------------------------------------
            
            /* "n_next_temps" should >= 2 */
            int n_next_temps=sysVar.get_n_regime_temps()[current_regime+1];
            if (n_next_temps<2) {n_next_temps=2;}
            
            /* get new temperatures of next regime (stored in vector<double>) */
            //------------------------------------------------------------------
            vector<double> hightolow=
            get_interpolated_Ts(sysVar,n_trl,n_sys,n_next_temps,Tl,Th);
            //------------------------------------------------------------------
            
            /* replace current temperatures with temperatures of next regime */
            //------------------------------------------------------------------
            if (sysVar.get_is_updatetinfo()) {
                sysVar.get_temperaturesInfo()[index]=hightolow;}
            //------------------------------------------------------------------
            
            /* clear content of container for new values */
            //******************************************************************
            sysVar.get_quenchingTs()[index].clear();
            //******************************************************************
            
            if (get_is_avgFitParams()) {
                
                if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg))) {
                    
                    for (int i=0; i<n_trial; ++i) {
                        for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) {
                            
                            //--------------------------------------------------
                            o.clear();
                            o.append(return_AnalysisFolderPath(sysVar));
                            o.append("/fit_data");
                            o.append("/tempsInfo_");
                            o.append(sysVar.get_usic());
                            o.append("_00"+to_string((long long int)i));
                            o.append("_"+sysVar.get_nameString(ii));
                            o.append(".dat");
                            ofstream write_tempsinfo(output.c_str(),ofstream::app);
                            //--------------------------------------------------
                            // NOTE:
                            // temperatures in this file has been normalized to the "expanded" format
                            // (ex. 1234567 for both real (1234.567) and lj (1.234567) units)
                            // so "setprecision" needs to be set to zero and with fixed format
                            //--------------------------------------------------
                            write_tempsinfo << fixed << setprecision(0);
                            
                            int index_i=i*n_system+ii;
                            
                            if (sysVar.get_is_updatetinfo()) {
                                sysVar.get_temperaturesInfo()[index_i]=hightolow;}
                            
                            size_t equTSize
                            = sysVar.get_equilibratedTs()[index_i].size();
                            double equT0
                            = sysVar.get_equilibratedTs()[index_i][equTSize-1];
                            //--------------------------------------------------
                            // NOTE:
                            // quenchingT's should be < lowest equilibrated T
                            //--------------------------------------------------
                            int count=0;
                            do {
                                double Ti
                                = sysVar.get_temperaturesInfo()[index_i][count];
                                
                                if (Ti<equT0) {
                                    sysVar.get_quenchingTs()[index_i].push_back(Ti);
                                }
                                ++count;
                            } while (count<sysVar.get_temperaturesInfo()[index_i].size());
                            //--------------------------------------------------
                        }
                    }
                }
                
            }
            else {
                
                size_t equTSize = sysVar.get_equilibratedTs()[index].size();
                double equT0 = sysVar.get_equilibratedTs()[index][equTSize-1];
                //--------------------------------------------------------------
                // NOTE:
                // quenchingT's should be < lowest equilibrated T
                //--------------------------------------------------------------
                int count=0;
                do {
                    double Ti=sysVar.get_temperaturesInfo()[index][count];
                    
                    if (Ti<equT0) {
                        sysVar.get_quenchingTs()[index].push_back(Ti);
                    }
                    ++count;
                } while (count<sysVar.get_temperaturesInfo()[index].size());
                //--------------------------------------------------------------
            }
        }
        
        // write temperaturesInfo
        //----------------------------------------------------------------------
        if (sysVar.times_write_tempsinfo<(sysVar.max_times_retry+1)) {
            
            if (get_is_avgFitParams()) {
                
                if (index==((n_trial-1)*n_system+(n_sys_end-n_sys_beg))) {
                    
                    for (int i=0; i<n_trial; ++i) {
                        for (int ii=0; ii<=(n_sys-n_sys_beg); ++ii) {
                            
                            int index_i=i*n_system+ii;
                            for (size_t iii=0; iii<sysVar.get_temperaturesInfo()[index_i].size(); ++iii) {
                                write_tempsinfo << sysVar.get_temperaturesInfo()[index_i][iii] << "\n";
                            }
                        }
                    }
                }
            }
            else {
                for (size_t i=0; i<sysVar.get_temperaturesInfo()[index].size(); ++i) {
                    write_tempsinfo << sysVar.get_temperaturesInfo()[index][i] << "\n";
                }
            }
        }
        //----------------------------------------------------------------------
        
        write_tempsinfo.close();
    }
    else {
        cout << "shootForNewTemperatures: 'tempsInfo.dat' cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
    }
}





/*==( public setters )==*/

/* string */
void FitData::set_analysispart(const std::string& str){analysispart=str;}
void FitData::set_relaxation_target(const std::string& str){relaxation_target=str;}

/* bool */
void FitData::set_is_fit_sExp(const bool b){is_fit_sExp=b;}
void FitData::set_is_find_DWF(const bool b){is_find_DWF=b;}
void FitData::set_is_fit_Arrhenius(const bool b){is_fit_Arrhenius=b;}
void FitData::set_is_fit_tauFit(const bool b){is_fit_tauFit=b;}
void FitData::set_is_shootForNext(const bool b){is_shootForNext=b;}
void FitData::set_is_avgFitParams(const bool b){is_avgFitParams=b;}
void FitData::set_is_fit_by_TA(const bool b){is_fit_by_TA=b;}
void FitData::set_is_fitByEveryPoint(const bool b){is_fitByEveryPoint=b;}
void FitData::set_is_fit_by_Tc(const bool b){is_fit_by_Tc=b;}
void FitData::set_is_applytauFitcut(const bool b){is_applytauFitcut=b;}
void FitData::set_is_fit_full_alpha(const bool b){is_fit_full_alpha=b;}
void FitData::set_is_use_gammafunc(const bool b){is_use_gammafunc=b;}
void FitData::set_is_imposeStrictTeq(const bool b){is_imposeStrictTeq=b;}
void FitData::set_is_calc_thermoData(const bool b){is_calc_thermoData=b;}

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
void FitData::set_deviateTA(const double d){deviateTA=d;}
void FitData::set_Tc_MCT(const double d){Tc_MCT=d;}
void FitData::set_Tc_SOU(const double d){Tc_SOU=d;}
void FitData::set_tauFit_cut(const double d){tauFit_cut=d;}
