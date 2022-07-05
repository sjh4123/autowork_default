//
//  theorytest.cpp
//  cppWork
//
//  Created by SJH on 9/30/15.
//  Copyright (c) 2015 Sean Jh H. All rights reserved.
//

#include "theorytest.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

TheoryTest::TheoryTest(const StructureClass& sysVar,
                       const WorkScripts& ws,
                       const AmdatAnalysis& aa,
                       const vector<vector<double>>& tinfo_org,
                       const string& model_org):

/* Base Class (Fitdata) Initialization */
FitData(sysVar,ws,aa),

/* bool */
is_use_ngp_peak_frame(true),
is_use_counting1(true),
is_data_smoothing(true),

/* int */
ngp_peak_frame(0),
ngp_block_frame(0),
peak_frame(0),

/* double */
frame_time(0),
mean_strings(0),
mean_length(0),
mean_length_count1(0),
order_parameter(0),
stringlen(0),
peak_time(0),
fast_threshold(0.15),
fit_tau_threshold(3.0),
stringlen_threshold(0.0),
n_distr_cutoff(1e-2),
p1_extrp(0),

delmiu(0),
delHa(0),
delSa(0),

/* string */
model_stringlen_distr("exp_string"), // exp_string, pwr_exp_string
model(model_org),

/* temperatures of current regime */
tinfo(tinfo_org)
{
    
    /* Vibrational relaxation time (approx. 0.1 ps) */
    if (sysVar.get_systemUnit()=="real")
    {
        tau_hiT=100; /* time: fs */
    }
    else if (sysVar.get_systemUnit()=="lj")
    {
        tau_hiT=0.1; /* time: tau */
    }
    else if (sysVar.get_systemUnit()=="metal")
    {
        tau_hiT=0.1; /* time: ps */
    }
}





void TheoryTest::amdatstringanalysis(const StructureClass& sysVar,
                                     const WorkScripts& ws,
                                     AmdatAnalysis& aa)
{
    
    aa.set_is_strFac(false);
    aa.set_is_msd(false);
    aa.set_is_ngp(false);
    aa.set_is_isfs(false);
    aa.set_is_composition(false);
    aa.set_is_u2dist(false);
    aa.set_is_stiffness_dist(false);
    aa.set_is_isf(false);
    
    aa.set_is_strings(true); // NOTE!
    
    string o=test_AnalysisFolder(sysVar,aa);
    //call_system_bash(o);
    
    vector<string> amdat_targets;
    
    /* find L(T) at peak frame of ngp */
    //==========================================================================
    if (is_use_ngp_peak_frame) {
        
        aa.set_is_peak_frame(true);
        aa.make_amdatInputFile(sysVar,ws);
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    
                    /* find highest peak of NGP */
                    find_ngp_peak_frame(sysVar,n_trl,n_sys,tinfo[index][T]);
                    
                    /* AMDAT submission files */
                    aa.make_amdatSubFile
                    (sysVar,ws,n_trl,n_sys,tinfo[index][T],ngp_peak_frame);
                    
                    /* Specify target files to be watched for */
                    amdat_targets.push_back
                    (aa.amdat_target(sysVar,n_trl,n_sys,tinfo[index][T],"strings"));
                }
            }
        }
        /* AMDAT Jobs Submissions */
        aa.make_amdatSubScript(sysVar,tinfo);
        
        /* Watch and Hold */
        if (sysVar.get_is_watch_hold()) {
            autoWork::watch_hold
            (sysVar,amdat_targets,"analysis",current_regime);}
        amdat_targets.clear();
        
        
        /* write out stringlen to file (all Ts) */
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    find_mean_stringlen
                    (sysVar,n_trl,n_sys,tinfo[index][T]);
                    //find_string_length
                    //(sysVar,n_trl,n_sys,tinfo[index][T]);
                    //find_fastfrac
                    //(sysVar,n_trl,n_sys,tinfo[index][T]);
                }
            }
        }
        
    } /* find L(T) at peak frame of ngp */
    //==========================================================================
    
    
    
    
    
    /* find L(T) at in-block frames */
    //==========================================================================
    else {
        
        /* find block frame index (same for same regime temperatures) */
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                // locate position of current system
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                find_ngp_block_frame(sysVar,n_trl,n_sys,tinfo[index][0]);
            }
        }
        
        aa.set_is_peak_frame(false);
        //aa.make_amdatInputFile_strings_loop(sysVar,ws,ngp_block_frame); // NOTE
        
        int beg_frame=0;
        int end_frame=0;
        
        vector<vector<int>> calc_frames; // NOTE: 2D int vector
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    
                    /* find highest peak of NGP */
                    find_ngp_peak_frame(sysVar,n_trl,n_sys,tinfo[index][T]);
                    
                    vector<int> tmp=aa.make_amdatInputFile_strings_loop
                    (sysVar,ws,n_trl,n_sys,tinfo[index][T],ngp_peak_frame,ngp_block_frame);
                    
                    beg_frame=tmp[0];
                    end_frame=tmp[1];
                    calc_frames.push_back({beg_frame,end_frame});
                    
                    /* AMDAT submission files */
                    aa.make_amdatSubFile
                    (sysVar,ws,n_trl,n_sys,tinfo[index][T],ngp_block_frame);
                    
                    /* Specify target files to be watched for */
                    amdat_targets.push_back
                    (aa.amdat_target(sysVar,n_trl,n_sys,tinfo[index][T],"strings",end_frame));
                }
            }
        }
        /* AMDAT Jobs Submissions */
        aa.make_amdatSubScript(sysVar,tinfo);
        
        /* Watch and Hold */
        if (sysVar.get_is_watch_hold()) {
            autoWork::watch_hold
            (sysVar,amdat_targets,"analysis",current_regime);}
        amdat_targets.clear();
        
        
        /* write out stringlen to file (all Ts) */
        int index1=0;
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    beg_frame=calc_frames[index1][0];
                    end_frame=calc_frames[index1][1];
                    ++index1;
                    for (int frame_index=beg_frame; frame_index<=end_frame; ++frame_index) {
                        find_time_stringlen
                        (sysVar,n_trl,n_sys,tinfo[index][T],frame_index);
                        //find_time_fastfrac
                        //(sysVar,n_trl,n_sys,tinfo[index][T],frame_index);
                    }
                    write_peak_stringlen
                    (sysVar,n_trl,n_sys,tinfo[index][T]);
                    //write_peak_fastfrac
                    //(sysVar,n_trl,n_sys,tinfo[index][T]);
                }
            }
        }
        
    } /* find L(T) at in-block frames */
    //==========================================================================
    
}





void TheoryTest::AGtest(const StructureClass& sysVar)
{
    
    //write_fastfrac_equ(sysVar);
    //write_fastfrac_equ_avg(sysVar);
    
    write_stringlen_equ(sysVar);
    write_stringlen_equ_avg(sysVar);
    
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/AGtest_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_AGtestraw(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/AGtest_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_AGtest(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/AGtest_log_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_AGtest_log(o.c_str(),ofstream::app); // 'append'
        
        bool is_cutbyTA=false;
        //----------------------------------------------------------------------
        double TA_d=0;
        if (is_cutbyTA) {
            /* find TA */
            read_all_taueq_data(sysVar,n_sys);
            find_TA(sysVar,n_sys,taueq_sortdecreasing);
            TA_d=TA_avg;
        }
        else {
            read_all_taueq_data(sysVar,n_sys);
            TA_d=taueq_sortdecreasing[0][0];
        }
        //----------------------------------------------------------------------
        
        int n_taueq=0;
        int strSize=0;
        
        if (is_data_smoothing) {
            fit_taueq(sysVar,n_sys,TA_d);
            fit_stringlen(sysVar,n_sys,TA_d);
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            strSize=(int)stringlen_sortdecreasing_fit.size();
        }
        else {
            read_all_taueq_data(sysVar,n_sys,TA_d);
            read_all_equ_stringlen(sysVar,n_sys,TA_d);
            n_taueq=(int)taueq_sortdecreasing.size();
            strSize=(int)stringlen_sortdecreasing.size();
        }
        
        if (n_taueq!=strSize) {
            cout << "n_taueq .NOT. equal to strSize!" << "\n";
            //n_taueq=strSize;
        }
        
        cout
        << "\n"
        << "n_taueq " << n_taueq << "\n"
        << "strSize " << strSize << "\n";
        
        for (int outter=0; outter<2; ++outter) {
            
            for (int inner=0; inner<n_taueq; ++inner) {
                
                int count_data=0;
                
                /* tau_alpha(T), L(T) */
                //-----------------------------------------------------
                if (is_data_smoothing) {
                    T_actual=taueq_sortdecreasing_fit[inner][0];
                    tau_T   =taueq_sortdecreasing_fit[inner][1];
                    for (int i=0; i<strSize; ++i) {
                        if (T_actual==stringlen_sortdecreasing_fit[i][0]) {
                            ++count_data;
                            stringlen=stringlen_sortdecreasing_fit[i][1];
                            break;
                        }
                    }
                }
                else {
                    read_all_taueq_data(sysVar,n_sys,TA_d);
                    read_all_equ_stringlen(sysVar,n_sys,TA_d);
                    T_actual=taueq_sortdecreasing[inner][0];
                    tau_T   =taueq_sortdecreasing[inner][1];
                    for (int i=0; i<strSize; ++i) {
                        if (T_actual==stringlen_sortdecreasing[i][0]) {
                            ++count_data;
                            stringlen=stringlen_sortdecreasing[i][1];
                            break;
                        }
                    }
                }
                
                /* find tau0 (from model fit) */
                //-----------------------------------------------------
                find_all_tau0(sysVar,n_sys);
                
                
                /* Activation energy from COOP */
                //-----------------------------------------------------
                if (model=="COOP") {
                    double E_inf=c_COOP[1];
                    double u    =c_COOP[2];
                    double b    =c_COOP[3];
                    Ea_COOP=E_inf*(1+exp(-u*((T_actual/E_inf)-b)));
                }
                
                
                /* find L(TA) */
                //-----------------------------------------------------
                find_LA(sysVar,n_sys);
                
                
                /* find Ea/k (from Arrhenius fit) */
                //-----------------------------------------------------
                find_all_Arrhenius(sysVar,n_sys);
                
                
                /* find delta miu (high T activation free energy) */
                //-----------------------------------------------------
                double delmiu_TA=0;
                //double delmiu_use=0;
                
                delmiu_TA = TA*log(tauA/tau_hiT); // activation free energy at TA
                delHa     = Ea;                   // actiavtion enthalpy
                delSa     = (delHa-delmiu_TA)/TA; // activation entropy
                
                delG      = T_actual*log(pow(10,tau_T)/tau_hiT); // delG(T)
                delmiu    = delHa-T_actual*delSa;                // delmiu(T)
                
                /* the actual delmiu used in string-AG test */
                //delmiu_use = delmiu;
                
                
                /* in first loop, collect raw data for AG test */
                //-----------------------------------------------------
                if(outter==0) {
                    
                    if (count_data==0) continue;
                    
                    if ((tau_T<fit_tau_threshold)&&
                        (T_actual<=TA)&&
                        (stringlen>stringlen_threshold)) { // NOTE!
                        
                        write_AGtestraw
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << tau_T*log(10) // NOTE: natural log form
                        << "\n";
                    }
                }
                
                /* in second loop, write out AG test results */
                //-----------------------------------------------------
                else if (outter==1) {
                    
                    if (inner==0) {
                        
                        write_AGtestraw.close();
                        /* NOTE: need to close before it can be read */
                        
                        //------------------------------------------------------
                        write_AGtest
                        << "tau0("<<model<<") = " << tau0_model     << "\n"
                        << "tau0(Arrhenius) = "   << tau0_Arrhenius << "\n"
                        << "delHa   = "           << delHa          << "\n"
                        << "delSa   = "           << delSa          << "\n"
                        << "TA      = "           << TA             << "\n"
                        << "tauA    = "           << tauA           << "\n"
                        << "L(TA)   = "           << LA             << "\n\n";
                        
                        write_AGtest_log
                        << "tau0("<<model<<") = " << tau0_model     << "\n"
                        << "tau0(Arrhenius) = "   << tau0_Arrhenius << "\n"
                        << "delHa   = "           << delHa          << "\n"
                        << "delSa   = "           << delSa          << "\n"
                        << "TA      = "           << TA             << "\n"
                        << "tauA    = "           << tauA           << "\n"
                        << "L(TA)   = "           << LA             << "\n\n";
                        //------------------------------------------------------
                        
                        find_AG_tau0_fit(sysVar,n_sys);
                        
                        //------------------------------------------------------
                        write_AGtest
                        << "tau0_fit= "           << exp(tau0_fit)  << "\n"
                        << "slope   = "           << slope_fit      << "\n"
                        << "R^2     = "           << rep.r2         << "\n\n"
                        << "T "
                        << (1000.0/(corrFacforT/precision))
                        << "/T "
                        << "tau(T) L(T) TA/T L/LA delG/delmiu delmiu delG Ea_COOP T(-delS) (L/LA)(delmiu/kT) Ln(tau/tau0_fit)"
                        << "\n";
                        
                        write_AGtest_log
                        << "tau0_fit= "           << exp(tau0_fit)  << "\n"
                        << "slope   = "           << slope_fit      << "\n"
                        << "R^2     = "           << rep.r2         << "\n\n"
                        << "T "
                        << (1000.0/(corrFacforT/precision))
                        << "/T "
                        << "tau(T) L(T) TA/T L/LA delG/delmiu delmiu delG Ea_COOP T(-delS) (L/LA)(delmiu/kT) Log(tau/tau0_fit)"
                        << "\n";
                        //------------------------------------------------------
                    }
                    
                    if (count_data==0) continue;
                    
                    if (stringlen>stringlen_threshold) {
                        
                        write_AGtest
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delG/delmiu << " "
                        << delmiu << " "
                        << delG << " "
                        << Ea_COOP << " "
                        << delG-Ea_COOP << " "
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << (tau_T*log(10))-tau0_fit
                        << "\n";
                        
                        write_AGtest_log
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delG/delmiu << " "
                        << delmiu << " "
                        << delG << " "
                        << Ea_COOP << " "
                        << delG-Ea_COOP << " "
                        << ((stringlen/LA)*(delmiu/T_actual))*log10(exp(1)) << " "
                        << ((tau_T*log(10))-tau0_fit)*log10(exp(1))
                        << "\n";
                    }
                    
                }
            }
        }
        write_AGtest.close();
        write_AGtest_log.close();
    }
}





void TheoryTest::GLMtest(const StructureClass& sysVar)
{
    
    //string model="GLM"; // 3 free parameters
    string model="GLM1"; // 1 free parameter
    
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/GLMtest_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtestraw(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/GLMtest_fit_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMfit(o.c_str());
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/GLMtest_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtest(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/GLMtest_log_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtest_log(o.c_str(),ofstream::app); // 'append'
        
        /* find <u2(TA)> */
        //-----------------------------------------------------
        find_uA2(sysVar,n_sys);
        const double TA_d=TA;
        
        // Curve-fitting of GLM
        //----------------------------------------------------------------------
        set_fitParams(model);
        read_all_equ_DWF(sysVar,n_sys,TA_d);
        fitdata_processing_lsfit(dwf_tau_data); // NOTE: dwf_tau_data: (DWF,tau)
        alglib_lsfit(model);
        
        const double tau0 = c[0]; // tau0 in GLM, tauA in GLM1
        const double u02  = c[1]; // u02 in GLM, uA2 in GLM1
        const double alpha= c[2];
        
        write_fitModel(write_GLMfit,model);
        write_fitarrays(write_GLMfit);
        write_stopcond(write_GLMfit);
        write_fitinfo(write_GLMfit);
        write_errcurve_actual(write_GLMfit,model,dwf_sortdecreasing);
        
        write_GLMfit.close();
        //----------------------------------------------------------------------
        
        int n_taueq=0;
        int dwfSize=0;
        
        if (is_data_smoothing) {
            fit_taueq(sysVar,n_sys,TA_d);
            fit_DWF(sysVar,n_sys,TA_d);
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            dwfSize=(int)dwf_sortdecreasing_fit.size();
        }
        else {
            read_all_taueq_data(sysVar,n_sys,TA_d);
            read_all_equ_DWF(sysVar,n_sys,TA_d);
            n_taueq=(int)taueq_sortdecreasing.size();
            dwfSize=(int)dwf_sortdecreasing.size();
        }
        
        if (n_taueq!=dwfSize) {
            cout << "n_taueq .NOT. equal to dwfSize!" << "\n";
            //n_taueq=strSize;
        }
        
        cout
        << "\n"
        << "n_taueq " << n_taueq << "\n"
        << "dwfSize " << dwfSize << "\n";
        
        for (int outter=0; outter<2; ++outter) {
            
            for (int inner=0; inner<n_taueq; ++inner) {
                
                int count_data=0;
                
                /* tau_alpha(T), DWF */
                //-----------------------------------------------------
                if (is_data_smoothing) {
                    T_actual=taueq_sortdecreasing_fit[inner][0];
                    tau_T   =taueq_sortdecreasing_fit[inner][1];
                    for (int i=0; i<dwfSize; ++i) {
                        if (T_actual==dwf_sortdecreasing_fit[i][0]) {
                            ++count_data;
                            DWF=dwf_sortdecreasing_fit[i][1];
                            break;
                        }
                    }
                }
                else {
                    read_all_taueq_data(sysVar,n_sys,TA_d);
                    read_all_equ_DWF(sysVar,n_sys,TA_d);
                    T_actual=taueq_sortdecreasing[inner][0];
                    tau_T   =taueq_sortdecreasing[inner][1];
                    for (int i=0; i<dwfSize; ++i) {
                        if (T_actual==dwf_sortdecreasing[i][0]) {
                            ++count_data;
                            DWF=dwf_sortdecreasing[i][1];
                            break;
                        }
                    }
                }
                
                /* in first loop, collect raw data for GLM test */
                //-----------------------------------------------------
                if(outter==0) {
                    
                    if (count_data==0) continue;
                    
                    if ((tau_T<fit_tau_threshold)&&
                        (T_actual<=TA_d)) { // NOTE!
                        
                        if (model=="GLM") {
                            write_GLMtestraw
                            << pow(u02/DWF,alpha/2) << " "
                            << tau_T*log(10) // NOTE: natural log form
                            << "\n";
                        }
                        else if (model=="GLM1") {
                            write_GLMtestraw
                            << pow(u02/DWF,alpha/2)-1.0 << " "
                            << tau_T*log(10) // NOTE: natural log form
                            << "\n";
                        }
                        
                    }
                }
                
                /* in second loop, write out AG test results */
                //-----------------------------------------------------
                else if (outter==1) {
                    
                    if (inner==0) {
                        
                        write_GLMtestraw.close();
                        /* NOTE: need to close before it can be read */
                        
                        //------------------------------------------------------
                        if (model=="GLM") {
                            
                            write_GLMtest
                            << "tau0("<<model<<") = " << tau0  << "\n"
                            << "u0^2    = "           << u02   << "\n"
                            << "alpha   = "           << alpha << "\n\n";
                            
                            write_GLMtest_log
                            << "tau0("<<model<<") = " << tau0  << "\n"
                            << "u0^2    = "           << u02   << "\n"
                            << "alpha   = "           << alpha << "\n\n";
                        }
                        else if (model=="GLM1") {
                            
                            write_GLMtest
                            << "TA    = "             << TA_d  << "\n"
                            << "tauA  = "             << tau0  << "\n"
                            << "uA^2  = "             << u02   << "\n"
                            << "alpha = "             << alpha << "\n\n";
                            
                            write_GLMtest_log
                            << "TA    = "             << TA_d  << "\n"
                            << "tauA  = "             << tau0  << "\n"
                            << "uA^2  = "             << u02   << "\n"
                            << "alpha = "             << alpha << "\n\n";
                        }
                        //------------------------------------------------------
                        
                        find_GLM_tau0_fit(sysVar,n_sys);
                        
                        //------------------------------------------------------
                        if (model=="GLM") {
                            
                            write_GLMtest
                            << "tau0_fit = " << exp(tau0_fit)   << "\n"
                            << "slope    = " << slope_fit      << "\n"
                            << "R^2      = " << rep.r2         << "\n\n"
                            << "T "
                            << (1000.0/(corrFacforT/precision))
                            << "/T "
                            << "tau(T) DWF(T) (u02/u2)^(a/2) Ln(tau/tau0)"
                            << "\n";
                            
                            write_GLMtest_log
                            << "tau0_fit = " << exp(tau0_fit)   << "\n"
                            << "slope    = " << slope_fit      << "\n"
                            << "R^2      = " << rep.r2         << "\n\n"
                            << "T "
                            << (1000.0/(corrFacforT/precision))
                            << "/T "
                            << "tau(T) DWF(T) (u02/u2)^(a/2) Log(tau/tau0)"
                            << "\n";
                        }
                        else if (model=="GLM1") {
                            
                            write_GLMtest
                            << "tauA_fit = " << exp(tau0_fit)   << "\n"
                            << "slope    = " << slope_fit      << "\n"
                            << "R^2      = " << rep.r2         << "\n\n"
                            << "T "
                            << (1000.0/(corrFacforT/precision))
                            << "/T "
                            << "tau(T) DWF(T) (uA2/u2)^(a/2)-1 Ln(tau/tauA)"
                            << "\n";
                            
                            write_GLMtest_log
                            << "tauA_fit = " << exp(tau0_fit)   << "\n"
                            << "slope    = " << slope_fit      << "\n"
                            << "R^2      = " << rep.r2         << "\n\n"
                            << "T "
                            << (1000.0/(corrFacforT/precision))
                            << "/T "
                            << "tau(T) DWF(T) (uA2/u2)^(a/2)-1 Log(tau/tauA)"
                            << "\n";
                        }
                        //------------------------------------------------------
                    }
                    
                    if (count_data==0) continue;
                    
                    if (model=="GLM") {
                        
                        write_GLMtest
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << DWF << " "
                        << pow(u02/DWF,alpha/2) << " "
                        << (tau_T*log(10))-tau0_fit
                        //<< (tau_T*log(10))-tau0
                        << "\n";
                        
                        write_GLMtest_log
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << DWF << " "
                        << (pow(u02/DWF,alpha/2))*log10(exp(1)) << " "
                        << ((tau_T*log(10))-tau0_fit)*log10(exp(1))
                        //<< ((tau_T*log(10))-tau0)*log10(exp(1))
                        << "\n";
                    }
                    else if(model=="GLM1") {
                        
                        write_GLMtest
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << DWF << " "
                        << pow(u02/DWF,alpha/2)-1.0 << " "
                        << (tau_T*log(10))-tau0_fit
                        //<< (tau_T*log(10))-tau0
                        << "\n";
                        
                        write_GLMtest_log
                        << T_actual << " "
                        << (1000.0/(corrFacforT/precision))/(T_actual) << " "
                        << tau_T << " "
                        << DWF << " "
                        << (pow(u02/DWF,alpha/2)-1.0)*log10(exp(1)) << " "
                        << ((tau_T*log(10))-tau0_fit)*log10(exp(1))
                        //<< ((tau_T*log(10))-tau0)*log10(exp(1))
                        << "\n";
                    }
                }
            }
        }
        write_GLMtest.close();
        write_GLMtest_log.close();
    }
}





void TheoryTest::find_ngp_peak_frame(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d)
{
    ngp_data_sortdecreasing.clear();
    vector<vector<double>> ngp_data;
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/statistics/ngp/");
    o.append(get_analysispart());
    o.append("/peakframe_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append("."+get_analysispart()+".dat");
    //cout << o << "\n";
    ofstream write_peak(o.c_str());
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/ngp/");
    i.append(get_analysispart());
    i.append("/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream read_ngp(i.c_str());
    if (read_ngp.is_open()) {
        
        string lineContent;
        string amdat_version;
        string strData;
        int    frame_index=0; // zero-based index
        double time=0,ngp=0,time_collected=0;
        
        /* read 1st line */
        getline(read_ngp,lineContent);
        amdat_version = lineContent;
        
        while (getline(read_ngp,lineContent)) {
            
            istringstream iss(lineContent);
            iss >> time; // time
            
            // only collect "in-block" data
            // i.e. 0 < time <= time(first block)
            
            if ((time>0)&&(time<=tauFit_cut)) {
                time_collected=time;
                iss >> ngp;  // non-Gaussian parameter, alpha2(time)
                ngp_data.push_back({time,ngp,(double)frame_index});
                
                write_peak
                << frame_index << " "
                << time_collected << " "
                << ngp
                << "\n";
            }
            
            ++frame_index; // NOTE: zero-based
        }
        
        write_peak
        << "\n"
        << "final collected time = " << time_collected << "\n"
        << "final time           = " << time << "\n";
        
        read_ngp.close();
    }
    else {
        cout
        << "read_individual_ngp_data: ngp file cannot open."
        << "\n";
    }
    
    try {
        if (ngp_data.size()==0) throw 0;
        std::sort(ngp_data.begin(),ngp_data.end(),sortDecreasing1);
    }
    catch (int i) {
        cout
        << "in read_individual_ngp_data:" << "\n"
        << "ngp size = " << i << "\n";
        ngp_data.push_back({0,0});
        std::sort(ngp_data.begin(),ngp_data.end(),sortDecreasing1);
    }
    
    ngp_data_sortdecreasing = ngp_data;
    
    // (time,ngp,frame_index)
    ngp_peak_frame = ngp_data[0][2];
    
    write_peak
    << "\n"
    << "time  ngp  frame" << "\n"
    << ngp_data[0][0] << " "
    << ngp_data[0][1] << " "
    << ngp_data[0][2] << "\n";
    
    write_peak.close();
}





void TheoryTest::find_ngp_block_frame(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    vector<vector<double>> ngp_data;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/ngp/");
    i.append(get_analysispart());
    i.append("/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream read_ngp(i.c_str());
    if (read_ngp.is_open()) {
        
        string lineContent;
        string amdat_version;
        string strData;
        int    frame_index=0; // zero-based index
        double time=0,ngp=0,time_collected=0;
        
        /* read 1st line */
        getline(read_ngp,lineContent);
        amdat_version = lineContent;
        
        while (getline(read_ngp,lineContent)) {
            
            istringstream iss(lineContent);
            iss >> time; // time
            
            // only collect "in-block" data
            // i.e. 0 < time <= time(first block)
            
            if ((time>0)&&(time<=tauFit_cut)) {
                time_collected=time;
                iss >> ngp;
                ngp_data.push_back({time,ngp,(double)frame_index});
            }
            
            ++frame_index; // NOTE: zero-based
        }
        
        read_ngp.close();
    }
    else {
        cout
        << "read_individual_ngp_data: ngp file cannot open."
        << "\n";
    }
    
    // (time,ngp,frame_index)
    size_t s=ngp_data.size();
    ngp_block_frame = ngp_data[s-1][2];
}





void TheoryTest::find_mean_stringlen(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d)
{
    // Read-in stringlen distr. data & do curve-fitting & extrapolation at n=1
    //--------------------------------------------------------------------------
    vector<vector<double>> n_distr; // (n,P(n))
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/strings_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_string(o.c_str());
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strings/");
    i.append(get_analysispart());
    i.append("/strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream read_strings(i.c_str());
    if (read_strings.is_open()) {
        
        string lineContent;
        string strData;
        vector<string> vec_strData;
        double dubData=0;
        vector<double> vec_dubData;
        int intData=0;
        vector<int> vec_intData;
        
        /* 1st line */
        // amdat version
        getline(read_strings,lineContent);
        write_string << lineContent << "\n";
        
        /* 2nd line */
        getline(read_strings,lineContent);
        istringstream iss(lineContent);
        // first five are statistics numbers
        for (int i=0; i<5; ++i) {
            iss >> strData;
            vec_strData.push_back(strData);
        }
        // the rest are string length distribution
        while (iss>>intData) {
            vec_intData.push_back(intData);
        }
        
        /* 3rd line */
        getline(read_strings,lineContent);
        istringstream iss1(lineContent);
        for (int i=0; i<5; ++i) {
            iss1 >> dubData;
            write_string
            << vec_strData[i] << " " << dubData << "\n";
            
            if      (i==0) frame_time         = dubData;
            else if (i==1) mean_strings       = dubData;
            else if (i==2) mean_length        = dubData;
            else if (i==3) mean_length_count1 = dubData;
            else if (i==4) order_parameter    = dubData;
        }
        // the rest are string length distribution
        
        write_string
        << "\n"
        << "Original Distribution of string length"
        << "\n";
        int count_index=0;
        while (iss1>>dubData) {
            write_string
            << vec_intData[count_index] << " " << dubData << "\n";
            n_distr.push_back({(double)count_index,dubData});
            ++count_index;
        }
        write_string << "\n";
        
        if (is_use_counting1) stringlen=mean_length_count1;
        else stringlen=mean_length;
        
        read_strings.close();
        
    }
    else {
        cout
        << "find_mean_stringlen: file cannot open."
        << "\n";
    }
    
    
    
    if (true) { // do the extrapolation to n=1 to get P(n=1)
        
        /* impose cutoff on P(n) */
        vector<vector<double>> fit_n_distr;
        for (int i=0; i<(int)n_distr.size(); ++i) {
            if (n_distr[i][1]>n_distr_cutoff) {
                fit_n_distr.push_back({n_distr[i][0],n_distr[i][1]});
            }
        }
        try {
            if(fit_n_distr.size()==0) throw 0;
            std::sort(fit_n_distr.begin(),fit_n_distr.end(),sortDecreasing);
        }
        catch (int i) {
            
            cout
            <<
            sysVar.get_usic()+"_00"+to_string((long long int)n_trl)+
            "_"+sysVar.get_nameString(n_sys)+"_T"+to_string((long long int)Temp_d)
            << "\n"
            << "in find_mean_stringlen:" << "\n"
            << "fit_n_distr size = " << i << "\n";
            fit_n_distr.push_back({0,0});
            std::sort(fit_n_distr.begin(),fit_n_distr.end(),sortDecreasing);
        }
        
        /* fit P(n) to a chosen model */
        set_fitParams(model_stringlen_distr);
        fitdata_processing_lsfit(fit_n_distr);
        alglib_lsfit(model_stringlen_distr);
        write_fitModel(write_string,model_stringlen_distr);
        write_fitarrays(write_string);
        write_stopcond(write_string);
        write_fitinfo(write_string);
        
        p1_extrp=calc_time_given_temp(c,1.0,model_stringlen_distr); // P(n=1)_extrapolated
        
        n_distr[1][1]=p1_extrp; // set value to P(n=1)
        
        stringlen=0; // NOTE
        
        double sum_distr=0;
        for (int i=1; i<(int)n_distr.size(); ++i) { // from n=1
            sum_distr += n_distr[i][1];
        }
        for (int i=1; i<(int)n_distr.size(); ++i) { // from n=1
            stringlen += (double)i*(n_distr[i][1]/sum_distr);
        }
        
        write_string
        << "mean_length "
        << stringlen
        << "\n";
        
        write_string
        << "\n"
        << "Fit Distribution of string length" << "\n"
        << "0 0" << "\n";
        for (int i=1; i<(int)n_distr.size(); ++i) { // from n=1
            write_string
            << i
            << " "
            << n_distr[i][1]/sum_distr
            << "\n";
        }
        
    }
    //--------------------------------------------------------------------------
    
    
    
    // Write-out processed stringlen distr. data to file
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T stringlen frame_index time
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringlen
    << stringlen
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_stringleninv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << stringlen
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_string.close();
    write_stringlen.close();
    write_stringleninv.close();
    //--------------------------------------------------------------------------
}





void TheoryTest::find_string_length(const StructureClass& sysVar,
                                    const int n_trl,
                                    const int n_sys,
                                    const double Temp_d)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/strings_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_string(o.c_str());
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strings/");
    i.append(get_analysispart());
    i.append("/strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream read_strings(i.c_str());
    if (read_strings.is_open()) {
        
        string lineContent;
        string strData;
        vector<string> vec_strData;
        double dubData=0;
        vector<double> vec_dubData;
        int intData=0;
        vector<int> vec_intData;
        
        /* 1st line */
        // amdat version
        getline(read_strings,lineContent);
        write_string << lineContent << "\n";
        
        /* 2nd line */
        getline(read_strings,lineContent);
        istringstream iss(lineContent);
        // first five are statistics numbers
        for (int i=0; i<5; ++i) {
            iss >> strData;
            vec_strData.push_back(strData);
        }
        // the rest are string length distribution
        while (iss>>intData) {
            vec_intData.push_back(intData);
        }
        
        /* 3rd line */
        getline(read_strings,lineContent);
        istringstream iss1(lineContent);
        for (int i=0; i<5; ++i) {
            iss1 >> dubData;
            write_string
            << vec_strData[i] << " " << dubData << "\n";
            
            if      (i==0) frame_time         = dubData;
            else if (i==1) mean_strings       = dubData;
            else if (i==2) mean_length        = dubData;
            else if (i==3) mean_length_count1 = dubData;
            else if (i==4) order_parameter    = dubData;
        }
        // the rest are string length distribution
        write_string
        << "\n"
        << "Distribution of string length" << "\n";
        int count_index=0;
        while (iss1>>dubData) {
            write_string
            << vec_intData[count_index] << " " << dubData << "\n";
            ++count_index;
        }
        
        if (is_use_counting1) stringlen=mean_length_count1;
        else stringlen=mean_length;
        
        read_strings.close();
        write_string.close();
    }
    else {
        cout
        << "find_string_length: file cannot open."
        << "\n";
    }
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T stringlen frame_index time
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_string
    << stringlen
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_stringleninv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << stringlen
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_stringlen.close();
    write_stringleninv.close();
}





void TheoryTest::find_time_stringlen(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const int frame)
{
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strings/");
    i.append(get_analysispart());
    i.append("/strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("_frame"+to_string((long long int)frame)); // NOTE!
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    
    ifstream read_strings(i.c_str());
    if (read_strings.is_open()) {
        
        string lineContent;
        string strData;
        vector<string> vec_strData;
        double dubData=0;
        vector<double> vec_dubData;
        int intData=0;
        vector<int> vec_intData;
        
        /* 1st line */
        // amdat version
        getline(read_strings,lineContent);
        
        
        /* 2nd line */
        getline(read_strings,lineContent);
        istringstream iss(lineContent);
        // first five are statistics numbers
        for (int i=0; i<5; ++i) {
            iss >> strData;
            vec_strData.push_back(strData);
        }
        // the rest are string length distribution
        while (iss>>intData) {
            vec_intData.push_back(intData);
        }
        
        
        /* 3rd line */
        getline(read_strings,lineContent);
        istringstream iss1(lineContent);
        for (int i=0; i<5; ++i) {
            iss1 >> dubData;
            
            if      (i==0) frame_time         = dubData;
            else if (i==1) mean_strings       = dubData;
            else if (i==2) mean_length        = dubData;
            else if (i==3) mean_length_count1 = dubData;
            else if (i==4) order_parameter    = dubData;
        }
        
        if (is_use_counting1) stringlen=mean_length_count1;
        else stringlen=mean_length;
        
        read_strings.close();
    }
    else {
        cout
        << "find_time_stringlen: file cannot open."
        << "\n";
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/strings_time_stringlen_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_string(o.c_str(),ofstream::app); // 'append'
    
    // (frame_index,stringlen,time)
    
    write_string
    << frame
    << " "
    << stringlen
    << " "
    << frame_time
    << "\n";
    
    write_string.close();
}





void TheoryTest::find_fastfrac(const StructureClass& sysVar,
                               const int n_trl,
                               const int n_sys,
                               const double Temp_d)
{
    double fast_fraction=0;
    double slow_fraction=0;
    vector<vector<double>> vec_fraction;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strings/");
    i.append(get_analysispart());
    i.append("/gaussian_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    ifstream readgaussian(i.c_str());
    
    if (readgaussian.is_open()) {
        string lineContent;
        double dubData=0;
        string strData;
        
        /* 1-6 lines */
        for (int i=1; i<=6; ++i) {
            getline(readgaussian,lineContent);
            if (i==3) {
                istringstream iss(lineContent);
                iss >> dubData;
                frame_time=dubData;
            }
        }
        // Gaussian slow_frac fast_frac
        // Actual   slow_frac fast_frac
        // Excess   slow_frac fast_frac
        for (int i=0; i<3; ++i) {
            getline(readgaussian,lineContent);
            istringstream iss(lineContent);
            iss >> strData;
            iss >> dubData; slow_fraction=dubData;
            iss >> dubData; fast_fraction=dubData;
            vec_fraction.push_back({slow_fraction,fast_fraction});
        }
        readgaussian.close();
    }
    else {
        cout
        << "find_fastfrac: gaussian file cannot open."
        << "\n";
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/fastfrac_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fastfrac(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/fastfrac_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fastfracinv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T fast_frac frame_index time
    streamsize ss=write_fastfrac.precision();
    
    write_fastfrac
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_fastfrac.precision(ss);
    write_fastfrac << resetiosflags(ios::fixed|ios::showpoint);
    
    write_fastfrac
    << vec_fraction[1][1]
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_fastfracinv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << vec_fraction[1][1]
    << " "
    << ngp_peak_frame
    << " "
    << frame_time
    << "\n";
    
    write_fastfrac.close();
    write_fastfracinv.close();
}





void TheoryTest::find_time_fastfrac(const StructureClass& sysVar,
                                    const int n_trl,
                                    const int n_sys,
                                    const double Temp_d,
                                    const int frame)
{
    double fast_fraction=0;
    double slow_fraction=0;
    vector<vector<double>> vec_fraction;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strings/");
    i.append(get_analysispart());
    i.append("/gaussian_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("_frame"+to_string((long long int)frame)); // NOTE
    i.append("."+get_analysispart()+".dat");
    //cout << i << "\n";
    ifstream readgaussian(i.c_str());
    
    if (readgaussian.is_open()) {
        string lineContent;
        double dubData=0;
        string strData;
        
        /* 1-6 lines */
        for (int i=1; i<=6; ++i) {
            getline(readgaussian,lineContent);
            if (i==3) {
                istringstream iss(lineContent);
                iss >> dubData;
                frame_time=dubData;
            }
        }
        // Gaussian slow_frac fast_frac
        // Actual   slow_frac fast_frac
        // Excess   slow_frac fast_frac
        for (int i=0; i<3; ++i) {
            getline(readgaussian,lineContent);
            istringstream iss(lineContent);
            iss >> strData;
            iss >> dubData; slow_fraction=dubData;
            iss >> dubData; fast_fraction=dubData;
            vec_fraction.push_back({slow_fraction,fast_fraction});
        }
        readgaussian.close();
    }
    else {
        cout
        << "find_time_fastfrac: gaussian file cannot open."
        << "\n";
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/strings_time_fastfrac_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d)); // NOTE
    o.append(".dat");
    ofstream write_fastfrac(o.c_str(),ofstream::app); // 'append'
    
    // (frame_index,fastfrac,time)
    
    write_fastfrac
    << frame
    << " "
    << vec_fraction[1][1]
    << " "
    << frame_time
    << "\n";
    
    write_fastfrac.close();
}





void TheoryTest::fit_stringlen(const StructureClass& sysVar,
                               const int n_sys)
{
    string model_string="COOP_string";
    vector<vector<double>> stringlen_fit;
    
    set_fitParams(model_string);
    read_all_equ_stringlen(sysVar,n_sys); // NOTE: read equilibrium strings
    
    double T_fit=0,str_fit=0;
    for (int i=0; i<stringlen_sortdecreasing.size(); ++i) {
        // NOTE:
        // use normalized T for fitting stringlen vs T
        T_fit  =stringlen_sortdecreasing[i][0]*pow(1000.0/(corrFacforT/precision),-1);
        str_fit=stringlen_sortdecreasing[i][1];
        stringlen_fit.push_back({T_fit,str_fit});
    }
    fitdata_processing_lsfit(stringlen_fit);
    alglib_lsfit(model_string);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/AGtest_fit_stringlen_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model_string);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_log(write_tau_inf,model_string,stringlen_fit);
    
    vector<vector<double>> vec2d;
    
    for (int i=0; i<(int)stringlen_sortdecreasing.size(); ++i) {
        
        T_actual=stringlen_sortdecreasing[i][0];
        T_fit   =stringlen_sortdecreasing[i][0]*pow(1000.0/(corrFacforT/precision),-1);
        str_fit =log10(calc_time_given_temp(c,T_fit,model_string));
        
        vec2d.push_back({T_actual,str_fit}); // NOTE: use {T_actual,str_fit}
    }
    stringlen_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_stringlen(const StructureClass& sysVar,
                               const int n_sys,
                               const double TA_d)
{
    string model_string="COOP_string";
    vector<vector<double>> stringlen_fit;
    
    set_fitParams(model_string);
    read_all_equ_stringlen(sysVar,n_sys,TA_d); // NOTE: fit T<=TA
    
    double T_fit=0,str_fit=0,T_actual=0;
    for (int i=0; i<stringlen_sortdecreasing.size(); ++i) {
        // NOTE:
        // use normalized T for fitting stringlen vs T
        T_fit  =stringlen_sortdecreasing[i][0]*pow(1000.0/(corrFacforT/precision),-1);
        str_fit=stringlen_sortdecreasing[i][1];
        stringlen_fit.push_back({T_fit,str_fit});
    }
    fitdata_processing_lsfit(stringlen_fit);
    alglib_lsfit(model_string);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/AGtest_fit_stringlen_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model_string);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_log(write_tau_inf,model_string,stringlen_fit);
    
    vector<vector<double>> vec2d;
    
    for (int i=0; i<(int)stringlen_sortdecreasing.size(); ++i) {
        
        T_actual=stringlen_sortdecreasing[i][0];
        T_fit   =stringlen_sortdecreasing[i][0]*pow(1000.0/(corrFacforT/precision),-1);
        str_fit =log10(calc_time_given_temp(c,T_fit,model_string));
        
        vec2d.push_back({T_actual,str_fit}); // NOTE: use {T_actual,str_fit}
    }
    stringlen_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_DWF(const StructureClass& sysVar,
                         const int n_sys)
{
    string model_DWF="COOP_DWF";
    
    set_fitParams(model_DWF);
    read_all_equ_DWF(sysVar,n_sys);
    fitdata_processing_lsfit(dwf_invT_sortdecreasing); // NOTE: use "dwf_invT"
    alglib_lsfit(model_DWF);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/GLMtest_fit_DWF_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fit_DWF(o.c_str());
    
    write_fitModel(write_fit_DWF,model_DWF);
    write_fitarrays(write_fit_DWF);
    write_stopcond(write_fit_DWF);
    write_fitinfo(write_fit_DWF);
    write_errcurve_log(write_fit_DWF,model_DWF,dwf_invT_sortdecreasing);
    
    vector<vector<double>> vec2d;
    double temp=0,dwfr=0;
    for (int i=0; i<(int)dwf_invT_sortdecreasing.size(); ++i) {
        temp=dwf_sortdecreasing[i][0]; // NOTE: is "dwf" .NOT. dwf_invT
        dwfr=log10(calc_time_given_temp(c,temp,model_DWF));
        vec2d.push_back({temp,dwfr});
    }
    dwf_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_DWF(const StructureClass& sysVar,
                         const int n_sys,
                         const double TA_d)
{
    string model_DWF="COOP_DWF";
    
    set_fitParams(model_DWF);
    read_all_equ_DWF(sysVar,n_sys,TA_d); // NOTE: fit T<=TA
    fitdata_processing_lsfit(dwf_invT_sortdecreasing); // NOTE: use "dwf_invT"
    alglib_lsfit(model_DWF);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/GLMtest_fit_DWF_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fit_DWF(o.c_str());
    
    write_fitModel(write_fit_DWF,model_DWF);
    write_fitarrays(write_fit_DWF);
    write_stopcond(write_fit_DWF);
    write_fitinfo(write_fit_DWF);
    write_errcurve_log(write_fit_DWF,model_DWF,dwf_invT_sortdecreasing);
    
    vector<vector<double>> vec2d;
    double temp=0,dwfr=0;
    for (int i=0; i<(int)dwf_invT_sortdecreasing.size(); ++i) {
        temp=dwf_sortdecreasing[i][0]; // NOTE: is "dwf" .NOT. dwf_invT
        dwfr=log10(calc_time_given_temp(c,temp,model_DWF));
        vec2d.push_back({temp,dwfr});
    }
    dwf_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_taueq(const StructureClass& sysVar,
                           const int n_sys)
{
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys); // all
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    
    vector<vector<double>> vec2d;
    double temp=0,taue=0;
    for (int i=0; i<(int)taueq_sortdecreasing.size(); ++i) {
        temp=taueq_sortdecreasing[i][0];
        taue=log10(calc_time_given_temp(c,temp,model));
        vec2d.push_back({temp,taue});
    }
    taueq_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_taueq(const StructureClass& sysVar,
                           const int n_sys,
                           const double TA_d)
{
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys,TA_d); // NOTE: fit T<=TA
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    
    vector<vector<double>> vec2d;
    double temp=0,taue=0;
    for (int i=0; i<(int)taueq_sortdecreasing.size(); ++i) {
        temp=taueq_sortdecreasing[i][0];
        taue=log10(calc_time_given_temp(c,temp,model));
        vec2d.push_back({temp,taue});
    }
    taueq_sortdecreasing_fit=vec2d;
}





void TheoryTest::write_peak_stringlen(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    vector<vector<double>> stringlen;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data/");
    i.append("/strings_time_stringlen_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readstringlendata(i.c_str());
    
    if (readstringlendata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double time=0,stringlenData=0;
        while (getline(readstringlendata,lineContent)) {
            istringstream iss(lineContent); // (frame_index,stringlen,time)
            iss >> frame_index;
            iss >> stringlenData;
            iss >> time;
            stringlen.push_back({(double)frame_index,stringlenData,time});
        }
        readstringlendata.close();
    }
    else {
        cout
        << "write_peak_stringlen: strings_time_stringlen file cannot open."
        << "\n";
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing1);
        // NOTE: sort 2nd column
    }
    catch (int i) {
        cout
        << "in write_peak_stringlen:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing1);
    }
    
    
    peak_frame=stringlen[0][0];
    peak_time =stringlen[0][2];
    
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // stringlen(frame_index,stringlen,time)
    
    // T stringlen frame_index time
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringlen
    << stringlen[0][1]
    << " "
    << stringlen[0][0]
    << " "
    << stringlen[0][2]
    << "\n";
    
    write_stringleninv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << stringlen[0][1]
    << " "
    << stringlen[0][0]
    << " "
    << stringlen[0][2]
    << "\n";
    
    write_stringlen.close();
    write_stringleninv.close();
}





void TheoryTest::write_peak_fastfrac(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d)
{
    vector<vector<double>> fastfrac;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data/");
    i.append("/strings_time_fastfrac_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readfastfracdata(i.c_str());
    
    if (readfastfracdata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double time=0,fastfracData=0;
        while (getline(readfastfracdata,lineContent)) {
            istringstream iss(lineContent); // (frame_index,fastfrac,time)
            iss >> frame_index;
            iss >> fastfracData;
            iss >> time;
            fastfrac.push_back({(double)frame_index,fastfracData,time});
        }
        readfastfracdata.close();
    }
    else {
        cout
        << "write_peak_fastfrac: strings_time_fastfrac file cannot open."
        << "\n";
    }
    
    // NOTE:
    //--------------------------------------------------------------------------
    // No need to sort the data,
    // but need to know the frame where the peak string length is present.
    //--------------------------------------------------------------------------
    int peak_index=0;
    for (int i=0; i<(int)fastfrac.size(); ++i) {
        if (fastfrac[i][0]==(double)peak_frame) {
            peak_index=i;
            break;
        }
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/fastfrac_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fastfrac(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data/");
    o.append("/fastfrac_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_fastfracinv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // fastfrac(frame_index,fastfrac,time)
    
    // T fast_frac frame_index time
    streamsize ss=write_fastfrac.precision();
    
    write_fastfrac
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_fastfrac.precision(ss);
    write_fastfrac << resetiosflags(ios::fixed|ios::showpoint);
    
    write_fastfrac
    << fastfrac[peak_index][1]
    << " "
    << fastfrac[peak_index][0]
    << " "
    << fastfrac[peak_index][2]
    << "\n";
    
    write_fastfracinv
    << (1000.0/(corrFacforT/precision))/(T_actual)
    << " "
    << fastfrac[peak_index][1]
    << " "
    << fastfrac[peak_index][0]
    << " "
    << fastfrac[peak_index][2]
    << "\n";
    
    write_fastfrac.close();
    write_fastfracinv.close();
}





void TheoryTest::read_individual_fastfrac(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data/");
    i.append("/fastfrac_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readfastfracdata(i.c_str());
    
    if (readfastfracdata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double Tdata=0,fastfracData=0,time=0;
        while (getline(readfastfracdata,lineContent)) {
            istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
            iss >> Tdata;
            iss >> fastfracData;
            iss >> frame_index;
            iss >> time;
            fastfrac.push_back({Tdata,fastfracData,(double)frame_index,time});
        }
        readfastfracdata.close();
    }
    else {
        cout
        << "read_individual_fastfrac_data: fastfrac datafile cannot open."
        << "\n";
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_individual_fastfrac_data:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
}





void TheoryTest::read_individual_equ_fastfrac(const StructureClass& sysVar,
                                              const int n_trl,
                                              const int n_sys)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/fastfrac_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readfastfracdata(i.c_str());
    
    if (readfastfracdata.is_open()) {
        string lineContent;
        int    frame_inedx=0;
        double Tdata=0,fastfracData=0,time=0;
        while (getline(readfastfracdata,lineContent)) {
            istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
            iss >> Tdata;
            iss >> fastfracData;
            iss >> frame_inedx;
            iss >> time;
            if (1) {
                fastfrac.push_back({Tdata,fastfracData,(double)frame_inedx,time});
            }
        }
        readfastfracdata.close();
    }
    else {
        cout
        << "read_individual_equ_fastfrac: fastfrac datafile cannot open."
        << "\n";
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_individual_equ_fastfrac:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    
    /* NOTE: individual */
    //-----------------------------------
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    
    int n_taueq =(int)taueq_sortdecreasing.size();
    int fastSize=(int)fastfrac_sortdecreasing.size();
    double temp=0;
    double fast=0;
    
    fastfrac.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<fastSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==fastfrac_sortdecreasing[ii][0]) {
                if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                    temp=fastfrac_sortdecreasing[ii][0];
                    fast=fastfrac_sortdecreasing[ii][1];
                    fastfrac.push_back({temp,fast});
                }
                break;
            }
        }
    }
    fastfrac_sortdecreasing=fastfrac;
}





void TheoryTest::read_all_fastfrac(const StructureClass& sysVar,
                                   const int n_sys)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data/");
        i.append("/fastfrac_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readfastfracdata(i.c_str());
        
        if (readfastfracdata.is_open()) {
            string lineContent;
            int    frame_index=0;
            double Tdata=0,fastfracData=0,time=0;
            while (getline(readfastfracdata,lineContent)) {
                istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
                iss >> Tdata;
                iss >> fastfracData;
                iss >> frame_index;
                iss >> time;
                fastfrac.push_back({Tdata,fastfracData,(double)frame_index,time});
            }
            readfastfracdata.close();
        }
        else {
            cout
            << "read_all_fastfrac_data: fastfrac datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_fastfrac_data:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
}





void TheoryTest::read_all_fastfrac(const StructureClass& sysVar,
                                   const int n_sys,
                                   const double TA_d)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data/");
        i.append("/fastfrac_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readfastfracdata(i.c_str());
        
        if (readfastfracdata.is_open()) {
            string lineContent;
            int    frame_index=0;
            double Tdata=0,fastfracData=0,time=0;
            while (getline(readfastfracdata,lineContent)) {
                istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
                iss >> Tdata;
                iss >> fastfracData;
                iss >> frame_index;
                iss >> time;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    fastfrac.push_back({Tdata,fastfracData,(double)frame_index,time});}
            }
            readfastfracdata.close();
        }
        else {
            cout
            << "read_all_fastfrac_data: fastfrac datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_fastfrac_data:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
}





void TheoryTest::read_all_equ_fastfrac(const StructureClass& sysVar,
                                       const int n_sys)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/fastfrac_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readfastfracdata(i.c_str());
        
        if (readfastfracdata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,fastfracData=0,time=0;
            while (getline(readfastfracdata,lineContent)) {
                istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
                iss >> Tdata;
                iss >> fastfracData;
                iss >> frame_inedx;
                iss >> time;
                if (1) {
                    fastfrac.push_back({Tdata,fastfracData,(double)frame_inedx,time});
                }
            }
            readfastfracdata.close();
        }
        else {
            cout
            << "read_all_string_equ_length: fastfrac datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_string_equ_length:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys);
    
    int n_taueq =(int)taueq_sortdecreasing.size();
    int fastSize=(int)fastfrac_sortdecreasing.size();
    double temp=0;
    double fast=0;
    
    fastfrac.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<fastSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==fastfrac_sortdecreasing[ii][0]) {
                if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                    temp=fastfrac_sortdecreasing[ii][0];
                    fast=fastfrac_sortdecreasing[ii][1];
                    fastfrac.push_back({temp,fast});
                }
                break;
            }
        }
    }
    fastfrac_sortdecreasing=fastfrac;
}





void TheoryTest::read_all_equ_fastfrac(const StructureClass& sysVar,
                                       const int n_sys,
                                       const double TA_d)
{
    fastfrac_sortdecreasing.clear();
    vector<vector<double>> fastfrac;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/fastfrac_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readfastfracdata(i.c_str());
        
        if (readfastfracdata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,fastfracData=0,time=0;
            while (getline(readfastfracdata,lineContent)) {
                istringstream iss(lineContent); // (T,fastfrac,frame_index,time)
                iss >> Tdata;
                iss >> fastfracData;
                iss >> frame_inedx;
                iss >> time;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    fastfrac.push_back({Tdata,fastfracData,(double)frame_inedx,time});
                }
            }
            readfastfracdata.close();
        }
        else {
            cout
            << "read_all_string_equ_length: fastfrac datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(fastfrac.size()==0) throw 0;
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_string_equ_length:" << "\n"
        << "fastfrac size = " << i << "\n";
        fastfrac.push_back({0,0});
        std::sort(fastfrac.begin(),fastfrac.end(),sortDecreasing);
        averaging_sorted_data(fastfrac,fastfrac_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys,TA_d);
    
    int n_taueq =(int)taueq_sortdecreasing.size();
    int fastSize=(int)fastfrac_sortdecreasing.size();
    double temp=0;
    double fast=0;
    
    fastfrac.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<fastSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==fastfrac_sortdecreasing[ii][0]) {
                if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                    temp=fastfrac_sortdecreasing[ii][0];
                    fast=fastfrac_sortdecreasing[ii][1];
                    fastfrac.push_back({temp,fast});
                }
                break;
            }
        }
    }
    fastfrac_sortdecreasing=fastfrac;
}





void TheoryTest::write_fastfrac_equ(const StructureClass& sysVar)
{
    string o;
    
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/fastfrac_equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            //cout << o << "\n";
            ofstream write_fastfrac(o.c_str(),ofstream::app); // 'append'
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/fastfrac_equ_inverseT_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_fastfracinv(o.c_str(),ofstream::app); // 'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            read_individual_fastfrac(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_taueq=(int)taueq_sortdecreasing.size();
            int fastfracSize=(int)fastfrac_sortdecreasing.size();
            
            for (int i=0; i<n_taueq; ++i) {
                for (int ii=0; ii<fastfracSize; ++ii) {
                    if (taueq_sortdecreasing[i][0]==fastfrac_sortdecreasing[ii][0]) {
                        if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                            
                            write_fastfrac
                            << fastfrac_sortdecreasing[ii][0] << " "
                            << fastfrac_sortdecreasing[ii][1] << "\n";
                            write_fastfracinv
                            << (1000.0/(corrFacforT/precision))/fastfrac_sortdecreasing[ii][0] << " "
                            << fastfrac_sortdecreasing[ii][1] << "\n";
                        }
                        break;
                    }
                }
            }
            write_fastfrac.close();
            write_fastfracinv.close();
        }
    }
}





void TheoryTest::write_fastfrac_equ_avg(const StructureClass& sysVar)
{
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fastfrac_equ_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        //cout << o << "\n";
        ofstream write_fastfrac(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/fastfrac_equ_avg_inverseT_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_fastfracinv(o.c_str(),ofstream::app); // 'append'
        
        /* NOTE: all */
        //---------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        read_all_fastfrac(sysVar,n_sys);
        //---------------------------------------------
        
        int n_taueq=(int)taueq_sortdecreasing.size();
        int fastfracSize=(int)fastfrac_sortdecreasing.size();
        
        for (int i=0; i<n_taueq; ++i) {
            for (int ii=0; ii<fastfracSize; ++ii) {
                if (taueq_sortdecreasing[i][0]==fastfrac_sortdecreasing[ii][0]) {
                    if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                        
                        write_fastfrac
                        << fastfrac_sortdecreasing[ii][0] << " "
                        << fastfrac_sortdecreasing[ii][1] << "\n";
                        write_fastfracinv
                        << (1000.0/(corrFacforT/precision))/fastfrac_sortdecreasing[ii][0] << " "
                        << fastfrac_sortdecreasing[ii][1] << "\n";
                    }
                    break;
                }
            }
        }
        write_fastfrac.close();
        write_fastfracinv.close();
    }
}





void TheoryTest::read_individual_stringlen(const StructureClass& sysVar,
                                           const int n_trl,
                                           const int n_sys)
{
    stringlen_sortdecreasing.clear();
    vector<vector<double>> stringlen;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/stringlen_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readstringlendata(i.c_str());
    
    if (readstringlendata.is_open()) {
        string lineContent;
        int    frame_inedx=0;
        double Tdata=0,stringlenData=0,time=0;
        while (getline(readstringlendata,lineContent)) {
            istringstream iss(lineContent); // (T,stringlen,frame_index,time)
            iss >> Tdata;
            iss >> stringlenData;
            iss >> frame_inedx;
            iss >> time;
            if (1) { //stringlenData!=0
                stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
            }
        }
        readstringlendata.close();
    }
    else {
        cout
        << "read_individual_stringlen_data: stringlen datafile cannot open."
        << "\n";
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_individual_stringlen_data:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
}





void TheoryTest::read_individual_equ_stringlen(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys)
{
    stringlen_sortdecreasing.clear();
    vector<vector<double>> stringlen;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/stringlen_all_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append(".dat");
    //cout << i << "\n";
    ifstream readstringlendata(i.c_str());
    
    if (readstringlendata.is_open()) {
        string lineContent;
        int    frame_inedx=0;
        double Tdata=0,stringlenData=0,time=0;
        while (getline(readstringlendata,lineContent)) {
            istringstream iss(lineContent); // (T,stringlen,frame_index,time)
            iss >> Tdata;
            iss >> stringlenData;
            iss >> frame_inedx;
            iss >> time;
            if (1) { //stringlenData!=0
                stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
            }
        }
        readstringlendata.close();
    }
    else {
        cout
        << "read_individual_equ_stringlen: stringlen datafile cannot open."
        << "\n";
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_individual_equ_stringlen:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    
    /* NOTE: individual */
    //-----------------------------------
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    //read_individual_fastfrac(sysVar,n_trl,n_sys); // NOTE: read "all"
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    double temp=0;
    double strl=0;
    
    stringlen.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                //if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                //}
                break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
}





void TheoryTest::read_all_stringlen(const StructureClass& sysVar,
                                    const int n_sys)
{
    stringlen_sortdecreasing.clear();
    vector<vector<double>> stringlen;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/stringlen_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readstringlendata(i.c_str());
        
        if (readstringlendata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringlenData=0,time=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent); // (T,stringlen,frame_index,time)
                iss >> Tdata;
                iss >> stringlenData;
                iss >> frame_inedx;
                iss >> time;
                if (1) { //stringlenData!=0
                    stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
                }
            }
            readstringlendata.close();
        }
        else {
            cout
            << "read_all_stringlen_data: stringlen datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_stringlen_data:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
}





void TheoryTest::read_all_equ_stringlen(const StructureClass& sysVar,
                                        const int n_sys)
{
    stringlen_sortdecreasing.clear();
    vector<vector<double>> stringlen;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/stringlen_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readstringlendata(i.c_str());
        
        if (readstringlendata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringlenData=0,time=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent); // (T,stringlen,frame_index,time)
                iss >> Tdata;
                iss >> stringlenData;
                iss >> frame_inedx;
                iss >> time;
                if (1) { //stringlenData!=0
                    stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
                }
            }
            readstringlendata.close();
        }
        else {
            cout
            << "read_all_string_equ_length: stringlen datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_string_equ_length:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys);
    //read_all_fastfrac(sysVar,n_sys); // NOTE: read "all"
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    double temp=0;
    double strl=0;
    
    stringlen.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                //if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                //}
                break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
}





void TheoryTest::read_all_equ_stringlen(const StructureClass& sysVar,
                                        const int n_sys,
                                        const double TA_d)
{
    stringlen_sortdecreasing.clear();
    vector<vector<double>> stringlen;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/stringlen_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readstringlendata(i.c_str());
        
        if (readstringlendata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringlenData=0,time=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent); // (T,stringlen,frame_index,time)
                iss >> Tdata;
                iss >> stringlenData;
                iss >> frame_inedx;
                iss >> time;
                if (Tdata<=TA_d) { // only collect when T<=TA
                    stringlen.push_back
                    ({Tdata,stringlenData,(double)frame_inedx,time});
                }
                
            }
            readstringlendata.close();
        }
        else {
            cout
            << "read_all_string_equ_length: stringlen datafile cannot open."
            << "\n";
        }
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_all_string_equ_length:" << "\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
    }
    
    /* NOTE: all */
    //-----------------------------------
    read_all_taueq_data(sysVar,n_sys,TA_d);
    //read_all_fastfrac(sysVar,n_sys,TA_d);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    double temp=0;
    double strl=0;
    
    stringlen.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                //if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                //}
                break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
}





void TheoryTest::write_stringlen_equ(const StructureClass& sysVar)
{
    string o;
    
    /* write out "equilibrium" stringlen to file (similar to taueq file) */
    //----------------------------------------------------------------------
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/stringlen_equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            //cout << o << "\n";
            ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/stringlen_equ_inverseT_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            read_individual_stringlen(sysVar,n_trl,n_sys);
            //read_individual_fastfrac(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_taueq=(int)taueq_sortdecreasing.size();
            int strSize=(int)stringlen_sortdecreasing.size();
            
            for (int i=0; i<n_taueq; ++i) {
                for (int ii=0; ii<strSize; ++ii) {
                    if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                        //if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                        
                        write_stringlen
                        << stringlen_sortdecreasing[ii][0] << " "
                        << stringlen_sortdecreasing[ii][1] << "\n";
                        write_stringleninv
                        << (1000.0/(corrFacforT/precision))/stringlen_sortdecreasing[ii][0] << " "
                        << stringlen_sortdecreasing[ii][1] << "\n";
                        //}
                        break;
                    }
                }
            }
            write_stringlen.close();
            write_stringleninv.close();
        }
    }
}





void TheoryTest::write_stringlen_equ_avg(const StructureClass& sysVar)
{
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/stringlen_equ_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        //cout << o << "\n";
        ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/stringlen_equ_avg_inverseT_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
        
        /* NOTE: all */
        //-----------------------------------
        read_all_taueq_data(sysVar,n_sys);
        read_all_stringlen(sysVar,n_sys);
        //read_all_fastfrac(sysVar,n_sys);
        //-----------------------------------
        
        int n_taueq=(int)taueq_sortdecreasing.size();
        int strSize=(int)stringlen_sortdecreasing.size();
        
        for (int i=0; i<n_taueq; ++i) {
            for (int ii=0; ii<strSize; ++ii) {
                if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                    //if (fastfrac_sortdecreasing[ii][1]<fast_threshold) { // NOTE
                    
                    write_stringlen
                    << stringlen_sortdecreasing[ii][0] << " "
                    << stringlen_sortdecreasing[ii][1] << "\n";
                    write_stringleninv
                    << (1000.0/(corrFacforT/precision))/stringlen_sortdecreasing[ii][0] << " "
                    << stringlen_sortdecreasing[ii][1] << "\n";
                    //}
                    break;
                }
            }
        }
        write_stringlen.close();
        write_stringleninv.close();
    }
}





void TheoryTest::find_LA(const StructureClass& sysVar,
                         const int n_sys)
{
    // find TA
    read_all_taueq_data(sysVar,n_sys);
    find_TA(sysVar,n_sys,taueq_sortdecreasing);
    
    TA   = TA_avg;
    tauA = tauA_avg;
    
    /* interpolate L(T) to get L(TA) */
    double T_pre=0,T_pos=0;
    double stringlen_pre=0,stringlen_pos=0;
    double slope=0;
    
    read_all_equ_stringlen(sysVar,n_sys); // NOTE: use "equ" data
    
    for (int i=0; i<(int)(stringlen_sortdecreasing.size()-1); ++i) {
        
        T_pre = stringlen_sortdecreasing[i][0];
        T_pos = stringlen_sortdecreasing[i+1][0];
        
        if (((TA-T_pre)*(TA-T_pos))<0) {
            
            stringlen_pre=stringlen_sortdecreasing[i][1];
            stringlen_pos=stringlen_sortdecreasing[i+1][1];
            
            slope = (stringlen_pos-stringlen_pre)/(T_pos-T_pre);
            LA    = stringlen_pre+slope*(TA-T_pre);
            
            return;
        }
    }
}





void TheoryTest::find_uA2(const StructureClass& sysVar,
                          const int n_sys)
{
    // find TA
    read_all_taueq_data(sysVar,n_sys);
    find_TA(sysVar,n_sys,taueq_sortdecreasing);
    
    TA   = TA_avg;
    tauA = tauA_avg;
    
    /* interpolate <u2(T)> to get <u2(TA)> */
    double T_pre=0,T_pos=0;
    double DWF_pre=0,DWF_pos=0;
    double slope=0;
    
    read_all_equ_DWF(sysVar,n_sys); // NOTE: use "equ" data
    
    for (int i=0; i<(int)(dwf_sortdecreasing.size()-1); ++i) {
        
        T_pre = dwf_sortdecreasing[i][0];
        T_pos = dwf_sortdecreasing[i+1][0];
        
        if (((TA-T_pre)*(TA-T_pos))<0) {
            
            DWF_pre=dwf_sortdecreasing[i][1];
            DWF_pos=dwf_sortdecreasing[i+1][1];
            
            slope = (DWF_pos-DWF_pre)/(T_pos-T_pre);
            uA2   = DWF_pre+slope*(TA-T_pre);
            
            return;
        }
    }
}





void TheoryTest::find_individual_Arrhenius(const StructureClass& sysVar,
                                           const int n_trl,
                                           const int n_sys)
{
    set_fitParams("Arrhenius");
    read_individual_taueq_data(sysVar,n_trl,n_sys,"Arrhenius"); // individual
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit("Arrhenius");
    tau0_Arrhenius=c[0];
    Ea=c[1];
}





void TheoryTest::find_all_Arrhenius(const StructureClass& sysVar,
                                    const int n_sys)
{
    set_fitParams("Arrhenius");
    read_all_taueq_data(sysVar,n_sys,"Arrhenius"); // all
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit("Arrhenius");
    tau0_Arrhenius=c[0];
    Ea=c[1];
}





void TheoryTest::find_individual_tau0(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys)
{
    set_fitParams(model);
    read_individual_taueq_data(sysVar,n_trl,n_sys); // individual
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    tau0_model=c[0];
    
    if (model=="COOP") c_COOP = c;
}





void TheoryTest::find_all_tau0(const StructureClass& sysVar,
                               const int n_sys)
{
    set_fitParams(model);
    read_all_taueq_data(sysVar,n_sys); // all
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model);
    tau0_model=c[0];
    
    if (model=="COOP") c_COOP = c;
}




void TheoryTest::find_AG_tau0_fit(const StructureClass& sysVar,
                                  const int n_sys)
{
    string model="linear";
    
    set_fitParams(model);
    read_AG_raw(sysVar,n_sys);
    fitdata_processing_lsfit(dataxy_sortdecreasing); // NOTE: dataxy_sortdecreasing
    alglib_lsfit(model);
    slope_fit=c[0];
    tau0_fit=c[1]; // NOTE: natural log form
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/AGtest_tau0_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_actual(write_tau_inf,model,dataxy_sortdecreasing);
}





void TheoryTest::read_AG_raw(const StructureClass& sysVar,
                             const int n_sys)
{
    dataxy_sortdecreasing.clear();
    vector<vector<double>> dataxy;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data/AGtest_raw_");
        i.append(sysVar.get_usic());
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readdataxydata(i.c_str());
        
        if (readdataxydata.is_open()) {
            string lineContent;
            double datax=0,datay=0;
            while (getline(readdataxydata,lineContent)) {
                istringstream iss(lineContent); // (x,y)
                iss >> datax;
                iss >> datay;
                dataxy.push_back({datax,datay});
            }
            readdataxydata.close();
        }
        else {
            cout
            << "read_AG_raw: AGtest_raw file cannot open." << "\n";
        }
    }
    
    try {
        if(dataxy.size()==0) throw 0;
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_AG_raw:" << "\n"
        << "dataxy size = " << i << "\n";
        dataxy.push_back({0,0});
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    }
}





void TheoryTest::find_GLM_tau0_fit(const StructureClass& sysVar,
                                   const int n_sys)
{
    string model="linear";
    
    set_fitParams(model);
    read_GLM_raw(sysVar,n_sys);
    fitdata_processing_lsfit(dataxy_sortdecreasing); // NOTE: dataxy_sortdecreasing
    alglib_lsfit(model);
    slope_fit=c[0];
    tau0_fit=c[1]; // NOTE: natural log form
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/GLMtest_tau0_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_actual(write_tau_inf,model,dataxy_sortdecreasing);
}





void TheoryTest::read_GLM_raw(const StructureClass& sysVar,
                              const int n_sys)
{
    dataxy_sortdecreasing.clear();
    vector<vector<double>> dataxy;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data/GLMtest_raw_");
        i.append(sysVar.get_usic());
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        //cout << i << "\n";
        ifstream readdataxydata(i.c_str());
        
        if (readdataxydata.is_open()) {
            string lineContent;
            double datax=0,datay=0;
            while (getline(readdataxydata,lineContent)) {
                istringstream iss(lineContent); // (x,y)
                iss >> datax;
                iss >> datay;
                dataxy.push_back({datax,datay});
            }
            readdataxydata.close();
        }
        else {
            cout
            << "read_GLM_raw: GLMtest_raw file cannot open." << "\n";
        }
    }
    
    try {
        if(dataxy.size()==0) throw 0;
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in read_GLM_raw:" << "\n"
        << "dataxy size = " << i << "\n";
        dataxy.push_back({0,0});
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    }
}





/*==( public setters )==*/
void TheoryTest::set_is_use_ngp_peak_frame(const bool b){is_use_ngp_peak_frame=b;}
void TheoryTest::set_is_use_counting1(const bool b){is_use_counting1=b;}
void TheoryTest::set_is_data_smoothing(const bool b){is_data_smoothing=b;}

