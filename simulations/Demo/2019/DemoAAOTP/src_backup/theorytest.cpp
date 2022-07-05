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

#include "theorytest.h"
#include "functions.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

TheoryTest::TheoryTest(const StructureClass& sysVar,
                       const WorkScripts& ws,
                       const AmdatAnalysis& aa,
                       const FitData& fd):
/* Base class */
FitData(fd),

/* bool */
is_AGtest(false),
is_RFOTtest(false),
is_GLMtest(false),
is_fitu2T(false),
is_Leporinitest(false),
is_HWtest(false),
is_reltaur(false),
is_check_amdat_version(aa.get_is_check_amdat_version()),
is_use_stringLen_smoothed(false),
is_use_fast_frac_smoothed(false),
is_use_counting1(true),
is_backwardExtrp(false),
is_data_smoothing(false),
is_find_fastfrac(false),
is_use_peak_stringlen(true),
is_stringmb(aa.get_is_stringmb()),
is_use_displacement_list(aa.get_is_use_displacement_list()),
is_fit_univILP(false),
is_use_Lbeta(false),

/* int */
peak_frame(0),
n_trajs(0),
n_total_frames(0),
n_moments(aa.get_n_moments()),

/* double */
time_stringlen(0),
DWF(0),
mean_strings(0),
mean_length(0),
mean_length_count1(0),
order_parameter(0),
stringlen(0),
stringsRg(0),
stringsRg_cutlo(2),
stringsRg_cuthi(10),
peak_time(0),
fast_threshold(0.15),
fit_tau_threshold(10.0),//log10-based
stringlen_threshold(0.0),
time_stringlen_threshold(3),//log10-based
n_distr_cutoff(0.01),
p1_extrp(0),
delG(0),
delmiu(0),
delmiuTA(0),
Ea_COOP(0),
delHa(0),
delSa(0),
tau0_vib(0),
tau0_shift(0),
strings_threshold(aa.get_strings_threshold()),
fixedfastpercentage(100.0-aa.get_threshold_percentile()),

/* string */
analysispart(aa.get_analysispart()),
model_stringlen_distr("exp_string"),
stringtype("stringlen"),
cusstrfolder(aa.get_cusstrfolder()),

/** STL containers **/
//------------------------------------------------------------------------------
dataxy_sortdecreasing(),   // (x,y)
stringlen_sortdecreasing(),// (T,stringlen)
fast_frac_sortdecreasing(),// (T,fast_fraction)
taueq_sortdecreasing_fit(),
stringlen_sortdecreasing_fit(),
dwf_sortdecreasing_fit(),
sigmaMatrix(aa.get_sigmaMatrix()),

/** Penalizaed spline interpolation in ALGLIB **/
//------------------------------------------------------------------------------
rho_Lt(0.0),

/** GLM **/
//------------------------------------------------------------------------------
tau0_GLM(0),
u02(0),
alpha(0),
lntau0(0),
dfit(0),

/** Leporini **/
//------------------------------------------------------------------------------
logtauref(2.0),//fixed ref timescale in ps
Tref(0),
u2ref(0),
Xref(0),
u2g(0),
Yg(0),
mg(0),

c_COOP()
{
    if (sysVar.get_systemUnit()=="real") {
        tau0_vib=pow(10,2);//fs
        tau0_shift=tau0_vib;
        time_stringlen_threshold=6;
        logtauref+=3.0;//ps
        fit_tau_threshold+=3;//ps
    } else if (sysVar.get_systemUnit()=="lj") {
        tau0_vib=pow(10,-1);//tau
        tau0_shift=tau0_vib;
        time_stringlen_threshold=2;
    } else if (sysVar.get_systemUnit()=="metal") {
        tau0_vib=pow(10,-1);//ps
        tau0_shift=tau0_vib;
        time_stringlen_threshold=3;
    }
}





void TheoryTest::amdatstringanalysis(StructureClass& sysVar,
                                     const WorkScripts& ws,
                                     AmdatAnalysis& aa,
                                     const vector<vector<double>>& tinfo)
{
    aa.set_is_strFac(false);
    aa.set_is_msd(false);
    aa.set_is_ngp(false);
    aa.set_is_isfs(false);
    aa.set_is_strings(true);//NOTE
    aa.set_is_use_voroNeighbors(false);//NOTE
    
    aa.test_AnalysisFolder(sysVar);
    //call_system_bash(o);
    
    vector<string> amdat_targets;
    
    
    /* find L(T) at the peak frame of NGP */
    ////////////////////////////////////////////////////////////////////////////
    if (is_use_ngp_peak_frame)
    {
        /* find L(T) at the peak frame of NGP */
        //----------------------------------------------------------------------
        aa.set_is_peak_frame(true);
        aa.make_amdatInputFile(sysVar,ws);
        
        vector<vector<int>> peakNGP_frames,frame_times;
        vector<int> peakNGP_tmp,frame_tmp;
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    
                    /** find the peak frame of NGP **/
                    find_ngp_peak_frame(sysVar,n_trl,n_sys,tinfo.at(index).at(T));
                    peakNGP_tmp.push_back(ngp_peak_frame);
                    frame_tmp.push_back(frame_time);
                    
                    /* AMDAT submission files */
                    aa.make_amdatSubFile
                    (sysVar,ws,n_trl,n_sys,tinfo.at(index).at(T),ngp_peak_frame);
                    
                    /* Specify target files to be watched for */
                    amdat_targets.push_back
                    (aa.amdat_target(sysVar,n_trl,n_sys,tinfo.at(index).at(T),"strings"));
                }
                peakNGP_frames.push_back(peakNGP_tmp);
                frame_times.push_back(frame_tmp);
            }
        }
        /* AMDAT Jobs Submissions */
        aa.make_amdatSubScript(sysVar,tinfo);
        
        /* Watch and Hold */
        if (sysVar.get_is_watch_hold()) {
            autoWork::watch_hold(sysVar,amdat_targets,"analysis",current_regime);
        } amdat_targets.clear();
        
        /* write out stringlen to file (all Ts) */
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                int n_Temp = (int)(tinfo[index].size());
                for (int T=0; T<n_Temp; ++T) {
                    int frame = peakNGP_frames.at(index).at(T);
                    double frametime = frame_times.at(index).at(T);
                    find_mean_stringlen
                    (sysVar,n_trl,n_sys,tinfo.at(index).at(T),frametime,frame);
                    if (is_find_fastfrac) {
                        find_mean_fast_frac
                        (sysVar,n_trl,n_sys,tinfo.at(index).at(T),frametime,frame);
                    }
                }
            }
        }
    } /* find L(T) at the peak frame of NGP */
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    /* find L(T) at in-block frames */
    ////////////////////////////////////////////////////////////////////////////
    else
    {
        /** find block frame **/
        //----------------------------------------------------------------------
        /** find block frame index (same for same regime temperatures) **/
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
        
        bool is_divide=true;//whether divide string analysis into multiple runs
        
        if (sysVar.get_is_directSub()) is_divide=true;
        else is_divide=false;
        
        /** make and run string analysis **/
        //----------------------------------------------------------------------
        vector<vector<int>> calc_frames;//NOTE: 2D int vector
        
        if (!is_divide) //keep string analysis in one run
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    int n_Temp = (int)(tinfo[index].size());
                    for (int T=0; T<n_Temp; ++T) {
                        
                        /** find the peak frame of NGP **/
                        //find_ngp_peak_frame(sysVar,n_trl,n_sys,tinfo.at(index).at(T));
                        
                        vector<int> tmp=aa.make_amdatInputFile_strings_loop
                        (sysVar,ws,n_trl,n_sys,tinfo.at(index).at(T),
                         ngp_peak_frame,ngp_block_frame);
                        
                        beg_frame=tmp[0];
                        end_frame=tmp[1];
                        calc_frames.push_back({beg_frame,end_frame});
                        
                        /* AMDAT submission files */
                        aa.make_amdatSubFile
                        (sysVar,ws,n_trl,n_sys,tinfo.at(index).at(T),ngp_block_frame);
                        
                        /* Specify target files to be watched for */
                        amdat_targets.push_back
                        (aa.amdat_target
                         (sysVar,n_trl,n_sys,tinfo.at(index).at(T),
                          cusstrfolder,"Rg_strings",end_frame));//NOTE: end_frame
                    }
                }
            }
            /* AMDAT Jobs Submissions */
            aa.make_amdatSubScript(sysVar,tinfo);
            
            /* Watch and Hold */
            if (sysVar.get_is_watch_hold()) {
                autoWork::watch_hold(sysVar,amdat_targets,"analysis",current_regime);
            } amdat_targets.clear();
            
            
        } else { //divide string analysis into multiple runs to avoid memory runout
            
            int framelimit=10;//KG8000,BinLJ8000
            
            int times_run=1;
            int tmpblockframe=ngp_block_frame;
            if (tmpblockframe>framelimit) {
                do { ++times_run;
                } while ((double)tmpblockframe/(double)times_run>framelimit);
            } int res=tmpblockframe%framelimit;
            
            if (tmpblockframe>0&&res==0) {
                res=framelimit;
            } int beg=0,end=0;
            
            for (int run=0;run<times_run;++run)
            {
                if (run!=times_run-1) {
                    beg=framelimit*run;
                    end=beg+framelimit-1;//because in amdat input: run [beg,end+1)
                } else {
                    if (times_run==1) {
                        beg=1;
                        end=res;
                    } else {
                        beg=framelimit*run;
                        end=beg+res;//NOTE: actual frame used
                    }
                }
                /* wait a moment for memory to clearup */
                if (run>0) {
                    system_wait(10);
                }
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        int n_Temp = (int)(tinfo[index].size());
                        for (int T=0; T<n_Temp; ++T) {
                            
                            vector<int> tmp=aa.make_amdatInputFile_strings_loop
                            (sysVar,ws,n_trl,n_sys,tinfo.at(index).at(T),
                             ngp_peak_frame,ngp_block_frame,beg,end);
                            
                            /* AMDAT submission files */
                            aa.make_amdatSubFile
                            (sysVar,ws,n_trl,n_sys,tinfo.at(index).at(T),
                             ngp_block_frame);
                            
                            /* Specify target files to be watched for */
                            amdat_targets.push_back
                            (aa.amdat_target
                             (sysVar,n_trl,n_sys,tinfo.at(index).at(T),
                              cusstrfolder,"Rg_strings",end));//NOTE: end
                        }
                    }
                }
                /* AMDAT Jobs Submissions */
                aa.make_amdatSubScript(sysVar,tinfo);
                
                /* Watch and Hold */
                if (sysVar.get_is_watch_hold()) {
                    autoWork::watch_hold(sysVar,amdat_targets,"analysis",current_regime);
                } amdat_targets.clear();
            }
        }
        
        // Collect strings data generated by AMDAT
        //----------------------------------------------------------------------
        ////////////////////////////////////////////////////////////////////////
        vector<vector<int>> which_frame;//NOTE:array of frame indices
        if (true)
        {
            /* write out stringlen to file (all Ts) */
            int index1=0;
            for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                    
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    int n_Temp = (int)(tinfo[index].size());
                    
                    vector<int> tmp;
                    
                    for (int T=0; T<n_Temp; ++T)
                    {
                        // Note the actual frames used in amdat run
                        //------------------------------------------------------
                        if (aa.get_species()=="KGPolymer"||aa.get_species()=="binaryLJ") {
                            beg_frame=1;
                            end_frame=ngp_block_frame-1;//NOTE: actual frame used
                        } else {
                            beg_frame=1;
                            end_frame=ngp_block_frame;//NOTE: actual frame used
                        }
                        //------------------------------------------------------
                        
                        ++index1;
                        for (int frame_index=beg_frame;frame_index<=end_frame;++frame_index)
                        {
                            find_time_stringlen
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),frame_index);
                            
                            find_time_stringsRg
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),frame_index);
                            
                            if (is_find_fastfrac) {
                                find_time_fast_frac
                                (sysVar,n_trl,n_sys,tinfo.at(index).at(T),frame_index);
                            }
                        }
                        int chosen_frame=0;
                        if (is_use_peak_stringlen) {
                            /** use frame of peak L **/
                            write_peak_stringlen(sysVar,n_trl,n_sys,tinfo.at(index).at(T));
                            tmp.push_back(peak_frame);
                            chosen_frame=peak_frame;
                            write_frame_stringsRg
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),peak_frame);
                        } else {
                            /** use frame at taualpha **/
                            write_time_stringlen(sysVar,n_trl,n_sys,tinfo.at(index).at(T));
                            tmp.push_back(taualpha_frame);
                            chosen_frame=taualpha_frame;
                            write_frame_stringsRg
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),taualpha_frame);
                        }
                        if (is_find_fastfrac) {
                            write_peak_fast_frac
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        fit_all_strings_rgvsm
                        (sysVar,n_trl,n_sys,tinfo.at(index).at(T),chosen_frame);
                        
                    } which_frame.push_back(tmp);
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////
        //----------------------------------------------------------------------
        
        if (false)
        {
            // for mobile particles
            //------------------------------------------------------------------
            if (!sysVar.get_is_fitData())
            {
                FitData fd(sysVar,ws,aa);
                
                fd.set_is_use_FG(true);
                fd.set_is_fit_Fs_by_spline(true);//NOTE
                fd.set_is_fit_full_alpha(true);
                fd.set_is_use_gammafunc(false);
                
                fd.set_relaxation_target("isfs");
                fd.set_is_fit_sExp(true);
                
                fd.set_is_fit_Arrhenius(true);
                fd.set_is_fit_tauFit(true);
                
                fd.set_is_applytauFitcut(false);//NOTE
                fd.set_is_use_KWWassist(true);
                fd.set_is_find_DWF(true);
                fd.set_is_calc_thermoData(false);
                
                fd.set_definitionofTA("ArrNor");
                
                /** get equ tau_alpha **/
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                    {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        sysVar.set_n_Temp((int)(tinfo.at(index).size()));
                        
                        /** clear content of contianers for new values **/
                        //******************************************************
                        sysVar.get_equilibratedTs().at(index).clear();
                        sysVar.get_tauFit().at(index).clear();
                        sysVar.get_tauEqu().at(index).clear();
                        //******************************************************
                        for (int T=0; T<sysVar.get_n_Temp(); ++T) {
                            /** fit relaxation profile to KWW function **/
                            fd.fit_sExp
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),"KWW",
                             which_frame.at(index).at(T));
                        }
                        /** write temperatures to file **/
                        fd.write_qchTs(sysVar,n_trl,n_sys);
                        fd.write_equTs(sysVar,n_trl,n_sys);
                    }
                }
                /** get eqi DWF **/
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                    {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        sysVar.set_n_Temp((int)(tinfo.at(index).size()));
                        for (int T=0; T<sysVar.get_n_Temp(); ++T) {
                            fd.write_DWF
                            (sysVar,n_trl,n_sys,tinfo.at(index).at(T),
                             which_frame.at(index).at(T));
                        }
                    }
                }
                if (current_regime==sysVar.get_regime_end())
                {
                    fd.write_DWF_equ(sysVar);
                    fd.write_DWF_equ_avg(sysVar);
                }
            }
        }
        
    } /* find L(T) at in-block frames */
    ////////////////////////////////////////////////////////////////////////////
}





void TheoryTest::AGtest(const StructureClass& sysVar)
{
    is_use_Lbeta=false;
    
    stringtype="stringlen";//use "strings" (mass scale)
    
    write_stringlen_equ(sysVar);
    write_stringlen_equ_avg(sysVar);
    
    vector<vector<double>> diagonal,loglog_diagonal;
    
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_diagonal_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_diag_raw(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_AG_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_AGtest(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_AG_log_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_AGtest_log(o.c_str(),ofstream::app);//'append'
        
        bool is_write_allsys=true;
        o.clear();
        o.append(sysVar.get_Path());
        o.append("/simulations");
        o.append("/ATheoryTest_AG_log_"+sysVar.get_usic());
        o.append(".dat");
        ofstream write_AGtest_log_allsys(o.c_str(),ofstream::app);//'append'
        
        bool is_cutbyTA=true;
        
        bool is_fit_stringmodel=true;
        
        string stringmodel="transitState";
        /** available models:
         ** "string2param" -- [delH,delS]
         ** "string3param" -- [tau0,delH,delS]
         ** "transitState" -- [delH,delS] **/
        
        /** find TA and L(TA) **/
        //----------------------------------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);
        find_LA(sysVar,n_sys,TA);
        
        if (is_fit_stringmodel) {
            fit_stringmodel(sysVar,n_sys,TA,stringmodel);
        } else {
            delmiuTA = TA*log(tauA/tau0_vib);
            delHa    = Ea;
            delSa    = (delHa-delmiuTA)*pow(TA,-1);
        }
        
        /** smooth data or use raw **/
        //----------------------------------------------------------------------
        int n_taueq=0;
        int strSize=0;
        if (is_data_smoothing) {
            if (is_cutbyTA) {
                fit_taueq(sysVar,n_sys,TA,"COOP");
                fit_stringlen(sysVar,n_sys,TA,"COOP_string");
            } else {
                fit_taueq(sysVar,n_sys,"COOP");
                fit_stringlen(sysVar,n_sys,"COOP_string");
            }
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            strSize=(int)stringlen_sortdecreasing_fit.size();
        } else {
            if (is_cutbyTA) {
                read_all_taueq_data(sysVar,n_sys,TA);
                read_all_equ_stringlen(sysVar,n_sys,TA);
            } else {
                read_all_taueq_data(sysVar,n_sys);
                read_all_equ_stringlen(sysVar,n_sys);
            }
            n_taueq=(int)taueq_sortdecreasing.size();
            strSize=(int)stringlen_sortdecreasing.size();
        }
        if (n_taueq!=strSize) {
            cout << "n_taueq .NOT. equal to strSize!" << "\n";
        }
        cout << "\n"
        << "n_taueq " << n_taueq << "\n"
        << "strSize " << strSize << "\n";
        
        int count_data=0;
        static int count_access=0;
        
        /** write fit data to file **/
        //----------------------------------------------------------------------
        for (int outter=0; outter<2; ++outter) {        // control processes
            for (int inner=0; inner<n_taueq; ++inner) { // control temperatures
                
                count_data=0;
                ++count_access;
                
                /** At each T, get tau_alpha(T) and L(T) **/
                //--------------------------------------------------------------
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
                } else {
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
                
                /** activation "free" energy **/
                //--------------------------------------------------------------
                delmiu = delHa-T_actual*delSa;
                delG   = T_actual*log(pow(10,tau_T)/tau0_vib);
                
                
                /** in first outter loop, collect raw data for AG test **/
                //--------------------------------------------------------------
                if(outter==0) {
                    
                    if (count_data==0) continue;
                    
                    if ((tau_T<fit_tau_threshold)&&
                        (stringlen>stringlen_threshold)) { // NOTE!
                        write_diag_raw
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << tau_T*log(10)
                        << "\n";
                    }
                }
                
                /** in second outter loop, write out AG test results **/
                //--------------------------------------------------------------
                else if (outter==1) {
                    
                    if (inner==0) {
                        
                        write_diag_raw.close();
                        /* NOTE: need to close before it can be read */
                        
                        //------------------------------------------------------
                        write_AGtest
                        << "delu(TA) = " << delHa-TA*delSa << "\n"
                        << "delHa    = " << delHa          << "\n"
                        << "delSa    = " << delSa          << "\n"
                        << "TA       = " << TA             << "\n"
                        << "tauA     = " << tauA           << "\n"
                        << "L(TA)    = " << LA             << "\n"
                        << "a_param  = " << strings_threshold << "\n"
                        << "sigma    = ";
                        write_AGtest << "{";
                        for (size_t i=0;i<sigmaMatrix.size();++i) {
                            for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                if (ii==0) write_AGtest << "{";
                                write_AGtest << sigmaMatrix.at(i).at(ii);
                                if (ii!=sigmaMatrix.at(i).size()-1) {
                                    write_AGtest << ",";
                                } else {
                                    write_AGtest << "}";
                                }
                            } if (i!=sigmaMatrix.size()-1) write_AGtest << ",";
                        } write_AGtest << "}\n";
                        if (is_use_displacement_list) {
                            write_AGtest
                            << "fast%    = "<<fixedfastpercentage<<"\n";
                        } write_AGtest << "\n";
                        
                        write_AGtest_log
                        << "delu(TA) = " << delHa-TA*delSa << "\n"
                        << "delHa    = " << delHa          << "\n"
                        << "delSa    = " << delSa          << "\n"
                        << "TA       = " << TA             << "\n"
                        << "tauA     = " << tauA           << "\n"
                        << "L(TA)    = " << LA             << "\n"
                        << "a_param  = " << strings_threshold << "\n"
                        << "sigma    = ";
                        write_AGtest_log << "{";
                        for (size_t i=0;i<sigmaMatrix.size();++i) {
                            for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                if (ii==0) write_AGtest_log << "{";
                                write_AGtest_log << sigmaMatrix.at(i).at(ii);
                                if (ii!=sigmaMatrix.at(i).size()-1) {
                                    write_AGtest_log << ",";
                                } else {
                                    write_AGtest_log << "}";
                                }
                            } if (i!=sigmaMatrix.size()-1) write_AGtest_log << ",";
                        } write_AGtest_log << "}\n";
                        if (is_use_displacement_list) {
                            write_AGtest_log
                            << "fast%    = "<<fixedfastpercentage<<"\n";
                        } write_AGtest_log << "\n";
                        
                        if (is_write_allsys)
                        {
                            write_AGtest_log_allsys
                            << "\n"
                            << "delu(TA) = " << delHa-TA*delSa << "\n"
                            << "delHa    = " << delHa          << "\n"
                            << "delSa    = " << delSa          << "\n"
                            << "TA       = " << TA             << "\n"
                            << "tauA     = " << tauA           << "\n"
                            << "L(TA)    = " << LA             << "\n"
                            << "a_param  = " << strings_threshold << "\n"
                            << "sigma    = ";
                            write_AGtest_log_allsys << "{";
                            for (size_t i=0;i<sigmaMatrix.size();++i) {
                                for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                    if (ii==0) write_AGtest_log_allsys << "{";
                                    write_AGtest_log_allsys << sigmaMatrix.at(i).at(ii);
                                    if (ii!=sigmaMatrix.at(i).size()-1) {
                                        write_AGtest_log_allsys << ",";
                                    } else {
                                        write_AGtest_log_allsys << "}";
                                    }
                                } if (i!=sigmaMatrix.size()-1) write_AGtest_log_allsys << ",";
                            } write_AGtest_log_allsys << "}\n";
                            if (is_use_displacement_list) {
                                write_AGtest_log_allsys
                                << "fast%    = "<<fixedfastpercentage<<"\n";
                            } write_AGtest_log_allsys << "\n";
                        }
                        //------------------------------------------------------
                        
                        /** fit to diagonal **/
                        //------------------------------------------------------
                        find_diagonal(sysVar,n_sys,"AGdiagonal");
                        
                        if (is_fit_stringmodel) {
                            tau0_shift=tauA*exp(delSa-(delHa/TA));
                            //tau0_shift=exp(lntau0);
                        } else {
                            tau0_shift=exp(lntau0);
                        }
                        
                        //------------------------------------------------------
                        write_AGtest
                        << "tau0  = "  << tau0_shift << "\n"
                        << "slope = "  << slope_fit  << "\n"
                        << "R^2   = "  << rep.r2     << "\n\n"
                        << "T "
                        << convert << "/T "
                        << "Log.tau(T) Z(T) TA/T Z/ZA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                        << "(Z/ZA)(delmiu/kT) Ln(tau/tau0) "
                        << "Z/ZA (KT/delmiu)Ln(tau/tau0) "
                        << "\n";
                        
                        write_AGtest_log
                        << "tau0  = "  << tau0_shift  << "\n"
                        << "slope = "  << slope_fit   << "\n"
                        << "R^2   = "  << rep.r2      << "\n\n"
                        << "T "
                        << convert << "/T "
                        << "Log.tau(T) Z(T) TA/T Z/ZA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                        << "(Z/ZA)(delmiu/kT) Log(tau/tau0) "
                        << "Z/ZA (KT/delmiu)Ln(tau/tau0) "
                        << "\n";
                        if (is_write_allsys)
                        {
                            write_AGtest_log_allsys
                            << "tau0  = "  << tau0_shift  << "\n"
                            << "slope = "  << slope_fit   << "\n"
                            << "R^2   = "  << rep.r2      << "\n\n"
                            << "T "
                            << convert << "/T "
                            << "Log.tau(T) Z(T) TA/T Z/ZA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                            << "(Z/ZA)(delmiu/kT) Log(tau/tau0) "
                            << "Z/ZA (KT/delmiu)Ln(tau/tau0) "
                            << "\n";
                        }
                        //------------------------------------------------------
                    }
                    
                    if (count_data==0) continue;
                    
                    if (stringlen>stringlen_threshold)
                    {
                        write_AGtest
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delmiu << " "
                        << delG/delmiu << " "
                        << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea << " "
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << (tau_T*log(10))-log(tau0_shift)  << " "
                        << (stringlen/LA) << " "
                        << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        
                        write_AGtest_log
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delmiu << " "
                        << delG/delmiu << " "
                        << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea   << " "
                        << (stringlen/LA)*(delmiu/T_actual)*log10(exp(1))  << " "
                        << ((tau_T*log(10))-log(tau0_shift))*log10(exp(1)) << " "
                        << (stringlen/LA) << " "
                        << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        
                        if (is_write_allsys)
                        {
                            write_AGtest_log_allsys
                            << T_actual << " "
                            << convert/T_actual << " "
                            << tau_T << " "
                            << stringlen << " "
                            << TA/T_actual << " "
                            << stringlen/LA << " "
                            << delmiu << " "
                            << delG/delmiu << " "
                            << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea   << " "
                            << (stringlen/LA)*(delmiu/T_actual)*log10(exp(1))  << " "
                            << ((tau_T*log(10))-log(tau0_shift))*log10(exp(1)) << " "
                            << (stringlen/LA) << " "
                            << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        }
                        
                        double x=stringlen/LA;
                        double y=(T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift));
                        diagonal.push_back({x,y});
                        loglog_diagonal.push_back({log10(x),log10(y)});
                    }
                }
            }
        }
        write_AGtest.close();
        write_AGtest_log.close();
    }
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_AG_slope_");
    o.append(sysVar.get_usic());
    o.append(".dat");
    ofstream write_loglog_slope(o.c_str(),ofstream::app);//'append'
    
    string a=to_string(strings_threshold);
    write_message(write_loglog_slope,"string threshold (a) value = "+a);
    
    set_fitParams("linear");
    fitdata_processing_lsfit(loglog_diagonal);
    alglib_lsfit("linear");
    
    write_fitModel(write_loglog_slope,"linear");
    write_fitarrays(write_loglog_slope);
    //
    write_message(write_loglog_slope,"Show original diagonal values");
    fitdata_processing_lsfit(diagonal);//set new x & y data
    write_fitarrays(write_loglog_slope);
    //
    write_stopcond(write_loglog_slope);
    write_fitinfo(write_loglog_slope);
    //write_yvsxby2dvecdub(write_loglog_slope,loglog_diagonal);
    write_errcurve_actual(write_loglog_slope,"linear",loglog_diagonal);
    write_loglog_slope.close();
    
    o.clear();
    o.append(sysVar.get_Path());
    o.append("/simulations");
    o.append("/ATheoryTest_AG_allsys.dat");
    ofstream write_allsys(o.c_str(),ofstream::app);//'append'
    write_allsys
    << sysVar.get_usic() <<" "
    << strings_threshold <<" "
    << c[0] <<" "
    << rep.errpar[0] <<" "
    << "\n";
}





void TheoryTest::RFOTtest(const StructureClass& sysVar)
{
    is_use_Lbeta=false;
    
    stringtype="stringsRg";//use "Rg of strings" (length scale)
    
    write_stringlen_equ(sysVar);
    write_stringlen_equ_avg(sysVar);
    
    vector<vector<double>> diagonal,loglog_diagonal;
    
    string o;
    
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_diagonal_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_diag_raw(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_RFOT_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_RFOTtest(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_RFOT_log_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_RFOTtest_log(o.c_str(),ofstream::app);//'append'
        
        bool is_write_allsys=true;
        o.clear();
        o.append(sysVar.get_Path());
        o.append("/simulations");
        o.append("/ATheoryTest_RFOT_log_"+sysVar.get_usic());
        o.append(".dat");
        ofstream write_RFOTtest_log_allsys(o.c_str(),ofstream::app);//'append'
        
        bool is_cutbyTA=true;
        
        bool is_fit_stringmodel=true;
        
        string stringmodel="transitState";
        /** available models:
         ** "string2param" -- [delH,delS]
         ** "string3param" -- [tau0,delH,delS]
         ** "transitState" -- [delH,delS] **/
        
        /** find TA and L(TA) **/
        //----------------------------------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);
        find_LA(sysVar,n_sys,TA);
        
        if (is_fit_stringmodel) {
            fit_stringmodel(sysVar,n_sys,TA,stringmodel);
        } else {
            delmiuTA = TA*log(tauA/tau0_vib);
            delHa    = Ea;
            delSa    = (delHa-delmiuTA)*pow(TA,-1);
        }
        
        /** smooth data or use raw **/
        //----------------------------------------------------------------------
        int n_taueq=0;
        int strSize=0;
        if (is_data_smoothing) {
            if (is_cutbyTA) {
                fit_taueq(sysVar,n_sys,TA,"COOP");
                fit_stringlen(sysVar,n_sys,TA,"COOP_string");
            } else {
                fit_taueq(sysVar,n_sys,"COOP");
                fit_stringlen(sysVar,n_sys,"COOP_string");
            }
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            strSize=(int)stringlen_sortdecreasing_fit.size();
        } else {
            if (is_cutbyTA) {
                read_all_taueq_data(sysVar,n_sys,TA);
                read_all_equ_stringlen(sysVar,n_sys,TA);
            } else {
                read_all_taueq_data(sysVar,n_sys);
                read_all_equ_stringlen(sysVar,n_sys);
            }
            n_taueq=(int)taueq_sortdecreasing.size();
            strSize=(int)stringlen_sortdecreasing.size();
        }
        if (n_taueq!=strSize) {
            cout << "n_taueq .NOT. equal to strSize!" << "\n";
        }
        cout << "\n"
        << "n_taueq " << n_taueq << "\n"
        << "strSize " << strSize << "\n";
        
        int count_data=0;
        static int count_access=0;
        
        /** write fit data to file **/
        //----------------------------------------------------------------------
        for (int outter=0; outter<2; ++outter) {        // control processes
            for (int inner=0; inner<n_taueq; ++inner) { // control temperatures
                
                count_data=0;
                ++count_access;
                
                /** At each T, get tau_alpha(T) and L(T) **/
                //--------------------------------------------------------------
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
                } else {
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
                
                /** activation "free" energy **/
                //--------------------------------------------------------------
                delmiu = delHa-T_actual*delSa;
                delG   = T_actual*log(pow(10,tau_T)/tau0_vib);
                
                
                /** in first outter loop, collect raw data for RFOT test **/
                //--------------------------------------------------------------
                if(outter==0) {
                    
                    if (count_data==0) continue;
                    
                    if ((tau_T<fit_tau_threshold)&&
                        (stringlen>stringlen_threshold)) { // NOTE!
                        write_diag_raw
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << tau_T*log(10)
                        << "\n";
                    }
                }
                
                /** in second outter loop, write out RFOT test results **/
                //--------------------------------------------------------------
                else if (outter==1) {
                    
                    if (inner==0) {
                        
                        write_diag_raw.close();
                        /* NOTE: need to close before it can be read */
                        
                        //------------------------------------------------------
                        write_RFOTtest
                        << "delu(TA) = " << delHa-TA*delSa << "\n"
                        << "delHa    = " << delHa          << "\n"
                        << "delSa    = " << delSa          << "\n"
                        << "TA       = " << TA             << "\n"
                        << "tauA     = " << tauA           << "\n"
                        << "L(TA)    = " << LA             << "\n"
                        << "a_param  = " << strings_threshold << "\n"
                        << "sigma    = ";
                        write_RFOTtest << "{";
                        for (size_t i=0;i<sigmaMatrix.size();++i) {
                            for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                if (ii==0) write_RFOTtest << "{";
                                write_RFOTtest << sigmaMatrix.at(i).at(ii);
                                if (ii!=sigmaMatrix.at(i).size()-1) {
                                    write_RFOTtest << ",";
                                } else {
                                    write_RFOTtest << "}";
                                }
                            } if (i!=sigmaMatrix.size()-1) write_RFOTtest << ",";
                        } write_RFOTtest << "}\n";
                        if (is_use_displacement_list) {
                            write_RFOTtest
                            << "fast%    = "<<fixedfastpercentage<<"\n";
                        } write_RFOTtest << "\n";
                        
                        write_RFOTtest_log
                        << "delu(TA) = " << delHa-TA*delSa << "\n"
                        << "delHa    = " << delHa          << "\n"
                        << "delSa    = " << delSa          << "\n"
                        << "TA       = " << TA             << "\n"
                        << "tauA     = " << tauA           << "\n"
                        << "L(TA)    = " << LA             << "\n"
                        << "a_param  = " << strings_threshold << "\n"
                        << "sigma    = ";
                        write_RFOTtest_log << "{";
                        for (size_t i=0;i<sigmaMatrix.size();++i) {
                            for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                if (ii==0) write_RFOTtest_log << "{";
                                write_RFOTtest_log << sigmaMatrix.at(i).at(ii);
                                if (ii!=sigmaMatrix.at(i).size()-1) {
                                    write_RFOTtest_log << ",";
                                } else {
                                    write_RFOTtest_log << "}";
                                }
                            } if (i!=sigmaMatrix.size()-1) write_RFOTtest_log << ",";
                        } write_RFOTtest_log << "}\n";
                        if (is_use_displacement_list) {
                            write_RFOTtest_log
                            << "fast%    = "<<fixedfastpercentage<<"\n";
                        } write_RFOTtest_log << "\n";
                        
                        if (is_write_allsys)
                        {
                            write_RFOTtest_log_allsys
                            << "\n"
                            << "delu(TA) = " << delHa-TA*delSa << "\n"
                            << "delHa    = " << delHa          << "\n"
                            << "delSa    = " << delSa          << "\n"
                            << "TA       = " << TA             << "\n"
                            << "tauA     = " << tauA           << "\n"
                            << "L(TA)    = " << LA             << "\n"
                            << "a_param  = " << strings_threshold << "\n"
                            << "sigma    = ";
                            write_RFOTtest_log_allsys << "{";
                            for (size_t i=0;i<sigmaMatrix.size();++i) {
                                for (size_t ii=0;ii<sigmaMatrix.at(i).size();++ii) {
                                    if (ii==0) write_RFOTtest_log_allsys << "{";
                                    write_RFOTtest_log_allsys << sigmaMatrix.at(i).at(ii);
                                    if (ii!=sigmaMatrix.at(i).size()-1) {
                                        write_RFOTtest_log_allsys << ",";
                                    } else {
                                        write_RFOTtest_log_allsys << "}";
                                    }
                                } if (i!=sigmaMatrix.size()-1) write_RFOTtest_log_allsys << ",";
                            } write_RFOTtest_log_allsys << "}\n";
                            if (is_use_displacement_list) {
                                write_RFOTtest_log_allsys
                                << "fast%    = "<<fixedfastpercentage<<"\n";
                            } write_RFOTtest_log_allsys << "\n";
                        }
                        //------------------------------------------------------
                        
                        /** fit to diagonal **/
                        //------------------------------------------------------
                        find_diagonal(sysVar,n_sys,"RFOTdiagonal");
                        
                        if (is_fit_stringmodel) {
                            tau0_shift=tauA*exp(delSa-(delHa/TA));
                            //tau0_shift=exp(lntau0);
                        } else {
                            tau0_shift=exp(lntau0);
                        }
                        
                        //------------------------------------------------------
                        write_RFOTtest
                        << "tau0  = "  << tau0_shift << "\n"
                        << "slope = "  << slope_fit  << "\n"
                        << "R^2   = "  << rep.r2     << "\n\n"
                        
                        << "T "
                        << convert << "/T "
                        << "Log.tau(T) l(T) TA/T l/lA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                        << "(l/lA)(delmiu/kT) Ln(tau/tau0) "
                        << "l/lA (KT/delmiu)Ln(tau/tau0) "
                        << "\n";
                        
                        write_RFOTtest_log
                        << "tau0  = "  << tau0_shift  << "\n"
                        << "slope = "  << slope_fit   << "\n"
                        << "R^2   = "  << rep.r2      << "\n\n"
                        
                        << "T "
                        << convert << "/T "
                        << "Log.tau(T) l(T) TA/T l/lA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                        << "(l/lA)(delmiu/kT) Log(tau/tau0) "
                        << "l/lA (KT/delmiu)Ln(tau/tau0) "
                        << "\n";
                        if (is_write_allsys)
                        {
                            write_RFOTtest_log_allsys
                            << "tau0  = "  << tau0_shift  << "\n"
                            << "slope = "  << slope_fit   << "\n"
                            << "R^2   = "  << rep.r2      << "\n\n"
                            << "T "
                            << convert << "/T "
                            << "Log.tau(T) Z(T) TA/T Z/ZA delmiu delG/delmiu TLn(tau/tauArr)/Ea "
                            << "(Z/ZA)(delmiu/kT) Log(tau/tau0) "
                            << "Z/ZA (KT/delmiu)Ln(tau/tau0) "
                            << "\n";
                        }
                        //------------------------------------------------------
                    }
                    
                    if (count_data==0) continue;
                    
                    if (stringlen>stringlen_threshold)
                    {
                        write_RFOTtest
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delmiu << " "
                        << delG/delmiu << " "
                        << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea << " "
                        << (stringlen/LA)*(delmiu/T_actual) << " "
                        << (tau_T*log(10))-log(tau0_shift)  << " "
                        << (stringlen/LA) << " "
                        << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        
                        write_RFOTtest_log
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << stringlen << " "
                        << TA/T_actual << " "
                        << stringlen/LA << " "
                        << delmiu << " "
                        << delG/delmiu << " "
                        << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea   << " "
                        << (stringlen/LA)*(delmiu/T_actual)*log10(exp(1))  << " "
                        << ((tau_T*log(10))-log(tau0_shift))*log10(exp(1)) << " "
                        << (stringlen/LA) << " "
                        << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        
                        if (is_write_allsys)
                        {
                            write_RFOTtest_log_allsys
                            << T_actual << " "
                            << convert/T_actual << " "
                            << tau_T << " "
                            << stringlen << " "
                            << TA/T_actual << " "
                            << stringlen/LA << " "
                            << delmiu << " "
                            << delG/delmiu << " "
                            << T_actual*log(pow(10,tau_T)/tau0_Arrhenius)/Ea   << " "
                            << (stringlen/LA)*(delmiu/T_actual)*log10(exp(1))  << " "
                            << ((tau_T*log(10))-log(tau0_shift))*log10(exp(1)) << " "
                            << (stringlen/LA) << " "
                            << (T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift)) << "\n";
                        }
                        
                        double x=stringlen/LA;
                        double y=(T_actual/delmiu)*((tau_T*log(10))-log(tau0_shift));
                        diagonal.push_back({x,y});
                        loglog_diagonal.push_back({log10(x),log10(y)});
                    }
                }
            }
        }
        write_RFOTtest.close();
        write_RFOTtest_log.close();
    }
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_RFOT_slope_");
    o.append(sysVar.get_usic());
    o.append(".dat");
    ofstream write_loglog_slope(o.c_str(),ofstream::app);//'append'
    
    string a=to_string(strings_threshold);
    write_message(write_loglog_slope,"string threshold (a) value = "+a);
    
    set_fitParams("linear");
    fitdata_processing_lsfit(loglog_diagonal);
    alglib_lsfit("linear");
    
    write_fitModel(write_loglog_slope,"linear");
    write_fitarrays(write_loglog_slope);
    //
    write_message(write_loglog_slope,"Show original diagonal values");
    fitdata_processing_lsfit(diagonal);//set new x & y data
    write_fitarrays(write_loglog_slope);
    //
    write_stopcond(write_loglog_slope);
    write_fitinfo(write_loglog_slope);
    //write_yvsxby2dvecdub(write_loglog_slope,loglog_diagonal);
    write_errcurve_actual(write_loglog_slope,"linear",loglog_diagonal);
    write_loglog_slope.close();
    
    o.clear();
    o.append(sysVar.get_Path());
    o.append("/simulations");
    o.append("/ATheoryTest_RFOT_allsys.dat");
    ofstream write_allsys(o.c_str(),ofstream::app);//'append'
    write_allsys
    << sysVar.get_usic() <<" "
    << strings_threshold <<" "
    << c[0] <<" "
    << rep.errpar[0] <<" "
    << "\n";
}





void TheoryTest::GLMtest(const StructureClass& sysVar)
{
    string model;
    //model="GLM";//3-param GLM
    model="GLM1";//1-param GLM (tau vs <u2>)
    //model="univGLM";//1-param alternate-GLM form (tau vs T)
    //model="univGT";
    //model="GLM1_hallwolynes";
    
    is_fitu2T=false;
    string u2Tmodel;
    //u2Tmodel="univu2T";//univu2T,linear2
    
    string o;
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_GLM_diagonal_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtestraw(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_u2Tdiagonal_raw_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        string path_u2T=o;
        ofstream write_u2Traw(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_GLM_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtest(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_GLM_log_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_GLMtest_log(o.c_str(),ofstream::app);//'append'
        
        bool is_cutbyTA=true;
        
        /** find TA and <u2(TA)> **/
        //----------------------------------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);//T*
        find_u2A(sysVar,n_sys,TA);//u2*
        find_rhoA(sysVar,n_sys,TA);//density*
        
        /** get highT activation free energy **/
        //----------------------------------------------------------------------
        //delmiuTA = TA*log(tauA/tau0_vib);
        //delHa    = Ea;
        //delSa    = (delHa-delmiuTA)*pow(TA,-1);
        
        /** get a new u2A from fitting <u2> to T/TA **/
        //----------------------------------------------------------------------
        //if (is_fitu2T) {
        //    if (u2Tmodel=="univu2T") {
        //        fit_u2_Tr(sysVar,n_sys,write_u2Traw,is_cutbyTA,u2Tmodel,"u2Tdiagonal");
        //    }
        //}
        
        if (is_reltaur) logtauref=log10(xtauA*tauA);
        
        if (false) logtauref=log10(tauA);//if true, enforce tuaref=tauA
        
        
        /** fit DWF by GLM model **/
        //----------------------------------------------------------------------
        alglib::real_1d_array c_linear;
        vector<alglib::real_1d_array> c_model;
        fit_GLMmodel(sysVar,n_sys,is_cutbyTA,model);
        c_model.push_back(c);
        
        is_fit_univILP=false;//20170703 tmp test
        if (is_fit_univILP) {
            fit_GLMmodel(sysVar,n_sys,is_cutbyTA,"univILP");
            c_model.push_back(c);
        }
        
        /** GLM model with reduced DWF **/
        //----------------------------------------------------------------------
        if (model=="GLM"||model=="GLM1"||model=="GLM1_hallwolynes") {
            write_GLMref(sysVar,n_sys,model);
        }
        
        /** write reference scaling info to file **/
        //----------------------------------------------------------------------
        //if (is_reltaur) {
        //    write_ref_scaling_u2
        //    (sysVar,xtauA,alpha,{TA,log10(tauA),u2A,Tref,logtauref,u2ref});
        //}
        
        
        /** calc Tg and fragility **/
        //----------------------------------------------------------------------
        alglib::real_1d_array c_org=c;;
        set_fitParams("COOP");
        read_all_taueq_data(sysVar,n_sys);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit("COOP");
        double Tg=calc_x_given_y(c,extrp_time,"COOP")[0];
        double mg=calc_x_given_y(c,extrp_time,"COOP")[1];
        c=c_org;//resume original c value
        
        
        /** smooth data or use raw **/
        //----------------------------------------------------------------------
        vector<vector<double>> r2teq;
        std::vector<std::vector<double>> RHO;
        std::vector<std::vector<double>> VOL,VOL2;
        int n_taueq=0;
        int n_dwfac=0;
        int n_therm=0;
        double rho=0,rhostd=0;
        double vol=0,volstd=0;
        double vol2=0,vol2std=0;
        if (is_data_smoothing) {
            if (is_cutbyTA) {
                fit_taueq(sysVar,n_sys,TA,"COOP");
                fit_DWF(sysVar,n_sys,TA,"COOP_DWF");
            } else {
                fit_taueq(sysVar,n_sys,"COOP");
                fit_DWF(sysVar,n_sys,"COOP_DWF");
            }
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            n_dwfac=(int)dwf_sortdecreasing_fit.size();
        } else {
            if (is_cutbyTA) {
                read_all_r2teq_data(sysVar,n_sys,TA);
                r2teq=taueq_sortdecreasing;
                read_all_taueq_data(sysVar,n_sys,TA);
                read_all_equ_DWF(sysVar,n_sys,TA);
                /** thermo data from production log **/
                //-----------------------------------------------
                read_all_thermo_data(sysVar,n_sys,"Density",TA);
                RHO=thermo_sortdecreasing;
                read_all_thermo_data(sysVar,n_sys,"Volume",TA);
                VOL=thermo_sortdecreasing;
                read_all_thermo_data(sysVar,n_sys,"Volume2",TA);
                VOL2=thermo_sortdecreasing;
                //-----------------------------------------------
            } else {
                read_all_r2teq_data(sysVar,n_sys);
                r2teq=taueq_sortdecreasing;
                read_all_taueq_data(sysVar,n_sys);
                read_all_equ_DWF(sysVar,n_sys);
                /** thermo data from production log **/
                //-----------------------------------------------
                read_all_thermo_data(sysVar,n_sys,"Density");
                RHO=thermo_sortdecreasing;
                read_all_thermo_data(sysVar,n_sys,"Volume");
                VOL=thermo_sortdecreasing;
                read_all_thermo_data(sysVar,n_sys,"Volume2");
                VOL2=thermo_sortdecreasing;
                //-----------------------------------------------
            }
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
            n_therm=(int)thermo_sortdecreasing.size();
            /** NOTE:
             ** n_therm may be more than n_taueq & n_dwfac becuase it counts
             ** all data in thermo files, only avoiding repeated T but accepting
             ** data that are not in-equilibrium.
             **/
        }
        if (n_taueq!=n_dwfac)
        {
            cout
            << "in TheoryTest::GLMtest():\n"
            << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
            make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
            cout
            << "in TheoryTest::GLMtest():\n"
            << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
            system_wait(1);
        } else {
            cout
            << "\n"
            << "n_taueq " << n_taueq << "\n"
            << "n_dwfac " << n_dwfac << "\n"
            << "n_therm " << n_therm << "\n";
        }
        
        /** write fit data to file **/
        //----------------------------------------------------------------------
        // ** NOTE: make sure output path is what's expected **
        string o;
        if (false) {
            o.append(sysVar.get_Path());
            o.append("/simulations/");
            o.append(sysVar.get_simType());
            o.append("/");
            o.append(sysVar.get_year());
            o.append("/scaling_universal.dat");
        } else {
            o.append(sysVar.get_Path());
            o.append("/simulations/AGtest/2015");
            o.append("/scaling_universal.dat");
        } ofstream writeUniversal(o.c_str(),ofstream::app);//'append'
        //----------------------------------------------------------------------
        
        
        /** interpolate log(tau/t*) vs log(tau) and return data container **/
        /**------------------------------------------------------------------**/
        double logr2t=0;
        double logSEx=0,logSEy=0,logSEyf=0;
        double logSErx=0,logSEry=0;
        double logSEmx=0,logSEmy=0;
        
        /** interpolate log(tau/t*) vs log(tau) and return data container **/
        //----------------------------------------------------------------------
        vector<vector<double>> logSEvvd,logSEtauAvvd,logSEminvvd;
        if (r2teq.size()!=taueq_sortdecreasing.size()) {
            cout
            << "in TheoryTest::GLMtest(): "
            << "r2teq.size()!=taueq_sortdecreasing.size()\n";
            exit(EXIT_FAILURE);
        }
        for (int i=0;i<(int)r2teq.size();++i) {
            double logtau=taueq_sortdecreasing.at(i).at(1);
            double logSEy=logtau-r2teq.at(i).at(1);
            logSEvvd.push_back({logtau,logSEy});
        }
        build_penalizedspline_interpolant(logSEvvd,0,1);
        alglib::spline1dinterpolant intrp_logSE=interpolant;
        
        /** fit SE decoupling exponent **/
        //----------------------------------------------------------------------
        double epsSEfit=0;
        double r2SEfit=0;
        vector<vector<double>> logSEfit;
        for (int i=0;i<(int)logSEvvd.size();++i) {
            double logSExfit=logSEvvd.at(i).at(0);
            double logSEyfit=logSEvvd.at(i).at(1);
            if (sysVar.get_systemUnit()=="real") logSExfit-=3.0;//use 'ps'
            if (logSExfit>=2.0) { //only take data above log.tau threshold
                logSEfit.push_back({logSExfit,logSEyfit});
            }
        }
        if (logSEfit.size()>=5) { //only fit when # of points >= 5
            set_fitParams("linear");
            fitdata_processing_lsfit(logSEfit);
            alglib_lsfit("linear");
            epsSEfit=(double)c[0];
            r2SEfit =(double)rep.r2;
        }
        
        /** normalize x and y axes by corrsponding values at TA **/
        //----------------------------------------------------------------------
        double logSEatTA=spline1dcalc(intrp_logSE,log10(tauA));//logSE at tauA
        double xr=0,yr=0;
        for (int i=0;i<(int)logSEvvd.size();++i) {
            xr=logSEvvd.at(i).at(0)-log10(tauA);
            yr=logSEvvd.at(i).at(1)-logSEatTA;
            logSEtauAvvd.push_back({xr,yr});
        }
        /** normalize x and y axes by corrsponding values at minimum **/
        //----------------------------------------------------------------------
        vector<vector<double>> data;
        int n_intervals=10000;
        double xlo=logSEvvd.at(0).at(0);
        double xhi=logSEvvd.at(logSEvvd.size()-1).at(0);
        double dx=(xhi-xlo)/double(n_intervals);
        double xin=xlo;
        for (int i=0;i<n_intervals;++i) {
            xin=xlo+dx*(double)i;
            data.push_back({xin,spline1dcalc(intrp_logSE,xin)});
        } sort(data.begin(),data.end(),sortIncreasing1);//sort logSEy in ascending order
        double logSEmin=spline1dcalc(intrp_logSE,data[0][0]);//logSE minimum
        double xrmin=0,yrmin=0;
        for (int i=0;i<(int)logSEvvd.size();++i) {
            xrmin=logSEvvd.at(i).at(0)-data[0][0];
            yrmin=logSEvvd.at(i).at(1)-logSEmin;
            logSEminvvd.push_back({xrmin,yrmin});
        }
        /**------------------------------------------------------------------**/
        
        
        /** interpolate ln(tau) vs 1/T data for (apparent) activation E **/
        //----------------------------------------------------------------------
        string mode_intrp="pspline";//pspline,twopoint
        
        double Eapp_TA=0;
        vector<vector<double>> lntauinvT;
        vector<vector<double>> lntauinvT_allT;
        
        vector<vector<double>> tmpvvd;
        tmpvvd=taueq_sortdecreasing;//to avoid data contamination
        
        /** all T **/
        read_all_taueq_data(sysVar,n_sys);//NOTE: get data at all T
        for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
            //double invT=convert/taueq_sortdecreasing[i][0];
            double invT=pow(taueq_sortdecreasing[i][0],-1);//NOTE: use 1/T
            double lntu=taueq_sortdecreasing[i][1]*log(10);
            lntauinvT_allT.push_back({invT,lntu});
        } taueq_sortdecreasing=tmpvvd;//NOTE: resume data from backed up container
        
        for (int i=0;i<(int)taueq_sortdecreasing.size();++i) {
            //double invT=convert/taueq_sortdecreasing[i][0];
            double invT=pow(taueq_sortdecreasing[i][0],-1);//NOTE: use 1/T
            double lntu=taueq_sortdecreasing[i][1]*log(10);
            lntauinvT.push_back({invT,lntu});
        }
        
        alglib::spline1dinterpolant intrp_lntauinvT;//pspline interpolant
        vector<vector<double>> Eapp_twopoint;//apparent Ea by two-point method
        
        if (mode_intrp=="pspline") {
            size_t s1=lntauinvT_allT.size();
            build_penalizedspline_interpolant(lntauinvT_allT,0,1,s1);//NOTE: use all T
            //intrp_lntauinvT=interpolant;
            double f=0,df=0,d2f=0;
            vector<vector<double>> EavsinvT;
            
            /* append 10 additional points */
            vector<vector<double>> tmpvvd=lntauinvT_allT;
            double x=lntauinvT_allT[lntauinvT_allT.size()-1][0];
            double dx=fabs(lntauinvT_allT[1][0]-lntauinvT_allT[0][0]);
            for (int i=1;i<=10;++i) {
                tmpvvd.push_back({x+dx*i,spline1dcalc(interpolant,x+dx*i)});
            } s1=tmpvvd.size();
            build_penalizedspline_interpolant(tmpvvd,0,1,s1);
            intrp_lntauinvT=interpolant;
            
            for (int i=0;i<(int)lntauinvT_allT.size();++i) {
                spline1ddiff(intrp_lntauinvT,lntauinvT_allT[i][0],f,df,d2f);
                EavsinvT.push_back({lntauinvT_allT[i][0],df});
                cout
                <<lntauinvT_allT[i][0]<<" "// 1/T
                <<lntauinvT_allT[i][1]<<" "// ln(tau)
                <<f<<" "<<df<<"\n";
            }
            /* get Eapp at TA by interpolation of all T's data */
            size_t s2=EavsinvT.size();
            build_penalizedspline_interpolant(EavsinvT,0,1,s2/2);
            Eapp_TA=spline1dcalc(interpolant,1.0/TA);
            
            /* rebuild pspline interpolant for T<TA */
            if (false) {
                size_t s3=lntauinvT.size();
                build_penalizedspline_interpolant(lntauinvT,0,1,s3,3.0);//NOTE: use all T
                intrp_lntauinvT=interpolant;
                for (int i=0;i<(int)lntauinvT.size();++i) {
                    spline1ddiff(intrp_lntauinvT,lntauinvT[i][0],f,df,d2f);
                    cout
                    <<lntauinvT[i][0]<<" "// 1/T
                    <<lntauinvT[i][1]<<" "// ln(tau)
                    <<f<<" "<<df<<"\n";
                }
            }
        }
        else if (mode_intrp=="twopoint") {
            double slope=0; vector<double> slopevd;
            for (int i=0;i<(int)lntauinvT.size()-1;++i) {//NOTE:'lntauinvT_ltTA.size()-1'
                slope=
                (lntauinvT[i+1][1]-lntauinvT[i][1])/
                (lntauinvT[i+1][0]-lntauinvT[i][0]);
                slopevd.push_back(slope);
            } slopevd.push_back(slope);//slope at last point = that of second to last
            
            for (int i=0;i<(int)lntauinvT.size();++i) {
                double T=pow(lntauinvT[i][0],-1);
                if (T<TA) {
                    Eapp_twopoint.push_back({lntauinvT[i][0],slopevd[i]});
                }
            }
            /* get Eapp at TA */
            //NOTE: fit a straight line to data "below TA"
            set_fitParams("linear");
            fitdata_processing_lsfit(Eapp_twopoint);
            alglib_lsfit("linear");
            Eapp_TA=calc_y_given_x(c,1.0/TA,"linear");
        }
        else {
            cout
            << "in TheoryTest::GLMtest():"
            << "mode_intrp should be 'pspline' or 'twopoint'.\n";
            exit(EXIT_FAILURE);
        }
        //----------------------------------------------------------------------
        
        
        for (int outter=0; outter<2; ++outter) {        //control process
            for (int inner=0; inner<n_taueq; ++inner) { //control temperature
                int count_data=0;
                
                logr2t =r2teq[inner][1];
                logSEx =logSEvvd[inner][0];
                logSEy =logSEvvd[inner][1];
                logSEyf=spline1dcalc(intrp_logSE,logSEx);
                logSErx=logSEtauAvvd[inner][0];
                logSEry=logSEtauAvvd[inner][1];
                logSEmx=logSEminvvd[inner][0];
                logSEmy=logSEminvvd[inner][1];
                
                /** At each T, get tau_alpha(T) and DWF(T) **/
                //--------------------------------------------------------------
                if (is_data_smoothing) {
                    T_actual=taueq_sortdecreasing_fit[inner][0];
                    tau_T   =taueq_sortdecreasing_fit[inner][1];
                    for (int i=0; i<n_dwfac; ++i) {
                        if (T_actual==dwf_sortdecreasing_fit[i][0]) {
                            ++count_data;
                            DWF=dwf_sortdecreasing_fit[i][1];
                            break;
                        }
                    }
                } else {
                    T_actual=taueq_sortdecreasing[inner][0];
                    tau_T   =taueq_sortdecreasing[inner][1];
                    for (int i=0; i<n_dwfac; ++i) {
                        if (T_actual==dwf_sortdecreasing[i][0]) {
                            ++count_data;
                            DWF=dwf_sortdecreasing[i][1];
                            break;
                        }
                    }
                    for (int i=0; i<n_therm; ++i) {
                        if (T_actual==RHO[i][0]) {
                            //++count_data;
                            rho=RHO[i][1];
                            rhostd=RHO[i][2];
                            vol=VOL[i][1];
                            volstd=VOL[i][2];
                            vol2=VOL2[i][1];
                            vol2std=VOL2[i][2];
                            break;
                        }
                    }
                }
                
                /** in first outter loop, collect raw data for GLM test **/
                //--------------------------------------------------------------
                if(outter==0) {
                    
                    if (count_data==0) continue;
                    
                    if (tau_T<fit_tau_threshold) {
                        
                        if (model=="GLM")
                        {
                            //NOTE:
                            //this section is originally designed for 3-param
                            //GLM, where tau0 is going to be determined from
                            //the diagonal fit as a shift parameter
                            
                            double xinput=0,yinput=0;
                            xinput=pow(u02/DWF,alpha/2.0);
                            yinput=tau_T*log(10);//no tau0, to be replaced by shift param
                            write_GLMtestraw
                            << xinput << " "
                            << yinput << "\n";
                        }
                        else if (model=="GLM1"||model=="GLM1_hallwolynes")
                        {
                            //NOTE:
                            //for GLM1, there is no y-axis shift parameter, and
                            //this section is reserved only for showing r2 on
                            //a fit to the diagonal line where slope=1
                            //therefore the y and x vlaues are directly from
                            //original GLM1 fit result
                            
                            double xinput=0,yinput=0;
                            if (GLM1mode=="mode1") {
                                xinput=pow(u2A/DWF,alpha/2.0)-1.0;
                                yinput=tau_T*log(10)-log(tauA);
                            } else if (GLM1mode=="mode2") {
                                xinput=pow(c_model[0][1],alpha/2)*(pow(u2A/DWF,alpha/2.0)-1.0);
                                yinput=tau_T*log(10)-log(tauA);
                            } else if (GLM1mode=="mode3") {
                                xinput=(Ea*pow(alpha*TA,-1))*(pow(u2A/DWF,alpha/2.0)-1.0);
                                yinput=tau_T*log(10)-log(tauA);
                            } else {
                                cout
                                << "in TheoryTest::GLMtest(): "
                                << "mode used for GLM1 ("+GLM1mode+") is not found.\n";
                                exit(EXIT_FAILURE);
                            }
                            write_GLMtestraw
                            << xinput << " "
                            << yinput << "\n";
                            
                            //if (is_fitu2T) {
                            //    if (u2Tmodel=="linear2") {
                            //        write_u2Traw
                            //        << T_actual/TA << " "
                            //        << DWF/u2A     << "\n";
                            //    }
                            //}
                        }
                        else if (model=="univGLM")
                        {
                            double shift=c_model[0][1];
                            double preft=c_model[0][2];
                            double Tr=T_actual/TA;
                            double taur=pow(10,tau_T)/tauA;
                            write_GLMtestraw
                            << log(taur) << " " //measured ln(tau/tauA)
                            << pow(preft*Tr-shift,-alpha/2.0)-1.0 //predicted
                            << "\n";
                        }
                    }
                }
                
                /** in second outter loop, write out GLM test results **/
                //--------------------------------------------------------------
                else if (outter==1)
                {
                    if (inner==0)
                    {
                        /* NOTE: need to close before it can be read */
                        write_GLMtestraw.close();
                        write_u2Traw.close();
                        string rm="rm "+path_u2T;
                        if (!is_fitu2T) system(rm.c_str());
                        
                        //------------------------------------------------------
                        if (model=="GLM")
                        {
                            write_GLMtest
                            << "tau0("<<model<<") = " << tau0_GLM << "\n"
                            << "u0^2       = "        << u02      << "\n"
                            << "alpha      = "        << alpha    << "\n\n";
                            write_GLMtest_log
                            << "tau0("<<model<<") = " << tau0_GLM << "\n"
                            << "u0^2       = "        << u02      << "\n"
                            << "alpha      = "        << alpha    << "\n\n";
                        }
                        else if (model=="GLM1"||model=="GLM1_hallwolynes")
                        {
                            write_GLMtest
                            << "Ea       = "          << Ea         << "\n"
                            << "Ea/TA    = "          << Ea/TA      << "\n"
                            << "TA       = "          << TA         << "\n"
                            << "Log.tauA = "          << log10(tauA)<< "\n"
                            << "uA^2     = "          << u2A        << "\n"
                            << "alpha    = "          << alpha      << "\n\n";
                            write_GLMtest_log
                            << "Ea       = "          << Ea         << "\n"
                            << "Ea/TA    = "          << Ea/TA      << "\n"
                            << "TA       = "          << TA         << "\n"
                            << "Log.tauA = "          << log10(tauA)<< "\n"
                            << "uA^2     = "          << u2A        << "\n"
                            << "alpha    = "          << alpha      << "\n\n";
                        }
                        //------------------------------------------------------
                        
                        fit_GLM_diagonal(sysVar,n_sys,"GLMdiagonal");
                        double slope_GLMdiag=slope_fit;
                        alglib::lsfitreport rep_GLMdiag=rep;
                        
                        //double slope_u2Tdiag=0;
                        //alglib::lsfitreport rep_u2Tdiag;
                        //if (is_fitu2T) {
                        //    if (u2Tmodel=="linear2") {
                        //        fit_u2T_linear(sysVar,n_sys,u2Tmodel,"u2Tdiagonal");
                        //        c_linear=c;
                        //        slope_u2Tdiag=slope_fit;
                        //        rep_u2Tdiag=rep;
                        //    }
                        //}
                        
                        //------------------------------------------------------
                        if (model=="GLM")
                        {
                            write_GLMtest
                            << "tau0_shift = " << exp(lntau0)    << "\n"
                            << "slope      = " << slope_GLMdiag  << "\n"
                            << "R^2        = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "Log(tau) Logtau.GLM u2(T) (u02/u2)^(a/2) Ln(tau/tau0)"
                            << "\n";
                            
                            write_GLMtest_log
                            << "tau0_shift = " << exp(lntau0)    << "\n"
                            << "slope      = " << slope_GLMdiag  << "\n"
                            << "R^2        = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "Log(tau) Logtau.GLM u2(T) (u02/u2)^(a/2) Log(tau/tau0)"
                            << "\n";
                        }
                        else if (model=="GLM1"||model=="GLM1_hallwolynes")
                        {
                            write_GLMtest
                            << "slope = " << slope_GLMdiag  << "\n"
                            << "R^2   = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "Log(tau) Logtau.GLM u2(T) ";
                            
                            if (GLM1mode=="mode1") {
                                write_GLMtest
                                <<"(u2A/u2)^(a/2)-1 Ln(tau/tauA)";
                            } else if (GLM1mode=="mode2") {
                                write_GLMtest
                                <<"A^(a/2)*[(u2A/u2)^(a/2)-1] Ln(tau/tauA)";
                            } else if (GLM1mode=="mode3") {
                                write_GLMtest
                                <<"(Ea/a/TA)*[(u2A/u2)^(a/2)-1] Ln(tau/tauA)";
                            } write_GLMtest << "\n";
                            
                            write_GLMtest_log
                            << "slope = " << slope_GLMdiag  << "\n"
                            << "R^2   = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "Log(tau) Logtau.GLM u2(T) ";
                            
                            if (GLM1mode=="mode1") {
                                write_GLMtest_log
                                <<"(u2A/u2)^(a/2)-1 Log(tau/tauA)";
                            } else if (GLM1mode=="mode2") {
                                write_GLMtest_log
                                <<"A^(a/2)*[(u2A/u2)^(a/2)-1] Log(tau/tauA)";
                            } else if (GLM1mode=="mode3") {
                                write_GLMtest_log
                                <<"(Ea/a/TA)*[(u2A/u2)^(a/2)-1] Log(tau/tauA)";
                            } write_GLMtest_log << "\n";
                        }
                        else if (model=="univGLM")
                        {
                            write_GLMtest
                            << "slope = " << slope_GLMdiag  << "\n"
                            << "R^2   = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "T/TA tau/tauA Logtau.m Logtau.p Ln(tau/tauA).m Ln(tau/tauA).p"
                            << "\n";
                            
                            write_GLMtest_log
                            << "slope = " << slope_GLMdiag  << "\n"
                            << "R^2   = " << rep_GLMdiag.r2 << "\n\n"
                            
                            << "T "
                            << convert << "/T "
                            << "T/TA tau/tauA Logtau.m Logtau.p Log(tau/tauA).m Log(tau/tauA).p"
                            << "\n";
                        }
                        //------------------------------------------------------
                    }
                    
                    if (count_data==0) continue;
                    
                    if (model=="GLM")
                    {
                        vector<double> inpVec={logtauref,u2ref,DWF,T_actual,tau_T,tauA,u2A};
                        double tauGLM=calc_time_given_vecinput(c_model[0],inpVec,model);
                        
                        write_GLMtest
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << log10(tauGLM) << " "
                        << DWF << " "
                        << pow(u02/DWF,alpha/2.0) << " "
                        << (tau_T*log(10))-lntau0
                        << "\n";
                        
                        write_GLMtest_log
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << log10(tauGLM) << " "
                        << DWF << " "
                        << (pow(u02/DWF,alpha/2.0))*log10(exp(1)) << " "
                        << ((tau_T*log(10))-lntau0)*log10(exp(1))
                        << "\n";
                    }
                    else if(model=="GLM1"||model=="GLM1_hallwolynes")
                    {
                        vector<double> inpVec={logtauref,u2ref,DWF,T_actual,tau_T,tauA,u2A};
                        double tauGLM=calc_time_given_vecinput(c_model[0],inpVec,model);
                        
                        double lntaurGLM=0;
                        if (GLM1mode=="mode1") {
                            lntaurGLM = pow(u2A/DWF,alpha/2)-1.0;
                        } else if (GLM1mode=="mode2") {
                            lntaurGLM = pow(c_model[0][1],alpha/2.0)*(pow(u2A/DWF,alpha/2.0)-1.0);
                        } else if (GLM1mode=="mode3") {
                            lntaurGLM = (Ea*pow(alpha*TA,-1))*(pow(u2A/DWF,alpha/2.0)-1.0);
                        }
                        
                        write_GLMtest
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << log10(tauGLM) << " "
                        << DWF << " "
                        << lntaurGLM << " "
                        << tau_T*log(10)-log(tauA)
                        << "\n";
                        
                        write_GLMtest_log
                        << T_actual << " "
                        << convert/T_actual << " "
                        << tau_T << " "
                        << log10(tauGLM) << " "
                        << DWF << " "
                        << lntaurGLM*log10(exp(1)) << " "
                        << (tau_T*log(10)-log(tauA))*log10(exp(1))
                        << "\n";
                    }
                    else if (model=="univGLM")
                    {
                        double shift=c_model[0][1];
                        double preft=c_model[0][2];
                        double Tr=T_actual/TA;
                        double taur=pow(10,tau_T)/tauA;
                        double lntauA=log(tauA);
                        double lntaur_p=pow(preft*Tr-shift,-alpha/2.0)-1.0;
                        double lntau_p=lntauA+lntaur_p;
                        
                        write_GLMtest
                        << T_actual << " "
                        << convert/T_actual << " "
                        << Tr << " "
                        << taur << " "
                        << tau_T << " "
                        << lntau_p*log10(exp(1)) << " "
                        << log(taur) << " "
                        << lntaur_p
                        << "\n";
                        
                        write_GLMtest_log
                        << T_actual << " "
                        << convert/T_actual << " "
                        << Tr << " "
                        << taur << " "
                        << tau_T << " "
                        << lntau_p*log10(exp(1)) << " "
                        << log10(taur) << " "
                        << lntaur_p*log10(exp(1))
                        << "\n";
                    }
                    /** write scaling info to file to a same file **/
                    //----------------------------------------------------------
                    
                    /* Apparent Ea from ln(tau) vs 1/T interpolation */
                    //----------------------------------------------------------
                    //double invT=convert/T_actual;
                    double invT=1.0/T_actual;
                    double Eapp=0,Eapp_r=0;
                    if (mode_intrp=="pspline") {
                        double f=0,df=0,d2f=0;
                        spline1ddiff(intrp_lntauinvT,invT,f,df,d2f);
                        Eapp=df;
                    }
                    else if (mode_intrp=="twopoint") {
                        if (Eapp_twopoint.size()!=n_taueq) {
                            cout
                            << "in TheoryTest::GLMtest():\n"
                            << "Eapp_twopoint.size()!=n_taueq \n"
                            << "Eapp_twopoint.size() = "<<Eapp_twopoint.size()<<"\n"
                            << "n_taueq              = "<<n_taueq<<"\n";
                            exit(EXIT_FAILURE);
                        } Eapp=Eapp_twopoint[inner][1];
                    } Eapp_r=Eapp/Eapp_TA;
                    //----------------------------------------------------------
                    
                    vector<vector<double>> import;
                    
                    vector<double> measured = {
                        TA,log10(tauA),u2A,T_actual,tau_T,DWF,Tg,mg,rhoA,Ea,Eapp,Eapp_r
                    };  import.push_back(measured);
                    
                    vector<double> thermo = {
                        rho,rhostd,vol,volstd,vol2,vol2std
                    };  import.push_back(thermo);
                    
                    vector<double> logSE = {
                        logr2t,
                        logSEx,logSEy,logSEyf,
                        logSErx,logSEry,logSEmx,logSEmy,
                        epsSEfit,r2SEfit
                    };  import.push_back(logSE);
                    
                    vector<string> models = {
                        model,u2Tmodel
                    };
                    
                    if (model=="GLM1"||model=="univGLM"||model=="GLM1_hallwolynes") {
                        write_universal_scaling_u2
                        (sysVar,models,c_model,c_linear,import,writeUniversal);
                    } else if (model=="univGT") {
                        write_universal_scaling_dfit(sysVar,models,c_model[0],c_linear,measured);
                    }
                }
            }
        } writeUniversal << "\n";
        write_GLMtest.close();
        write_GLMtest_log.close();
    }
}





void TheoryTest::HWtest(const StructureClass& sysVar)
{
    string model;
    //model="hallwolynes";
    model="hallwolynes1param";
    
    string o;
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_hallwolynes_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_hw(o.c_str(),ofstream::app); // 'append'
        
        string o;
        o.append(sysVar.get_Path());
        o.append("/simulations/");
        o.append(sysVar.get_simType());
        o.append("/");
        o.append(sysVar.get_year());
        o.append("/HallWolynes_test_all.dat");
        ofstream writeFile(o.c_str(),ofstream::app);//'append'
        
        bool is_cutbyTA=true;
        
        /** find TA and <u2(TA)> **/
        //----------------------------------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);
        find_u2A(sysVar,n_sys,TA);
        
        /** fit tau(u2) **/
        //----------------------------------------------------------------------
        fit_HWmodel(sysVar,n_sys,is_cutbyTA,model);
        alglib::real_1d_array c_model=c;
        
        /** smooth data or use raw **/
        //----------------------------------------------------------------------
        int n_taueq=0;
        int n_dwfac=0;
        if (is_data_smoothing) {
            if (is_cutbyTA) {
                fit_taueq(sysVar,n_sys,TA,"COOP");
                fit_DWF(sysVar,n_sys,TA,"COOP_DWF");
            } else {
                fit_taueq(sysVar,n_sys,"COOP");
                fit_DWF(sysVar,n_sys,"COOP_DWF");
            }
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            n_dwfac=(int)dwf_sortdecreasing_fit.size();
        } else {
            if (is_cutbyTA) {
                read_all_taueq_data(sysVar,n_sys,TA);
                read_all_equ_DWF(sysVar,n_sys,TA);
            } else {
                read_all_taueq_data(sysVar,n_sys);
                read_all_equ_DWF(sysVar,n_sys);
            }
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
        }
        if (n_taueq!=n_dwfac)
        {
            cout
            << "in TheoryTest::HWtest():\n"
            << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
            make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
            cout
            << "in TheoryTest::HWtest():\n"
            << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
            system_wait(1);
        } else {
            cout
            << "\n"
            << "n_taueq " << n_taueq << "\n"
            << "n_dwfac " << n_dwfac << "\n";
        }
        
        static int count_access=0;
        
        for (int inner=0; inner<n_taueq; ++inner) { // control temperatures
            
            ++count_access;
            
            /** At each T, get tau_alpha(T) and DWF(T) **/
            //--------------------------------------------------------------
            if (is_data_smoothing) {
                T_actual=taueq_sortdecreasing_fit[inner][0];
                tau_T   =taueq_sortdecreasing_fit[inner][1];
                for (int i=0; i<n_dwfac; ++i) {
                    if (T_actual==dwf_sortdecreasing_fit[i][0]) {
                        DWF=dwf_sortdecreasing_fit[i][1];
                        break;
                    }
                }
            } else {
                T_actual=taueq_sortdecreasing[inner][0];
                tau_T   =taueq_sortdecreasing[inner][1];
                for (int i=0; i<n_dwfac; ++i) {
                    if (T_actual==dwf_sortdecreasing[i][0]) {
                        DWF=dwf_sortdecreasing[i][1];
                        break;
                    }
                }
            }
            
            if (count_access==1)
            {
                if (model=="hallwolynes")
                {
                    write_hw
                    << "TA         = " << TA          << "\n"
                    << "Log.tauA   = " << log10(tauA) << "\n"
                    << "R^2        = " << rep.r2      << "\n\n"
                    
                    << "T "
                    << convert << "/T "
                    << "u2(T) Log(tau) Logtau.HW u20/u2 Log(tau/tau0)"
                    << "\n";
                }
                else if (model=="hallwolynes1param")
                {
                    write_hw
                    << "TA         = " << TA          << "\n"
                    << "Log.tauA   = " << log10(tauA) << "\n"
                    << "R^2        = " << rep.r2      << "\n\n"
                    
                    << "T "
                    << convert << "/T "
                    << "u2(T) Log(tau) Logtau.HW1p u20r*(1/u2r-1) Log(taur)"
                    << "\n";
                }
            }
            
            if (model=="hallwolynes")
            {
                double logtauHW=
                log10(calc_time_given_vecinput(c_model,{DWF,T_actual,tau_T},model));
                
                if (sysVar.get_systemUnit()=="real") {
                    write_hw
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T-3.0 << " "
                    << logtauHW-3.0 << " "
                    << (c_model[1]/DWF)*log10(exp(1)) << " "
                    << tau_T-log10(c_model[0])
                    << "\n";
                    writeFile
                    << sysVar.get_usic() <<" "
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T-3.0 << " "
                    << logtauHW-3.0 << " "
                    << (c_model[1]/DWF)*log10(exp(1)) << " "
                    << tau_T-log10(c_model[0])
                    << "\n";
                } else {
                    write_hw
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T << " "
                    << logtauHW << " "
                    << (c_model[1]/DWF)*log10(exp(1)) << " "
                    << tau_T-log10(c_model[0])
                    << "\n";
                    writeFile
                    << sysVar.get_usic() <<" "
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T << " "
                    << logtauHW << " "
                    << (c_model[1]/DWF)*log10(exp(1)) << " "
                    << tau_T-log10(c_model[0])
                    << "\n";
                }
            }
            else if (model=="hallwolynes1param")
            {
                double logtauHW1p=
                log10(calc_time_given_vecinput(c_model,{DWF,u2A,tauA,T_actual,tau_T},model));
                
                if (sysVar.get_systemUnit()=="real") {
                    write_hw
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T-3.0 << " "
                    << logtauHW1p-3.0 << " "
                    << ((c_model[0]/u2A)*((u2A/DWF)-1))*log10(exp(1)) << " "
                    << tau_T-log10(tauA)
                    << "\n";
                    writeFile
                    << sysVar.get_usic() <<" "
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T-3.0 << " "
                    << logtauHW1p-3.0 << " "
                    << ((c_model[0]/u2A)*((u2A/DWF)-1))*log10(exp(1)) << " "
                    << tau_T-log10(tauA)
                    << "\n";
                } else {
                    write_hw
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T << " "
                    << logtauHW1p << " "
                    << ((c_model[0]/u2A)*((u2A/DWF)-1))*log10(exp(1)) << " "
                    << tau_T-log10(tauA)
                    << "\n";
                    writeFile
                    << sysVar.get_usic() <<" "
                    << T_actual << " "
                    << convert/T_actual << " "
                    << DWF << " "
                    << tau_T << " "
                    << logtauHW1p << " "
                    << ((c_model[0]/u2A)*((u2A/DWF)-1))*log10(exp(1)) << " "
                    << tau_T-log10(tauA)
                    << "\n";
                }
            }
        }
        write_hw.close();
        writeFile.close();
    }
}





void TheoryTest::Leporinitest(const StructureClass& sysVar)
{
    string model="leporiniref";
    
    string o;
    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/ATheoryTest_Leporini_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_leporini(o.c_str(),ofstream::app); // 'append'
        
        o.clear();
        o.append(sysVar.get_Path());
        o.append("/simulations/");
        o.append(sysVar.get_simType());
        o.append("/");
        o.append(sysVar.get_year());
        o.append("/Leporini_test_all.dat");
        ofstream writeFile(o.c_str(),ofstream::app);//'append'
        
        bool is_cutbyTA=true;
        
        /** find TA and <u2(TA)> **/
        //----------------------------------------------------------------------
        read_all_taueq_data(sysVar,n_sys);
        find_TA(sysVar,n_sys,taueq_sortdecreasing);
        find_u2A(sysVar,n_sys,TA);
        
        /** fit reduced DWF by Leporini universal scaling **/
        //----------------------------------------------------------------------
        is_1paramLeporini=true;
        fit_Leporinimodel(sysVar,n_sys,is_cutbyTA,model);
        alglib::real_1d_array c_model=c;
        
        /** smooth data or use raw **/
        //----------------------------------------------------------------------
        int n_taueq=0;
        int n_dwfac=0;
        if (is_data_smoothing) {
            if (is_cutbyTA) {
                fit_taueq(sysVar,n_sys,TA,"COOP");
                fit_DWF(sysVar,n_sys,TA,"COOP_DWF");
            } else {
                fit_taueq(sysVar,n_sys,"COOP");
                fit_DWF(sysVar,n_sys,"COOP_DWF");
            }
            n_taueq=(int)taueq_sortdecreasing_fit.size();
            n_dwfac=(int)dwf_sortdecreasing_fit.size();
        } else {
            if (is_cutbyTA) {
                read_all_taueq_data(sysVar,n_sys,TA);
                read_all_equ_DWF(sysVar,n_sys,TA);
            } else {
                read_all_taueq_data(sysVar,n_sys);
                read_all_equ_DWF(sysVar,n_sys);
            }
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
        }
        if (n_taueq!=n_dwfac)
        {
            cout
            << "in TheoryTest::Leporinitest():\n"
            << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
            make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
            n_taueq=(int)taueq_sortdecreasing.size();
            n_dwfac=(int)dwf_sortdecreasing.size();
            cout
            << "in TheoryTest::Leporinitest():\n"
            << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
            system_wait(1);
        } else {
            cout
            << "\n"
            << "n_taueq " << n_taueq << "\n"
            << "n_dwfac " << n_dwfac << "\n";
        }
        
        static int count_access=0;
        
        for (int inner=0; inner<n_taueq; ++inner) { // control temperatures
            
            ++count_access;
            
            /** At each T, get tau_alpha(T) and DWF(T) **/
            //--------------------------------------------------------------
            if (is_data_smoothing) {
                T_actual=taueq_sortdecreasing_fit[inner][0];
                tau_T   =taueq_sortdecreasing_fit[inner][1];
                for (int i=0; i<n_dwfac; ++i) {
                    if (T_actual==dwf_sortdecreasing_fit[i][0]) {
                        DWF=dwf_sortdecreasing_fit[i][1];
                        break;
                    }
                }
            } else {
                T_actual=taueq_sortdecreasing[inner][0];
                tau_T   =taueq_sortdecreasing[inner][1];
                for (int i=0; i<n_dwfac; ++i) {
                    if (T_actual==dwf_sortdecreasing[i][0]) {
                        DWF=dwf_sortdecreasing[i][1];
                        break;
                    }
                }
            }
            
            if (count_access==1)
            {
                write_leporini
                << "TA         = " << TA          << "\n"
                << "Log.tauA   = " << log10(tauA) << "\n"
                << "Tref       = " << Tref        << "\n"
                << "Log.tauref = " << logtauref   << "\n"
                << "u2ref      = " << u2ref       << "\n"
                << "R^2        = " << rep.r2      << "\n\n"
                
                << "Yg         = " << Yg          << "\n"
                << "mg         = " << mg          << "\n"
                << "u2g        = " << u2g         << "\n\n"
                
                << "T "
                << convert << "/T "
                << "u2(T) u2ref/u2 Log(tau) Logtau.Leporini"
                << "\n";
            }
            
            if (sysVar.get_systemUnit()=="real") {
                write_leporini
                << T_actual << " "
                << convert/T_actual << " "
                << DWF << " "
                << u2ref*pow(DWF,-1) << " "
                << tau_T-3.0 << " "
                << log10(calc_time_given_vecinput(c_model,{logtauref,u2ref,DWF},model))-3.0
                << "\n";
                writeFile
                << sysVar.get_usic() <<" "
                << T_actual << " "
                << convert/T_actual << " "
                << DWF << " "
                << u2ref*pow(DWF,-1) << " "
                << tau_T-3.0 << " "
                << log10(calc_time_given_vecinput(c_model,{logtauref,u2ref,DWF},model))-3.0
                << "\n";
            } else {
                write_leporini
                << T_actual << " "
                << convert/T_actual << " "
                << DWF << " "
                << u2ref*pow(DWF,-1) << " "
                << tau_T << " "
                << log10(calc_time_given_vecinput(c_model,{logtauref,u2ref,DWF},model))
                << "\n";
                writeFile
                << sysVar.get_usic() <<" "
                << T_actual << " "
                << convert/T_actual << " "
                << DWF << " "
                << u2ref*pow(DWF,-1) << " "
                << tau_T << " "
                << log10(calc_time_given_vecinput(c_model,{logtauref,u2ref,DWF},model))
                << "\n";
            }
        }
        write_leporini.close();
        writeFile.close();
        
        /** fit the universal scaling form from Leporini et al. **/
        fit_Leporinimodel(sysVar,n_sys,is_cutbyTA,"leporini_universal");
    }
}





void TheoryTest::find_mean_stringlen(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const double frametime,
                                     const int frame)
{
    vector<vector<double>> n_distr; // (n,P(n))
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/strings_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_string(o.c_str());
    
    string i;
    if (is_stringmb)
    {
        cout
        << "in TheoryTest::find_mean_stringlen():\n"
        << "method for stringlen from string multibodies not implemented\n";
        exit(EXIT_FAILURE);
    }
    
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/"+cusstrfolder+"/strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        int intData=0;
        vector<int> vec_intData;
        double dubData=0;
        vector<double> vec_dubData;
        string lineContent;
        string strData;
        vector<string> vec_strData;
        
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in TheoryTest::find_mean_stringlen():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        } write_string << lineContent << "\n";
        
        getline(readFile,lineContent);// 2nd line
        istringstream iss(lineContent);
        for (int i=0; i<5; ++i) {
            // first five are statistics numbers
            iss >> strData;
            vec_strData.push_back(strData);
        } while (iss>>intData) {
            // the rest are string length distribution
            vec_intData.push_back(intData);
        }
        getline(readFile,lineContent); // 3rd line
        istringstream iss1(lineContent);
        for (int i=0; i<5; ++i) {
            iss1 >> dubData;
            write_string
            << vec_strData[i] << " "
            << dubData << "\n";
            if      (i==0) frame_time         = dubData;
            else if (i==1) mean_strings       = dubData;
            else if (i==2) mean_length        = dubData;
            else if (i==3) mean_length_count1 = dubData;
            else if (i==4) order_parameter    = dubData;
        } readFile.close();
        if (is_use_counting1) stringlen=mean_length_count1;
        else stringlen=mean_length;
        
        write_string << "\nL  P(L)_raw\n";
        int count_index=0;
        while (iss1>>dubData) {
            write_string
            << vec_intData.at(count_index) << " "
            << dubData << "\n";
            n_distr.push_back({(double)count_index,dubData});
            ++count_index;
        } write_string << "\n";
    } else {
        cout
        << "in TheoryTest::find_mean_stringlen():\n"
        << i << " cannot open.\n";
    }
    if (is_backwardExtrp)
    {
        /** impose cutoff on P(n) **/
        vector<vector<double>> fit_n_distr;
        for (int i=0; i<(int)n_distr.size(); ++i) {
            if (n_distr.at(i).at(1)>n_distr_cutoff) {
                fit_n_distr.push_back({n_distr.at(i).at(0),n_distr.at(i).at(1)});
            }
        }
        try {
            if(fit_n_distr.size()==0) throw 0;
            std::sort(fit_n_distr.begin(),fit_n_distr.end(),sortDecreasing0);
        } catch (int i) {
            cout <<
            sysVar.get_usic()+"_00"+to_string((long long int)n_trl)+
            "_"+sysVar.get_nameString(n_sys)+"_T"+to_string((long long int)Temp_d)
            << "\n"
            << "in TheoryTest::find_mean_stringlen():\n"
            << "fit_n_distr size = " << i << "\n";
            fit_n_distr.push_back({0,0});
        }
        /** fit P(n) to a chosen model **/
        set_fitParams(model_stringlen_distr);
        fitdata_processing_lsfit(fit_n_distr);
        alglib_lsfit(model_stringlen_distr);
        write_fitModel(write_string,model_stringlen_distr);
        write_fitarrays(write_string);
        write_stopcond(write_string);
        write_fitinfo(write_string);
        
        //p1_extrp=calc_y_given_x(c,1.0,model_stringlen_distr); // P(n=1)_extrapolated
        //n_distr[1][1]=p1_extrp; // set value to P(n=1)
        stringlen=1+c[1]; // analytic 1st moment: 1+B, Re[1/B]>0
        write_string
        << "mean_length "
        << stringlen
        << "\n";
    }
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T frame_index time stringlen
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringlen
    << frame
    << " "
    << frametime
    << " "
    << stringlen
    << "\n";
    
    write_stringleninv
    << convert/T_actual
    << " "
    << frame
    << " "
    << frametime
    << " "
    << stringlen
    << "\n";
    
    write_string.close();
    write_stringlen.close();
    write_stringleninv.close();
}





double TheoryTest::get_frame_time(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d,
                                  const int frame)
{
    vector<vector<double>> ngp_data;
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/ngp/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        string strData;
        int    frame_index=0;//zero-based index
        double time=0;
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in TheoryTest::find_time_stringlen():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        } //write_string << lineContent << "\n";
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> time;//time
            if (frame_index==frame) {
                frame_time=time;
                break;
            } ++frame_index;//NOTE: zero-based
        } readFile.close();
    } else {
        cout
        << "in TheoryTest::get_frame_time():\n"
        << i << " cannot open.\n";
        exit(EXIT_FAILURE);
    } return frame_time;
}





void TheoryTest::find_time_stringlen(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const int frame)
{
    vector<vector<double>> n_distr; // (n,P(n))
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlenDistr_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_string(o.c_str(),ofstream::app); // 'append'
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/"+cusstrfolder+"/strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("_frame"+to_string((long long int)frame)); // NOTE!
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    
    if (is_stringmb)
    { /* Front is_stringmb=true */
        
        if (readFile.is_open())
        {
            int intData=0;
            vector<int> vec_intData;
            double dubData=0;
            vector<double> size,freq;
            vector<double> ith_moment,moment;
            string lineContent;
            string strData;
            vector<string> keyWords;
            
            /** 1st:amdat_version **/
            getline(readFile,lineContent);//amdat_version
            if (extract_amdat_svn(lineContent)!=amdat_svn) {
                if (is_check_amdat_version) {
                    cout
                    << "in TheoryTest::find_time_stringlen():\n"
                    << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                    << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
                }
            }
            /** 2nd:String Size **/
            getline(readFile,lineContent);
            istringstream iss(lineContent);
            iss >> strData; keyWords.push_back(strData);//0:Size
            while (iss>>dubData) size.push_back(dubData);
            /** 3rd:String Size Frequency **/
            getline(readFile,lineContent);iss.clear();iss.str(lineContent);
            iss >> strData; keyWords.push_back(strData);//1:Frequency
            while (iss>>dubData) freq.push_back(dubData);
            /** 4th: blank **/
            getline(readFile,lineContent);
            /** 5th :String Moment **/
            getline(readFile,lineContent);iss.clear();iss.str(lineContent);
            iss >> strData; keyWords.push_back(strData);//2:Moment
            while (iss>>dubData) ith_moment.push_back(dubData);
            if (ith_moment.size()!=n_moments) {
                cout
                << "TheoryTest::find_time_stringlen():\n"
                << "n_moment from file ("<<ith_moment.size()<<") != "
                << "n_moment from AMDAT ("<<n_moments<<") \n";
                exit(EXIT_FAILURE);
            }
            /** 6th: Moment Value **/
            getline(readFile,lineContent);iss.clear();iss.str(lineContent);
            iss >> strData; keyWords.push_back(strData);//3:Value
            while (iss>>dubData) moment.push_back(dubData);
            /** 7th: blank **/
            getline(readFile,lineContent);
            /** 8th: Total_Multibodies **/
            getline(readFile,lineContent);iss.clear();iss.str(lineContent);
            iss >> strData; keyWords.push_back(strData);//4:Total_Multibodies
            iss >> intData; //int Total_Multibodies=intData;
            /** 9th: Mean_Multibodies **/
            getline(readFile,lineContent);iss.clear();iss.str(lineContent);
            iss >> strData; keyWords.push_back(strData);//5:Mean_Multibodies
            iss >> intData; //int Mean_Multibodies=intData;
            
            frame_time=get_frame_time(sysVar,n_trl,n_sys,Temp_d,frame);
            
            bool is_n_distr_cutoff=false;//NOTE
            
            vector<vector<double>> n_distr,fit_n_distr;
            n_distr.push_back(size);
            n_distr.push_back(freq);
            transpose(n_distr);
            
            if (is_n_distr_cutoff) {
                for (int i=0;i<n_distr.size();++i) {
                    if (n_distr.at(i).at(1)>n_distr_cutoff) {
                        fit_n_distr.push_back({n_distr.at(i).at(0),n_distr.at(i).at(1)});
                    } else {
                        fit_n_distr.push_back({n_distr.at(i).at(0),0.0});//NOTE:P(n)=0
                    }
                }
                double sum=0;
                for (int i=0;i<fit_n_distr.size();++i) {
                    sum+=fit_n_distr.at(i).at(1);
                }
                for (int i=0;i<fit_n_distr.size();++i) {
                    fit_n_distr.at(i).at(1)*=(1.0/sum);
                } sum=0;
                for (int i=0;i<fit_n_distr.size();++i) {
                    sum+=(fit_n_distr.at(i).at(0)*fit_n_distr.at(i).at(1));
                } stringlen=sum;
            } else {
                fit_n_distr=n_distr;
                for (int i=0;i<ith_moment.size();++i) {
                    if (ith_moment.at(i)==1) {
                        stringlen=moment.at(i);
                        break;
                    }
                }
            }
            
            write_string
            << "frame_"<<frame<<"; "<< "time_"<<frame_time<<"\n"
            << keyWords.at(0)<<" "<<keyWords.at(1)<<"\n";
            for (int i=0;i<fit_n_distr.size();++i) {
                write_string
                << fit_n_distr.at(i).at(0) << " "
                << fit_n_distr.at(i).at(1) << "\n";
            } write_string << "\n1st_moment "<<stringlen<<"\n\n";
            
        } else {
            cout
            << "in TheoryTest::find_time_stringlen():\n"
            << i << " cannot open.\n";
            exit(EXIT_FAILURE);
        }
    } /* End is_stringmb=true */
    else
    { /* Front is_stringmb=false */
        
        if (readFile.is_open())
        {
            int intData=0;
            vector<int> vec_intData;
            double dubData=0;
            vector<double> vec_dubData;
            string lineContent;
            string strData;
            vector<string> vec_strData;
            
            getline(readFile,lineContent);//amdat_version
            if (extract_amdat_svn(lineContent)!=amdat_svn) {
                if (is_check_amdat_version) {
                    cout
                    << "in TheoryTest::find_time_stringlen():\n"
                    << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                    << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
                }
            } //write_string << lineContent << "\n";
            
            getline(readFile,lineContent); // 2nd line
            istringstream iss(lineContent);
            for (int i=0; i<5; ++i) {
                // first five are statistics numbers
                iss >> strData;
                vec_strData.push_back(strData);
            } while (iss>>intData) {
                // the rest are string length distribution
                vec_intData.push_back(intData);
            }
            getline(readFile,lineContent); // 3rd line
            istringstream iss1(lineContent);
            write_string << "Frame_" << frame << "\n";
            for (int i=0; i<5; ++i) {
                iss1 >> dubData;
                write_string
                << vec_strData[i] << " "
                << dubData << "\n";
                if      (i==0) frame_time         = dubData;
                else if (i==1) mean_strings       = dubData;
                else if (i==2) mean_length        = dubData;
                else if (i==3) mean_length_count1 = dubData;
                else if (i==4) order_parameter    = dubData;
            } readFile.close();
            if (is_use_counting1) stringlen=mean_length_count1;
            else stringlen=mean_length;
            
            write_string << "\nL  P(L)_raw\n";
            int count_index=0;
            while (iss1>>dubData) {
                write_string
                << vec_intData.at(count_index) << " "
                << dubData << "\n";
                n_distr.push_back({(double)count_index,dubData});
                ++count_index;
            } write_string << "\n";
        } else {
            cout
            << "in TheoryTest::find_time_stringlen():\n"
            << i << " cannot open.\n";
            exit(EXIT_FAILURE);
        }
        if (is_backwardExtrp)
        {
            /** impose cutoff on P(n) **/
            vector<vector<double>> fit_n_distr;
            for (int i=0; i<(int)n_distr.size(); ++i) {
                if (n_distr.at(i).at(1)>n_distr_cutoff) {
                    fit_n_distr.push_back
                    ({n_distr.at(i).at(0),n_distr.at(i).at(1)});
                }
            }
            try {
                if(fit_n_distr.size()==0) throw 0;
                std::sort(fit_n_distr.begin(),fit_n_distr.end(),sortDecreasing0);
            } catch (int i) {
                cout <<
                sysVar.get_usic()+"_00"+to_string((long long int)n_trl)+
                "_"+sysVar.get_nameString(n_sys)+"_T"+to_string((long long int)Temp_d)
                +"_frame"+to_string((long long int)frame)
                << "\n"
                << "in TheoryTest::find_time_stringlen():" << "\n"
                << "fit_n_distr size = " << i << "\n";
                fit_n_distr.push_back({0,0});
            }
            if (fit_n_distr.size()>=3) {
                /** fit P(n) to a chosen model **/
                set_fitParams(model_stringlen_distr);
                fitdata_processing_lsfit(fit_n_distr);
                alglib_lsfit(model_stringlen_distr);
                write_fitModel(write_string,model_stringlen_distr);
                write_fitarrays(write_string);
                write_stopcond(write_string);
                write_fitinfo(write_string);
            }
            if (rep.r2>0.8) {
                //p1_extrp=calc_y_given_x(c,1.0,model_stringlen_distr); // P(n=1)_extrapolated
                //n_distr[1][1]=p1_extrp; // set value to P(n=1)
                stringlen=1+c[1]; // analytic 1st moment: 1+B, Re[1/B]>0
                write_string
                << "r2(frame_"<<frame<<")=" << rep.r2 << "\n"
                << "mean_length " << stringlen << "\n";
                
                write_string << "\nL  P(L)_fit\n";
                double slen=0;
                for (int i=0; i<(int)n_distr.size(); ++i) {
                    slen=n_distr.at(i).at(0);
                    write_string
                    << slen << " "
                    << calc_y_given_x(c,slen,model_stringlen_distr)
                    << "\n";
                } write_string << "\n\n";
            }
        } /* End is_backwardExtrp=true */
    } /* End is_stringmb=false */
    
    
    // Write time and <L(time)> to file
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/strings_time_stringlen_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    write_stringlen
    << frame
    << " "
    << frame_time
    << " "
    << stringlen
    << "\n";
    //--------------------------------------------------------------------------
    write_string.close();
    write_stringlen.close();
}





void TheoryTest::find_time_stringsRg(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const int frame)
{
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/"+cusstrfolder+"/Rg_strings_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("_frame"+to_string((long long int)frame));//NOTE!
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringsRgDistr_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_stringsRg(o.c_str(),ofstream::app);//'append'
    
    if (readFile.is_open())
    {
        double dubData=0;
        vector<int> vec_intData;
        vector<string> body,frac,rg;
        vector<double> mean;
        string lineContent;
        string strData;
        vector<string> keyWords;
        
        /** 1st: amdat_version **/
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in TheoryTest::find_time_stringsRg():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        /** 2nd: bodies **/
        getline(readFile,lineContent);
        istringstream iss(lineContent);
        iss >> strData;
        keyWords.push_back(strData);//0:bodies
        while (iss>>strData) body.push_back(strData);
        /** 3rd: fraction **/
        getline(readFile,lineContent);iss.clear();iss.str(lineContent);
        iss >> strData;
        keyWords.push_back(strData);//1:fraction
        while (iss>>strData) frac.push_back(strData);
        /** 4th: Rg **/
        getline(readFile,lineContent);iss.clear();iss.str(lineContent);
        iss >> strData;
        keyWords.push_back(strData);//2:rg
        while (iss>>strData) rg.push_back(strData);
        
        /** 5th: blank **/
        getline(readFile,lineContent);
        /** 6th: overall mean **/
        getline(readFile,lineContent);iss.clear();iss.str(lineContent);
        iss >> strData;
        iss >> strData;
        iss >> dubData;
        mean.push_back(dubData);
        
        stringsRg=dubData;
        
        if (body.size()!=frac.size()) {
            cout
            << "body.size() = "<<body.size()<<"\n"
            << "frac.size() = "<<frac.size()<<"\n\n";
            cout
            << "in TheoryTest::find_time_stringsRg():\n"
            << "body.size()!=frac.size()\n"; exit(EXIT_FAILURE);
        }
        if (body.size()!=rg.size()) {
            cout
            << "body.size() = "<<body.size()<<"\n"
            << "rg.size()   = "<<rg.size()  <<"\n\n";
            cout
            << "in TheoryTest::find_time_stringsRg():\n"
            << "body.size()!=rg.size()\n"; exit(EXIT_FAILURE);
        }
        
        frame_time=get_frame_time(sysVar,n_trl,n_sys,Temp_d,frame);
        
        vector<vector<double>> strings_rg;
        for (int i=0;i<(int)body.size();++i)
        {
            double bodyval=0,fracval=0,rgval=0;
            if (rg[i]=="nan") {
                continue;
            } else {
                bodyval=atof(body[i].c_str());
                fracval=atof(frac[i].c_str());
                rgval  =atof(rg[i].c_str());
                strings_rg.push_back({bodyval,fracval,rgval});
            }
        }
        if (strings_rg.size()==0) {
            cout
            << "in TheoryTest::find_time_stringsRg():\n"
            << "strings_rg.size()==0\n"; exit(EXIT_FAILURE);
        }
        
        write_stringsRg
        << "frame_"<<frame<<"; "<< "time_"<<frame_time<<"\n"
        << keyWords[0] << " "
        << keyWords[1] << " "
        << keyWords[2] << "\n";
        if (strings_rg.size()>0) {
            for (int i=0;i<strings_rg.size();++i) {
                write_stringsRg
                << strings_rg[i][0] << " "  //n
                << strings_rg[i][1] << " "  //frac
                << strings_rg[i][2] << "\n";//rg
            }
        } write_stringsRg << "Overall_Mean "<<stringsRg<<"\n\n";
        
    } else {
        cout
        << "in TheoryTest::find_time_stringsRg():\n"
        << i << " cannot open.\n"; exit(EXIT_FAILURE);
    }
    
    //--------------------------------------------------------------------------
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/strings_time_stringsRg_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_time_stringsRg(o.c_str(),ofstream::app);//'append'
    write_time_stringsRg
    << frame
    << " "
    << frame_time
    << " "
    << stringsRg
    << "\n";
    //--------------------------------------------------------------------------
    write_stringsRg.close();
    write_time_stringsRg.close();
}





void TheoryTest::find_time_fast_frac(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const int frame)
{
    double fast_frac=0;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/"+cusstrfolder+"/fast_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("_frame"+to_string((long long int)frame));//NOTE
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        double dubData=0;
        string strData,lineContent;
        getline(readFile,lineContent);//Average_trajectories: n_traj_avg
        istringstream iss(lineContent);
        iss >> strData;
        iss >> dubData;
        find_n_trajs(sysVar,n_trl,n_sys,Temp_d);
        if (n_trajs>0) {
            fast_frac=dubData/double(n_trajs);
        } else {
            cout
            << "in TheoryTest::find_time_fast_frac():\n"
            << "n_trajs="<<n_trajs<<"; "
            << "it should be eqaul to number of trajectories in the system.\n";
            exit(EXIT_FAILURE);
        } readFile.close();
    } else {
        cout
        << "in TheoryTest::find_time_fast_frac():\n"
        << i << " cannot open.\n";
    }
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fastFrac_time_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream writeFile(o.c_str(),ofstream::app); // 'append'
    writeFile
    << frame
    << " "
    << frame_time
    << " "
    << fast_frac
    << "\n";
    writeFile.close();
}





void TheoryTest::find_n_trajs(const StructureClass& sysVar,
                              const int n_trl,
                              const int n_sys,
                              const double Temp_d)
{
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/screen/");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".amdat.screen");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        string strData,lineContent;
        getline(readFile,lineContent);//amdat_version
        if (extract_amdat_svn(lineContent)!=amdat_svn) {
            if (is_check_amdat_version) {
                cout
                << "in TheoryTest::find_n_trajs():\n"
                << "AMDAT version from file ("+extract_amdat_svn(lineContent)+") != "
                << "AMDAT version in use ("+amdat_svn+")\n"; exit(EXIT_FAILURE);
            }
        }
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData;
            if (strData=="Reading") {
                iss >> strData;
                if (strData=="a") {
                    iss >> n_total_frames;
                    while (iss>>strData) {
                        if (strData=="of") {
                            iss >> n_trajs;
                            break;
                        }
                    }
                }
            }
        }
    } else {
        cout
        << "in TheoryTest::find_n_trajs():\n"
        << i << " cannot open.\n";
    }
}





void TheoryTest::find_mean_fast_frac(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d,
                                     const double frametime,
                                     const int frame)
{
    double fast_frac=0;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/"+cusstrfolder+"/fast_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        double dubData=0;
        string strData,lineContent;
        getline(readFile,lineContent);//Average_trajectories: n_traj_avg
        istringstream iss(lineContent);
        iss >> strData;
        iss >> dubData;
        find_n_trajs(sysVar,n_trl,n_sys,Temp_d);
        if (n_trajs>0) {
            fast_frac=dubData/double(n_trajs);
        } else {
            cout
            << "in TheoryTest::find_mean_fast_frac():\n"
            << "n_trajs="<<n_trajs<<"; "
            << "it should be eqaul to number of trajectories in the system.\n";
            exit(EXIT_FAILURE);
        } readFile.close();
    } else {
        cout
        << "in TheoryTest::find_mean_fast_frac():\n"
        << i << " cannot open.\n";
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fastFrac_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writeFile(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T fast_frac frame_index time
    streamsize ss=writeFile.precision();
    
    writeFile
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    writeFile.precision(ss);
    writeFile << resetiosflags(ios::fixed|ios::showpoint);
    
    writeFile
    << frame
    << " "
    << frametime
    << "\n";
    writeFile.close();
}





void TheoryTest::fit_stringlen(const StructureClass& sysVar,
                               const int n_sys,
                               const std::string& model_d)
{
    vector<vector<double>> stringlen_fit;
    
    set_fitParams(model_d);
    read_all_equ_stringlen(sysVar,n_sys); // NOTE: read equilibrium strings
    
    double T_fit=0,str_fit=0;
    for (int i=0; i<stringlen_sortdecreasing.size(); ++i) {
        // NOTE:
        // use normalized T for fitting stringlen vs T
        T_fit  =stringlen_sortdecreasing[i][0]*pow(convert,-1);
        str_fit=stringlen_sortdecreasing[i][1];
        stringlen_fit.push_back({T_fit,str_fit});
    }
    fitdata_processing_lsfit(stringlen_fit);
    alglib_lsfit(model_d);
    
    cout
    << "fit_stringlen ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    //cout << o << "\n";
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model_d);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_log(write_tau_inf,model_d,stringlen_fit);
    write_tau_inf.close();
    
    vector<vector<double>> vec2d;
    for (int i=0; i<(int)stringlen_sortdecreasing.size(); ++i) {
        T_actual=stringlen_sortdecreasing[i][0];
        T_fit   =stringlen_sortdecreasing[i][0]*pow(convert,-1);
        str_fit =log10(calc_y_given_x(c,T_fit,model_d));
        vec2d.push_back({T_actual,str_fit}); // NOTE: use {T_actual,str_fit}
    } stringlen_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_stringlen(const StructureClass& sysVar,
                               const int n_sys,
                               const double TA_d,
                               const std::string& model_d)
{
    vector<vector<double>> stringlen_fit;
    
    set_fitParams(model_d);
    read_all_equ_stringlen(sysVar,n_sys,TA_d); // NOTE: fit T<=TA
    
    double T_fit=0,str_fit=0,T_actual=0;
    for (int i=0; i<stringlen_sortdecreasing.size(); ++i) {
        // NOTE:
        // use normalized T for fitting stringlen vs T
        T_fit  =stringlen_sortdecreasing[i][0]*pow(convert,-1);
        str_fit=stringlen_sortdecreasing[i][1];
        stringlen_fit.push_back({T_fit,str_fit});
    }
    fitdata_processing_lsfit(stringlen_fit);
    alglib_lsfit(model_d);
    
    cout
    << "fit_stringlen ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tau_inf(o.c_str());
    
    write_fitModel(write_tau_inf,model_d);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_log(write_tau_inf,model_d,stringlen_fit);
    write_tau_inf.close();
    
    vector<vector<double>> vec2d;
    for (int i=0; i<(int)stringlen_sortdecreasing.size(); ++i) {
        T_actual=stringlen_sortdecreasing[i][0];
        T_fit   =stringlen_sortdecreasing[i][0]*pow(convert,-1);
        str_fit =log10(calc_y_given_x(c,T_fit,model_d));
        vec2d.push_back({T_actual,str_fit}); // NOTE: use {T_actual,str_fit}
    } stringlen_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_DWF(const StructureClass& sysVar,
                         const int n_sys,
                         const std::string& model_d)
{
    set_fitParams(model_d);
    read_all_equ_DWF(sysVar,n_sys);
    fitdata_processing_lsfit(dwf_invT_sortdecreasing);//NOTE: "dwf_invT"
    alglib_lsfit(model_d);
    
    cout
    << "fit_DWF ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fit_DWF(o.c_str());
    write_fitModel(write_fit_DWF,model_d);
    write_fitarrays(write_fit_DWF);
    write_stopcond(write_fit_DWF);
    write_fitinfo(write_fit_DWF);
    if (model_d=="secondpoly") {
        write_errcurve_actual(write_fit_DWF,model_d,dwf_invT_sortdecreasing);
    } else {
        write_errcurve_log(write_fit_DWF,model_d,dwf_invT_sortdecreasing);
    } write_fit_DWF.close();
    
    vector<vector<double>> vec2d;
    double temp=0,invT=0,dwfr=0;
    for (int i=0; i<(int)dwf_invT_sortdecreasing.size(); ++i) {
        temp=dwf_sortdecreasing[i][0]; // NOTE: is "dwf" .NOT. dwf_invT
        invT=dwf_invT_sortdecreasing[i][0];//NOTE
        if (model_d=="secondpoly") {
            dwfr=calc_y_given_x(c,invT,model_d);
        } else {
            dwfr=log10(calc_y_given_x(c,invT,model_d));
        } vec2d.push_back({temp,dwfr});
    } dwf_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_DWF(const StructureClass& sysVar,
                         const int n_sys,
                         const double TA_d,
                         const std::string& model_d)
{
    set_fitParams(model_d);
    read_all_equ_DWF(sysVar,n_sys,TA_d);//NOTE: fit T<=TA
    fitdata_processing_lsfit(dwf_invT_sortdecreasing);//NOTE: use "dwf_invT"
    alglib_lsfit(model_d);
    
    cout
    << "fit_DWF ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fit_DWF(o.c_str());
    write_fitModel(write_fit_DWF,model_d);
    write_fitarrays(write_fit_DWF);
    write_stopcond(write_fit_DWF);
    write_fitinfo(write_fit_DWF);
    if (model_d=="secondpoly") {
        write_errcurve_actual(write_fit_DWF,model_d,dwf_invT_sortdecreasing);
    } else {
        write_errcurve_log(write_fit_DWF,model_d,dwf_invT_sortdecreasing);
    } write_fit_DWF.close();
    
    vector<vector<double>> vec2d;
    double temp=0,invT=0,dwfr=0;
    for (int i=0; i<(int)dwf_invT_sortdecreasing.size(); ++i) {
        temp=dwf_sortdecreasing[i][0]; // NOTE: is "dwf" .NOT. dwf_invT
        invT=dwf_invT_sortdecreasing[i][0];//NOTE
        if (model_d=="secondpoly") {
            dwfr=calc_y_given_x(c,invT,model_d);
        } else {
            dwfr=log10(calc_y_given_x(c,invT,model_d));
        } vec2d.push_back({temp,dwfr});
    } dwf_sortdecreasing_fit=vec2d;
    
}





void TheoryTest::fit_taueq(const StructureClass& sysVar,
                           const int n_sys,
                           const std::string& model_d)
{
    set_fitParams(model_d);
    read_all_taueq_data(sysVar,n_sys);//all
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model_d);
    
    cout
    << "fit_taueq ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fit_taueq(o.c_str());
    write_fitModel(write_fit_taueq,model_d);
    write_fitarrays(write_fit_taueq);
    write_stopcond(write_fit_taueq);
    write_fitinfo(write_fit_taueq);
    write_errcurve_log(write_fit_taueq,model_d,taueq_sortdecreasing);
    write_fit_taueq.close();
    
    vector<vector<double>> vec2d;
    double temp=0,taue=0;
    for (int i=0; i<(int)taueq_sortdecreasing.size(); ++i) {
        temp=taueq_sortdecreasing[i][0];
        taue=log10(calc_y_given_x(c,temp,model_d));
        vec2d.push_back({temp,taue});
    } taueq_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_taueq(const StructureClass& sysVar,
                           const int n_sys,
                           const double TA_d,
                           const std::string& model_d)
{
    set_fitParams(model_d);
    read_all_taueq_data(sysVar,n_sys,TA_d); // NOTE: fit T<=TA
    fitdata_processing_lsfit(taueq_sortdecreasing);
    alglib_lsfit(model_d);
    
    cout
    << "fit_taueq ("<<model_d<<")\n"
    << "R^2="<< rep.r2 << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_fit_"+model_d+"_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fit_taueq(o.c_str());
    write_fitModel(write_fit_taueq,model_d);
    write_fitarrays(write_fit_taueq);
    write_stopcond(write_fit_taueq);
    write_fitinfo(write_fit_taueq);
    write_errcurve_log(write_fit_taueq,model_d,taueq_sortdecreasing);
    write_fit_taueq.close();
    
    vector<vector<double>> vec2d;
    double temp=0,taue=0;
    for (int i=0; i<(int)taueq_sortdecreasing.size(); ++i) {
        temp=taueq_sortdecreasing[i][0];
        taue=log10(calc_y_given_x(c,temp,model_d));
        vec2d.push_back({temp,taue});
    } taueq_sortdecreasing_fit=vec2d;
}





void TheoryTest::fit_all_strings_rgvsm(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d,
                                       const int frame)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringsRg_Rgvsn_powerlawfit_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".dat");
    ofstream write_stringsRgdata(o.c_str(),ofstream::app);//'append'
    
    //1. read strings rg vs n (in log-log form)
    read_frame_strings_rgvsn(sysVar,n_trl,n_sys,Temp_d,frame);
    
    double df=0,r2=0;
    
    //if (strings_rgvsnloglog_sortincreasing.size()==0) return;
    if (strings_rgvsnloglog_sortincreasing.size()>0)
    {
        //2. get df from the slope of a linear fit
        set_fitParams("linear");
        fitdata_processing_lsfit(strings_rgvsnloglog_sortincreasing);
        alglib_lsfit("linear");
        
        df=c[0];
        r2=rep.r2;
        
        string msg="T"+to_string(std::round(Temp_d))+"_frame"+to_string(frame);
        write_message(write_stringsRgdata,msg);
        write_fitModel(write_stringsRgdata,"linear");
        write_fitarrays(write_stringsRgdata);
        //to show original Rg vs n values
        write_message(write_stringsRgdata,"Show original Rg vs n");
        fitdata_processing_lsfit(strings_rgvsn_sortincreasing);//set new x & y data
        write_fitarrays(write_stringsRgdata);
        fitdata_processing_lsfit(strings_rgvsnloglog_sortincreasing);//resume x & y data
        //
        write_stopcond(write_stringsRgdata);
        write_fitinfo(write_stringsRgdata);
        write_yvsxby2dvecdub(write_stringsRgdata,strings_rgvsnloglog_sortincreasing);
    }
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringsRg_Rgvsn_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringsRg(o.c_str(),ofstream::app);//'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1);//NOTE: convert!
    
    // T frame time df r2
    streamsize ss=write_stringsRg.precision();
    
    write_stringsRg
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringsRg.precision(ss);
    write_stringsRg << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringsRg
    << frame
    << " "
    << get_frame_time(sysVar,n_trl,n_sys,Temp_d,frame)
    << " "
    << df
    << " "
    << r2
    << "\n";
}





void TheoryTest::fit_stringmodel(const StructureClass& sysVar,
                                 const int n_sys,
                                 const double TA_d,
                                 const string& model_d)
{
    /** fit high Ts (T>TA) **/
    read_all_taueq_data(sysVar,n_sys,TA_d,false);    //NOTE: for T>TA
    read_all_equ_stringlen(sysVar,n_sys,TA_d,false); //NOTE: for T>TA
    
    /** fit low Ts (T<=TA) **/
    //read_all_taueq_data(sysVar,n_sys,TA_d);        //NOTE: for T<=TA
    //read_all_equ_stringlen(sysVar,n_sys,TA_d);     //NOTE: for T<=TA
    
    /** fit all Ts **/
    //read_all_taueq_data(sysVar,n_sys);
    //read_all_equ_stringlen(sysVar,n_sys);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    //cout << n_taueq << "\n";
    //cout << strSize << "\n";
    if (n_taueq==strSize) {
        double Teqtau=0,TeqL=0,LT=0,taueq=0;
        taueq2d_sortdecreasing.clear();
        for (indexi=0; indexi<n_taueq; ++indexi) {
            Teqtau = taueq_sortdecreasing.at(indexi).at(0);
            taueq  = taueq_sortdecreasing.at(indexi).at(1);
            TeqL   = stringlen_sortdecreasing.at(indexi).at(0);
            LT     = stringlen_sortdecreasing.at(indexi).at(1);
            if (Teqtau==TeqL) {
                if (model_d=="string2param") {
                    taueq2d_sortdecreasing.push_back({Teqtau,LT,LA,tauA,TA,taueq});
                    //input: {T,L(T),LA,tauA,TA,taueq}
                } else if (model_d=="string3param") {
                    taueq2d_sortdecreasing.push_back({Teqtau,LT,LA,taueq});
                    //input: {T,L(T),LA,taueq}
                } else if (model_d=="transitState") {
                    taueq2d_sortdecreasing.push_back({Teqtau,tau0_vib,taueq});
                    //input: {T,tau0_vib,taueq}
                }
            }
        }
    } else {
        cout
        << "in TheoryTest::fit_stringmodel()" << "\n"
        << "n_taueq!=strSize (1st part)" << "\n";
        exit(EXIT_FAILURE);
    }
    fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    set_fitParams(model_d);
    alglib_lsfit(model_d);
    
    if (model_d=="string2param") {
        delHa=c[0];
        delSa=c[1];
    } else if (model_d=="string3param") {
        tau0_vib=c[0];
        delHa=c[1];
        delSa=c[2];
    } else if (model_d=="transitState") {
        delHa=c[0];
        delSa=c[1];
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_AGstringmodel_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_stopcond(writefile);
    write_fitinfo(writefile);
    write_errcurve2d_log(writefile,model_d,taueq2d_sortdecreasing);
    writefile.close();
}





void TheoryTest::fit_GLMmodel(const StructureClass& sysVar,
                              const int n_sys,
                              const bool is_cutbyTA,
                              const string& model_d)
{
    string fit="pSpline";//COOP,pSpline
    
    alglib::real_1d_array c_taufit,c_dwffit;
    alglib::spline1dinterpolant tauip,dwfip;
    
    //NOTE: use entire T range to determine Tref
    if (fit=="COOP") {
        fit_taueq(sysVar,n_sys,"COOP");      c_taufit=c;
        fit_DWF(sysVar,n_sys,"COOP_DWF");    c_dwffit=c;
    } else if (fit=="pSpline") {
        read_all_taueq_data(sysVar,n_sys);//all
        build_penalizedspline_interpolant(taueqinvT_sortdecreasing,1,0);//NOTE:(y,x)
        tauip=interpolant;
        read_all_equ_DWF(sysVar,n_sys);//all
        build_penalizedspline_interpolant(dwf_invT_sortdecreasing,0,1);
        dwfip=interpolant;
    }
    
    if (fit=="COOP") {
        Tref =calc_x_given_y(c_taufit,pow(10,logtauref),"COOP").at(0);
        u2ref=log10(calc_y_given_x(c_dwffit,convert/Tref,"COOP_DWF"));
    } else if (fit=="pSpline") {
        Tref =convert/spline1dcalc(tauip,logtauref);
        u2ref=spline1dcalc(dwfip,convert/Tref);
    } else {
        cout
        << "in TheoryTest::fit_GLMmodel(): "
        << "fitting model ("<<fit<<") not found.\n";
        exit(EXIT_FAILURE);
    }
    
    
    //fit GLM form
    //--------------------------------------------------------------------------
    bool is_use_translational_time=false;
    
    if (is_cutbyTA) {
        read_all_equ_DWF(sysVar,n_sys,TA);
        if (is_use_translational_time) read_all_r2teq_data(sysVar,n_sys,TA);
        else read_all_taueq_data(sysVar,n_sys,TA);
    } else {
        read_all_equ_DWF(sysVar,n_sys);
        if (is_use_translational_time) read_all_r2teq_data(sysVar,n_sys);
        else read_all_taueq_data(sysVar,n_sys);
    }
    int n_taueq=(int)taueq_sortdecreasing.size();
    int n_dwfac=(int)dwf_sortdecreasing.size();
    if (n_taueq!=n_dwfac)
    {
        cout
        << "in TheoryTest::fit_GLMmodel():\n"
        << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
        make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
        n_taueq=(int)taueq_sortdecreasing.size();
        n_dwfac=(int)dwf_sortdecreasing.size();
        cout
        << "in TheoryTest::fit_GLMmodel():\n"
        << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
        system_wait(1);
    }
    
    double Teqtau=0,Teqdwf=0,dwfac=0,taueq=0;
    taueq2d_sortdecreasing.clear();
    standard2dcontainer.clear();
    for (indexi=0; indexi<n_taueq; ++indexi)
    {
        Teqtau = taueq_sortdecreasing.at(indexi).at(0);
        taueq  = taueq_sortdecreasing.at(indexi).at(1);
        Teqdwf = dwf_sortdecreasing.at(indexi).at(0);
        dwfac  = dwf_sortdecreasing.at(indexi).at(1);
        if (model_d=="univGLM")
        {
            //v=({tauA,TA,T,tau}
            standard2dcontainer.push_back({tauA,TA,Teqtau,taueq});
            
            /** NOTE: final param in input must be taueq **/
            //input: {tauA,TA,T,tau}
            taueq2d_sortdecreasing.push_back({tauA,TA,Teqtau,taueq});
        }
        else if (model_d=="univGT")
        {
            //v={T/TA,u2/u2A,tauA,tau0_Arrhenius}
            standard2dcontainer.push_back({Teqtau/TA,dwfac/u2A,tauA,tau0_Arrhenius});
            
            /** NOTE: final param in input must be taueq **/
            //input: {T/TA,u2/u2A,tauA,tau0_Arrhenius,T,tau}
            taueq2d_sortdecreasing.push_back
            ({Teqtau/TA,dwfac/u2A,tauA,tau0_Arrhenius,Teqtau,taueq});
        }
        else if (model_d=="GLM"||model_d=="GLM1"||model_d=="GLM1_hallwolynes")//GLM,GLM1
        {
            //v={logtauref,u2ref,u2,T,tau,tauA,u2A,EA,TA}
            standard2dcontainer.push_back({logtauref,u2ref,dwfac,Teqtau,taueq,tauA,u2A,Ea,TA});
            
            /** NOTE: final param in input must be taueq **/
            //input: {tauA,u2A,u2,Ea,TA,T,tau}
            taueq2d_sortdecreasing.push_back({tauA,u2A,dwfac,Ea,TA,Teqtau,taueq});
        }
        else if (model_d=="univILP")
        {
            taueq2d_sortdecreasing.push_back({dwfac/u2A,taueq});
        }
        else
        {
            cout
            << "in TheoryTest::fit_GLMmodel():\n"
            << "model ("<<model_d<<") not found.\n";
            exit(EXIT_FAILURE);
        }
    }
    
    if (model_d=="GLM"||model_d=="GLM1"||model_d=="GLM1_hallwolynes") {
        fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    } else if (model_d=="univGLM") {
        fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    } else if (model_d=="univGT") {
        fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    } else if (model_d=="univILP") {
        fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    }
    
    //NOTE: use F-mode for GLM fitting
    //--------------------------------------------------------------------------
    bool is_use_FG_org=is_use_FG; is_use_FG=false;
    
    set_fitParams(model_d);
    alglib_lsfit(model_d);
    
    if (model_d=="GLM") {
        tau0_GLM = c[0];
        u02      = c[1];
        alpha    = c[2];
    } else if (model_d=="GLM1"||model_d=="GLM1_hallwolynes") {
        alpha    = c[0];
    } else if (model_d=="univGLM") {
        alpha    = c[0];
    } else if (model_d=="univGT") {
        dfit     = c[0];
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+model_d+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_stopcond(writefile);
    write_fitinfo(writefile);
    if (model_d=="GLM") {
        write_errcurve_log(writefile,model_d,dwf_tau_data);
    } else if (model_d=="GLM1"||model_d=="GLM1_hallwolynes") {
        write_errcurve2d_log(writefile,model_d,standard2dcontainer);
    } else if (model_d=="univGLM") {
        write_errcurve2d_log(writefile,model_d,standard2dcontainer);
    } else if (model_d=="univGT") {
        write_errcurve2d_log(writefile,model_d,standard2dcontainer);
    } writefile.close();
    is_use_FG=is_use_FG_org;
    //--------------------------------------------------------------------------
}





void TheoryTest::write_GLMref(const StructureClass& sysVar,
                              const int n_sys,
                              const string& model_d)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_GLMref_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    
    vector<double> ref={logtauref,u2ref,Tref};
    
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_stopcond(writefile);
    write_Yg_fragility(writefile,c,ref,model_d);
    write_fitinfo(writefile);
    write_tauFitfit_vs_reduceddwf(writefile,c,standard2dcontainer,ref,0.0,model_d);
    writefile.close();
}





void TheoryTest::fit_HWmodel(const StructureClass& sysVar,
                             const int n_sys,
                             const bool is_cutbyTA,
                             const string& model_d)
{
    string fit="pSpline";//COOP,pSpline
    
    alglib::real_1d_array c_taufit,c_dwffit;
    alglib::spline1dinterpolant tauip,dwfip;
    
    //NOTE: use entire T range to determine Tref
    fit_taueq(sysVar,n_sys,"COOP");      c_taufit=c;
    fit_DWF(sysVar,n_sys,"COOP_DWF");    c_dwffit=c;
    build_penalizedspline_interpolant(taueqinvT_sortdecreasing,1,0);//NOTE:(y,x)
    tauip=interpolant;
    build_penalizedspline_interpolant(dwf_invT_sortdecreasing,0,1);
    dwfip=interpolant;
    
    if (fit=="COOP") {
        Tref =calc_x_given_y(c_taufit,pow(10,logtauref),"COOP").at(0);
        u2ref=log10(calc_y_given_x(c_dwffit,convert/Tref,"COOP_DWF"));
    } else if (fit=="pSpline") {
        Tref =convert/spline1dcalc(tauip,logtauref);
        u2ref=spline1dcalc(dwfip,convert/Tref);
    } else {
        cout
        << "in TheoryTest::fit_HWmodel(): "
        << "fitting model ("<<fit<<") not found.\n";
        exit(EXIT_FAILURE);
    }
    
    if (is_cutbyTA) {
        read_all_taueq_data(sysVar,n_sys,TA);
        read_all_equ_DWF(sysVar,n_sys,TA);
    } else {
        read_all_taueq_data(sysVar,n_sys);
        read_all_equ_DWF(sysVar,n_sys);
    }
    int n_taueq=(int)taueq_sortdecreasing.size();
    int n_dwfac=(int)dwf_sortdecreasing.size();
    if (n_taueq!=n_dwfac)
    {
        cout
        << "in TheoryTest::fit_HWmodel():\n"
        << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
        make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
        n_taueq=(int)taueq_sortdecreasing.size();
        n_dwfac=(int)dwf_sortdecreasing.size();
        cout
        << "in TheoryTest::fit_HWmodel():\n"
        << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
        system_wait(1);
    }
    
    double Teqtau=0,Teqdwf=0,dwfac=0,taueq=0;
    taueq2d_sortdecreasing.clear();
    standard2dcontainer.clear();
    for (indexi=0; indexi<n_taueq; ++indexi)
    {
        Teqtau = taueq_sortdecreasing.at(indexi).at(0);
        taueq  = taueq_sortdecreasing.at(indexi).at(1);
        Teqdwf = dwf_sortdecreasing.at(indexi).at(0);
        dwfac  = dwf_sortdecreasing.at(indexi).at(1);
        if (model_d=="hallwolynes") {
            standard2dcontainer.push_back({dwfac,Teqtau,taueq});
            //v={u2,T,tau}
            taueq2d_sortdecreasing.push_back({dwfac,Teqtau,taueq});
            //input: {u2,T,tau}
        } else if (model_d=="hallwolynes1param") {
            standard2dcontainer.push_back({dwfac,u2A,tauA,Teqtau,taueq});
            //v={u2,u2A,tauA,T,tau}
            taueq2d_sortdecreasing.push_back({dwfac,u2A,tauA,Teqtau,taueq});
            //input: {u2,u2A,tauA,T,tau}
        } else {
            cout
            << "in TheoryTest::fit_HWmodel():\n"
            << "model ("<<model_d<<") not found.\n";
            exit(EXIT_FAILURE);
        }
    }
    
    //NOTE: use F-mode for Leporini fitting
    //--------------------------------------------------------------------------
    bool is_use_FG_org=is_use_FG; is_use_FG=false;
    
    fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    set_fitParams(model_d);
    alglib_lsfit(model_d);
    
    vector<double> ref={logtauref,u2ref,Tref};
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+model_d+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_stopcond(writefile);
    write_fitinfo(writefile);
    write_errcurve2d_log(writefile,model_d,standard2dcontainer);
    writefile.close();
    is_use_FG=is_use_FG_org;
    //--------------------------------------------------------------------------
}





void TheoryTest::fit_Leporinimodel(const StructureClass& sysVar,
                                   const int n_sys,
                                   const bool is_cutbyTA,
                                   const string& model_d)
{
    string fit="pSpline";//COOP,pSpline
    
    alglib::real_1d_array c_taufit,c_dwffit;
    alglib::spline1dinterpolant tauip,dwfip;
    
    //NOTE: use entire T range to determine Tref
    fit_taueq(sysVar,n_sys,"COOP");      c_taufit=c;
    fit_DWF(sysVar,n_sys,"COOP_DWF");    c_dwffit=c;
    build_penalizedspline_interpolant(taueqinvT_sortdecreasing,1,0);//NOTE:(y,x)
    tauip=interpolant;
    build_penalizedspline_interpolant(dwf_invT_sortdecreasing,0,1);
    dwfip=interpolant;
    
    if (fit=="COOP") {
        Tref =calc_x_given_y(c_taufit,pow(10,logtauref),"COOP").at(0);
        u2ref=log10(calc_y_given_x(c_dwffit,convert/Tref,"COOP_DWF"));
    } else if (fit=="pSpline") {
        Tref =convert/spline1dcalc(tauip,logtauref);
        u2ref=spline1dcalc(dwfip,convert/Tref);
    } else {
        cout
        << "in TheoryTest::fit_Leporinimodel(): "
        << "fitting model ("<<fit<<") not found.\n";
        exit(EXIT_FAILURE);
    }
    
    if (is_cutbyTA) {
        read_all_taueq_data(sysVar,n_sys,TA);
        read_all_equ_DWF(sysVar,n_sys,TA);
    } else {
        read_all_taueq_data(sysVar,n_sys);
        read_all_equ_DWF(sysVar,n_sys);
    }
    int n_taueq=(int)taueq_sortdecreasing.size();
    int n_dwfac=(int)dwf_sortdecreasing.size();
    if (n_taueq!=n_dwfac)
    {
        cout
        << "in TheoryTest::fit_Leporinimodel():\n"
        << "n_taueq("<<n_taueq<<")!=n_dwfac("<<n_dwfac<<")\n";
        make1dconsistent(taueq_sortdecreasing,dwf_sortdecreasing);
        n_taueq=(int)taueq_sortdecreasing.size();
        n_dwfac=(int)dwf_sortdecreasing.size();
        cout
        << "in TheoryTest::fit_Leporinimodel():\n"
        << "After make 1d consistent: n_taueq("<<n_taueq<<")=n_dwfac("<<n_dwfac<<")\n";
        system_wait(1);
    }
    
    double Teqtau=0,Teqdwf=0,dwfac=0,taueq=0;
    taueq2d_sortdecreasing.clear();
    standard2dcontainer.clear();
    for (indexi=0; indexi<n_taueq; ++indexi)
    {
        Teqtau = taueq_sortdecreasing.at(indexi).at(0);
        taueq  = taueq_sortdecreasing.at(indexi).at(1);
        Teqdwf = dwf_sortdecreasing.at(indexi).at(0);
        dwfac  = dwf_sortdecreasing.at(indexi).at(1);
        standard2dcontainer.push_back({logtauref,u2ref,dwfac,Teqtau,taueq});
        //v={logtauref,u2ref,u2,T,tau}
        taueq2d_sortdecreasing.push_back({logtauref,u2ref,dwfac,Teqtau,taueq});
        //input: {logtauref,u2ref,u2,T,tau}
    }
    
    //NOTE: use F-mode for Leporini fitting
    //--------------------------------------------------------------------------
    bool is_use_FG_org=is_use_FG; is_use_FG=false;
    
    fitdata_processing_lsfit_2dx(taueq2d_sortdecreasing);
    set_fitParams(model_d);
    alglib_lsfit(model_d);
    
    vector<double> ref={logtauref,u2ref,Tref};
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+model_d+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream writefile(o.c_str());
    write_fitModel(writefile,model_d);
    write_fitarrays(writefile);
    write_stopcond(writefile);
    vector<double> glass_extrp=
    write_Yg_fragility(writefile,c,ref,model_d);
    Yg  = glass_extrp.at(0);
    mg  = glass_extrp.at(1);
    u2g = glass_extrp.at(2);
    write_fitinfo(writefile);
    write_errcurve2d_log(writefile,model_d,standard2dcontainer);
    write_tauFitfit_vs_reduceddwf
    (writefile,c,standard2dcontainer,ref,0.0,model_d);
    writefile.close();
    is_use_FG=is_use_FG_org;
    //--------------------------------------------------------------------------
}





void TheoryTest::write_peak_stringlen(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    // find stringlen at the peak of stringlen as a function of time
    
    vector<vector<double>> stringlen;
    vector<vector<double>> stringlen_logtime;
    vector<vector<double>> stringlen_smoothed;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/strings_time_stringlen_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readstringlendata(i.c_str());
    if (readstringlendata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double time=0,stringlenData=0;
        while (getline(readstringlendata,lineContent)) {
            istringstream iss(lineContent);//(frame_index,time,stringlen)
            iss >> frame_index;
            iss >> time;
            iss >> stringlenData;
            stringlen.push_back({(double)frame_index,time,stringlenData});
            stringlen_logtime.push_back({(double)frame_index,log10(time),stringlenData});
        } readstringlendata.close();
    } else {
        cout
        << "in TheoryTest::write_peak_stringlen():\n"
        << i << " cannot open.\n";
    }
    
    /** noise reduction on L(time) **/
    //--------------------------------------------------------------------------
    /** NOTE: append smoothed data to original file **/
    ofstream write_string(i.c_str(),ofstream::app);//'append'
    int size=(int)stringlen_logtime.size();
    build_penalizedspline_interpolant(stringlen_logtime,1,2,size,rho_Lt);
    double logtime=0,f=0;
    write_string << "\nframe  time  stringlen_smoothed\n";
    for (indexi=0; indexi<(int)stringlen.size(); ++indexi) {
        logtime=stringlen_logtime.at(indexi).at(1);//log10(time)
        f=spline1dcalc(interpolant,logtime);
        write_string
        << stringlen.at(indexi).at(0) << " " //frame
        << stringlen.at(indexi).at(1) << " " //time
        << f << "\n";
        stringlen_smoothed.push_back
        ({stringlen.at(indexi).at(0),stringlen.at(indexi).at(1),f});
    } write_string.close();
    if (is_use_stringLen_smoothed) stringlen=stringlen_smoothed;
    //--------------------------------------------------------------------------
    
    //sort by L in descending order
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing2);
    } catch (int i) {
        cout
        << "in TheoryTest::write_peak_stringlen():\n"
        << "stringlen size = " << i << "\n";
        stringlen.push_back({0,0});
    }
    
    peak_frame=stringlen.at(0).at(0);
    peak_time =stringlen.at(0).at(1);
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringlen(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringleninv(o.c_str(),ofstream::app);//'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1);//NOTE: convert!
    
    // T frame_index time stringlen
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringlen
    << stringlen.at(0).at(0) //frame
    << " "
    << stringlen.at(0).at(1) //time
    << " "
    << stringlen.at(0).at(2) //stringlen
    << "\n";
    
    write_stringleninv
    << convert/T_actual
    << " "
    << stringlen.at(0).at(0)
    << " "
    << stringlen.at(0).at(1)
    << " "
    << stringlen.at(0).at(2)
    << "\n";
    
    write_stringlen.close();
    write_stringleninv.close();
}





void TheoryTest::write_frame_stringsRg(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d,
                                       const int frame)
{
    vector<vector<double>> stringsRg;
    vector<vector<double>> stringsRg_logtime;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/strings_time_stringsRg_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readstringsRgdata(i.c_str());
    if (readstringsRgdata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double time=0,stringsRgData=0;
        while (getline(readstringsRgdata,lineContent)) {
            istringstream iss(lineContent);//(frame_index,time,stringsRg)
            iss >> frame_index;
            iss >> time;
            iss >> stringsRgData;
            stringsRg.push_back({(double)frame_index,time,stringsRgData});
            stringsRg_logtime.push_back({(double)frame_index,log10(time),stringsRgData});
        } readstringsRgdata.close();
    } else {
        cout
        << "in TheoryTest::write_frame_stringsRg():\n"
        << i << " cannot open.\n";
    }
    
    double time=0,stringrg=0;
    for (int i=0;i<(int)stringsRg.size();++i) {
        if (frame==stringsRg[i][0]) {
            time=stringsRg[i][1];
            stringrg=stringsRg[i][2];
            break;
        }
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringsRg_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringsRg(o.c_str(),ofstream::app);//'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringsRg_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringsRginv(o.c_str(),ofstream::app);//'append'
    
    double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T frame_index time stringsRg
    streamsize ss=write_stringsRg.precision();
    
    write_stringsRg
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringsRg.precision(ss);
    write_stringsRg << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringsRg
    << frame
    << " "
    << time
    << " "
    << stringrg
    << "\n";
    
    write_stringsRginv
    << convert/T_actual
    << " "
    << frame
    << " "
    << time
    << " "
    << stringrg
    << "\n";
    
    write_stringsRg.close();
    write_stringsRginv.close();
}





void TheoryTest::write_time_stringlen(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    // find stringlen at tau_alpha
    
    /** find tau_alpha frame **/
    find_taualpha_frame(sysVar,n_trl,n_sys,Temp_d);
    
    vector<vector<double>> stringlen;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/strings_time_stringlen_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readstringlendata(i.c_str());
    if (readstringlendata.is_open()) {
        string lineContent;
        int    frame_index=0;
        double time=0,stringlenData=0;
        while (getline(readstringlendata,lineContent)) {
            istringstream iss(lineContent);//(frame_index,time,stringlen)
            iss >> frame_index;
            iss >> time;
            iss >> stringlenData;
            stringlen.push_back({(double)frame_index,time,stringlenData});
        } readstringlendata.close();
    } else {
        cout
        << "in TheoryTest::write_time_stringlen():\n"
        << i << " cannot open.\n";
    }
    
    /** find stringlen at tau_alpha frame **/
    double stringlen_taualpha=0;
    for (size_t i=0;i<stringlen.size();++i) {
        if (stringlen.at(i).at(0)==taualpha_frame) {
            stringlen_taualpha=stringlen.at(i).at(2);
        }
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringlen(o.c_str(),ofstream::app); // 'append'
    
    o.clear();
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/stringlen_all_inverseT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_stringleninv(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T frame_index time stringlen
    streamsize ss=write_stringlen.precision();
    
    write_stringlen
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_stringlen.precision(ss);
    write_stringlen << resetiosflags(ios::fixed|ios::showpoint);
    
    write_stringlen
    << taualpha_frame
    << " "
    << taualpha_frame_time
    << " "
    << stringlen_taualpha
    << "\n";
    
    write_stringleninv
    << convert/T_actual
    << " "
    << taualpha_frame
    << " "
    << taualpha_frame_time
    << " "
    << stringlen_taualpha
    << "\n";
    
    write_stringlen.close();
    write_stringleninv.close();
}





void TheoryTest::write_peak_fast_frac(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    vector<vector<double>> fast_frac;
    vector<vector<double>> fast_frac_logtime;
    vector<vector<double>> fast_frac_smoothed;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/fastFrac_time_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        string lineContent;
        int    frame=0;
        double time=0,fast=0;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);//(frame_index,time,fast)
            iss >> frame;
            iss >> time;
            iss >> fast;
            fast_frac.push_back({(double)frame,time,fast});
            fast_frac_logtime.push_back({(double)frame,log10(time),fast});
        } readFile.close();
    } else {
        cout
        << "in TheoryTest::write_peak_fast_frac():\n"
        << i << " cannot open.\n";
    }
    
    /** noise reduction on L(time) **/
    //--------------------------------------------------------------------------
    /** NOTE: append smoothed data to original file **/
    ofstream writeFile(i.c_str(),ofstream::app);//'append'
    build_penalizedspline_interpolant(fast_frac_logtime,1,2,50,rho_Lt);
    double logtime=0,f=0;
    writeFile << "\nframe  time  fastFrac_smoothed\n";
    for (indexi=0; indexi<(int)fast_frac.size(); ++indexi) {
        logtime=fast_frac_logtime.at(indexi).at(1);//log10(time)
        f=spline1dcalc(interpolant,logtime);
        writeFile
        << fast_frac.at(indexi).at(0) << " " //frame
        << fast_frac.at(indexi).at(1) << " " //time
        << f << "\n";
        fast_frac_smoothed.push_back
        ({fast_frac.at(indexi).at(0),fast_frac.at(indexi).at(1),f});
    } writeFile.close();
    if (is_use_fast_frac_smoothed) fast_frac=fast_frac_smoothed;
    //--------------------------------------------------------------------------
    
    try {
        if(fast_frac.size()==0) throw 0;
        std::sort(fast_frac.begin(),fast_frac.end(),sortDecreasing2);
    } catch (int i) {
        cout
        << "in TheoryTest::write_peak_fast_frac():\n"
        << "fast_frac size = " << i << "\n";
        fast_frac.push_back({0,0});
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/fastFrac_all_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_fast_frac(o.c_str(),ofstream::app); // 'append'
    
    double T_actual = Temp_d*pow(corrFacforT,-1); // NOTE: convert!
    
    // T frame_index time fast_frac
    streamsize ss=write_fast_frac.precision();
    
    write_fast_frac
    << fixed << setprecision(10) // NOTE: be consistent in temperature format
    << T_actual
    << " ";
    
    write_fast_frac.precision(ss);
    write_fast_frac << resetiosflags(ios::fixed|ios::showpoint);
    
    write_fast_frac
    << fast_frac.at(0).at(0)
    << " "
    << fast_frac.at(0).at(1)
    << " "
    << fast_frac.at(0).at(2)
    << "\n";
    
    write_fast_frac.close();
}





void TheoryTest::read_individual_stringlen(const StructureClass& sysVar,
                                           const int n_trl,
                                           const int n_sys)
{
    stringlen_sortdecreasing.clear();
    stringsdf_sortdecreasing.clear();
    vector<vector<double>> stringlen,stringsdf;
    
    string i;
    
    /** L extracted from beta_kww **/
    //--------------------------------------------------------------------------
    if (is_use_Lbeta)
    {
        i.clear();
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/qtest_bKWW_min_vs_T_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("."+analysispart+".dat");
        ifstream readstringlendata(i.c_str());
        if (readstringlendata.is_open()) {
            string lineContent;
            double Tdata=0,stringlenData=0;
            double qmin=0,bmin=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent);
                iss >> Tdata;
                iss >> qmin;
                iss >> bmin;
                stringlenData=(2*pi())/qmin;
                stringlen.push_back({Tdata,stringlenData});
            } readstringlendata.close();
        } else {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << i << " cannot open.\n";
        }
    }
    /** L from string calculations by AMADT **/
    //--------------------------------------------------------------------------
    else
    {
        i.clear();
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/"+stringtype+"_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readstringlendata(i.c_str());
        if (readstringlendata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringlenData=0,time=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent);
                iss >> Tdata;
                iss >> frame_inedx;
                iss >> time;
                iss >> stringlenData;
                stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
            } readstringlendata.close();
        } else {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << i << " cannot open.\n";
        }
        /** fractal dimension (df) from Rg vs n **/
        //----------------------------------------------------------------------
        if (stringtype=="stringsRg")
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/stringsRg_Rgvsn_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readstringRgdata(i.c_str());
            if (readstringRgdata.is_open()) {
                string lineContent;
                int    frame_inedx=0;
                double Tdata=0,time=0,df=0,r2=0;
                while (getline(readstringRgdata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> frame_inedx;
                    iss >> time;
                    iss >> df;
                    iss >> r2;
                    stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                } readstringRgdata.close();
            } else {
                cout
                << "in TheoryTest::read_individual_stringlen():\n"
                << i << " cannot open.\n";
            }
        }
        //----------------------------------------------------------------------
        /** fractal dimension (df) from Rg vs n **/
    }
    //--------------------------------------------------------------------------
    /** L from string calculations by AMADT **/
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing0);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
        if (stringtype=="stringsRg") {
            if(stringsdf.size()==0) throw 1;
            std::sort(stringsdf.begin(),stringsdf.end(),sortDecreasing0);
            averaging_sorted_data(stringsdf,stringsdf_sortdecreasing);
        }
    } catch (int i) {
        if (i==0) {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << "stringlen size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
        if (i==1) {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << "stringsdf size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
    }
}





void TheoryTest::read_individual_equ_stringlen(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys)
{
    stringlen_sortdecreasing.clear();
    stringsdf_sortdecreasing.clear();
    vector<vector<double>> stringlen,stringsdf;
    
    string i;
    
    /** L extracted from beta_kww **/
    //--------------------------------------------------------------------------
    if (is_use_Lbeta)
    {
        i.clear();
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/qtest_bKWW_min_vs_T_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("."+analysispart+".dat");
        ifstream readstringlendata(i.c_str());
        if (readstringlendata.is_open()) {
            string lineContent;
            double Tdata=0,stringlenData=0;
            double qmin=0,bmin=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent);
                iss >> Tdata;
                iss >> qmin;
                iss >> bmin;
                stringlenData=(2*pi())/qmin;
                stringlen.push_back({Tdata,stringlenData});
            } readstringlendata.close();
        } else {
            cout
            << "in TheoryTest::read_individual_equ_stringlen():\n"
            << i << " cannot open.\n";
        }
    }
    /** L from string calculations by AMADT **/
    //--------------------------------------------------------------------------
    else
    {
        i.clear();
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/"+stringtype+"_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readstringlendata(i.c_str());
        if (readstringlendata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringlenData=0,time=0;
            while (getline(readstringlendata,lineContent)) {
                istringstream iss(lineContent);
                iss >> Tdata;
                iss >> frame_inedx;
                iss >> time;
                iss >> stringlenData;
                stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
            } readstringlendata.close();
        } else {
            cout
            << "in TheoryTest::read_individual_equ_stringlen():\n"
            << i << " cannot open.\n";
        }
        /** fractal dimension (df) from Rg vs n **/
        //----------------------------------------------------------------------
        if (stringtype=="stringsRg")
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/stringsRg_Rgvsn_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readstringRgdata(i.c_str());
            if (readstringRgdata.is_open()) {
                string lineContent;
                int    frame_inedx=0;
                double Tdata=0,time=0,df=0,r2=0;
                while (getline(readstringRgdata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> frame_inedx;
                    iss >> time;
                    iss >> df;
                    iss >> r2;
                    stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                } readstringRgdata.close();
            } else {
                cout
                << "in TheoryTest::read_individual_equ_stringlen():\n"
                << i << " cannot open.\n";
            }
        }
        //----------------------------------------------------------------------
        /** fractal dimension (df) from Rg vs n **/
    }
    //--------------------------------------------------------------------------
    /** L from string calculations by AMADT **/
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing0);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
        if (stringtype=="stringsRg") {
            if(stringsdf.size()==0) throw 1;
            std::sort(stringsdf.begin(),stringsdf.end(),sortDecreasing0);
            averaging_sorted_data(stringsdf,stringsdf_sortdecreasing);
        }
    } catch (int i) {
        if (i==0) {
            cout
            << "in TheoryTest::read_individual_equ_stringlen():\n"
            << "stringlen size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
        if (i==1) {
            cout
            << "in TheoryTest::read_individual_equ_stringlen():\n"
            << "stringsdf size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
    }
    
    /* NOTE: read individual */
    read_individual_taueq_data(sysVar,n_trl,n_sys);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    
    double temp=0;
    double strl=0;
    double stdf=0;
    
    if (stringtype=="stringsRg") {
        if (stringlen.size()!=stringsdf.size()) {
            cout
            << "in TheoryTest::read_individual_equ_stringlen():\n"
            << "stringlen.size() "<<stringlen.size() <<"\n"
            << "stringsdf.size() "<<stringsdf.size() <<"\n"
            << "stringlen.size()!=stringsdf.size()\n"; exit(EXIT_FAILURE);
        }
    }
    stringlen.clear();
    stringsdf.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                if (stringtype=="stringsRg") {
                    stdf=stringsdf_sortdecreasing[ii][1];
                    stringsdf.push_back({temp,stdf});
                } break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
    if (stringtype=="stringsRg") {
        stringsdf_sortdecreasing=stringsdf;
    }
}





void TheoryTest::read_all_stringlen(const StructureClass& sysVar,
                                    const int n_sys,
                                    const bool is_cuttime,
                                    const double logcuttime)
{
    stringlen_sortdecreasing.clear();
    stringsdf_sortdecreasing.clear();
    vector<vector<double>> stringlen,stringsdf;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        
        /** L extracted from beta_kww **/
        //----------------------------------------------------------------------
        if (is_use_Lbeta)
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_bKWW_min_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("."+analysispart+".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                double Tdata=0,stringlenData=0;
                double qmin=0,bmin=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> qmin;
                    iss >> bmin;
                    stringlenData=(2*pi())/qmin;
                    stringlen.push_back({Tdata,stringlenData});
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_stringlen():\n"
                << i << " cannot open.\n";
            }
        }
        /** L from string calculations by AMADT **/
        //----------------------------------------------------------------------
        else
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/"+stringtype+"_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                int    frame_inedx=0;
                double Tdata=0,stringlenData=0,time=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> frame_inedx;
                    iss >> time;
                    iss >> stringlenData;
                    if (is_cuttime) {
                        if (log10(time)<logcuttime) {
                            stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
                        }
                    } else {
                        stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
                    }
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_stringlen():\n"
                << i << " cannot open.\n";
            }
            /** fractal dimension (df) from Rg vs n **/
            //------------------------------------------------------------------
            if (stringtype=="stringsRg")
            {
                i.clear();
                i.append(return_AnalysisFolderPath(sysVar));
                i.append("/fit_data");
                i.append("/stringsRg_Rgvsn_all_");
                i.append(sysVar.get_usic());
                i.append("_00"+to_string((long long int)n_trl));
                i.append("_"+sysVar.get_nameString(n_sys));
                i.append(".dat");
                ifstream readstringRgdata(i.c_str());
                if (readstringRgdata.is_open()) {
                    string lineContent;
                    int    frame_inedx=0;
                    double Tdata=0,time=0,df=0,r2=0;
                    while (getline(readstringRgdata,lineContent)) {
                        istringstream iss(lineContent);
                        iss >> Tdata;
                        iss >> frame_inedx;
                        iss >> time;
                        iss >> df;
                        iss >> r2;
                        stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                    } readstringRgdata.close();
                } else {
                    cout
                    << "in TheoryTest::read_all_stringlen():\n"
                    << i << " cannot open.\n";
                }
            }
            //------------------------------------------------------------------
            /** fractal dimension (df) from Rg vs n **/
        }
        //----------------------------------------------------------------------
        /** L from string calculations by AMADT **/
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing0);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
        if (stringtype=="stringsRg") {
            if(stringsdf.size()==0) throw 1;
            std::sort(stringsdf.begin(),stringsdf.end(),sortDecreasing0);
            averaging_sorted_data(stringsdf,stringsdf_sortdecreasing);
        }
    } catch (int i) {
        if (i==0) {
            cout
            << "in TheoryTest::read_all_stringlen():\n"
            << "stringlen size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
        if (i==1) {
            cout
            << "in TheoryTest::read_all_stringlen():\n"
            << "stringsdf size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
    }
}





void TheoryTest::read_all_equ_stringlen(const StructureClass& sysVar,
                                        const int n_sys)
{
    stringlen_sortdecreasing.clear();
    stringsdf_sortdecreasing.clear();
    vector<vector<double>> stringlen,stringsdf;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        
        /** L extracted from beta_kww **/
        //----------------------------------------------------------------------
        if (is_use_Lbeta)
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_bKWW_min_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("."+analysispart+".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                double Tdata=0,stringlenData=0;
                double qmin=0,bmin=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> qmin;
                    iss >> bmin;
                    stringlenData=(2*pi())/qmin;
                    stringlen.push_back({Tdata,stringlenData});
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_equ_stringlen():\n"
                << i << " cannot open.\n";
            }
        }
        /** L from string calculations by AMADT **/
        //----------------------------------------------------------------------
        else
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/"+stringtype+"_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                int    frame_inedx=0;
                double Tdata=0,stringlenData=0,time=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> frame_inedx;
                    iss >> time;
                    iss >> stringlenData;
                    stringlen.push_back({Tdata,stringlenData,(double)frame_inedx,time});
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_equ_stringlen():\n"
                << i << " cannot open.\n";
            }
            /** fractal dimension (df) from Rg vs n **/
            //------------------------------------------------------------------
            if (stringtype=="stringsRg")
            {
                i.clear();
                i.append(return_AnalysisFolderPath(sysVar));
                i.append("/fit_data");
                i.append("/stringsRg_Rgvsn_all_");
                i.append(sysVar.get_usic());
                i.append("_00"+to_string((long long int)n_trl));
                i.append("_"+sysVar.get_nameString(n_sys));
                i.append(".dat");
                ifstream readstringRgdata(i.c_str());
                if (readstringRgdata.is_open()) {
                    string lineContent;
                    int    frame_inedx=0;
                    double Tdata=0,time=0,df=0,r2=0;
                    while (getline(readstringRgdata,lineContent)) {
                        istringstream iss(lineContent);
                        iss >> Tdata;
                        iss >> frame_inedx;
                        iss >> time;
                        iss >> df;
                        iss >> r2;
                        stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                    } readstringRgdata.close();
                } else {
                    cout
                    << "in TheoryTest::read_all_equ_stringlen():\n"
                    << i << " cannot open.\n";
                }
            }
            //------------------------------------------------------------------
            /** fractal dimension (df) from Rg vs n **/
        }
        //----------------------------------------------------------------------
        /** L from string calculations by AMADT **/
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing0);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
        if (stringtype=="stringsRg") {
            if(stringsdf.size()==0) throw 1;
            std::sort(stringsdf.begin(),stringsdf.end(),sortDecreasing0);
            averaging_sorted_data(stringsdf,stringsdf_sortdecreasing);
        }
    } catch (int i) {
        if (i==0) {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << "stringlen size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
        if (i==1) {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << "stringsdf size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
    }
    
    /* NOTE: read all */
    read_all_taueq_data(sysVar,n_sys);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    
    double temp=0;
    double strl=0;
    double stdf=0;
    
    if (stringtype=="stringsRg") {
        if (stringlen.size()!=stringsdf.size()) {
            cout
            << "in TheoryTest::read_individual_stringlen():\n"
            << "stringlen.size() "<<stringlen.size() <<"\n"
            << "stringsdf.size() "<<stringsdf.size() <<"\n"
            << "stringlen.size()!=stringsdf.size()\n"; exit(EXIT_FAILURE);
        }
    }
    stringlen.clear();
    stringsdf.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                if (stringtype=="stringsRg") {
                    stdf=stringsdf_sortdecreasing[ii][1];
                    stringsdf.push_back({temp,stdf});
                } break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
    if (stringtype=="stringsRg") {
        stringsdf_sortdecreasing=stringsdf;
    }
}





void TheoryTest::read_all_equ_strings_rg(const StructureClass& sysVar,
                                         const int n_sys)
{
    stringsRg_sortdecreasing.clear();
    vector<vector<double>> stringsRg;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data");
        i.append("/stringsRg_all_");
        i.append(sysVar.get_usic());
        i.append("_00"+to_string((long long int)n_trl));
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readstringsRgdata(i.c_str());
        if (readstringsRgdata.is_open()) {
            string lineContent;
            int    frame_inedx=0;
            double Tdata=0,stringsRgData=0,time=0;
            while (getline(readstringsRgdata,lineContent)) {
                istringstream iss(lineContent);
                iss >> Tdata;
                iss >> frame_inedx;
                iss >> time;
                iss >> stringsRgData;
                stringsRg.push_back
                ({Tdata,stringsRgData,(double)frame_inedx,time});
            } readstringsRgdata.close();
        } else {
            cout
            << "in TheoryTest::read_all_equ_strings_rg():\n"
            << i << " cannot open.\n";
        }
    }
    try {
        if(stringsRg.size()==0) throw 0;
        std::sort(stringsRg.begin(),stringsRg.end(),sortDecreasing0);
        averaging_sorted_data(stringsRg,stringsRg_sortdecreasing);
    } catch (int i) {
        cout
        << "in TheoryTest::read_all_equ_strings_rg()\n"
        << "stringsRg size = " << i << "\n";
        stringsRg.push_back({0,0});
    }
    
    /* NOTE: read all */
    read_all_taueq_data(sysVar,n_sys);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringsRg_sortdecreasing.size();
    
    double temp=0;
    double strl=0;
    
    stringsRg.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringsRg_sortdecreasing[ii][0]) {
                temp=stringsRg_sortdecreasing[ii][0];
                strl=stringsRg_sortdecreasing[ii][1];
                stringsRg.push_back({temp,strl});
                break;
            }
        }
    }
    stringsRg_sortdecreasing=stringsRg;
}





void TheoryTest::read_frame_strings_rgvsn(const StructureClass& sysVar,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d,
                                          const int frame)
{
    strings_rgvsn_sortincreasing.clear();
    strings_rgvsnloglog_sortincreasing.clear();
    
    vector<vector<double>> rgvsn,rgvsnloglog;
    
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/fit_data");
    i.append("/stringsRgDistr_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append(".dat");
    ifstream readFile(i.c_str());
    if (readFile.is_open())
    {
        bool is_finishreaddata=false;
        double dubvar=0,n=0,rg=0;
        string target="frame_"+to_string(frame)+";";//ex. "frame_53;"
        string str;
        string lineContent;
        while(getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> str;
            if (str==target) //ex. "frame_53;"
            {
                iss >> str;  //ex. "time_65.52"
                getline(readFile,lineContent);//"bodies fraction Rg"
                while(true) {
                    getline(readFile,lineContent);
                    iss.clear();iss.str(lineContent);
                    iss >> str;
                    if (str=="Overall_Mean") {
                        is_finishreaddata=true;
                        break;
                    } else {
                        n=atof(str.c_str());//bodies
                    }
                    iss >> dubvar;          //fraction
                    iss >> rg;              //Rg
                    if (n>=stringsRg_cutlo&&n<=stringsRg_cuthi) {
                        rgvsn.push_back({n,rg});
                        rgvsnloglog.push_back({log10(n),log10(rg)});
                    }
                } if (is_finishreaddata) break;
            }
        }
    } else {
        cout
        << "in TheoryTest::read_equ_strings_rgvsn():\n"
        << i << " cannot open.\n"; exit(EXIT_FAILURE);
    }
    try {
        if (rgvsn.size()==0) throw 0;
        std::sort(rgvsn.begin(),rgvsn.end(),sortIncreasing0);
        std::sort(rgvsnloglog.begin(),rgvsnloglog.end(),sortIncreasing0);
        strings_rgvsn_sortincreasing=rgvsn;
        strings_rgvsnloglog_sortincreasing=rgvsnloglog;
    } catch (int i) {
        cout
        << "in TheoryTest::read_frame_strings_rgvsn():\n"
        << "rgvsn size = " << i << "\n";
        strings_rgvsn_sortincreasing.clear();
        strings_rgvsnloglog_sortincreasing.clear();
        //exit(EXIT_FAILURE);
    }
}





void TheoryTest::read_all_equ_stringlen(const StructureClass& sysVar,
                                        const int n_sys,
                                        const double TA_d,
                                        const bool is_ltTA)
{
    stringlen_sortdecreasing.clear();
    stringsdf_sortdecreasing.clear();
    vector<vector<double>> stringlen,stringsdf;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        
        /** L extracted from beta_kww **/
        //----------------------------------------------------------------------
        if (is_use_Lbeta)
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/qtest_bKWW_min_vs_T_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("."+analysispart+".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                double Tdata=0,stringlenData=0;//,time=0;
                double qmin=0,bmin=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> qmin;
                    iss >> bmin;
                    stringlenData=(2*pi())/qmin;
                    if (is_ltTA) {
                        if (Tdata<=TA_d) /** collect T<=TA **/
                        {
                            stringlen.push_back({Tdata,stringlenData});
                        }
                    } else {
                        if (Tdata>TA_d) /** collect T>TA **/
                        {
                            stringlen.push_back({Tdata,stringlenData});
                        }
                    }
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_equ_stringlen():\n"
                << i << " cannot open.\n";
            }
        }
        /** L from string calculations by AMADT **/
        //----------------------------------------------------------------------
        else
        {
            i.clear();
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/fit_data");
            i.append("/"+stringtype+"_all_");
            i.append(sysVar.get_usic());
            i.append("_00"+to_string((long long int)n_trl));
            i.append("_"+sysVar.get_nameString(n_sys));
            i.append(".dat");
            ifstream readstringlendata(i.c_str());
            if (readstringlendata.is_open()) {
                string lineContent;
                int    frame_inedx=0;
                double Tdata=0,stringlenData=0,time=0;
                while (getline(readstringlendata,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> Tdata;
                    iss >> frame_inedx;
                    iss >> time;
                    iss >> stringlenData;
                    if (is_ltTA) {
                        if (Tdata<=TA_d) /** collect T<=TA **/
                        {
                            stringlen.push_back
                            ({Tdata,stringlenData,(double)frame_inedx,time});
                        }
                    } else {
                        if (Tdata>TA_d) /** collect T>TA **/
                        {
                            stringlen.push_back
                            ({Tdata,stringlenData,(double)frame_inedx,time});
                        }
                    }
                } readstringlendata.close();
            } else {
                cout
                << "in TheoryTest::read_all_equ_stringlen():\n"
                << i << " cannot open.\n";
            }
            /** fractal dimension (df) from Rg vs n **/
            //------------------------------------------------------------------
            if (stringtype=="stringsRg")
            {
                i.clear();
                i.append(return_AnalysisFolderPath(sysVar));
                i.append("/fit_data");
                i.append("/stringsRg_Rgvsn_all_");
                i.append(sysVar.get_usic());
                i.append("_00"+to_string((long long int)n_trl));
                i.append("_"+sysVar.get_nameString(n_sys));
                i.append(".dat");
                ifstream readstringRgdata(i.c_str());
                if (readstringRgdata.is_open()) {
                    string lineContent;
                    int    frame_inedx=0;
                    double Tdata=0,time=0,df=0,r2=0;
                    while (getline(readstringRgdata,lineContent)) {
                        istringstream iss(lineContent);
                        iss >> Tdata;
                        iss >> frame_inedx;
                        iss >> time;
                        iss >> df;
                        iss >> r2;
                        if (is_ltTA) {
                            if (Tdata<=TA_d) /** collect T<=TA **/
                            {
                                stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                            }
                        } else {
                            if (Tdata>TA_d) /** collect T>TA **/
                            {
                                stringsdf.push_back({Tdata,df,r2,(double)frame_inedx,time});
                            }
                        }
                    } readstringRgdata.close();
                } else {
                    cout
                    << "in TheoryTest::read_all_equ_stringlen():\n"
                    << i << " cannot open.\n";
                }
            }
            //------------------------------------------------------------------
            /** fractal dimension (df) from Rg vs n **/
        }
        //------------------------------------------------------------------
        /** L from string calculations by AMADT **/
    }
    
    try {
        if(stringlen.size()==0) throw 0;
        std::sort(stringlen.begin(),stringlen.end(),sortDecreasing0);
        averaging_sorted_data(stringlen,stringlen_sortdecreasing);
        if (stringtype=="stringsRg") {
            if(stringsdf.size()==0) throw 1;
            std::sort(stringsdf.begin(),stringsdf.end(),sortDecreasing0);
            averaging_sorted_data(stringsdf,stringsdf_sortdecreasing);
        }
    } catch (int i) {
        if (i==0) {
            cout
            << "in TheoryTest::read_all_equ_stringlen():\n"
            << "stringlen size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
        if (i==1) {
            cout
            << "in TheoryTest::read_all_equ_stringlen():\n"
            << "stringsdf size = " << i << "\n";
            exit(EXIT_FAILURE);
        }
    }
    
    /* NOTE: read all */
    read_all_taueq_data(sysVar,n_sys,TA_d,is_ltTA);
    
    int n_taueq=(int)taueq_sortdecreasing.size();
    int strSize=(int)stringlen_sortdecreasing.size();
    
    double temp=0;
    double strl=0;
    double stdf=0;
    
    if (stringtype=="stringsRg") {
        if (stringlen.size()!=stringsdf.size()) {
            cout
            << "in TheoryTest::read_all_equ_stringlen():\n"
            << "stringlen.size() "<<stringlen.size() <<"\n"
            << "stringsdf.size() "<<stringsdf.size() <<"\n"
            << "stringlen.size()!=stringsdf.size()\n"; exit(EXIT_FAILURE);
        }
    }
    stringlen.clear();
    stringsdf.clear();
    
    for (int i=0; i<n_taueq; ++i) {
        for (int ii=0; ii<strSize; ++ii) {
            if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                temp=stringlen_sortdecreasing[ii][0];
                strl=stringlen_sortdecreasing[ii][1];
                stringlen.push_back({temp,strl});
                if (stringtype=="stringsRg") {
                    stdf=stringsdf_sortdecreasing[ii][1];
                    stringsdf.push_back({temp,stdf});
                } break;
            }
        }
    }
    stringlen_sortdecreasing=stringlen;
    if (stringtype=="stringsRg") {
        stringsdf_sortdecreasing=stringsdf;
    }
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
            o.append("/"+stringtype+"_equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_stringlen(o.c_str(),ofstream::app);//'append'
            
            o.clear();
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/"+stringtype+"_equ_inverseT_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
            ofstream write_stringleninv(o.c_str(),ofstream::app);//'append'
            
            o.clear();
            if (stringtype=="stringsRg") {
                o.append(return_AnalysisFolderPath(sysVar));
                o.append("/fit_data");
                o.append("/stringsRg_fractaldf_equ_");
                o.append(sysVar.get_usic());
                o.append("_00"+to_string((long long int)n_trl));
                o.append("_"+sysVar.get_nameString(n_sys));
                o.append(".dat");
            } ofstream write_stringsdf(o.c_str(),ofstream::app);//'append'
            
            /* NOTE: individual */
            //---------------------------------------------
            read_individual_taueq_data(sysVar,n_trl,n_sys);
            read_individual_stringlen(sysVar,n_trl,n_sys);
            //---------------------------------------------
            
            int n_taueq=(int)taueq_sortdecreasing.size();
            int strSize=(int)stringlen_sortdecreasing.size();
            
            for (int i=0; i<n_taueq; ++i) {
                for (int ii=0; ii<strSize; ++ii) {
                    if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                        write_stringlen
                        << stringlen_sortdecreasing[ii][0] << " "
                        << stringlen_sortdecreasing[ii][1] << "\n";
                        if (stringtype=="stringsRg") {
                            write_stringsdf
                            << stringsdf_sortdecreasing[ii][0] << " "
                            << stringsdf_sortdecreasing[ii][1] << "\n";
                        }
                        write_stringleninv
                        << convert/stringlen_sortdecreasing[ii][0] << " "
                        << stringlen_sortdecreasing[ii][1] << "\n";
                        break;
                    }
                }
            }
            write_stringlen.close();
            write_stringsdf.close();
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
        o.append("/"+stringtype+"_equ_avg_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_stringlen(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/"+stringtype+"_equ_avg_inverseT_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream write_stringleninv(o.c_str(),ofstream::app);//'append'
        
        o.clear();
        if (stringtype=="stringsRg") {
            o.append(return_AnalysisFolderPath(sysVar));
            o.append("/fit_data");
            o.append("/stringsRg_fractaldf_equ_avg_");
            o.append(sysVar.get_usic());
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append(".dat");
        } ofstream write_stringsdf(o.c_str(),ofstream::app);//'append'
        
        /* NOTE: all */
        //-----------------------------------
        read_all_taueq_data(sysVar,n_sys);
        read_all_stringlen(sysVar,n_sys);
        //-----------------------------------
        
        int n_taueq=(int)taueq_sortdecreasing.size();
        int strSize=(int)stringlen_sortdecreasing.size();
        
        for (int i=0; i<n_taueq; ++i) {
            for (int ii=0; ii<strSize; ++ii) {
                if (taueq_sortdecreasing[i][0]==stringlen_sortdecreasing[ii][0]) {
                    write_stringlen
                    << stringlen_sortdecreasing[ii][0] << " "
                    << stringlen_sortdecreasing[ii][1] << "\n";
                    if (stringtype=="stringsRg") {
                        write_stringsdf
                        << stringsdf_sortdecreasing[ii][0] << " "
                        << stringsdf_sortdecreasing[ii][1] << "\n";
                    }
                    write_stringleninv
                    << convert/stringlen_sortdecreasing[ii][0] << " "
                    << stringlen_sortdecreasing[ii][1] << "\n";
                    break;
                }
            }
        }
        write_stringlen.close();
        write_stringsdf.close();
        write_stringleninv.close();
    }
}





void TheoryTest::find_LA(const StructureClass& sysVar,
                         const int n_sys,
                         const double TA_d)
{
    /** smooth L(T) data by a chosen function and get LA at TA **/
    read_all_taueq_data(sysVar,n_sys);
    read_all_equ_stringlen(sysVar,n_sys);
    std::vector<std::vector<double>> stringlen;
    for (size_t i=0;i<taueq_sortdecreasing.size();++i) {
        double T  =taueq_sortdecreasing.at(i).at(0);
        double tau=taueq_sortdecreasing.at(i).at(1);
        if (T<TA_d) {
            for (size_t ii=0;ii<stringlen_sortdecreasing.size();++ii) {
                if (T==stringlen_sortdecreasing.at(ii).at(0)) {
                    double L=stringlen_sortdecreasing.at(ii).at(1);
                    if (tau<time_stringlen_threshold) stringlen.push_back({T,L});
                }
            }
        }
    }
    if (stringlen.size()<3) {
        cout
        << "in TheoryTest::find_LA():\n"
        << "stringlen data < 3 \n"; exit(EXIT_FAILURE);
    }
    
    bool is_use_pspline=true;
    
    string function="Arrhenius_string";
    if (is_use_pspline) {
        int m=(int)stringlen.size();
        build_penalizedspline_interpolant(stringlen,0,1,m,2.0);
        LA=spline1dcalc(interpolant,TA_d);
    } else {
        build_function_interpolant(function,stringlen,0,1);
        LA=log10(calc_y_given_x(c,TA_d,function));
    }
    
    if (true)
    {
        string o;
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/stringlen_cuttime_");
        o.append(sysVar.get_usic());
        o.append("_"+sysVar.get_nameString(n_sys));
        o.append(".dat");
        ofstream writeFile(o.c_str());
        for (size_t i=0;i<stringlen.size();++i) {
            double T=stringlen.at(i).at(0);
            double L=stringlen.at(i).at(1);
            double L_fit=0;
            if (is_use_pspline) {
                L_fit=spline1dcalc(interpolant,T);
            } else {
                L_fit=log10(calc_y_given_x(c,T,function));
            } writeFile << T << " " << L << " " << L_fit << "\n";
        }
    }
    
    if (false) /** OLD algorithm **/
    {
        /* interpolate L(T) to get L(TA) */
        double T_pre=0,T_pos=0;
        double stringlen_pre=0,stringlen_pos=0;
        double slope=0;
        read_all_equ_stringlen(sysVar,n_sys); // NOTE: use "equ" data
        for (int i=0; i<(int)(stringlen_sortdecreasing.size()-1); ++i)
        {
            T_pre = stringlen_sortdecreasing[i][0];
            T_pos = stringlen_sortdecreasing[i+1][0];
            
            if (((TA_d-T_pre)*(TA_d-T_pos))<0)
            {
                stringlen_pre=stringlen_sortdecreasing[i][1];
                stringlen_pos=stringlen_sortdecreasing[i+1][1];
                slope = (stringlen_pos-stringlen_pre)/(T_pos-T_pre);
                LA    = stringlen_pre+slope*(TA_d-T_pre);
                return;
            }
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





/*
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
 */





/*
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
 */





void TheoryTest::find_diagonal(const StructureClass& sysVar,
                               const int n_sys,
                               const string& keyword)
{
    string model="diagonal";
    
    set_fitParams(model);
    read_rawxydata(sysVar,n_sys,keyword);
    fitdata_processing_lsfit(dataxy_sortdecreasing);//NOTE: dataxy_sortdecreasing
    alglib_lsfit(model);
    slope_fit=c[0];
    intcp_fit=c[1];
    lntau0=intcp_fit;//NOTE: natural log form
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+keyword+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tau_inf(o.c_str());
    write_fitModel(write_tau_inf,model);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_actual(write_tau_inf,model,dataxy_sortdecreasing);
}





void TheoryTest::fit_GLM_diagonal(const StructureClass& sysVar,
                                  const int n_sys,
                                  const string& keyword)
{
    string model_d="diagonal";
    
    set_fitParams(model_d);
    read_rawxydata(sysVar,n_sys,keyword);
    fitdata_processing_lsfit(dataxy_sortdecreasing);//NOTE:dataxy_sortdecreasing
    alglib_lsfit(model_d);
    slope_fit=c[0];
    intcp_fit=c[1];
    lntau0=intcp_fit;//NOTE: natural log form
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+keyword+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tau_inf(o.c_str());
    write_fitModel(write_tau_inf,model_d);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_actual(write_tau_inf,model_d,dataxy_sortdecreasing);
}





void TheoryTest::fit_u2_Tr(const StructureClass& sysVar,
                           const int n_sys,
                           ofstream& writeFile,
                           const bool is_cutbyTA,
                           const string& u2Tmodel,
                           const string& keyword)
{
    if (is_cutbyTA) {
        read_all_taueq_data(sysVar,n_sys,TA);
        read_all_equ_DWF(sysVar,n_sys,TA);
    } else {
        read_all_taueq_data(sysVar,n_sys);
        read_all_equ_DWF(sysVar,n_sys);
    }
    int n_taueq=(int)taueq_sortdecreasing.size();
    int dwfSize=(int)dwf_sortdecreasing.size();
    
    for (int inner=0;inner<n_taueq;++inner) {
        T_actual=taueq_sortdecreasing[inner][0];
        tau_T   =taueq_sortdecreasing[inner][1];
        for (int i=0; i<dwfSize; ++i) {
            if (T_actual==dwf_sortdecreasing[i][0]) {
                DWF=dwf_sortdecreasing[i][1];
                break;
            }
        }
        writeFile
        << T_actual/TA << " "
        << DWF         << "\n";
    } writeFile.close();
    
    /* fit u2 to T/TA */
    fit_u2T_linear(sysVar,n_sys,u2Tmodel,keyword);
}





void TheoryTest::fit_u2T_linear(const StructureClass& sysVar,
                                const int n_sys,
                                const string& model_d,
                                const string& keyword)
{
    //string model_d="univu2T";//univu2T,linear2
    
    bool is_use_FG_org=is_use_FG;
    is_use_FG=false;
    set_fitParams(model_d);
    read_rawxydata(sysVar,n_sys,keyword);
    fitdata_processing_lsfit(dataxy_sortdecreasing);//NOTE:dataxy_sortdecreasing
    alglib_lsfit(model_d);
    
    if (model_d=="univu2T") {
        u2A      =c[0];
        slope_fit=c[1];
        intcp_fit=c[2];
    } else if (model_d=="linear2") {
        slope_fit=c[0];
        intcp_fit=c[1];
    }
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/ATheoryTest_"+keyword+"_fit_");
    o.append(sysVar.get_usic());
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append(".dat");
    ofstream write_tau_inf(o.c_str());
    write_fitModel(write_tau_inf,model_d);
    write_fitarrays(write_tau_inf);
    write_stopcond(write_tau_inf);
    write_fitinfo(write_tau_inf);
    write_errcurve_actual(write_tau_inf,model_d,dataxy_sortdecreasing);
    //write_errcurve_log(write_tau_inf,model_d,dataxy_sortdecreasing);
    is_use_FG=is_use_FG_org;
}





void TheoryTest::read_rawxydata(const StructureClass& sysVar,
                                const int n_sys,
                                const string& keyword)
{
    dataxy_sortdecreasing.clear();
    vector<vector<double>> dataxy;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/fit_data/ATheoryTest_diagonal_raw_");
        i.append(sysVar.get_usic());
        i.append("_"+sysVar.get_nameString(n_sys));
        i.append(".dat");
        ifstream readFile(i.c_str());
        if (readFile.is_open())
        {
            string lineContent;
            double datax=0,datay=0;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);//(x,y)
                iss >> datax;
                iss >> datay;
                dataxy.push_back({datax,datay});
            } readFile.close();
        } else {
            cout
            << "in TheoryTest::read_rawxydata():\n "
            << i+" cannot open." << "\n";
        }
    }
    try {
        if(dataxy.size()==0) throw 0;
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing0);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    } catch (int i) {
        cout
        << "in TheoryTest::read_rawxydata():\n"
        << "dataxy size = " << i << "\n";
        dataxy.push_back({0,0});
        std::sort(dataxy.begin(),dataxy.end(),sortDecreasing0);
        averaging_sorted_data(dataxy,dataxy_sortdecreasing);
    }
}





void TheoryTest::write_ref_scaling_u2(const StructureClass& sysVar,
                                      const double xtauA,
                                      const double alpha,
                                      const vector<double>& ref)
{
    double TA        = ref.at(0);
    double logtauA   = ref.at(1);
    double u2A       = ref.at(2);
    double Tref      = ref.at(3);
    double logtauref = ref.at(4);
    double u2ref     = ref.at(5);
    string o;
    o.append(sysVar.get_Path());
    o.append("/simulations/");
    o.append(sysVar.get_simType());
    o.append("/");
    o.append(sysVar.get_year());
    o.append("/scaling_ref.dat");
    ofstream writeFile(o.c_str(),ofstream::app);//'append'
    writeFile
    << sysVar.get_usic() <<"  "
    << xtauA <<"  "
    << alpha <<"  "
    << TA    <<"  "<< logtauA   <<"  "<< u2A   <<"  "
    << Tref  <<"  "<< logtauref <<"  "<< u2ref <<"\n";
}





void TheoryTest::write_universal_scaling_u2(const StructureClass& sysVar,
                                            const vector<string>& models,
                                            const vector<alglib::real_1d_array>& c_model,
                                            const alglib::real_1d_array& c_linear,
                                            const vector<vector<double>>& import,
                                            ofstream& writeFile)
{
    double shift=0.5;
    double preft=1.5;
    double logtref=logtauref;
    double logtb=log10(DWF_time);
    
    string model_d = models.at(0);
    string u2Tmodel= models.at(1);
    
    /** Imported Parameters **/
    //--------------------------------------------------------------------------
    vector<double> measured=import.at(0);
    double TA      = measured.at(0);
    double logtauA = measured.at(1); double lntauA=logtauA*log(10);
    double u2A     = measured.at(2);
    double T       = measured.at(3);
    double logtau  = measured.at(4);
    double u2      = measured.at(5);
    double Tg      = measured.at(6);
    double mg      = measured.at(7);
    double rhoA    = measured.at(8);
    double Ea_Arr  = measured.at(9);
    double Eapp    = measured.at(10);
    double Eapp_r  = measured.at(11);
    
    vector<double> thermo=import.at(1);
    double rho     = thermo.at(0);
    double rhostd  = thermo.at(1);
    double vol     = thermo.at(2);
    double volstd  = thermo.at(3);
    double vol2    = thermo.at(4);
    double vol2std = thermo.at(5);
    
    double Kt=(vol2-pow(vol,2))/(T*vol);//isothermal compressibility
    
    vector<double> logSE=import.at(2);
    double logr2t  = logSE.at(0);
    double logSEx  = logSE.at(1);
    double logSEy  = logSE.at(2);
    double logSEyf = logSE.at(3);
    double logSErx = logSE.at(4);
    double logSEry = logSE.at(5);
    double logSEmx = logSE.at(6);
    double logSEmy = logSE.at(7);
    double epsSEfit= logSE.at(8);
    double r2SEfit = logSE.at(9);
    
    /** Reduced Parameters **/
    //--------------------------------------------------------------------------
    double Tr      = T/TA;
    double taur    = pow(10,logtau)/pow(10,logtauA);
    double u2r     = u2/u2A;
    double logtaur = log10(taur);
    
    if (model_d=="univGLM") {
        shift=c_model[0][1];
        preft=c_model[0][2];
    }
    if (is_fitu2T) {
        if (u2Tmodel=="linear2") {
            shift=c_linear[1];//intercept
            preft=c_linear[0];//slope
            u2r=calc_y_given_x(c_linear,Tr,u2Tmodel);
        }
    }
    
    /** Predicted Parameters **/
    //--------------------------------------------------------------------------
    double lntaur_p   = pow(preft*Tr-shift,-alpha/2.0)-1.0;
    double logtaur_p  = lntaur_p*log10(exp(1));
    double lntau_p    = lntauA+lntaur_p;
    double logtau_p   = lntau_p*log10(exp(1));
    double Tr_p       = (1.0/3.0)*(1.0+2.0*pow(log(taur)+1,-2.0/alpha));
    double T_p        = TA*Tr_p;
    double logtGLMVFT = (lntauA-1+pow((2*(TA/3))/(T-(TA/3)),alpha/2))*log10(exp(1));
    
    double logtaurGLM = 0;
    if (GLM1mode=="mode1") {
        logtaurGLM = (pow(u2A/u2,alpha/2)-1.0)*log10(exp(1));
    } else if (GLM1mode=="mode2") {
        logtaurGLM = pow(c_model[0][1],alpha/2.0)*(pow(u2A/u2,alpha/2.0)-1.0)*log10(exp(1));
    } else if (GLM1mode=="mode3") {
        logtaurGLM = (Ea*pow(alpha*TA,-1))*(pow(u2A/u2,alpha/2.0)-1.0)*log10(exp(1));
    }
    
    double logtauruniv = 0;
    if (is_fit_univILP) {
        logtauruniv = ((3/(2*u2r+1))-1)*log10(exp(1));
    }
    
    if (sysVar.get_systemUnit()=="real") {
        logr2t  -=3.0;
        logtau  -=3.0;
        logtauA -=3.0;
        logtau_p-=3.0;
        logtref -=3.0;
        logtb   -=3.0;
        logSEx  -=3.0;
    }
    
    string o;
    if (false) {
        /** file path according to usic **/
        o.append(sysVar.get_Path());
        o.append("/simulations/");
        o.append(sysVar.get_simType());
        o.append("/");
        o.append(sysVar.get_year());
        o.append("/scaling_universal.dat");
    } else {
        /** fixed file path **/
        o.append("/home/jh148/AUTOWORK");
        o.append("/simulations");
        o.append("/AGtest/2015");
        o.append("/scaling_universal.dat");
    }
    //ofstream writeFile(o.c_str(),ofstream::app);//'append'
    
    static int c_access=0;
    if (false) {
        //if (c_access==0) {
        writeFile
        << "USIC "
        << "T inverse_T "
        << "Tg_100s m|Tg_100s "
        << "log(t)|DWF log(tau)|ref "
        << "Ea_Arrhenius Eapp (Eapp)r "
        << "Compressibility(Kt) "
        << "Tg/T log(tau) log(MSD|t) <u2> "
        << "TA log(tauA) <u2>|TA "
        
        << "log(tref) T/Tref <u2>/<u2>ref "
        << "Tr <u2>r alpha|GLM "
        << "log(tau) logtau_p log(tau|GLMVFT) "
        << "T_p T_p/TA "
        << "log(taurGLM) logtaur logtaur_p "
        << "3*Tr 3*Tr_p "
        
        << "avg(Volume) stdev(Volume) avg(Volume^2) stdev(Volume^2) "
        << "avg(rho) stdev(rho) rho|TA rho/rhoA "
        
        << "logSEx logSEy logSEyf "
        << "epsSEfit r2SEfit "
        << "logSErx logSEry logSEmx logSEmy "
        
        << "\n"; ++c_access;
    }
    if (true) {
        writeFile
        << sysVar.get_usic() <<" "
        << T          <<" "<< convert/T <<" "
        << Tg         <<" "<< mg        <<" "
        << logtb      <<" "<< taur      <<" "
        << Ea_Arr     <<" "<< Eapp      <<" "<< Eapp_r     <<" "
        << Kt         <<" "
        << Tg/T       <<" "<< logtau    <<" "<< logr2t     <<" "<< u2 <<" "
        << TA         <<" "<< logtauA   <<" "<< u2A        <<" "
        
        << logtref    <<" "<< T/Tref    <<" "<< u2/u2ref   <<" "
        << Tr         <<" "<< u2r       <<" "<< alpha      <<" "
        << logtau     <<" "<< logtau_p  <<" "<< logtGLMVFT <<" "
        << T_p        <<" "<< T_p/TA    <<" "
        << logtaurGLM <<" "<< logtaur   <<" "<< logtaur_p  <<" "
        << 3*Tr       <<" "<< 3*Tr_p    <<" "
        
        << vol        <<" "<< volstd    <<" "<< vol2       <<" "<< vol2std  <<" "
        << rho        <<" "<< rhostd    <<" "<< rhoA       <<" "<< rho/rhoA <<" "
        
        << logSEx     <<" "<< logSEy    <<" "<< logSEyf    <<" "
        << epsSEfit   <<" "<< r2SEfit   <<" "
        << logSErx    <<" "<< logSEry   <<" "<< logSEmx    <<" "<< logSEmy  <<" "
        
        << "\n";
    } else {
        writeFile
        << sysVar.get_usic() <<" "
        << T          <<" "
        << convert/T  <<" "
        << Tg         <<" "
        << mg         <<" "
        << logtb      <<" "
        << taur       <<" "
        << alpha      <<" "
        << Tg/T       <<" "<< logtau    <<" "<< logr2t     <<" "<< u2       <<" "
        << TA         <<" "<< logtauA   <<" "<< u2A        <<" "
        << Tr         <<" "<< u2r       <<" "
        << logtau     <<" "<< logtau_p  <<" "<< logtGLMVFT <<" "
        << logtauruniv<<" "<< logtaurGLM<<" "<< logtaur    <<" "<< logtaur_p<<" "
        << "\n";
    }
}





void TheoryTest::write_universal_scaling_dfit(const StructureClass& sysVar,
                                              const vector<string>& models,
                                              const alglib::real_1d_array& c_model,
                                              const alglib::real_1d_array& c_linear,
                                              const vector<double>& measured)
{
    //double dfit=c_model[0];
    //double shift=0.5;
    //double preft=1.5;
    
    string model_d = models.at(0);
    string u2Tmodel= models.at(1);
    
    /** Measured Parameters **/
    double TA      = measured.at(0);
    double logtauA = measured.at(1);
    //double lntauA  = logtauA*log(10);
    double u2A     = measured.at(2);
    double T       = measured.at(3);
    double logtau  = measured.at(4);
    double u2      = measured.at(5);
    double Tg      = measured.at(6);
    double mg      = measured.at(7);
    double rho     = measured.at(8);
    double rhoA    = measured.at(9);
    
    /** Reduced Parameters **/
    double Tr     = T/TA;
    double taur   = pow(10,logtau)/pow(10,logtauA);
    double u2r    = u2/u2A;
    //double logtaur= log10(taur);
    
    /** Predicted Parameters **/
    //double lntaur_p  = pow(preft*Tr-shift,-alpha/2.0)-1.0;
    //double logtaur_p = lntaur_p*log10(exp(1));
    //double lntau_p   = lntauA+lntaur_p;
    //double logtau_p  = lntau_p*log10(exp(1));
    //double Tr_p      = (1.0/3.0)*(1.0+2.0*pow(log(taur)+1,-2.0/alpha));
    //double T_p       = TA*Tr_p;
    
    if (sysVar.get_systemUnit()=="real") {
        logtau  -=3.0;
        logtauA -=3.0;
        //logtau_p-=3.0;
    }
    
    string o;
    o.append(sysVar.get_Path());
    o.append("/simulations/");
    o.append(sysVar.get_simType());
    o.append("/");
    o.append(sysVar.get_year());
    o.append("/scaling_dfit.dat");
    ofstream writeFile(o.c_str(),ofstream::app);//'append'
    writeFile
    << sysVar.get_usic() <<" "
    << taur    <<" "
    << dfit    <<" "
    << T       <<" "<< Tg        <<" "<< mg       <<" "
    << Tg/T    <<" "<< logtau    <<" "<< u2       <<" "
    << TA      <<" "<< logtauA   <<" "<< u2A      <<" "
    << rho     <<" "<< rhoA      <<" "<< rho/rhoA <<" "
    << Tr      <<" "<< u2r       <<"\n";
}





/** public setters **/
//------------------------------------------------------------------------------
/* bool */
void TheoryTest::set_is_AGtest(const bool b){is_AGtest=b;}
void TheoryTest::set_is_RFOTtest(const bool b){is_RFOTtest=b;}
void TheoryTest::set_is_GLMtest(const bool b){is_GLMtest=b;}
void TheoryTest::set_is_Leporinitest(const bool b){is_Leporinitest=b;}
void TheoryTest::set_is_HWtest(const bool b){is_HWtest=b;}
void TheoryTest::set_is_reltaur(const bool b){is_reltaur=b;}
void TheoryTest::set_is_use_counting1(const bool b){is_use_counting1=b;}
void TheoryTest::set_is_backwardExtrp(const bool b){is_backwardExtrp=b;}
void TheoryTest::set_is_data_smoothing(const bool b){is_data_smoothing=b;}
/* string */
void TheoryTest::set_stringtype(const std::string&str){stringtype=str;}



