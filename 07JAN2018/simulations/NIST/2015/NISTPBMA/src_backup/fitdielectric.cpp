//
//  fitdielectric.cpp
//  cppWork
//
//  Created by SJH on 11/7/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#include "fitdielectric.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

FitDielectric::FitDielectric(const StructureClass& sysVar):

/* base */
AlglibFittingKernel(),

/* string */
form("NIST"),
func("HN_loss"),
model("COOP"),
func_HN("Kremer"), // "wiki" or "Kremer"

/* bool */
is_use_dc_free(false),
is_manual_cutlowT(false),
is_manualFitSingleDS(false),

/* int */
n_fit_HN(3),
n_used_HN(100),
n_check(3),
n_each_side(200),

/* double */
slope_deviation(0.25),
cutlowT(0),
cutlowT_manual(0),
tau_sectops(12.0),
r2_threshold(0.90),
r2_peak_threshold(0.90),
df_threshold(0.5),
dc_threshold(0.3),
frac_threshold(0.25)
{
    is_use_FG   = false;
    extrp_time  = 100.0;
    systemUnit  = sysVar.get_systemUnit();
    corrFacforT = sysVar.get_corrFacforT();
    precision   = sysVar.get_precision();
}





vector<vector<double>> FitDielectric::dielectricData_preprocess(const StructureClass& sysVar)
{
    vector<vector<double>> tinfo_dielectric;
    vector<double> tinfo_cooling;
    vector<double> tinfo_heating;
    
    string input;
    
    if (form=="NPIC") {
        //
    } else if (form=="NIST") {
        
        /* input file */
        input.append(return_AnalysisFolderPath(sysVar));
        input.append("/statistics/dielectric");
        input.append("/"+sysVar.get_usic());
        input.append(".dat");
        
        /* output file */
        string output;
        string overall;
        overall.append(return_AnalysisFolderPath(sysVar));
        overall.append("/fit_data");
        overall.append("/overall_dielectric_");
        overall.append(sysVar.get_usic());
        overall.append(".dat");
        ofstream writeOverall(overall.c_str(),ofstream::app); // append
        
        ifstream readfile(input.c_str());
        if (readfile.is_open()) {
            int    lineindex=0;
            string lineContent;
            double sfreq=0,epsp=0,epspp=0;
            double T_now=0,T_pre=0;
            double Temp_d=0;
            bool   is_in_cooling=false;
            vector<double> sfreq_vd;
            vector<double> tnow_vd;
            vector<double> epsp_vd;
            vector<double> epspp_vd;
            
            /* first 3 lines are useless */
            getline(readfile,lineContent);
            getline(readfile,lineContent);
            getline(readfile,lineContent);
            
            while (getline(readfile,lineContent)){
                ++lineindex;
                istringstream iss(lineContent);
                iss >> sfreq;
                iss >> T_now; T_now+=273;
                iss >> epsp;
                iss >> epspp;
                
                if (lineindex>1) {
                    /* store data at the end of each temperature block */
                    if (T_now!=T_pre) {
                        writeOverall << "\n";
                        is_in_cooling=false;
                        Temp_d=std::trunc(T_pre*corrFacforT);
                        /* cooling process */
                        if (T_now<T_pre) {
                            is_in_cooling=true;
                            tinfo_cooling.push_back(Temp_d);
                            output.clear();
                            output.append(return_AnalysisFolderPath(sysVar));
                            output.append("/statistics/dielectric");
                            output.append("/cooling_eps_");
                            output.append(sysVar.get_usic());
                            output.append("_T"+to_string((long long int)Temp_d));
                            output.append(".dat");
                        }
                        /* heating process */
                        else if (T_now>T_pre) {
                            is_in_cooling=false;
                            tinfo_heating.push_back(Temp_d);
                            output.clear();
                            output.append(return_AnalysisFolderPath(sysVar));
                            output.append("/statistics/dielectric");
                            output.append("/heating_eps_");
                            output.append(sysVar.get_usic());
                            output.append("_T"+to_string((long long int)Temp_d));
                            output.append(".dat");
                        }
                        ofstream writefile(output.c_str());
                        for (int i=0; i<tnow_vd.size(); ++i) {
                            writefile
                            << 2*pi()*sfreq_vd[i] << " "    // afreq
                            << sfreq_vd[i]        << " "    // sfreq
                            << epspp_vd[i]        << " "    // loss
                            << epsp_vd[i]         << "\n";  // storage
                        } writefile.close();
                        sfreq_vd.clear();
                        tnow_vd.clear();
                        epsp_vd.clear();
                        epspp_vd.clear();
                    }
                }
                T_pre=T_now;
                sfreq_vd.push_back(sfreq);
                tnow_vd.push_back(T_now);
                epsp_vd.push_back(epsp);
                epspp_vd.push_back(epspp);
                
                writeOverall
                << sfreq << " "
                << 2*pi()*sfreq << " "
                << epspp  << " "
                << epsp << "\n";
                
            } readfile.close(); writeOverall.close();
            /* store data of the last T block */
            Temp_d=std::trunc(T_now*corrFacforT);
            if (is_in_cooling) {
                tinfo_cooling.push_back(Temp_d);
                output.clear();
                output.append(return_AnalysisFolderPath(sysVar));
                output.append("/statistics/dielectric");
                output.append("/cooling_eps_");
                output.append(sysVar.get_usic());
                output.append("_T"+to_string((long long int)Temp_d));
                output.append(".dat");
            } else {
                tinfo_heating.push_back(Temp_d);
                output.clear();
                output.append(return_AnalysisFolderPath(sysVar));
                output.append("/statistics/dielectric");
                output.append("/heating_eps_");
                output.append(sysVar.get_usic());
                output.append("_T"+to_string((long long int)Temp_d));
                output.append(".dat");
            }
            ofstream writefile(output.c_str());
            for (int i=0; i<tnow_vd.size(); ++i) {
                writefile
                << 2*pi()*sfreq_vd[i] << " "    // afreq
                << sfreq_vd[i]        << " "    // sfreq
                << epspp_vd[i]        << " "    // loss
                << epsp_vd[i]         << "\n";  // storage
            } writefile.close();
            
            if (tinfo_cooling.size()>0) {
                tinfo_dielectric.push_back(tinfo_cooling);
            }
            if (tinfo_heating.size()>0) {
                tinfo_dielectric.push_back(tinfo_heating);
            }
        }
        else {
            cout
            << "in FitDielectric::dielectricData_preprocess: dielectric file cannot open."
            << "\n";
        }
    } else if (form=="manuallySet") {
        /* input file */
        input.append(return_AnalysisFolderPath(sysVar));
        input.append("/statistics/dielectric");
        input.append("/"+sysVar.get_usic());
        input.append(".dat");
        /* output file */
        string output;
        output.append(return_AnalysisFolderPath(sysVar));
        output.append("/statistics/dielectric");
        output.append("/loss_");
        output.append(sysVar.get_usic());
        output.append(".dat");
        vector<vector<double>> dielectric_freq;
        ifstream readfile(input.c_str());
        if (readfile.is_open()) {
            int    lineindex=0;
            string lineContent;
            while (getline(readfile,lineContent)){
                ++lineindex;
                istringstream iss(lineContent);
                iss >> sfreq; if (is_log_data[0]) sfreq=pow(10,sfreq);
                iss >> loss;  if (is_log_data[1]) loss =pow(10,loss);
                dielectric_freq.push_back({sfreq,loss});
            } readfile.close();
        }
        ofstream writefile(output.c_str());
        for (int i=0; i<dielectric_freq.size(); ++i) {
            writefile
            << dielectric_freq[i][0] << " "
            << dielectric_freq[i][1] << "\n";
        } writefile.close();
    }
    return tinfo_dielectric;
}





void FitDielectric::read_all_dielectric_data(const StructureClass& sysVar,
                                             int process,
                                             double Temp_d)
{
    dielectric_freq_overall.clear();
    dielectric_freq_sortincreasing.clear();
    dc_curve_freq_sortincreasing.clear();
    
    vector<vector<double>> dielectric_freq;
    
    string input;
    
    /* store dielectric data into data container */
    //**************************************************************************
    if (form=="NPIC") {
        //
    } else if (form=="NIST") {
        if (process==0) {
            input.append(return_AnalysisFolderPath(sysVar));
            input.append("/statistics/dielectric/");
            input.append("/cooling_eps_");
            input.append(sysVar.get_usic());
            input.append("_T"+to_string((long long int)Temp_d));
            input.append(".dat");
        } else if (process==1) {
            input.append(return_AnalysisFolderPath(sysVar));
            input.append("/statistics/dielectric/");
            input.append("/heating_eps_");
            input.append(sysVar.get_usic());
            input.append("_T"+to_string((long long int)Temp_d));
            input.append(".dat");
        }
        ifstream readFile(input.c_str());
        if (readFile.is_open()) {
            string lineContent,strData;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);
                iss >> afreq;
                iss >> sfreq;
                iss >> loss;
                iss >> storage;
                dielectric_freq.push_back({afreq,sfreq,loss,storage});
            } readFile.close();
        } else {
            cout
            << "FitDielectric::read_all_dielectric_data: " << "\n"
            << "NIST format: dielectric file cannot open." << "\n";
        }
    } else if (form=="manuallySet") {
        input.append(return_AnalysisFolderPath(sysVar));
        input.append("/statistics/dielectric/");
        input.append("/loss_");
        input.append(sysVar.get_usic());
        input.append(".dat");
        ifstream readFile(input.c_str());
        if (readFile.is_open()) {
            string lineContent,strData;
            while (getline(readFile,lineContent)) {
                istringstream iss(lineContent);
                iss >> sfreq;
                iss >> loss;
                afreq=2*pi()*sfreq;
                storage=loss;
                dielectric_freq.push_back({afreq,sfreq,loss,storage});
            } readFile.close();
        } else {
            cout
            << "FitDielectric::read_all_dielectric_data: "        << "\n"
            << "manuallySet format: dielectric file cannot open." << "\n";
        }
    }
    /* NOTE:
     original dielectric data are sorted angular frequency (low to high) */
    try {
        if (dielectric_freq.size()==0) throw 0;
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    catch (int i) {
        cout
        << "in FitDielectric::read_all_dielectric_data:" << "\n"
        << "NO dielectric data were read!"               << "\n"
        << "System Aborted.              "               << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    //**************************************************************************
    
    afreq_MIN=dielectric_freq[0][0];
    sfreq_MIN=dielectric_freq[0][1];
    afreq_MAX=dielectric_freq[dielectric_freq.size()-1][0];
    sfreq_MAX=dielectric_freq[dielectric_freq.size()-1][1];
    
    log_afreq_MIN=log10(afreq_MIN);
    log_afreq_MAX=log10(afreq_MAX);
    log_sfreq_MIN=log10(sfreq_MIN);
    log_sfreq_MAX=log10(sfreq_MAX);
    
    if (is_use_dc_free) {
        
        vector<vector<double>> dielectric_dc_free;
        
        /* establish cubic spline interpolation routine */
        vector<vector<double>> interpolantdata;
        for (int i=0; i<dielectric_freq.size(); ++i)
        {
            interpolantdata.push_back
            ({
                log(fabs(dielectric_freq[i][0])), // NOTE: ln(afreq)
                dielectric_freq[i][3]             // NOTE: storage
            });
        }
        fitdata_processing_1d(interpolantdata);
        real_1d_array x_in=fit_xData.c_str();
        real_1d_array y_in=fit_yData.c_str();
        alglib::spline1dinterpolant s;
        alglib::spline1dbuildcubic(x_in,y_in,s);
        
        double afreqi=0,sfreqi=0,f=0,df=0,d2f=0;
        double storage_dc=0,loss_dc=0;
        
        for (int i=0; i<dielectric_freq.size(); ++i) {
            afreqi     = dielectric_freq[i][0];
            sfreqi     = dielectric_freq[i][1];
            storage_dc = dielectric_freq[i][3];
            alglib::spline1ddiff(s,log(afreqi),f,df,d2f); // NOTE: ln(afreq)
            loss_dc    = -(pi()/2)*df;
            
            // NOTE:
            // sotre freq domain and loss part in log form
            // storage part in original form
            dielectric_dc_free.push_back
            ({afreqi,sfreqi,loss_dc,storage_dc});
        }
        dielectric_freq_overall=dielectric_dc_free;
    }
    else {
        dielectric_freq_overall=dielectric_freq;
    }
    
    
    /* calculate 1st and 2nd derivatives of loss spline curve */
    calc_loglog_derivatives(dielectric_freq_overall,0,2);
    interpolantLossCurve=interpolant;
    
    
    /* sample data in the proximity of the peaks */
    peak_sampling_bySlope(dielectric_freq_overall);
    
    
    /* sample and fit dc curve data if normal mode is chosen */
    if (!is_use_dc_free) {
        dc_curve_sampling(dielectric_freq_overall);
        if(dc_curve_freq_sortincreasing.size()>1) fit_dc_loss();
    }
    
    
    /*----------------------------------------------------------------------
     Final Product Containers:
     
     "dielectric_freq_overall"
     --------------------------------
     <I dimension>
     individual dielectric data
     
     <II dimension>
     dielectric data:
     0: afreq
     1: sfreq
     2: loss
     3: storage
     
     "dielectric_freq_sortincreasing"
     --------------------------------
     <I dimension>
     individual peaks
     
     <II dimension>
     individual dielectric data
     
     <III dimension>
     dielectric data:
     0: log(afreq)
     1: log(loss)
     2: log(sfreq)
     3: storage
     
     "dc_curve_freq_sortincreasing"
     --------------------------------
     <I dimension>
     individual dielectric data
     
     <II dimension>
     dielectric data:
     0: log(afreq)
     1: log(loss)
     2: log(sfreq)
     3: storage
     ---------------------------------------------------------------------*/
}





void FitDielectric::peak_sampling_byValue(const vector<vector<double>>& data)
{
    vector<vector<double>> dielectric_freq_overall=data; // NOTE
    
    vector<vector<double>> dielectric_freq;
    vector<vector<double>> dielectric_freq_peak;
    
    /* 1. find the lower cutoff (frequency) */
    //----------------------------------------------------------------------
    int cutoff_lo=0;
    double loss_cut_lo=0;
    
    for (int i=0; i<dielectric_freq_overall.size(); ++i)
    {
        if (i>=n_fit_HN) // NOTE: ">="
        {
            int n_count=0;
            for (int count=1; count<=n_fit_HN; ++count)
            {
                if(dielectric_freq_overall[i][2]>
                   dielectric_freq_overall[i-count][2])
                {
                    ++n_count;
                }
            }
            if (n_count==n_fit_HN)
            {
                //cutoff_lo=i;
                cutoff_lo=i-n_fit_HN;
                
                cutoff_lo_afreq.push_back(dielectric_freq_overall[cutoff_lo][0]);
                cutoff_lo_sfreq.push_back(dielectric_freq_overall[cutoff_lo][1]);
                loss_cut_lo=dielectric_freq_overall[cutoff_lo][2];
                break;
            }
        }
    }
    
    /* 2. shrink the original data to only >= cutoff_lo */
    //----------------------------------------------------------------------
    dielectric_freq.clear();
    for (int i=cutoff_lo; i<dielectric_freq_overall.size(); ++i)
    {
        dielectric_freq.push_back
        ({
            dielectric_freq_overall[i][0], // afreq
            dielectric_freq_overall[i][1], // sfreq
            dielectric_freq_overall[i][2], // loss
            dielectric_freq_overall[i][3]  // storage
        });
    }
    try
    {
        if (dielectric_freq.size()==0) throw 0;
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    catch (int i)
    {
        cout
        << "in FitDielectric::peak_sampling_byValue1:" << "\n"
        << "dielectric_freq size = " << i << "\n";
        dielectric_freq.push_back({0,0});
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    dielectric_freq_peak.clear();
    dielectric_freq_peak=dielectric_freq;
    
    /* 3. find the frequency at the loss peak */
    //----------------------------------------------------------------------
    std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortDecreasing2);
    loss_max.push_back(dielectric_freq[0][2]);
    afreq_max.push_back(dielectric_freq[0][0]);
    sfreq_max.push_back(dielectric_freq[0][1]);
    
    /* 4. find the loss peak index */
    //----------------------------------------------------------------------
    int peak_index=0;
    for (int i=0; i<dielectric_freq_peak.size(); ++i)
    {
        if (dielectric_freq_peak[i][2]==loss_max[0])
        {
            peak_index=i;
            break;
        }
    }
    
    /* 5. find the higher cutoff (frequency) */
    //----------------------------------------------------------------------
    int  count_n_used_HN=0;
    int  n_used_HN_tmp=n_used_HN;
    int  cutoff_hi=0;
    bool is_cutoff_hi=false;
    bool is_fixed_points=true;
    
    dielectric_freq.clear();
    for (int i=0; i<dielectric_freq_peak.size(); ++i)
    {
        dielectric_freq.push_back
        ({
            // NOTE:
            //------------------------------------------------------
            // order changed to conform to fitdata_procrssing format
            // first 2 columns are fitting data
            //------------------------------------------------------
            dielectric_freq_peak[i][0], // afreq
            dielectric_freq_peak[i][2], // loss
            dielectric_freq_peak[i][1], // sfreq
            dielectric_freq_peak[i][3]  // storage
        });
        ++count_n_used_HN;
        
        if (is_fixed_points) // 1st method: use fixed number of points
        {
            is_cutoff_hi=
            (count_n_used_HN>=n_used_HN_tmp)&&
            (dielectric_freq_peak[i][0]>afreq_max[0]);
        }
        else // 2nd method: use horozontal cut (loss_cut)
        {
            is_cutoff_hi=
            (dielectric_freq_peak[i][2]<loss_cut_lo)&& // NOTE: "<"
            (dielectric_freq_peak[i][0]>afreq_max[0]);
        }
        
        if (is_cutoff_hi)
        {
            if (is_fixed_points)
            {
                if (abs(i-peak_index)<=5)
                {
                    n_used_HN_tmp += 5;
                    continue;
                }
            }
            cutoff_hi=i;
            cutoff_hi_afreq.push_back(dielectric_freq_peak[cutoff_hi][0]);
            cutoff_hi_sfreq.push_back(dielectric_freq_peak[cutoff_hi][1]);
            break;
        }
    }
    if (!is_cutoff_hi) // reach the end of data
    {
        cutoff_hi_afreq.push_back(dielectric_freq_peak[count_n_used_HN-1][0]);
        cutoff_hi_sfreq.push_back(dielectric_freq_peak[count_n_used_HN-1][1]);
    }
    try
    {
        if (dielectric_freq.size()==0) throw 0;
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    catch (int i)
    {
        cout
        << "in FitDielectric::peak_sampling_byValue2:" << "\n"
        << "dielectric_freq size = " << i << "\n";
        dielectric_freq.push_back({0,0});
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    dielectric_freq_peak.clear();
    dielectric_freq_peak=dielectric_freq;
    
    
    dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
}





void FitDielectric::peak_sampling_bySlope(const vector<vector<double>>& data)
{
    // total # of data points
    n_points=data.size();
    
    if (peak_indices_sampled.size()>0) {
        interpolant = interpolantReducedLoss;
    }
    else {
        interpolant = interpolantLossCurve;
    }
    
    /* 1. Distinguish Local Peaks */
    //**************************************************************************
    std::vector<int> peak_indices;
    
    for (indexii=1; indexii<n_points; ++indexii) { // NOTE: start from 1
        
        afreq_now=log10(fabs(data[indexii][0]));
        afreq_pre=log10(fabs(data[indexii-1][0]));
        alglib::spline1ddiff(interpolant,afreq_now,f_now,df_now,d2f_now);
        alglib::spline1ddiff(interpolant,afreq_pre,f_pre,df_pre,d2f_pre);
        
        if (df_pre*df_now<0) {
            
            if (data[indexii][2]>data[indexii-1][2]) { // NOTE: >
                peak_indexii=indexii;
            }
            else {
                peak_indexii=indexii-1;
            }
            
            /* NOTE:
             current peak should NOT repeat or be too close to peaks already
             collected */
            if (peak_indices_sampled.size()>0) {
                if (is_peaksRepeat(peak_indexii,peak_indices_sampled)) continue;
                if (is_peaksTooClose(peak_indexii,peak_indices_sampled)) continue;
            }
            if (peak_indices.size()>0) {
                if (is_peaksRepeat(peak_indexii,peak_indices)) continue;
                if (is_peaksTooClose(peak_indexii,peak_indices)) continue;
            }
            
            /* validate the current found peak */
            peak_validation(data,peak_indexii);
            
            //cout << peak_indexii << " " << n_count << " " << n_total << "\n";
            
            /* collect validated peaks */
            if (n_count==n_total) {
                peak_indices.push_back(peak_indexii);
            }
        }
    }
    //**************************************************************************
    
    //cout << "peak_indices.size() " << peak_indices.size() << "\n";
    
    /* 2. Data Sampling Around Found Local Peaks */
    //**************************************************************************
    for (whichpeak=0; whichpeak<peak_indices.size(); ++whichpeak) {
        
        dielectric_freq.clear();
        dielectric_freq_peak.clear();
        
        peak_indexii=peak_indices.at(whichpeak);
        
        dielectric_freq.push_back
        ({
            data[peak_indexii][0],
            data[peak_indexii][1],
            data[peak_indexii][2],
            data[peak_indexii][3]
        });
        
        /* sample fit data wrt to local peaks */
        backward_sampling(data,peak_indexii);
        forward_sampling(data,peak_indexii);
        
        if (dielectric_freq.size()<3) { // # of data too few, not taken the peak
            continue;
        } else {
            
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
            
            for (int i=0; i<dielectric_freq.size(); ++i) {
                dielectric_freq_peak.push_back
                ({
                    log10(fabs(dielectric_freq[i][0])), // log(afreq)
                    log10(fabs(dielectric_freq[i][2])), // log(loss)
                    log10(fabs(dielectric_freq[i][1])), // log(sfreq)
                    dielectric_freq[i][3]               // storage
                });
            }
            /* NOTE:
             only take the current peak if the r2 of curve-fitting is greater
             than r2_peak_threshold  */
            fitdata_processing_lsfit(dielectric_freq_peak);
            fit_loss_peak();
            
            if (rep.r2>r2_peak_threshold) {
                afreq_max.push_back(data[peak_indexii][0]);
                sfreq_max.push_back(data[peak_indexii][1]);
                loss_max.push_back(data[peak_indexii][2]);
                cutoff_lo_afreq.push_back(dielectric_freq[0][0]);
                cutoff_lo_sfreq.push_back(dielectric_freq[0][1]);
                cutoff_hi_afreq.push_back(dielectric_freq[dielectric_freq.size()-1][0]);
                cutoff_hi_sfreq.push_back(dielectric_freq[dielectric_freq.size()-1][1]);
                peak_indices_sampled.push_back(peak_indexii);
                dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
            } else {
                continue;
            }
            
        }
    }
    //**************************************************************************
}





void FitDielectric::dc_curve_sampling(const vector<vector<double>>& data)
{
    // total # of data points
    n_points=data.size();
    
    int firstLossPeakii=0;
    if (peak_indices_sampled.size()>0) firstLossPeakii=peak_indices_sampled[0];
    
    if (firstLossPeakii>n_check) {
        
        /* 1. Distinguish Valid Base (lowest point in loss w/ dc effect ) */
        //**********************************************************************
        for (indexii=0; indexii<firstLossPeakii; ++indexii) { // NOTE: start from 0
            
            afreq_now=log10(fabs(data.at(indexii)[0]));
            alglib::spline1ddiff(interpolantLossCurve,afreq_now,f_now,df_now,d2f_now);
            
            bool is_valid = df_now<0 && fabs(df_now+1.0)<dc_threshold;
            
            if (is_valid) {
                
                base_indexii=indexii;
                
                /* validate the current found base */
                base_validation(data,base_indexii);
                
                /* collect validated base */
                //--------------------------------------------------------------
                if(n_count==n_total) {
                    afreq_base.push_back(data[base_indexii][0]);
                    sfreq_base.push_back(data[base_indexii][1]);
                    loss_base.push_back(data[base_indexii][2]);
                }
            }
        }
        //**********************************************************************
        
        /* 2. Data Sampling by the Found Base */
        //**********************************************************************
        dielectric_freq.clear();
        dielectric_freq_peak.clear();
        
        dielectric_freq.push_back
        ({
            data[base_indexii][0],
            data[base_indexii][1],
            data[base_indexii][2],
            data[base_indexii][3]
        });
        
        /* sample data for dc fitting */
        base_sampling(data,base_indexii);
        
        if (dielectric_freq.size()<3) {
            return;
        }
        else {
            
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
            
            for (int i=0; i<dielectric_freq.size(); ++i) {
                dielectric_freq_peak.push_back
                ({
                    log10(fabs(dielectric_freq[i][0])), // log(afreq)
                    log10(fabs(dielectric_freq[i][2])), // log(loss)
                    log10(fabs(dielectric_freq[i][1])), // log(sfreq)
                    dielectric_freq[i][3]               // storage
                });
            } dc_curve_freq_sortincreasing = dielectric_freq_peak; // NOTE!
        }
        //**********************************************************************
    }
    else {
        return;
    }
}





void FitDielectric::mastercurve_superposition(const vector<vector<double>>& data)
{
    vector<vector<double>> mastercurve;
    double afreq=0;
    double loss_afreq=0;
    for (int i=0; i<data.size(); ++i) {
        afreq=data[i][0];
        loss_afreq=0;
        /* dc contribution */
        if (dc_curve_freq_sortincreasing.size()>1) {
            loss_afreq += dc_loss_value(c_dc_loss,afreq);
        }
        /* HN fits around local peaks */
        for (int ii=0; ii<c_losspeaks.size(); ++ii) {
            loss_afreq += HN_func_value(c_losspeaks[ii],afreq,func);
        } mastercurve.push_back({log10(afreq),log10(fabs(loss_afreq))}); // log-log
    } masterloss=mastercurve;
}





void FitDielectric::fit_dc_loss()
{
    fitdata_processing_lsfit(dc_curve_freq_sortincreasing);
    set_fitParams("dc_loss");
    alglib_lsfit("dc_loss");
    c_dc_loss=c;
    rep_dc_loss=rep;
}





void FitDielectric::fit_loss_peak()
{
    set_fitParams(func);
    alglib_lsfit(func);
}





void FitDielectric::calc_derivatives(const vector<vector<double>>& data,
                                     const int index_x,
                                     const int index_y)
{
    vector<vector<double>> interpolantdata;
    for (int i=0; i<data.size(); ++i) {
        interpolantdata.push_back
        ({
            data[i][index_x],
            data[i][index_y]
        });
    }
    fitdata_processing_1d(interpolantdata);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    alglib::spline1dbuildcubic(x_in,y_in,interpolant);
}





void FitDielectric::calc_loglog_derivatives(const vector<vector<double>>& data,
                                            const int index_x,
                                            const int index_y)
{
    vector<vector<double>> interpolantdata;
    for (int i=0; i<data.size(); ++i) {
        interpolantdata.push_back
        ({
            log10(fabs(data[i][index_x])),
            log10(fabs(data[i][index_y]))
        });
    }
    fitdata_processing_1d(interpolantdata);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    alglib::spline1dbuildcubic(x_in,y_in,interpolant);
}





void FitDielectric::calc_taueq_derivatives()
{
    /* establish cubic spline interpolation routine */
    vector<vector<double>> interpolantdata;
    double invT=0;
    double logtau=0;
    for (int i=0; i<taueq_sortdecreasing.size(); ++i)
    {
        invT=(1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0]; // NOTE: inverse T
        logtau=taueq_sortdecreasing[i][1]; // NOTE: log(taueq)
        interpolantdata.push_back({invT,logtau});
    }
    fitdata_processing_1d(interpolantdata);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    alglib::spline1dbuildcubic(x_in,y_in,interpolanttaueqCurve);
}





void FitDielectric::peak_validation(const vector<vector<double>>& data,
                                    const int peak_index)
{
    int n_points=data.size();
    
    /* forward checking:
     strict monotonic "decreasing" */
    //------------------------------------------------------------------
    n_count=0; n_total=n_check;
    
    for (int count=1; count<=n_check; ++count) {
        if ((peak_index+count)<n_points) {
            bool is_valid =
            data[peak_index+count-1][2]>data[peak_index+count][2];
            if (is_valid) ++n_count;
        }
    }
    
    /* backward checking:
     strict monotonic "decreasing" */
    //------------------------------------------------------------------
    bool is_at_hifreq_end = abs(peak_index-(n_points-1))<=(int)ceil(frac_threshold*n_points);
    bool is_at_lofreq_end = peak_index<=ceil(frac_threshold*n_points);
    
    if (is_at_hifreq_end||is_at_lofreq_end) {
        n_total=n_check*2;
        for (int count=1; count<=n_check; ++count) {
            if ((peak_index-count+1)<n_points && (peak_index-count)>=0) {
                bool is_valid =
                data[peak_index-count+1][2]>data[peak_index-count][2];
                if (is_valid) ++n_count;
            }
        }
    }
}





void FitDielectric::base_validation(const vector<vector<double>>& data,
                                    const int base_index)
{
    int n_points=data.size();
    
    /* forward checking:
     strict monotonic "decreasing" */
    //------------------------------------------------------------------
    n_count=0; n_total=n_check;
    
    for (int count=1; count<=n_check; ++count) {
        if ((base_index+count)<n_points && base_index+count-1>=0) {
            bool is_valid =
            data[base_index+count-1][2]>data[base_index+count][2];
            if (is_valid) ++n_count;
        }
    }
    /* NOTE:
     for conductivity effect, only forward checking for monotonic decreasing
     should be sufficient */
}





void FitDielectric::base_sampling(const vector<vector<double>>& data,
                                  const int base_index)
{
    int  n_points=data.size();
    
    int  index_lo=0;
    bool is_valid=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        /* backward sampling from base_index-1 */
        index_lo = base_index-j;
        is_valid = index_lo<n_points && index_lo>=0;
        
        if (is_valid) {
            
            afreq_now=log10(fabs(data.at(index_lo)[0]));
            alglib::spline1ddiff(interpolantLossCurve,afreq_now,f_now,df_now,d2f_now);
            
            is_taken =
            df_now<0
            && fabs(df_now+1.0)<dc_threshold;
            
            if (is_taken) {
                dielectric_freq.push_back
                ({
                    data[index_lo][0],
                    data[index_lo][1],
                    data[index_lo][2],
                    data[index_lo][3]
                });
            }
            else {
                break;
            }
        }
    } return;
}





void FitDielectric::backward_sampling(const vector<vector<double>>& data,
                                      const int peak_index)
{
    int  n_points=data.size();
    
    int  index_lo=0;
    bool is_valid=false;
    bool is_at_lofreq_end=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        index_lo = peak_index-j;
        is_valid = index_lo<n_points && index_lo-2>=0;
        
        if (is_valid) {
            
            is_at_lofreq_end =
            index_lo<=ceil(frac_threshold*n_points);
            
            afreq_now=log10(fabs(data.at(index_lo)[0]));   // j   behind peak
            afreq_pre=log10(fabs(data.at(index_lo-1)[0])); // j+1 behind peak
            afreq_prr=log10(fabs(data.at(index_lo-2)[0])); // j+2 behind peak
            alglib::spline1ddiff(interpolant,afreq_now,f_now,df_now,d2f_now);
            alglib::spline1ddiff(interpolant,afreq_pre,f_pre,df_pre,d2f_pre);
            alglib::spline1ddiff(interpolant,afreq_prr,f_prr,df_prr,d2f_prr);
            
            if (is_at_lofreq_end) {
                is_taken =
                df_now*df_pre>0
                && df_now*df_prr>0
                //&& d2f_now*d2f_pre>0
                && fabs(df_now-df_pre)<df_threshold;
            }
            else {
                is_taken=
                df_now*df_pre>0
                && df_now*df_prr>0;
                //&& d2f_now*d2f_pre>0;
            }
            
            if (is_taken) {
                dielectric_freq.push_back
                ({
                    data[index_lo][0],
                    data[index_lo][1],
                    data[index_lo][2],
                    data[index_lo][3]
                });
            }
            else {
                break;
            }
        }
    } return;
}





void FitDielectric::forward_sampling(const vector<vector<double>>& data,
                                     const int peak_index)
{
    int  n_points=data.size();
    
    int  index_hi=0;
    bool is_valid=false;
    bool is_at_hifreq_end=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        index_hi = peak_index+j;
        is_valid = index_hi+2<n_points && index_hi>=0;
        
        if (is_valid) {
            
            is_at_hifreq_end =
            abs(index_hi-(n_points-1))<=(int)ceil(frac_threshold*n_points);
            
            afreq_now=log10(fabs(data.at(index_hi)[0]));   // j   ahead of peak
            afreq_pre=log10(fabs(data.at(index_hi+1)[0])); // j+1 ahead of peak
            afreq_prr=log10(fabs(data.at(index_hi+2)[0])); // j+2 ahead of peak
            alglib::spline1ddiff(interpolant,afreq_now,f_now,df_now,d2f_now);
            alglib::spline1ddiff(interpolant,afreq_pre,f_pre,df_pre,d2f_pre);
            alglib::spline1ddiff(interpolant,afreq_prr,f_prr,df_prr,d2f_prr);
            
            if (is_at_hifreq_end) {
                is_taken=
                df_now*df_pre>0
                && df_now*df_prr>0
                //&& d2f_now*d2f_pre>0
                && fabs(df_now-df_pre)<df_threshold;
            }
            else {
                is_taken=
                df_now*df_pre>0
                && df_now*df_prr>0;
                //&& d2f_now*d2f_pre>0;
            }
            
            if (is_taken) {
                dielectric_freq.push_back
                ({
                    data[index_hi][0],
                    data[index_hi][1],
                    data[index_hi][2],
                    data[index_hi][3]
                });
            }
            else{
                break;
            }
        }
    } return;
}





bool FitDielectric::is_peaksRepeat(const int peak_index,
                                   std::vector<int>& peakindices)
{
    bool is_repeat=false;
    //--------------------------------------------------------------
    // "find" returns iterators in the range: [first,last); it would
    // return a pointer to the "last" element if nothing is found.
    //--------------------------------------------------------------
    vector<int>::iterator itr;
    itr=find(peakindices.begin(),peakindices.end(),peak_index);
    if (itr!=peakindices.end()) is_repeat=true;
    else is_repeat=false;
    return is_repeat;
}





bool FitDielectric::is_peaksTooClose(const int peak_index,
                                     const vector<int>& peakindices)
{
    bool is_tooClose=false;
    for (int i=0; i<peakindices.size(); ++i) {
        if (abs(peak_index-peakindices[i])<=n_check) {
            is_tooClose=true; break;
        }
    } return is_tooClose;
}





void FitDielectric::decide_cutlowT(const vector<vector<double>>& data)
{
    double invT_now=0,invT_pre=0;
    double logtau_now=0,logtau_pre=0;
    double f_now=0,df_now=0,d2f_now=0;
    double f_pre=0,df_pre=0,d2f_pre=0;
    
    for (int i=1; i<data.size(); ++i) // NOTE: start from 1
    {
        invT_now=(1000.0/(corrFacforT/precision))/data[i][0];
        invT_pre=(1000.0/(corrFacforT/precision))/data[i-1][0];
        logtau_now=data[i][1];
        logtau_pre=data[i-1][1];
        alglib::spline1ddiff(interpolanttaueqCurve,invT_now,f_now,df_now,d2f_now);
        alglib::spline1ddiff(interpolanttaueqCurve,invT_pre,f_pre,df_pre,d2f_pre);
        
        if (df_now==0) df_pre=1e-6;
        if (((fabs(df_now-df_pre)/fabs(df_pre))>slope_deviation)&&(df_now<df_pre)) {
            cutlowT=data[i-1][0];
            //cout << "\n\n" << "cutlowT " << cutlowT << "\n\n";
            break;
        }
        if(false) {
            if(i==1) {
                cout
                << invT_pre << " " << logtau_pre << " "
                << f_pre << " " << df_pre << " " << d2f_pre << "\n";
            }
            cout
            << invT_now << " " << logtau_now << " "
            << f_now << " " << df_now << " " << d2f_now << "\n";
        }
    }
}





void FitDielectric::lossCurveReduction(const vector<vector<double>>& loss,
                                       const vector<vector<double>>& masterloss)
{
    /* NOTE:
     "dielectric_freq_overall_reduced" should preserve the data structure:
     {afreq,sfreq,lossReduced,storage} */
    double loss_reduced=0; vector<vector<double>> loss_reduced_vec;
    for (int i=0; i<loss.size(); ++i) {
        //loss_reduced=pow(10,masterloss[i][1])-loss[i][2];
        loss_reduced=loss[i][2]-pow(10,masterloss[i][1]);
        loss_reduced_vec.push_back
        ({loss[i][0],loss[i][1],loss_reduced,loss[i][3]});
    } dielectric_freq_overall_reduced=loss_reduced_vec;
}





void FitDielectric::clear_containers()
{
    // NOTE:
    // containers whose content are storing T dependent information
    // should be clear each time when used to store new info of another T's
    
    loss_max.clear();  loss_base.clear();
    afreq_max.clear(); afreq_base.clear();
    sfreq_max.clear(); sfreq_base.clear();
    
    cutoff_lo_afreq.clear(); cutoff_lo_sfreq.clear();
    cutoff_hi_afreq.clear(); cutoff_hi_sfreq.clear();
    
    peak_indices_sampled.clear();
    c_losspeaks.clear();
    
    masterloss.clear();
}





void FitDielectric::cout_peaksinfo()
{
    if (T_nominal!=0) cout << "T_nominal " << T_nominal << " ";
    cout
    << "n_peaks " << n_peaks << "\n";
    for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
        cout
        //<< "peak_" << whichpeak+1 << " "
        << log10(afreq_max[whichpeak]) << "\n";
    }
}
void FitDielectric::write_peaksinfo(const StructureClass& sysVar)
{
    string output;
    output.append(return_AnalysisFolderPath(sysVar));
    output.append("/fit_data");
    output.append("/peaksinfo_");
    output.append(sysVar.get_usic());
    output.append(".dat");
    ofstream write_peaks(output.c_str(),ofstream::app); // append
    
    if (T_nominal!=0) {
        write_peaks << "T_nominal " << T_nominal << " ";
    }
    write_peaks << "n_peaks " << n_peaks << "\n";
    for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
        write_peaks
        //<< "peak_" << whichpeak+1 << " "
        << log10(afreq_max[whichpeak]) << "\n";
    }
}




void FitDielectric::fit_HN(StructureClass& sysVar,
                           int process,
                           double Temp_d)
{
    is_use_FG=false;
    
    T_nominal=Temp_d;
    T_actual=T_nominal*pow(corrFacforT,-1); // NOTE: convert!
    
    string o1,o2,o3,o4,o5,o6;
    
    if (process==0) {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/cooling_");
        o1.append("taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/cooling_");
        o2.append("taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/cooling_");
        o3.append("fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/cooling_");
        o4.append("fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)T_nominal));
        o4.append(".dat");
    } else if (process==1) {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/heating_");
        o1.append("taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/heating_");
        o2.append("taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/heating_");
        o3.append("fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/heating_");
        o4.append("fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)T_nominal));
        o4.append(".dat");
    } else if (process==2) {
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append(".dat");
    }
    
    if (process!=2) {
        o5.clear();
        o5.append(return_AnalysisFolderPath(sysVar));
        o5.append("/fit_data");
        o5.append("/overall_");
        o5.append("taueq_HN_");
        o5.append(sysVar.get_usic());
        o5.append(".dat");
        
        o6.clear();
        o6.append(return_AnalysisFolderPath(sysVar));
        o6.append("/fit_data");
        o6.append("/overall_");
        o6.append("taueq_HN_inverseT_");
        o6.append(sysVar.get_usic());
        o6.append(".dat");
    }
    
    ofstream write_taueq(o1.c_str(),ofstream::app);       // 'append'
    ofstream write_taueqinv(o2.c_str(),ofstream::app);    // 'append'
    ofstream write_params(o3.c_str(),ofstream::app);      // 'append'
    ofstream write_HN(o4.c_str(),ofstream::app);          // 'append'
    ofstream write_alltaueq(o5.c_str(),ofstream::app);    // 'append'
    ofstream write_alltaueqinv(o6.c_str(),ofstream::app); // 'append'
    
    /* Clear containers before use */
    clear_containers();
    
    /*--------------------------------------------------------------------------
     read_all_dielectric_data would give the data containers as follows:
     dielectric_freq_overall        -- 2D double
     dielectric_freq_sortincreasing -- 3D double
     dc_curve_freq_sortincreasing   -- 2D double
     -------------------------------------------------------------------------*/
    read_all_dielectric_data(sysVar,process,T_nominal);
    
    n_peaks=dielectric_freq_sortincreasing.size();
    
    if (n_peaks==0) { /* if no peaks were found */
        write_fit_HN(write_HN);
    }
    else { /* if n_peaks were found */
        
        /*----------------------------------------------------------------------
         check if any further peaks (ex. not found in first sampling)
         can be found by reducing the loss data by the master fit curve
         ---------------------------------------------------------------------*/
        for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
            fitdata_processing_lsfit(dielectric_freq_sortincreasing[whichpeak]);
            fit_loss_peak(); c_losspeaks.push_back(c);
        } mastercurve_superposition(dielectric_freq_overall);
        lossCurveReduction(dielectric_freq_overall,masterloss);
        calc_loglog_derivatives(dielectric_freq_overall_reduced,0,2);
        interpolantReducedLoss=interpolant;
        peak_sampling_bySlope(dielectric_freq_overall_reduced);
        c_losspeaks.clear();
        n_peaks=dielectric_freq_sortincreasing.size();
        
        for (whichpeak=0; whichpeak<n_peaks; ++whichpeak){
            
            // fit rountine
            n_points=dielectric_freq_sortincreasing[whichpeak].size();
            fitdata_processing_lsfit(dielectric_freq_sortincreasing[whichpeak]);
            fit_loss_peak(); c_losspeaks.push_back(c);
            if (rep.r2<r2_threshold) write_badR2(write_HN,r2_threshold);
            if (n_points<n_fit_HN) write_insufficientDara(write_HN,n_points,n_fit_HN);
            
            // write fitinfo to file
            write_peak_info(write_HN);
            write_fitModel(write_HN,func);
            write_fitarrays(write_HN);
            write_stopcond(write_HN);
            write_fitinfo(write_HN);
            write_fit_HN(write_HN);
            write_HN_params(write_params);
            
            // master curve superposition
            mastercurve_superposition(dielectric_freq_overall);
        }
        
        // write master curve info to file
        write_dc_loss(write_HN);
        write_master_curve(write_HN);
        
        // write relaxation time to file
        write_taueqToFile(write_taueq,write_taueqinv,write_alltaueq,write_alltaueqinv);
    }
    
    cout_peaksinfo();
    write_peaksinfo(sysVar);
    
    write_taueq.close();
    write_taueqinv.close();
    write_params.close();
    write_HN.close();
    write_alltaueq.close();
    write_alltaueqinv.close();
}





void FitDielectric::fit_tauFit(StructureClass& sysVar)
{
    is_use_FG=false; // NOTE: grad[0] is too large if time unit is "sec"
    
    string output;
    output.append(return_AnalysisFolderPath(sysVar));
    output.append("/fit_data");
    output.append("/fit_taueq_overall_");
    output.append(sysVar.get_usic());
    output.append(".dat");
    ofstream write_fit_tauFit(output.c_str());
    
    if (write_fit_tauFit.is_open()) {
        
        get_Theq(sysVar);
        set_fitParams(model);
        
        /* cutlowT determined manually */
        if (is_manual_cutlowT) {
            cutlowT=cutlowT_manual;
        }
        /* determine cutlowT by turning point at taueq */
        else {
            cutlowT=0;
            read_individual_taueq_data(sysVar);
            calc_taueq_derivatives();
            decide_cutlowT(taueq_sortdecreasing);
        }
        
        read_individual_taueq_data(sysVar);
        fitdata_processing_lsfit(taueq_sortdecreasing);
        alglib_lsfit(model);
        real_1d_array c_org=c;
        
        write_fitModel(write_fit_tauFit,model);
        write_fitcutoff(write_fit_tauFit);
        write_fitarrays(write_fit_tauFit);
        write_stopcond(write_fit_tauFit);
        write_Tg_fragility(write_fit_tauFit,c_org,model);
        write_fitinfo(write_fit_tauFit);
        write_errcurve(write_fit_tauFit,model);
        write_tauFitfit_vs_T(write_fit_tauFit,c_org,model);
        write_tauFit_vs_T(write_fit_tauFit);
        write_tauFitfit_vs_invT(write_fit_tauFit,c_org,model);
        write_tauFit_vs_invT(write_fit_tauFit);
        
        write_fit_tauFit.close();
    }
    else {
        cout << "FitDielectric::fit_tauFit: 'fit_tauFit.dat' cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
    }
}





void FitDielectric::read_individual_taueq_data(const StructureClass& sysVar)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string o1;
    o1.append(return_AnalysisFolderPath(sysVar));
    o1.append("/fit_data");
    o1.append("/overall_");
    o1.append("taueq_HN_");
    o1.append(sysVar.get_usic());
    o1.append(".dat");
    
    ifstream readtaueqdata(o1.c_str());
    
    if (readtaueqdata.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readtaueqdata,lineContent)) {
            istringstream iss(lineContent); // (T,taueq)
            iss >> Tdata;
            if (is_manual_cutlowT) {
                if (Tdata>=cutlowT_manual) { // NOTE: >=
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            }
            else {
                if (Tdata>=cutlowT) { // NOTE: >=
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            }
        } readtaueqdata.close();
    }
    else {
        cout
        << "FitDielectric::read_individual_taueq_data: taueq datafile cannot open." << "\n";
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
        << "in FitDielectric::read_individual_taueq_data:" << "\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortDecreasing);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
        std::sort(taueq.begin(),taueq.end(),sortIncreasing);
        averaging_sorted_data(taueq,taueq_sortincreasing);
    }
}





void FitDielectric::write_fit_HN(std::ofstream& outputFile)
{
    if (outputFile.is_open()) {
        
        if (n_peaks>0) {
            
            outputFile
            << "cutoff_lo_afreq " << log10(cutoff_lo_afreq.at(whichpeak)) << "\n"
            << "cutoff_lo_sfreq " << log10(cutoff_lo_sfreq.at(whichpeak)) << "\n\n";
            
            outputFile
            << "cutoff_hi_afreq " << log10(cutoff_hi_afreq.at(whichpeak)) << "\n"
            << "cutoff_hi_sfreq " << log10(cutoff_hi_sfreq.at(whichpeak)) << "\n\n";
            
            outputFile
            << "afreq_max " << log10(afreq_max.at(whichpeak)) << "\n"
            << "sfreq_max " << log10(sfreq_max.at(whichpeak)) << "\n"
            << "loss_max  " << log10(loss_max.at(whichpeak))  << "\n\n";
            
            outputFile
            << "tauMax(afreq) " << HN_tauFit_calc(c,func)[0]        << "\n"
            << "log(tauMax_a) " << log10(HN_tauFit_calc(c,func)[0]) << "\n"
            << "tauMax(sfreq) " << HN_tauFit_calc(c,func)[1]        << "\n"
            << "log(tauMax_s) " << log10(HN_tauFit_calc(c,func)[1]) << "\n\n";
            
            
            if (whichpeak==0)
            {
                // Raw Data
                //--------------------------------------------------------------
                outputFile
                << "Raw Data (#.="<< dielectric_freq_overall.size() << ")\n"
                << "log.afreq  log.sfreq  log.loss  storage" << "\n"
                << "=======================================" << "\n";
                for (size_t i=0; i<dielectric_freq_overall.size(); ++i) {
                    outputFile
                    << log10(fabs(dielectric_freq_overall[i][0])) << " "
                    << log10(fabs(dielectric_freq_overall[i][1])) << " "
                    << log10(fabs(dielectric_freq_overall[i][2])) << " "
                    << dielectric_freq_overall[i][3]        << "\n";
                }
                outputFile << "\n\n";
                //--------------------------------------------------------------
            }
            
            
            
            // Fit Data (Fitting Range)
            //------------------------------------------------------------------
            outputFile
            << "Fit Data (Fitting Range)                  " << "\n"
            << "log.afreq  log.sfreq  log.loss  log.HN_fit" << "\n"
            << "==========================================" << "\n";
            double afreq_tmp=0,sfreq_tmp=0,loss_tmp=0;
            for (size_t i=0; i<dielectric_freq_sortincreasing[whichpeak].size(); ++i) {
                
                afreq_tmp=dielectric_freq_sortincreasing[whichpeak][i][0];
                sfreq_tmp=dielectric_freq_sortincreasing[whichpeak][i][2];
                loss_tmp =dielectric_freq_sortincreasing[whichpeak][i][1];
                
                outputFile
                << afreq_tmp << " "
                << sfreq_tmp << " "
                << loss_tmp  << " "
                << log10(fabs(HN_func_value(c,pow(10,afreq_tmp),func)))
                << "\n";
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            
            int    equal_pieces=1000;
            double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
            double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
            
            
            
            if (whichpeak==0)
            {
                // Interpolated Data (Full Range)
                //--------------------------------------------------------------
                outputFile
                << "Cubic Spline Interpolation (Full Range)" << "\n"
                << "log.afreq  log.sfreq  f  df  d2f       " << "\n"
                << "=======================================" << "\n";
                double f=0,df=0,d2f=0;
                for (int i=0; i<=equal_pieces; ++i) {
                    
                    afreqi=log_afreq_MIN+daf*i;
                    sfreqi=log_sfreq_MIN+dsf*i;
                    
                    alglib::spline1ddiff(interpolantLossCurve,afreqi,f,df,d2f);
                    
                    outputFile
                    << afreqi << " "
                    << sfreqi << " "
                    << f      << " "   // f is log(loss)
                    << df     << " "   // 1st derivative of log(loss) vs log(afreq)
                    << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
                }
                outputFile << "\n\n";
                //--------------------------------------------------------------
            }
            
            
            
            // Fit Data (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "Fit Data (Full Range)           " << "\n"
            << "log.afreq  log.sfreq  log.HN_fit" << "\n"
            << "================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << log10(fabs(HN_func_value(c,pow(10,afreqi),func))) << "\n";
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            
            double a=c[1],b=c[2],tau=pow(10,c[3]); // A=c[0]
            double t=0;
            double y=0;
            double theta=0;
            double omega=0;
            double f_t=0;
            
            
            
            // Distribution of tau (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "Distribution of tau (Full Range) " << "\n"
            << "log.time  log.f(time)            " << "\n"
            << "=================================" << "\n";
            for (int i=0; i<dielectric_freq_overall.size(); ++i)
            {
                t=pow(dielectric_freq_overall[i][1],-1);
                y=t/tau;
                theta=atan(sin(pi()*a)/(pow(y,a)+cos(pi()*a)));
                omega=pow(1+2*pow(y,a)*cos(pi()*a)+pow(y,2*a),-b/2);
                f_t=pow(pi(),-1)*pow(y,a*b)*sin(b*theta)*omega;
                
                outputFile
                << log10(t)      << " "
                << log(fabs(f_t)) << "\n";
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
            
        } else {
            
            // Raw Data
            //------------------------------------------------------------------
            outputFile
            << "Raw Data                               " << "\n"
            << "log.afreq  log.sfreq  log.loss  storage" << "\n"
            << "=======================================" << "\n";
            for (size_t i=0; i<dielectric_freq_overall.size(); ++i)
            {
                outputFile
                << log10(fabs(dielectric_freq_overall[i][0])) << " "
                << log10(fabs(dielectric_freq_overall[i][1])) << " "
                << log10(fabs(dielectric_freq_overall[i][2])) << " "
                << dielectric_freq_overall[i][3]        << "\n";
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            int    equal_pieces=1000;
            double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
            double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
            
            
            // Interpolated Data (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "Interpolated Data (Full Range)   " << "\n"
            << "log.afreq  log.sfreq  f  df  d2f " << "\n"
            << "=================================" << "\n";
            double f=0,df=0,d2f=0;
            for (int i=0; i<=equal_pieces; ++i)
            {
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                
                alglib::spline1ddiff(interpolantLossCurve,afreqi,f,df,d2f);
                
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << f      << " "   // f is log(loss)
                << df     << " "   // 1st derivative of log(loss) vs log(afreq)
                << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
        }
    }
    else {
        cout
        << "FitDielectric::write_fit_HN: file cannot open."
        << "\n";
    }
}





void FitDielectric::write_dc_loss(std::ofstream& outputFile)
{
    if (dc_curve_freq_sortincreasing.size()>0) {
        
        calc_derivatives(dc_curve_freq_sortincreasing,0,1);
        interpolantDCLoss=interpolant;
        
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        size_t dcpoints=dc_curve_freq_sortincreasing.size();
        
        outputFile
        << "DC Loss Curve Fit                      " << "\n"
        << "sig/eps  s  r2                         " << "\n"
        << "=======================================" << "\n";
        for (int i=0; i<c_dc_loss.length(); ++i) outputFile << c_dc_loss[i] << " ";
        outputFile << rep_dc_loss.r2 << "\n\n";
        
        outputFile
        << "Fit Range                              " << "\n"
        << "=======================================" << "\n";
        for (size_t i=0; i<dcpoints; ++i) {
            outputFile
            << dc_curve_freq_sortincreasing[i][0] << " "
            << dc_curve_freq_sortincreasing[i][1] << "\n";
        }
        outputFile << "\n"
        //<< "log.afreq  log.sfreq  f_fit  f  df     " << "\n"
        << "Full Range                        " << "\n"
        << "log.afreq  log.sfreq  f_fit       " << "\n"
        << "==================================" << "\n";
        double f=0,df=0,d2f=0;
        for (int i=0; i<=equal_pieces; ++i) {
            
            afreqi=log_afreq_MIN+daf*i;
            sfreqi=log_sfreq_MIN+dsf*i;
            
            alglib::spline1ddiff(interpolantDCLoss,afreqi,f,df,d2f);
            
            outputFile
            << afreqi << " "
            << sfreqi << " "
            << log10(dc_loss_value(c_dc_loss,pow(10,afreqi))) << "\n";
            //<< f      << " "   // f is log(loss)
            //<< df     << "\n"; // 1st derivative of log(loss) vs log(afreq)
        }
        outputFile << "\n\n";
    }
    else {
        outputFile << "No dc (conductivity) contribution detected." << "\n";
        outputFile << "\n\n";
    }
}





void FitDielectric::write_master_curve(std::ofstream& outputFile)
{
    if (dielectric_freq_overall_reduced.size()>0) {
        
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        outputFile
        << "Reduced Curve of Loss Part (Full Range)" << "\n"
        << "log.afreq  log.sfreq  f  df  d2f       " << "\n"
        << "=======================================" << "\n";
        double f=0,df=0,d2f=0;
        for (int i=0; i<=equal_pieces; ++i) {
            
            afreqi=log_afreq_MIN+daf*i;
            sfreqi=log_sfreq_MIN+dsf*i;
            
            alglib::spline1ddiff(interpolantReducedLoss,afreqi,f,df,d2f);
            
            outputFile
            << afreqi << " "
            << sfreqi << " "
            << f      << " "   // f is log(loss)
            << df     << " "   // 1st derivative of log(loss) vs log(afreq)
            << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
        }
        outputFile << "\n\n";
    }
    
    if (masterloss.size()>0) {
        
        calc_derivatives(masterloss,0,1);
        interpolantMasterCurve=interpolant;
        
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        outputFile
        << "Master Curve of Loss Part (Full Range) " << "\n"
        << "log.afreq  log.sfreq  f  df  d2f       " << "\n"
        << "=======================================" << "\n";
        double f=0,df=0,d2f=0;
        for (int i=0; i<=equal_pieces; ++i) {
            
            afreqi=log_afreq_MIN+daf*i;
            sfreqi=log_sfreq_MIN+dsf*i;
            
            alglib::spline1ddiff(interpolantMasterCurve,afreqi,f,df,d2f);
            
            outputFile
            << afreqi << " "
            << sfreqi << " "
            << f      << " "   // f is log(loss)
            << df     << " "   // 1st derivative of log(loss) vs log(afreq)
            << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
        }
        outputFile << "\n\n";
    }
}





void FitDielectric::write_HN_params(std::ofstream& outputFile)
{
    outputFile << T_actual << " ";
    for (int i=0; i<c.length(); ++i) outputFile << c[i] << " ";
    outputFile << rep.r2 << "\n";
}





void FitDielectric::write_peak_info(std::ofstream& outputFile)
{
    if (whichpeak==0) {
        outputFile
        << "total # of peaks = " << n_peaks
        << "\n\n"
        << "From low frequency, the #" << whichpeak+1 << " peak"
        << "\n\n";
    }
    else {
        outputFile
        << "\n\n"
        << "From low frequency, the #" << whichpeak+1 << " peak"
        << "\n\n";
    }
}





void FitDielectric::write_taueqToFile(std::ofstream& taueq,
                                      std::ofstream& taueqinv,
                                      std::ofstream& alltaueq,
                                      std::ofstream& alltaueqinv)
{
    /* write tau_eq data to separate files (cooling or heating) */
    taueq       << T_actual << " ";
    taueqinv    << (1000.0/(corrFacforT/precision))/(T_actual) << " ";
    alltaueq    << T_actual << " ";
    alltaueqinv << (1000.0/(corrFacforT/precision))/(T_actual) << " ";
    for (int i=0; i<c_losspeaks.size(); ++i) {
        taueq       << log10(HN_tauFit_calc(c_losspeaks[i],func)[1]) << " ";
        taueqinv    << log10(HN_tauFit_calc(c_losspeaks[i],func)[1]) << " ";
        alltaueq    << log10(HN_tauFit_calc(c_losspeaks[i],func)[1]) << " ";
        alltaueqinv << log10(HN_tauFit_calc(c_losspeaks[i],func)[1]) << " ";
    }
    taueq       << "\n";
    taueqinv    << "\n";
    alltaueq    << "\n";
    alltaueqinv << "\n";
}





void FitDielectric::write_fitcutoff(std::ofstream& outputFile)
{
    outputFile
    << "fit temperature range in [" << Theq << "," << cutlowT << "]"
    << "\n\n";
}





vector<double> FitDielectric::compute_Tg_fragility(StructureClass& sysVar,
                                                   const alglib::real_1d_array& c,
                                                   const std::string& model)
{
    double Tg_extrp=calc_temp_given_time(c,extrp_time,model)[0];
    double m_extrp =calc_temp_given_time(c,extrp_time,model)[1];
    return {Tg_extrp,m_extrp};
}





double FitDielectric::calc_time_given_temp(const real_1d_array& c,
                                           const double x,
                                           const string& model)
{
    if (model=="VFT") {
        // tau = tau0*exp((D*T0)/(T-T0))
        // param: {pow(10,c[0])=tau0, c[1]=D, c[2]=T0}
        return pow(10,c[0])*exp((c[1]*c[2])/(x-c[2]));
    }
    else if (model=="COOP") {
        // tau = tau0*exp(E_inf*(1+exp(-u*((T/E_inf)-b)))/T)
        // param:{pow(10,c[0])=tau0, c[1]=E_inf, c[2]=u, c[3]=b}
        return pow(10,c[0])*exp((c[1]*(1+exp(-c[2]*((x/c[1])-c[3]))))/x);
    }
    else
        return 0;
}





vector<double> FitDielectric::calc_temp_given_time(const real_1d_array& c,
                                                   const double time,
                                                   const string& model)
{
    vector<double> glass_data;
    
    double lnTime=log(time/pow(10,c[0])); // F=ln(time/tau0)
    
    if (model=="VFT")
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
    else if (model=="COOP")
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
    else
    {
        cout
        << "in FitDielectric::calc_temp_given_time, " << model << " model not found!"
        << "\n"
        << "Program aborted."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    return glass_data;
}





void FitDielectric::write_Tg_fragility(std::ofstream& outputFile,
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
        //double Tg_compu = calc_temp_given_time(c,compu_time,model)[0];
        //double m_compu  = calc_temp_given_time(c,compu_time,model)[1];
        //------------------------------------------------------------------
        
        outputFile
        << "-------------------------------------" << "\n"
        << "Tg_extrp [@("<<extrp_time<<")timeUnit] "
        << Tg_extrp
        << "\n"
        << "m_extrp [@Tg_extrp] "
        << m_extrp
        << "\n"
        << "-------------------------------------" << "\n"
        << "\n";
        
    }
    else {
        cout
        << "FitDielectric::write_Tg_fragility: file cannot open."
        << "\n";
    }
}





void FitDielectric::write_tauFitfit_vs_T(std::ofstream& outputFile,
                                         const alglib::real_1d_array& c,
                                         const std::string& model)
{
    if (outputFile.is_open()) {
        
        double T_initial=Theq;
        
        int equal_pieces=1000;
        double Ti=0;
        
        /* format: T log10(tauFit_fit) */
        //------------------------------------------------------------------
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
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitDielectric::write_tauFitfit_vs_T: file cannot open."
        << "\n";
    }
    
}





void FitDielectric::write_tauFitfit_vs_invT(std::ofstream& outputFile,
                                            const alglib::real_1d_array& c,
                                            const std::string& model)
{
    if (outputFile.is_open()) {
        
        double T_initial=Theq;
        
        int equal_pieces=1000;
        double Ti=0;
        
        /* format: 1/T log10(tauFit_fit) */
        //------------------------------------------------------------------
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
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "FitDielectric::write_tauFitfit_vs_invT: file cannot open."
        << "\n";
    }
}





void FitDielectric::get_Theq(const StructureClass& sysVar)
{
    read_individual_taueq_data(sysVar);
    Theq=taueq_sortdecreasing[0][0];
}





void FitDielectric::alglib_solverobj(const alglib::real_2d_array& x,
                                     const alglib::real_1d_array& y,
                                     const alglib::real_1d_array& c,
                                     alglib::lsfitstate& state)
{
    if (get_is_use_FG())
    {
        lsfitcreatefg(x, y, c, true, state);
    }
    else
    {
        lsfitcreatef(x, y, c, diffstep, state);
    }
}





void FitDielectric::alglib_lsfit_form(const std::string& model,
                                      alglib::lsfitstate& state)
{
    if (model=="HN_loss")
    {
        if (func_HN=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (func_HN=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
    else if (model=="CC_loss")
    {
        if (func_HN=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (func_HN=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
    else if (model=="CD_loss")
    {
        if (func_HN=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (func_HN=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
    else if (model=="dc_loss")
    {
        alglib::lsfitfit(state, dc_loss);
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
    else if (model=="COOP")
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
}





void FitDielectric::set_fitParams(const std::string& model)
{
    epsf     = 0;
    epsx     = 1e-10;
    maxits   = 1e+5;
    diffstep = 1e-5;
    
    if (model=="HN_loss") // [A,a,b,tau]
    {
        // HN: shape params are free
        
        // NOTE:
        // tau uses log scale
        
        set_coeffs({       1.0,  0.5,  0.5,  -8.0});
        set_coeffs_scale({ 1.0,  0.5,  0.5,  -8.0});
        set_coeffs_bndl({  0.0,  0.0,  0.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
    }
    else if (model=="CC_loss") // [A,a,b,tau]
    {
        // Cole-Cole: special case of HN when beta=1
        
        // NOTE:
        // tau uses log scale
        
        set_coeffs({       1.0,  0.5,  1.0,  -8.0});
        set_coeffs_scale({ 1.0,  0.5,  1.0,  -8.0});
        set_coeffs_bndl({  0.0,  0.0,  1.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
    }
    else if (model=="CD_loss") // [A,a,b,tau]
    {
        // Cole-Davidson: special case of HN when alpha=1
        
        // NOTE:
        // tau uses log scale
        
        set_coeffs({       1.0,  1.0,  0.5,  -8.0});
        set_coeffs_scale({ 1.0,  1.0,  0.5,  -8.0});
        set_coeffs_bndl({  0.0,  1.0,  0.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
    }
    else if (model=="dc_loss") // [A,s]
    {
        // NOTE:
        // A represents sig0/eps0
        
        set_coeffs({       1.0,  1.0});
        set_coeffs_scale({ 1.0,  1.0});
        set_coeffs_bndl({  0.0, -inf});
        set_coeffs_bndu({  inf,  inf});
    }
    else if (model=="VFT") // [tau0,D,T0]
    {
        // tau0 uses log scale
        
        set_coeffs({      -12.0,  50.0, 1e+2});
        set_coeffs_scale({-12.0,  10.0, 1e+2});
        set_coeffs_bndl({ -20.0,  0.0,   0.0});
        set_coeffs_bndu({  20.0, 1e+3,  1e+3});
    }
    else if (model=="COOP") // [tau0,E_inf,u,b]
    {
        // tau0 uses log scale
        
        set_coeffs({      -12.0, 1e+3,   1.0,   0.1});
        set_coeffs_scale({-12.0, 1e+3,   1.0,   0.1});
        set_coeffs_bndl({ -20.0,  0.0, -5e+2, -1e+3});
        set_coeffs_bndu({  20.0, 1e+6,  5e+2,  1e+3});
    }
    else
    {
        cout
        << "in FitDielectric::set_fitParams, " << model << " model not found!"
        << "\n"
        << "Program aborted."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}



void FitDielectric::fitParams_correction(const std::string& model)
{
    if ((model=="HN_loss")||(model=="CC_loss")||(model=="CD_loss"))
    {
        if (cyclecount>0)
        {
            coeffs_vD[3] += 0.1; // NOTE: log scale
            set_coeffs
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            set_coeffs_scale
            ({coeffs_vD[0],coeffs_vD[1],coeffs_vD[2],coeffs_vD[3]});
            c = get_coeffs().c_str();
        }
    }
}





void FitDielectric::alglib_lsfit(const std::string& model)
{
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
        
    } while((rep.r2<r2_threshold)&&(cyclecount<100));
    usedt.push_back(time(NULL));
    //======================================================================
    ////////////////////////////////////////////////////////////////////////
}





void FitDielectric::HN_loss_wiki(const alglib::real_1d_array &c,
                                 const alglib::real_1d_array &x,
                                 double &func,
                                 void *ptr)
{
    // NOTE: A represents the dielectric strength
    // source: wiki
    double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]),afreq=pow(10,x[0]);
    double phi=atan((pow(afreq*tau,a)*sin(pi()*a/2))/(1+pow(afreq*tau,a)*cos(pi()*a/2)));
    func = log10(fabs(A*pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2)*sin(b*phi))); // log form
}





void FitDielectric::HN_loss_Kremer(const alglib::real_1d_array &c,
                                   const alglib::real_1d_array &x,
                                   double &func,
                                   void *ptr)
{
    // NOTE: A represents the dielectric strength
    // source: Kremer (Broadband Dielectric Spectroscopy) textbook Table 3.1
    double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]),afreq=pow(10,x[0]);
    double r_w=pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2);
    double p_w=atan(sin(pi()*a/2)/(pow(afreq*tau,-a)+cos(pi()*a/2)));
    func = log10(fabs(A*r_w*sin(b*p_w))); // log form
}





double FitDielectric::HN_func_value(const real_1d_array& c,
                                    const double afreq,
                                    const string& func)
{
    double Fvalue=0;
    if ((func=="HN_loss")||(func=="CC_loss")||(func=="CD_loss"))
    {
        if (func_HN=="wiki")
        {
            double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]);
            double phi=atan((pow(afreq*tau,a)*sin(pi()*a/2))/(1+pow(afreq*tau,a)*cos(pi()*a/2)));
            Fvalue=A*pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2)*sin(b*phi);
        }
        else if (func_HN=="Kremer")
        {
            double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]);
            double r_w=pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2);
            double p_w=atan(sin(pi()*a/2)/(pow(afreq*tau,-a)+cos(pi()*a/2)));
            Fvalue=A*r_w*sin(b*p_w);
        }
    }
    else
    {
        cout
        << "in FitDielectric::HN_func_value, " << func << " model not found!"
        << "\n"
        << "Program aborted."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    return Fvalue;
}





vector<double> FitDielectric::HN_tauFit_calc(const alglib::real_1d_array& c,
                                             const std::string& func)
{
    double tauFit_afreq=0,tauFit_sfreq=0;
    if ((func=="HN_loss")||(func=="CC_loss")||(func=="CD_loss"))
    {
        double a=c[1],b=c[2],tau=pow(10,c[3]); // A=c[0]
        //tauFit=tau*pow(sin((pi()*a*b)/(2*(1+b)))/sin((pi()*a)/(2*(1+b))),1/a);
        double fmax=pow(2*pi()*tau,-1)*pow(sin((pi()*a)/(2*(1+b)))/sin((pi()*a*b)/(2*(1+b))),1/a);
        tauFit_afreq=pow(2*pi()*fmax,-1); // afreq-based relaxation time
        tauFit_sfreq=pow(fmax,-1);        // sfreq-based relaxation time
    }
    else
    {
        cout
        << "in FitDielectric::HN_tauFit_calc, " << func << " model not found!"
        << "\n"
        << "Program aborted."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    return {tauFit_afreq,tauFit_sfreq};
}





void FitDielectric::dc_loss(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr)
{
    double A=c[0],s=c[1],afreq=pow(10,x[0]);
    func = log10(fabs(A*pow(afreq,-s))); // log form
}





double FitDielectric::dc_loss_value(const alglib::real_1d_array& c,
                                    const double afreq)
{
    double A=c[0],s=c[1];
    return A*pow(afreq,-s);
}





void FitDielectric::VFT_func(const real_1d_array &c,
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
    double tau0=pow(10,c[0]),D=c[1],T0=c[2],T=x[0];
    func = (log(tau0)+((D*T0)/(T-T0)))*log10(exp(1)); // log form
}
void FitDielectric::VFT_grad(const real_1d_array &c,
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
    double tau0=pow(10,c[0]),D=c[1],T0=c[2],T=x[0];
    func = (log(tau0)+((D*T0)/(T-T0)))*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    grad[1]= (T0/(T-T0))*log10(exp(1));
    grad[2]= ((D*T)/pow(T-T0,2))*log10(exp(1));
}





void FitDielectric::COOP_func(const alglib::real_1d_array &c,
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
    double tau0=pow(10,c[0]),E_inf=c[1],u=c[2],b=c[3],T=x[0];
    double exponent=exp(-u*((T/E_inf)-b));
    double E_coop=E_inf*exponent;
    double E=E_inf+E_coop;
    func = (log(tau0)+E/T)*log10(exp(1)); // log form
}
void FitDielectric::COOP_grad(const alglib::real_1d_array &c,
                              const alglib::real_1d_array &x,
                              double &func,
                              alglib::real_1d_array &grad,
                              void *ptr)
{
    double tau0=pow(10,c[0]),E_inf=c[1],u=c[2],b=c[3],T=x[0];
    double exponent=exp(-u*((T/E_inf)-b));
    double E_coop=E_inf*exponent;
    double E=E_inf+E_coop;
    func = (log(tau0)+E/T)*log10(exp(1));
    grad[0]= (pow(tau0,-1))*log10(exp(1));
    cout << grad[0] << "\n";
    grad[1]= ((((T*u*exponent)/E_inf)+exponent+1)/T)*log10(exp(1));
    grad[2]= ((E_inf*(b-T/E_inf)*exponent)/T)*log10(exp(1));
    grad[3]= ((E_inf*u*exponent)/T)*log10(exp(1));
}




void FitDielectric::initialize()
{
    if(is_manualFitSingleDS)
    {
        form="manuallySet";
        set_is_log_data({true,false});
    }
}





void FitDielectric::make_dielectric_folder(const StructureClass& sysVar)
{
    static int countaccess=0;
    ++countaccess;
    
    if (countaccess==1) {
        
        /* create folder for dielectric data */
        string dir,mkdir;
        dir.append(return_AnalysisFolderPath(sysVar));
        dir.append("/statistics/dielectric");
        mkdir="mkdir "+dir;
        system(mkdir.c_str());
        
        /* copy data from source to destination folder */
        if (form=="NPIC") {
            //
        } else if (form=="NIST") {
            string path,cp;
            path.append(sysVar.get_Path());
            path.append("/dielectric");
            cp="cp "+path+"/* "+dir;
            system(cp.c_str());
        } else if (form=="manuallySet") {
            string path,cp;
            path.append(sysVar.get_Path());
            path.append("/dielectric");
            cp="cp "+path+"/* "+dir;
            system(cp.c_str());
        }
        
        /* NOTE:
         naming rule: {species}_*  */
        
        string initial,final;
        string bashfile="mv_bash.txt";
        ofstream mv_bash(bashfile.c_str());
        
        if (form=="NPIC") {
            //
        } else if (form=="NIST") {
            if (false) {
                mv_bash
                << "cd "+dir
                << "\n\n"
                << "j=0                     " << "\n"
                << "for i in *              " << "\n"
                << "do                      " << "\n"
                << "  j=$[$j+1]             " << "\n"
                << "  mv $i ${i%%_*}_$j.dat " << "\n"
                << "done                    " << "\n";
            } else {
                mv_bash
                << "cd "+dir
                << "\n\n"
                << "for i in *                          " << "\n"
                << "do                                  " << "\n"
                << "  mv $i "<< sysVar.get_usic()<<".dat" << "\n"
                << "done                                " << "\n";
            }
        } else if (form=="manuallySet") {
            mv_bash
            << "cd "+dir
            << "\n\n"
            << "for i in *                          " << "\n"
            << "do                                  " << "\n"
            << "  mv $i "<< sysVar.get_usic()<<".dat" << "\n"
            << "done                                " << "\n";
        } mv_bash.close();
        string usebash="bash ./"+bashfile;
        system(usebash.c_str());
        string rmbash="rm "+bashfile;
        system(rmbash.c_str());
    }
}





/* public setters */

/* string */
void FitDielectric::set_form(const std::string& str){form=str;}
void FitDielectric::set_func(const std::string& str){func=str;}
void FitDielectric::set_model(const std::string& str){model=str;};
/* bool */
void FitDielectric::set_is_use_dc_free(const bool b){is_use_dc_free=b;}
void FitDielectric::set_is_manual_cutlowT(const bool b){is_manual_cutlowT=b;}
void FitDielectric::set_is_manualFitSingleDS(const bool b){is_manualFitSingleDS=b;}
void FitDielectric::set_is_log_data(const std::vector<bool>& vb){is_log_data=vb;}
/* double */
void FitDielectric::set_cutlowT_manual(const double d){cutlowT_manual=d;}

