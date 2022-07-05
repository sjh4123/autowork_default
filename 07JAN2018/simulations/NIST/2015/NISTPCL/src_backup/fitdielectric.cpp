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
dielectric_format("NIST"),
source("Kremer"),

/* bool */
is_use_dc_free(true),

/* int */
n_fit_HN(3),
n_used_HN(100),
n_check(3),
n_each_side(100),

/* double */
r2_threshold(0.90),
frac_threshold(0.25)
{
    is_use_FG   = false;
    
    corrFacforT = sysVar.get_corrFacforT();
    precision   = sysVar.get_precision();
}





vector<vector<double>> FitDielectric::dielectricData_preprocess(const StructureClass& sysVar,
                                                                const string& form)
{
    vector<vector<double>> tinfo_dielectric;
    vector<double> tinfo_cooling;
    vector<double> tinfo_heating;
    
    if (form=="NPIC")
    {
        //
    }
    else if (form=="NIST")
    {
        /* input file */
        string i;
        i.append(return_AnalysisFolderPath(sysVar));
        i.append("/statistics/dielectric");
        i.append("/"+sysVar.get_usic());
        i.append(".dat");
        
        /* output file */
        string o;
        
        ifstream readfile(i.c_str());
        if (readfile.is_open())
        {
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
            
            while (getline(readfile,lineContent))
            {
                ++lineindex;
                istringstream iss(lineContent);
                
                iss >> sfreq;
                iss >> T_now; T_now+=273;
                iss >> epsp;
                iss >> epspp;
                
                if (lineindex>1)
                {
                    /* store data at the end of each temperature block */
                    if (T_now!=T_pre)
                    {
                        is_in_cooling=false;
                        
                        Temp_d=std::trunc(T_pre*corrFacforT);
                        
                        /* cooling process */
                        
                        if (T_now<T_pre)
                        {
                            is_in_cooling=true;
                            tinfo_cooling.push_back(Temp_d);
                            
                            o.clear();
                            o.append(return_AnalysisFolderPath(sysVar));
                            o.append("/statistics/dielectric");
                            o.append("/cooling_eps_");
                            o.append(sysVar.get_usic());
                            o.append("_T"+to_string((long long int)Temp_d));
                            o.append(".dat");
                        }
                        
                        /* heating process */
                        
                        else if (T_now>T_pre)
                        {
                            is_in_cooling=false;
                            tinfo_heating.push_back(Temp_d);
                            
                            o.clear();
                            o.append(return_AnalysisFolderPath(sysVar));
                            o.append("/statistics/dielectric");
                            o.append("/heating_eps_");
                            o.append(sysVar.get_usic());
                            o.append("_T"+to_string((long long int)Temp_d));
                            o.append(".dat");
                        }
                        
                        ofstream writefile(o.c_str());
                        for (int i=0; i<tnow_vd.size(); ++i)
                        {
                            writefile
                            << 2*pi()*sfreq_vd[i] << " "    // afreq
                            << sfreq_vd[i]        << " "    // sfreq
                            << epspp_vd[i]        << " "    // loss
                            << epsp_vd[i]         << "\n";  // storage
                        }
                        writefile.close();
                        
                        sfreq_vd.clear();
                        tnow_vd.clear();
                        epsp_vd.clear();
                        epspp_vd.clear();
                    }
                }
                T_pre=T_now;
                
                sfreq_vd.push_back(sfreq);
                epsp_vd.push_back(epsp);
                epspp_vd.push_back(epspp);
                tnow_vd.push_back(T_now);
            }
            
            /* store data of the last T block */
            Temp_d=std::trunc(T_now*corrFacforT);
            
            if (is_in_cooling)
            {
                tinfo_cooling.push_back(Temp_d);
                o.clear();
                o.append(return_AnalysisFolderPath(sysVar));
                o.append("/statistics/dielectric");
                o.append("/cooling_eps_");
                o.append(sysVar.get_usic());
                o.append("_T"+to_string((long long int)Temp_d));
                o.append(".dat");
            }
            else
            {
                tinfo_heating.push_back(Temp_d);
                o.clear();
                o.append(return_AnalysisFolderPath(sysVar));
                o.append("/statistics/dielectric");
                o.append("/heating_eps_");
                o.append(sysVar.get_usic());
                o.append("_T"+to_string((long long int)Temp_d));
                o.append(".dat");
            }
            
            ofstream writefile(o.c_str());
            for (int i=0; i<tnow_vd.size(); ++i)
            {
                writefile
                << 2*pi()*sfreq_vd[i] << " "    // afreq
                << sfreq_vd[i]        << " "    // sfreq
                << epspp_vd[i]        << " "    // loss
                << epsp_vd[i]         << "\n";  // storage
            }
            writefile.close();
            
            
            if (tinfo_cooling.size()!=0)
            {
                tinfo_dielectric.push_back(tinfo_cooling);
            }
            if (tinfo_heating.size()!=0)
            {
                tinfo_dielectric.push_back(tinfo_heating);
            }
        }
        else
        {
            cout
            << "in dielectricData_preprocess: dielectric file cannot open."
            << "\n";
        }
    }
    
    return tinfo_dielectric;
}





void FitDielectric::read_all_dielectric_data(const StructureClass& sysVar,
                                             const int process,
                                             const double Temp_d,
                                             const string& form)
{
    dielectric_freq_sortincreasing.clear();
    dielectric_freq_overall.clear();
    
    vector<vector<double>> dielectric_freq;
    
    if (form=="NPIC")
    {
        //
    }
    else if (form=="NIST")
    {
        string i;
        if (process==0)
        {
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/statistics/dielectric/");
            i.append("/cooling_eps_");
            i.append(sysVar.get_usic());
            i.append("_T"+to_string((long long int)Temp_d));
            i.append(".dat");
            //cout << i << "\n";
        }
        else if (process==1)
        {
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/statistics/dielectric/");
            i.append("/heating_eps_");
            i.append(sysVar.get_usic());
            i.append("_T"+to_string((long long int)Temp_d));
            i.append(".dat");
            //cout << i << "\n";
        }
        
        ifstream readFile(i.c_str());
        if (readFile.is_open())
        {
            string lineContent;
            string strData;
            
            while (getline(readFile,lineContent))
            {
                istringstream iss(lineContent);
                
                iss >> afreq;
                iss >> sfreq;
                iss >> loss;
                iss >> storage;
                
                dielectric_freq.push_back({afreq,sfreq,loss,storage});
            }
            readFile.close();
        }
        else
        {
            cout
            << "read_all_dielectric_data: dielectric file cannot open."
            << "\n";
        }
        try
        {
            if (dielectric_freq.size()==0) throw 0;
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
        }
        catch (int i)
        {
            cout
            << "in read_all_dielectric_data:" << "\n"
            << "dielectric_freq size = " << i << "\n";
            dielectric_freq.push_back({0,0});
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
        }
        
        
        afreq_MIN=dielectric_freq[0][0];
        sfreq_MIN=dielectric_freq[0][1];
        afreq_MAX=dielectric_freq[dielectric_freq.size()-1][0];
        sfreq_MAX=dielectric_freq[dielectric_freq.size()-1][1];
        
        log_afreq_MIN=log10(afreq_MIN);
        log_afreq_MAX=log10(afreq_MAX);
        log_sfreq_MIN=log10(sfreq_MIN);
        log_sfreq_MAX=log10(sfreq_MAX);
        
        
        /* put read-in raw data into proper containers */
        //----------------------------------------------------------------------
        if (is_use_dc_free)
        {
            vector<vector<double>> dielectric_dc_free;
            
            /* establish cubic spline interpolation routine */
            vector<vector<double>> interpolatingdata;
            for (int i=0; i<dielectric_freq.size(); ++i)
            {
                interpolatingdata.push_back
                ({
                    log(abs(dielectric_freq[i][0])), // NOTE: ln(afreq)
                    dielectric_freq[i][3]            // NOTE: storage
                });
            }
            fitdata_processing_1d(interpolatingdata);
            real_1d_array x_in=fit_xData.c_str();
            real_1d_array y_in=fit_yData.c_str();
            alglib::spline1dinterpolant s;
            alglib::spline1dbuildcubic(x_in,y_in,s);
            
            double afreqi=0,sfreqi=0,f=0,df=0,d2f=0;
            double storage_dc=0,loss_dc=0;
            
            for (int i=0; i<dielectric_freq.size(); ++i)
            {
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
        else
        {
            dielectric_freq_overall=dielectric_freq;
        }
        //----------------------------------------------------------------------
        // now should get data structures containing all raw data:
        // dielectric_freq_sortincreasing
        // dielectric_freq_overall
        
        
        /* calculate 1st and 2nd derivatives of interpolated eps" curve */
        calc_loss_derivatives();
        
        
        /* sample data to be fit in the proximity of the peaks */
        peak_sampling_bySlope();
        //peak_sampling_byValue();
        
        
        /*----------------------------------------------------------------------
         Final Product Containers:
         
         "dielectric_freq_overall"
         --------------------------------
         <I dimension>
         individual data points on the entire freq range
         
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
         individual data points around each peak
         
         <III dimension>
         dielectric data:
         0: log(afreq)
         1: log(loss)
         2: log(sfreq)
         3: storage
         ---------------------------------------------------------------------*/
    }
}





void FitDielectric::peak_sampling_byValue()
{
    
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
        << "in read_all_dielectric_data:" << "\n"
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
        << "in read_all_dielectric_data:" << "\n"
        << "dielectric_freq size = " << i << "\n";
        dielectric_freq.push_back({0,0});
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
    }
    dielectric_freq_peak.clear();
    dielectric_freq_peak=dielectric_freq;
    
    
    dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
}





void FitDielectric::peak_sampling_bySlope()
{
    
    vector<int> peak_indices;
    vector<vector<double>> dielectric_freq;
    vector<vector<double>> dielectric_freq_peak;
    
    
    /* 1. Distinguish Local Peaks */
    //--------------------------------------------------------------------------
    int    peak_index=0;
    int    n_total=n_check;
    
    double afreq_now=0,afreq_pre=0;
    double f_now=0,df_now=0,d2f_now=0;
    double f_pre=0,df_pre=0,d2f_pre=0;
    
    
    // total # of data points
    size_t n_points=dielectric_freq_overall.size();
    
    
    for (int index=1; index<n_points; ++index) // NOTE: start from 1
    {
        afreq_now=log10(abs(dielectric_freq_overall[index][0]));
        afreq_pre=log10(abs(dielectric_freq_overall[index-1][0]));
        alglib::spline1ddiff(losscurve,afreq_now,f_now,df_now,d2f_now);
        alglib::spline1ddiff(losscurve,afreq_pre,f_pre,df_pre,d2f_pre);
        
        if (df_pre*df_now<0)
        {
            
            if (dielectric_freq_overall[index][2]>
                dielectric_freq_overall[index-1][2])
            {
                peak_index=index;
            }
            else
            {
                peak_index=index-1;
            }
            
            /* peak index should not repeat */
            //------------------------------------------------------------------
            // "find" returns iterators in the range: [first,last); it would
            // return a pointer to the "last" element if nothing is found.
            //------------------------------------------------------------------
            vector<int>::iterator itr0;
            itr0=find(peak_indices.begin(),peak_indices.end(),peak_index);
            if (itr0!=peak_indices.end()) continue;
            
            /* peaks should not be too close */
            //------------------------------------------------------------------
            bool is_too_close=false;
            for (int ii=0; ii<peak_indices.size(); ++ii)
            {
                if (abs(peak_index-peak_indices[ii])<=n_check)
                {
                    is_too_close=true;
                    break;
                }
            }
            if (is_too_close) continue;
            //------------------------------------------------------------------
            
            int n_count=0;
            
            n_total=n_check;
            
            /* forward checking */
            for (int count=1; count<=n_check; ++count)
            {
                if ((peak_index+count)<n_points)
                {
                    // strict monotonic decreasing
                    bool is_valid=
                    dielectric_freq_overall[peak_index+count-1][2]>
                    dielectric_freq_overall[peak_index+count][2];
                    
                    if (is_valid)
                    {
                        ++n_count;
                    }
                }
            }
            
            /* backward checking */
            // NOTE:
            // if index is close to the ends, by fraction of "frac_threshold"
            // out of total # of data points, then do backward checking as well
            bool is_at_margin=
            abs(index-(int)(n_points-1))<=ceil(frac_threshold*n_points);
            
            //(abs(index-(int)(n_points-1))<=ceil(frac_threshold*n_points))&&
            //(index<=ceil(frac_threshold*n_points));
            
            if (is_at_margin)
            {
                n_total=n_check*2;
                
                /* backward checking */
                for (int count=1; count<=n_check; ++count)
                {
                    if ((peak_index-count+1)<n_points)
                    {
                        // strict monotonic decreasing
                        bool is_valid=
                        dielectric_freq_overall[peak_index-count+1][2]>
                        dielectric_freq_overall[peak_index-count][2];;
                        
                        if (is_valid)
                        {
                            ++n_count;
                        }
                    }
                }
            }
            
            if(n_count==n_total)
            {
                peak_indices.push_back(peak_index);
                
                // NOTE:
                // here store the original value
                afreq_max.push_back(dielectric_freq_overall[peak_index][0]);
                sfreq_max.push_back(dielectric_freq_overall[peak_index][1]);
                loss_max.push_back(dielectric_freq_overall[peak_index][2]);
            }
        }
    }
    
    
    /* 2. Data Sampling Around Local Peaks */
    //--------------------------------------------------------------------------
    for (int whichpeak=0; whichpeak<peak_indices.size(); ++whichpeak)
    {
        dielectric_freq.clear();
        dielectric_freq_peak.clear();
        
        /* the peak index */
        int peak_index=peak_indices[whichpeak];
        
        dielectric_freq.push_back
        ({
            dielectric_freq_overall[peak_index][0],
            dielectric_freq_overall[peak_index][1],
            dielectric_freq_overall[peak_index][2],
            dielectric_freq_overall[peak_index][3]
        });
        
        for (int side=0; side<2; ++side)
        {
            /* backward sampling */
            if (side==0)
            {
                int index_lo=0;
                for (int j=1; j<=n_each_side; ++j)
                {
                    bool is_valid=
                    ((peak_indices.at(whichpeak)-j)<n_points)&&
                    ((peak_indices.at(whichpeak)-j-1)>=0);
                    
                    if (is_valid)
                    {
                        index_lo=peak_indices.at(whichpeak)-j;
                        
                        bool is_at_margin=
                        index_lo<=ceil(frac_threshold*n_points);
                        
                        //(abs(index_lo-(int)(n_points-1))<=ceil(frac_threshold*n_points))&&
                        //((index_lo+1)<=ceil(frac_threshold*n_points));
                        
                        bool is_taking=false;
                        
                        afreq_now=log10(abs(dielectric_freq_overall.at(index_lo)[0]));
                        afreq_pre=log10(abs(dielectric_freq_overall.at(index_lo-1)[0]));
                        alglib::spline1ddiff(losscurve,afreq_now,f_now,df_now,d2f_now);
                        alglib::spline1ddiff(losscurve,afreq_pre,f_pre,df_pre,d2f_pre);
                        
                        if (is_at_margin)
                        {
                            //is_taking=(df_now<0.5)&&(df_pre<0.5);
                            
                            is_taking=
                            (is_taking=df_now*df_pre>0)&&
                            (abs(df_now-df_pre)<0.5);
                        }
                        else
                        {
                            is_taking=df_now*df_pre>0;
                        }
                        
                        if (is_taking)
                        {
                            dielectric_freq.push_back
                            ({
                                dielectric_freq_overall[index_lo][0],
                                dielectric_freq_overall[index_lo][1],
                                dielectric_freq_overall[index_lo][2],
                                dielectric_freq_overall[index_lo][3]
                            });
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
            
            /* forward sampling */
            else
            {
                int index_hi=0;
                for (int j=1; j<=n_each_side; ++j)
                {
                    bool is_valid=
                    ((peak_indices.at(whichpeak)+j+1)<n_points)&&
                    ((peak_indices.at(whichpeak)+j)>=0);
                    
                    if (is_valid)
                    {
                        index_hi=peak_indices.at(whichpeak)+j;
                        
                        afreq_now=log10(abs(dielectric_freq_overall.at(index_hi)[0]));
                        afreq_pre=log10(abs(dielectric_freq_overall.at(index_hi+1)[0]));
                        alglib::spline1ddiff(losscurve,afreq_now,f_now,df_now,d2f_now);
                        alglib::spline1ddiff(losscurve,afreq_pre,f_pre,df_pre,d2f_pre);
                        
                        bool is_turning=df_now*df_pre<0;
                        
                        if (!is_turning)
                        {
                            dielectric_freq.push_back
                            ({
                                dielectric_freq_overall[index_hi][0],
                                dielectric_freq_overall[index_hi][1],
                                dielectric_freq_overall[index_hi][2],
                                dielectric_freq_overall[index_hi][3]
                            });
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }
        try
        {
            if (dielectric_freq.size()==0) throw 0;
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
        }
        catch (int i)
        {
            cout
            << "in dc_free_peak_sampling:" << "\n"
            << "dielectric_freq size = " << i << "\n";
            dielectric_freq.push_back({0,0});
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing);
        }
        
        cutoff_lo_afreq.push_back(dielectric_freq[0][0]);
        cutoff_lo_sfreq.push_back(dielectric_freq[0][1]);
        cutoff_hi_afreq.push_back(dielectric_freq[dielectric_freq.size()-1][0]);
        cutoff_hi_sfreq.push_back(dielectric_freq[dielectric_freq.size()-1][1]);
        
        for (int i=0; i<dielectric_freq.size(); ++i)
        {
            dielectric_freq_peak.push_back
            ({
                log10(abs(dielectric_freq[i][0])), // log(afreq)
                log10(abs(dielectric_freq[i][2])), // log(loss)
                log10(abs(dielectric_freq[i][1])), // log(sfreq)
                dielectric_freq[i][3]              // storage
            });
        }
        
        dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
    }
}





void FitDielectric::calc_loss_derivatives()
{
    /* establish cubic spline interpolation routine */
    vector<vector<double>> interpolatingdata;
    for (int i=0; i<dielectric_freq_overall.size(); ++i)
    {
        interpolatingdata.push_back
        ({
            log10(abs(dielectric_freq_overall[i][0])), // NOTE: log(afreq)
            log10(abs(dielectric_freq_overall[i][2]))  // NOTE: log(loss)
        });
    }
    fitdata_processing_1d(interpolatingdata);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    alglib::spline1dbuildcubic(x_in,y_in,losscurve);
}





void FitDielectric::clear_containers()
{
    loss_max.clear();
    afreq_max.clear();
    sfreq_max.clear();
    cutoff_lo_afreq.clear();
    cutoff_lo_sfreq.clear();
    cutoff_hi_afreq.clear();
    cutoff_hi_sfreq.clear();
}





void FitDielectric::fit_HN(StructureClass& sysVar,
                           const int process,
                           const double Temp_d,
                           const std::string& func,
                           const std::string& form)
{
    dielectric_format=form;
    
    string o1,o2,o3,o4;
    
    if (process==0)
    {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/cooling_eps_");
        o1.append("taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        //cout << o1 << "\n";
        
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/cooling_eps_");
        o2.append("taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        //cout << o2 << "\n";
        
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/cooling_eps_");
        o3.append("fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        //cout << o << "\n";
        
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/cooling_eps_");
        o4.append("fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)Temp_d));
        o4.append(".dat");
        //cout << o4 << "\n";
    }
    else if (process==1)
    {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/heating_eps_");
        o1.append("taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        //cout << o1 << "\n";
        
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/heating_eps_");
        o2.append("taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        //cout << o2 << "\n";
        
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/heating_eps_");
        o3.append("fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        //cout << o3 << "\n";
        
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/heating_eps_");
        o4.append("fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)Temp_d));
        o4.append(".dat");
        //cout << o4 << "\n";
    }
    
    ofstream write_taueq(o1.c_str(),ofstream::app);    // 'append'
    ofstream write_taueqinv(o2.c_str(),ofstream::app); // 'append'
    ofstream write_params(o3.c_str(),ofstream::app);   // 'append'
    ofstream write_HN(o4.c_str(),ofstream::app);       // 'append'
    
    
    /* Clear containers before each use */
    //--------------------------------------------------------------------------
    clear_containers();
    
    
    /* Read in data to be fit at each temperature */
    //--------------------------------------------------------------------------
    read_all_dielectric_data(sysVar,process,Temp_d,form);
    
    
    if (true)
    {
        cout
        << Temp_d << " "
        << dielectric_freq_sortincreasing.size() << " ";
        for (int i=0; i<dielectric_freq_sortincreasing.size(); ++i)
        {
            cout << afreq_max[i] << " ";
        }
        cout << "\n";
    }
    
    
    for (int whichpeak=0; whichpeak<dielectric_freq_sortincreasing.size(); ++whichpeak)
    {
        /* Data processing to ALGLIB readable form */
        //----------------------------------------------------------------------
        fitdata_processing_lsfit(dielectric_freq_sortincreasing[whichpeak]);
        
        
        /* Curve-fitting Kernel */
        //----------------------------------------------------------------------
        set_fitParams(func);
        alglib_lsfit(func);
        write_peak_info(write_HN,dielectric_freq_sortincreasing.size(),whichpeak+1);
        write_fitModel(write_HN,func);
        write_fitarrays(write_HN);
        write_stopcond(write_HN);
        write_fitinfo(write_HN);
        write_fit_HN(write_HN,func,whichpeak);
        write_HN_params(write_params,Temp_d);
        
        
        if (rep.r2<r2_threshold)
        {
            write_badR2(write_HN,r2_threshold);
        }
        if (dielectric_freq_sortincreasing[whichpeak].size()<n_fit_HN)
        {
            write_insufficientDara
            (write_HN,dielectric_freq_sortincreasing[whichpeak].size(),n_fit_HN);
        }
        //if ((rep.r2<r2_threshold)||(dielectric_freq_sortincreasing.size()<n_fit_HN)) return;
        
        
        /* Write tau_alpha to file */
        //----------------------------------------------------------------------
        if (whichpeak==(dielectric_freq_sortincreasing.size()-1))
        {
            double T_actual=Temp_d*pow(corrFacforT,-1); // NOTE: convert!
            
            /* write tau_eq data to file */
            write_taueq
            << T_actual
            << " "
            << log10(HN_tauFit_calc(c,func))
            << "\n";
            
            write_taueqinv
            << (1000.0/(corrFacforT/precision))/(T_actual)
            << " "
            << log10(HN_tauFit_calc(c,func))
            << "\n";
        }
    }
    
    write_taueq.close();
    write_taueqinv.close();
    write_params.close();
    write_HN.close();
}





double FitDielectric::HN_func_value(const real_1d_array& c,
                                    const double afreq,
                                    const string& tcalc)
{
    double Fvalue=0;
    
    if ((tcalc=="HN_loss")||(tcalc=="CC_loss")||(tcalc=="CD_loss"))
    {
        if (source=="wiki")
        {
            double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]);
            double phi=atan((pow(afreq*tau,a)*sin(pi()*a/2))/(1+pow(afreq*tau,a)*cos(pi()*a/2)));
            Fvalue=log10(abs(A*pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2)*sin(b*phi)));
        }
        else if (source=="Kremer")
        {
            double A=c[0],a=c[1],b=c[2],tau=pow(10,c[3]);
            double r_w=pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2);
            double p_w=atan(sin(pi()*a/2)/(pow(afreq*tau,-a)+cos(pi()*a/2)));
            Fvalue=log10(abs(A*r_w*sin(b*p_w)));
        }
    }
    else
    {
        cout
        << "in HN_func_value, " << tcalc << " model not found!" << "\n"
        << "Program aborted." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    return Fvalue;
}





double FitDielectric::HN_tauFit_calc(const alglib::real_1d_array& c,
                                     const std::string& tcalc)
{
    double tauFit=0;
    
    if ((tcalc=="HN_loss")||(tcalc=="CC_loss")||(tcalc=="CD_loss"))
    {
        double a=c[1],b=c[2],tau=pow(10,c[3]); // A=c[0]
        //tauFit=tau*pow(sin((pi()*a*b)/(2*(1+b)))/sin((pi()*a)/(2*(1+b))),1/a);
        
        double fmax=
        pow(2*pi()*tau,-1)*pow(sin((pi()*a)/(2*(1+b)))/sin((pi()*a*b)/(2*(1+b))),1/a);
        tauFit=pow(fmax,-1);
        //tauFit=pow(2*pi()*fmax,-1);
    }
    else
    {
        cout
        << "in HN_tauFit_calc, " << tcalc << " model not found!" << "\n"
        << "Program aborted." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    return tauFit;
}





void FitDielectric::write_fit_HN(std::ofstream& outputFile,
                                 const std::string& tcalc,
                                 const int whichpeak)
{
    if (outputFile.is_open()) {
        
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
        << "tauMax    "    << log10(HN_tauFit_calc(c,tcalc)) << "\n\n";
        
        
        if (whichpeak==0)
        {
            // Raw Data
            //------------------------------------------------------------------
            outputFile
            << "Raw Data                               " << "\n"
            << "log.afreq  log.sfreq  log.loss  storage" << "\n"
            << "=======================================" << "\n";
            for (size_t i=0; i<dielectric_freq_overall.size(); ++i) {
                outputFile
                << log10(abs(dielectric_freq_overall[i][0])) << " "
                << log10(abs(dielectric_freq_overall[i][1])) << " "
                << log10(abs(dielectric_freq_overall[i][2])) << " "
                << dielectric_freq_overall[i][3]        << "\n";
            }
            outputFile << "\n\n";
            //------------------------------------------------------------------
        }
        
        
        
        // Fit Data (In Fitting Range)
        //----------------------------------------------------------------------
        outputFile
        << "Fit Data (In Fitting Range)               " << "\n"
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
            << HN_func_value(c,pow(10,afreq_tmp),tcalc)
            << "\n";
        }
        outputFile << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        int    equal_pieces=1000;
        double daf=abs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=abs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        
        
        if (whichpeak==0)
        {
            // Interpolated Data (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "Interpolated Data (Full Range)   " << "\n"
            << "log.afreq  log.sfreq  f  df  d2f " << "\n"
            << "=================================" << "\n";
            double f=0,df=0,d2f=0;
            for (int i=0; i<=equal_pieces; ++i) {
                
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                
                alglib::spline1ddiff(losscurve,afreqi,f,df,d2f);
                
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
        
        
        
        // Fit Data (Full Range)
        //----------------------------------------------------------------------
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
            << HN_func_value(c,pow(10,afreqi),tcalc) << "\n";
        }
        outputFile << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        double a=c[1],b=c[2],tau=pow(10,c[3]); // A=c[0]
        double t=0;
        double y=0;
        double theta=0;
        double omega=0;
        double f_t=0;
        
        
        
        // Distribution of tau (Full Range)
        //----------------------------------------------------------------------
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
            << log(abs(f_t)) << "\n";
        }
        //----------------------------------------------------------------------
    }
    else {
        cout
        << "write_fit_HN: file cannot open."
        << "\n";
    }
}





void FitDielectric::write_HN_params(std::ofstream& outputFile,
                                    const double Temp_d)
{
    outputFile << Temp_d*pow(corrFacforT,-1) << " ";
    for (int i=0; i<c.length(); ++i) outputFile << c[i] << " ";
    outputFile << rep.r2 << "\n";
}





void FitDielectric::write_peak_info(std::ofstream& outputFile,
                                    const int n_peaks,
                                    const int whichpeak)
{
    if (whichpeak==1) // NOTE: 1-based when output to file
    {
        outputFile
        << "total # of peaks = " << n_peaks
        << "\n\n"
        << "From low frequency, the #" << whichpeak << " peak"
        << "\n\n";
    }
    else
    {
        outputFile
        << "\n\n"
        << "From low frequency, the #" << whichpeak << " peak"
        << "\n\n";
    }
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
        if (source=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (source=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
    else if (model=="CC_loss")
    {
        if (source=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (source=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
    else if (model=="CD_loss")
    {
        if (source=="wiki")
        {
            alglib::lsfitfit(state, HN_loss_wiki);
        }
        else if (source=="Kremer")
        {
            alglib::lsfitfit(state, HN_loss_Kremer);
        }
    }
}





void FitDielectric::set_fitParams(const std::string& model)
{
    epsf     = 0;
    epsx     = 1e-10;
    maxits   = 1e+3;
    diffstep = 1e-5;
    
    if (model=="HN_loss") // [A,a,b,tau]
    {
        // HN: shape params are free
        
        // NOTE:
        // tau is in log scale
        
        set_coeffs({       1.0,  0.5,  0.5,  -8.0});
        set_coeffs_scale({ 1.0,  0.5,  0.5,  -8.0});
        set_coeffs_bndl({  0.0,  0.0,  0.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
    }
    else if (model=="CC_loss") // [A,a,b,tau]
    {
        // Cole-Cole: special case of HN when beta=1
        
        // NOTE:
        // tau is in log scale
        
        set_coeffs({       1.0,  0.5,  1.0,  -8.0});
        set_coeffs_scale({ 1.0,  0.5,  1.0,  -8.0});
        set_coeffs_bndl({  0.0,  0.0,  1.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
    }
    else if (model=="CD_loss") // [A,a,b,tau]
    {
        // Cole-Davidson: special case of HN when alpha=1
        
        // NOTE:
        // tau is in log scale
        
        set_coeffs({       1.0,  1.0,  0.5,  -8.0});
        set_coeffs_scale({ 1.0,  1.0,  0.5,  -8.0});
        set_coeffs_bndl({  0.0,  1.0,  0.0, -20.0});
        set_coeffs_bndu({ 10.0,  1.0,  1.0,  20.0});
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
    func = log10(abs(A*pow(1+2*pow(afreq*tau,a)*cos(pi()*a/2)+pow(afreq*tau,2*a),-b/2)*sin(b*phi)));
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
    func = log10(abs(A*r_w*sin(b*p_w)));
}





/* public setters */
void FitDielectric::set_is_use_dc_free(const bool b){is_use_dc_free=b;}



