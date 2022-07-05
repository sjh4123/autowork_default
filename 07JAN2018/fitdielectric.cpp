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

#include "fitdielectric.h"
#include "functions.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

FitDielectric::FitDielectric(const StructureClass& sysVar):

/* Base class */
AlglibFittingKernel(),

/* bool */
is_use_dc_free(false),
is_manual_cutlowT(false),
is_manualFitSingleDS(false),

/* int */
n_fit_HN(3),
n_check(3),
n_used_HN(100),
n_each_side(200),
n_points(0),
n_dcpoints(0),
n_peaks(0),
n_total(0),
n_count(0),
n_reduction(0),
counter(0),
process_id(0),
n_points_min(5),//NOTE

/* double */
r2_threshold(0.98),
r2_peak_threshold(0.90),
r2_combined(0),
r2_global_pre(0),
r2_global_now(0),
n_peaks_pre(0),
n_peaks_now(0),
lossfrac_cut(0.01),
npoints_frac(0.25),
df_threshold(0.5),
dc_threshold(0.3),
slope_deviation(0.25),
cutlowT(0),
cutlowT_manual(0),

/* range parameters */
afreq_MIN(0),
afreq_MAX(0),
sfreq_MIN(0),
sfreq_MAX(0),
log_afreq_MIN(0),
log_afreq_MAX(0),
log_sfreq_MIN(0),
log_sfreq_MAX(0),
Theq(0),
Tleq(0),
T_actual(0),
T_nominal(0),

/* others */
tau_sectops(12.0),

/* string */
form(""),
func("HN_loss"),
model("COOP"),
func_HN("Kremer"),//"wiki" or "Kremer"
process("cooling"),

/** STL containers **/
//------------------------------------------------------------------------------
is_log_data(),
spc_id(),
T_id(),
loss_max(),
loss_base(),
afreq_max(),
afreq_base(),
sfreq_max(),
sfreq_base(),
cutoff_lo_afreq(),
cutoff_lo_sfreq(),
cutoff_hi_afreq(),
cutoff_hi_sfreq(),
sfreq_range({0.0,inf}),

/** internal dummy variables **/
//------------------------------------------------------------------------------
indexii(0),
peak_indexii(0),
base_indexii(0),
whichpeak(0),
afreq(0),
afreqi(0),    // angular frequecy 2*pi()*sfreq
sfreq(0),
sfreqi(0),    // typical frequency 1/sec
storage(0),   // eps'
loss(0),      // eps''
afreq_now(0),
afreq_pre(0),
afreq_prr(0),
f_now(0),
df_now(0),
d2f_now(0),
f_pre(0),
df_pre(0),
d2f_pre(0),
f_prr(0),
df_prr(0),
d2f_prr(0),

/** ALGLIB fitting DS **/
//------------------------------------------------------------------------------
c_dc_loss(),
rep_dc_loss(),
interpolantRawLoss(),
interpolantDCLoss(),
interpolanttaueqCurve(),
c_combined(),
c_global(),
c_losspeaks(),
rep_combined(),
rep_global(),
rep_losspeaks(),
info_losspeaks(),
usedt_losspeaks(),
cyclecount_losspeaks(),

/** Dielectric data containers **/
//------------------------------------------------------------------------------
peak_indices_sampled(),
dielectric_freq_sortincreasing(),
dielectric_freq_domains(),
dielectric_freq_reduced_all(),
dielectric_freq_global(),
dc_curve_freq_sortincreasing(),
dielectric_freq_overall(),
dielectric_freq_dcfree(),
dielectric_freq_reduced(),
dielectric_freq(),
dielectric_freq_peak(),
masterloss()
{
    /** assign values to members in parent class (AlglibFittingKernel) **/
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
    
    if (form=="NPIC_ht") // NPIC high-throughput form
    {
        /** NOTE:
         ** File names from NPIC are first written to a temporary file by bash;
         ** C++ then reads the file to extract file names **/
        
        /* input file */
        input.append(return_AnalysisFolderPath(sysVar));
        input.append("/statistics/dielectric");
        string tmp="cd "+input+";for i in *;do echo ${i} >> "+input+"/tmp.txt;done";
        system(tmp.c_str());
        string filepath=input+"/tmp.txt";
        ifstream readFile(filepath.c_str());
        
        /* output file */
        string input,output;
        string overall;
        
        bool   is_newblock=false;
        int    lineindex=0;
        double temp=0,sfreq=0,epsp=0,epspp=0;
        string lineContent,spcstr,tempstr,extstr;
        
        vector<string> species;
        vector<string> fileext;
        vector<double> tempFile_vd;
        vector<double> temp_vd;
        //vector<vector<string>> species_2d;
        vector<vector<string>> fileext_2d;
        vector<vector<double>> tempFile_2d;
        vector<vector<double>> temp_2d;
        
        storedspecies.clear();
        
        while (getline(readFile,lineContent))
        {
            spcstr.clear();
            tempstr.clear();
            extstr.clear();
            
            int dlm=-1,unit=-1,ext=-1;
            is_newblock=false;
            
            /* extract species from file name */
            for (int i=0;i<lineContent.size();++i) {
                if (lineContent.at(i)=='-') {
                    dlm=i; break;
                } spcstr+=lineContent.at(i);
            }
            
            /* store species */
            if (lineindex>0)
            {
                if (spcstr==species.at(lineindex-1)) {
                    species.push_back(spcstr);
                } else {
                    is_newblock=true;
                    lineindex=0;
                    //species_2d.push_back(species);
                    storedspecies.push_back(spcstr);//2nd to last species
                    species.clear();
                    species.push_back(spcstr);
                }
            } else {
                species.push_back(spcstr);
                storedspecies.push_back(spcstr);//1st species
            }
            
            /* extract temperature from file name */
            for (int i=dlm+1;i<lineContent.size();++i) {
                if (lineContent.at(i)=='C') {
                    unit=i; break;
                } tempstr+=lineContent.at(i);
            } temp=atof(tempstr.c_str());
            
            /* store temperatures */
            if (is_newblock==false) {
                tempFile_vd.push_back(temp);
                temp_vd.push_back(temp+273);
            } else {
                if (process=="cooling") {
                    std::sort(temp_vd.begin(),temp_vd.end(),sortDecreasing1D);
                    std::sort(tempFile_vd.begin(),tempFile_vd.end(),sortDecreasing1D);
                } else if (process=="heating") {
                    std::sort(temp_vd.begin(),temp_vd.end(),sortIncreasing1D);
                    std::sort(tempFile_vd.begin(),tempFile_vd.end(),sortIncreasing1D);
                } else {
                    cout
                    << "in dielectricData_preprocess():\n"
                    << "process specified ("+process+") not found. "
                    << "Valid processes are 'cooling' or 'heating.'\n";
                    exit(EXIT_FAILURE);
                }
                tempFile_2d.push_back(tempFile_vd);
                temp_2d.push_back(temp_vd);
                tempFile_vd.clear();
                temp_vd.clear();
                tempFile_vd.push_back(temp);
                temp_vd.push_back(temp+273);
            }
            
            /* extract file extension from file name */
            for (int i=unit+1;i<lineContent.size();++i) {
                if (lineContent.at(i)=='.') {
                    ext=i; continue;
                } extstr+=lineContent.at(i);
            }
            
            /* store file extensions */
            if (is_newblock==false) {
                fileext.push_back(extstr);
            } else {
                fileext_2d.push_back(fileext);
                fileext.clear();
                fileext.push_back(extstr);
            }
            
            /* exceptions */
            if (dlm<=0||unit<0||ext<0) {
                cout
                << "in dielectricData_preprocess():\n"
                << "Dielectric file ("+lineContent+") has wrong file name format.\n"
                << "Naming rule: species-TinC.ext, ex. 1-85C.csv\n";
                exit(EXIT_FAILURE);
            } ++lineindex;
            
        } readFile.close();
        
        /* store final block info */
        //species_2d.push_back(species);
        
        if (process=="cooling") {
            std::sort(temp_vd.begin(),temp_vd.end(),sortDecreasing1D);
            std::sort(tempFile_vd.begin(),tempFile_vd.end(),sortDecreasing1D);
        } else if (process=="heating") {
            std::sort(temp_vd.begin(),temp_vd.end(),sortIncreasing1D);
            std::sort(tempFile_vd.begin(),tempFile_vd.end(),sortIncreasing1D);
        }
        tempFile_2d.push_back(tempFile_vd);
        temp_2d.push_back(temp_vd);
        fileext_2d.push_back(fileext);
        
        tinfo_dielectric=temp_2d;//NOTE
        
        /* remove tmp file */
        //string rm="rm "+filepath+";";system(rm.c_str());
        
        dielectric_freq_NPIC.clear();
        spc_Temps.clear();
        
        for (int species=0;species<(int)storedspecies.size();++species)
        {
            vector<double> Temps;
            vector<vector<vector<double>>> tmpvvvd;
            
            for (int T=0;T<temp_2d.at(species).size();++T)
            {
                vector<vector<double>> tmpvvd;
                
                Temps.push_back(temp_2d.at(species).at(T));
                
                /* write dielectric data into separate files */
                output.clear();
                output.append(return_AnalysisFolderPath(sysVar));
                output.append("/statistics/dielectric/");
                output.append("species_"+storedspecies.at(species)+"_");
                output.append(to_string((long long int)tempFile_2d.at(species).at(T))+"C");
                output.append(".dat");
                ofstream writeFile(output.c_str(),ofstream::app);//append
                
                /* group all Ts of one species into one file */
                overall.clear();
                overall.append(return_AnalysisFolderPath(sysVar));
                overall.append("/fit_data");
                overall.append("/DSdata_species_"+storedspecies.at(species)+"_");
                overall.append(".dat");
                ofstream writeOverall(overall.c_str(),ofstream::app);//append
                
                /* input file */
                input.clear();
                input.append(return_AnalysisFolderPath(sysVar));
                input.append("/statistics/dielectric/");
                input.append(storedspecies.at(species)+"-");//NOTE:which species to read
                input.append(to_string((long long int)tempFile_2d.at(species).at(T))+"C");
                input.append("."+fileext_2d.at(species).at(T));
                ifstream readFile(input.c_str());
                
                if (readFile.is_open())
                {
                    /* skip first four lines */
                    getline(readFile,lineContent);
                    getline(readFile,lineContent);
                    getline(readFile,lineContent);
                    getline(readFile,lineContent);
                    string str;
                    vector<string> dlmContent;
                    
                    while (getline(readFile,lineContent))
                    {
                        dlmContent.clear();
                        for (int i=0;i<(int)lineContent.size();++i) {
                            if (lineContent.at(i)==',') {//NOTE: dlm=','
                                dlmContent.push_back(str); str.clear();
                                continue;
                            } str+=lineContent.at(i);
                        }
                        sfreq=atof(dlmContent.at(4).c_str()); //frequency(1/s)   index=4
                        epspp=atof(dlmContent.at(11).c_str());//loss response    index=11
                        epsp =atof(dlmContent.at(10).c_str());//storage response index=10
                        
                        if (sfreq>=sfreq_range.at(0)&&sfreq<=sfreq_range.at(1))
                        {
                            /** avoid contamination by negative signal **/
                            //--------------------------------------------------
                            if (epspp>0&&epsp>0)
                            {
                                tmpvvd.push_back
                                ({2*pi()*sfreq,sfreq,epspp,epsp});
                                
                                writeFile
                                << sfreq        << " "
                                << 2*pi()*sfreq << " "
                                << epspp        << " "
                                << epsp         << "\n";
                                
                                writeOverall
                                << to_string((long long int)tempFile_2d.at(species).at(T)) << " "
                                << to_string((long long int)temp_2d.at(species).at(T)) << " "
                                << sfreq        << " "
                                << 2*pi()*sfreq << " "
                                << epspp        << " "
                                << epsp         << "\n";
                            }
                            //--------------------------------------------------
                        }
                    }
                    
                    try {
                        if (tmpvvd.size()>0) {
                            tmpvvvd.push_back(tmpvvd);
                        } else {
                            throw tmpvvd.size();
                        }
                    } catch (size_t s) {
                        cout
                        << "in FitDielectric::dielectricData_preprocess():\n"
                        << "form="<<form<<": "
                        << "species#"<<species<<" at T="<<temp_2d.at(species).at(T)<<"\n"
                        << "no dielectric data collected.\n"; exit(EXIT_FAILURE);
                    }
                    
                    writeOverall << "\n";
                    readFile.close();
                    writeFile.close();
                    writeOverall.close();
                    
                } else {
                    cout
                    << "in dielectricData_preprocess():\n"
                    << "file ("+input+") canoot open.\n";
                    exit(EXIT_FAILURE);
                }
            }
            
            spc_Temps.push_back(Temps);
            
            try {
                if (tmpvvvd.size()>0) {
                    dielectric_freq_NPIC.push_back(tmpvvvd);
                } else {
                    throw tmpvvvd.size();
                }
            } catch (size_t s) {
                cout
                << "in FitDielectric::dielectricData_preprocess():\n"
                << "form="<<form<<": species#"<<species<<"\n"
                << "no dielectric data collected.\n"; exit(EXIT_FAILURE);
            }
        }
    }
    else if (form=="NPIC_cs") // NPIC Cryostat form
    {
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
        overall.append("/DSdata_overall_");
        overall.append(sysVar.get_usic());
        overall.append(".dat");
        ofstream writeOverall(overall.c_str(),ofstream::app);//append
        
        ifstream readFile(input.c_str());
        if (readFile.is_open())
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
            
            string str;
            vector<string> dlmContent;
            
            if(true)
            {
                /* skip first four lines */
                getline(readFile,lineContent);
                getline(readFile,lineContent);
                getline(readFile,lineContent);
                getline(readFile,lineContent);
            }
            while (getline(readFile,lineContent))
            {
                ++lineindex;
                dlmContent.clear();
                for (int i=0;i<(int)lineContent.size();++i) {
                    if (lineContent.at(i)==',') {//NOTE: dlm=','
                        dlmContent.push_back(str); str.clear();
                        continue;
                    } str+=lineContent.at(i);
                }
                sfreq=atof(dlmContent.at(4).c_str()); //frequency(1/s)   index=4
                epspp=atof(dlmContent.at(11).c_str());//loss response    index=11
                epsp =atof(dlmContent.at(10).c_str());//storage response index=10
                T_now=atof(dlmContent.at(8).c_str()); //T in Celcius     index=8;
                T_now+=273;
                
                if (sfreq<0) {
                    cout
                    << "Reading dielectric data (Data format:"<<form<<"):\n"
                    << "@T="<<T_now<<"K, frequency < 0\n";
                    exit(EXIT_FAILURE);
                }
                if (lineindex>1) //need to have at least one T value
                {
                    /* write data to file at the end of each temperature */
                    if (T_now!=T_pre)
                    {
                        writeOverall << "\n";
                        is_in_cooling=false;
                        Temp_d=std::round(T_pre*corrFacforT);
                        
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
                        else { /* T_now>T_pre */
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
                        /* create file for one of the processes */
                        ofstream writeFile(output.c_str());
                        for (int i=0; i<(int)tnow_vd.size(); ++i) {
                            writeFile
                            << 2*pi()*sfreq_vd[i] << " "  //afreq
                            << sfreq_vd[i]        << " "  //sfreq
                            << epspp_vd[i]        << " "  //loss
                            << epsp_vd[i]         << "\n";//storage
                        } writeFile.close();
                        
                        /* clear containers */
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
                << epspp << " "
                << epsp  << "\n";
                
            } readFile.close(); writeOverall.close();
            
            /* store data of the final T block */
            Temp_d=std::round(T_now*corrFacforT);
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
            ofstream writeFile(output.c_str());
            for (int i=0; i<(int)tnow_vd.size(); ++i) {
                writeFile
                << 2*pi()*sfreq_vd[i] << " "  //afreq
                << sfreq_vd[i]        << " "  //sfreq
                << epspp_vd[i]        << " "  //loss
                << epsp_vd[i]         << "\n";//storage
            } writeFile.close();
            
            if (tinfo_cooling.size()>0) {
                tinfo_dielectric.push_back(tinfo_cooling);
            }
            if (tinfo_heating.size()>0) {
                tinfo_dielectric.push_back(tinfo_heating);
            }
        }
        else {
            cout
            << "in FitDielectric::dielectricData_preprocess():\n"
            << "dielectric file cannot open.\n";
        }
    }
    else if (form=="NIST")
    {
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
        overall.append("/DSdata_overall_");
        overall.append(sysVar.get_usic());
        overall.append(".dat");
        ofstream writeOverall(overall.c_str(),ofstream::app);//append
        
        ifstream readFile(input.c_str());
        if (readFile.is_open())
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
            
            if(true)
            {
                /* first 3 lines are useless */
                getline(readFile,lineContent);
                getline(readFile,lineContent);
                getline(readFile,lineContent);
            }
            while (getline(readFile,lineContent))
            {
                ++lineindex;
                istringstream iss(lineContent);
                
                iss >> sfreq;
                iss >> T_now; T_now+=273;
                iss >> epsp;
                iss >> epspp;
                
                if (sfreq<0) {
                    cout
                    << "Reading dielectric data (Data format:"<<form<<"):\n"
                    << "@T="<<T_now<<"K, frequency < 0\n";
                    exit(EXIT_FAILURE);
                }
                if (lineindex>1) //need to have at least one T value
                {
                    /* write data to file at the end of each temperature */
                    if (T_now!=T_pre)
                    {
                        writeOverall << "\n";
                        is_in_cooling=false;
                        Temp_d=std::round(T_pre*corrFacforT);
                        
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
                        else { /* T_now>T_pre */
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
                        /* create file for one of the processes */
                        ofstream writeFile(output.c_str());
                        for (int i=0; i<(int)tnow_vd.size(); ++i) {
                            writeFile
                            << 2*pi()*sfreq_vd[i] << " "  //afreq
                            << sfreq_vd[i]        << " "  //sfreq
                            << epspp_vd[i]        << " "  //loss
                            << epsp_vd[i]         << "\n";//storage
                        } writeFile.close();
                        
                        /* clear containers */
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
                << epspp << " "
                << epsp  << "\n";
                
            } readFile.close(); writeOverall.close();
            
            /* store data of the final T block */
            Temp_d=std::round(T_now*corrFacforT);
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
            ofstream writeFile(output.c_str());
            for (int i=0; i<(int)tnow_vd.size(); ++i) {
                writeFile
                << 2*pi()*sfreq_vd[i] << " "  //afreq
                << sfreq_vd[i]        << " "  //sfreq
                << epspp_vd[i]        << " "  //loss
                << epsp_vd[i]         << "\n";//storage
            } writeFile.close();
            
            if (tinfo_cooling.size()>0) {
                tinfo_dielectric.push_back(tinfo_cooling);
            }
            if (tinfo_heating.size()>0) {
                tinfo_dielectric.push_back(tinfo_heating);
            }
        }
        else {
            cout
            << "in FitDielectric::dielectricData_preprocess():\n"
            << "dielectric file cannot open.\n";
        }
    }
    else if (form=="manuallySet")
    {
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
        ifstream readFile(input.c_str());
        if (readFile.is_open()) {
            int    lineindex=0;
            string lineContent;
            while (getline(readFile,lineContent)){
                ++lineindex;
                istringstream iss(lineContent);
                iss >> sfreq; if (is_log_data[0]) sfreq=pow(10,sfreq);
                iss >> loss;  if (is_log_data[1]) loss =pow(10,loss);
                if (sfreq<0) {
                    cout
                    << "Reading dielectric data:\n"
                    << "frequency < 0           \n";
                    exit(EXIT_FAILURE);
                } dielectric_freq.push_back({sfreq,loss});
            } readFile.close();
        }
        ofstream writeFile(output.c_str());
        for (int i=0; i<(int)dielectric_freq.size(); ++i) {
            writeFile
            << dielectric_freq[i][0] << " "
            << dielectric_freq[i][1] << "\n";
        } writeFile.close();
    }
    return tinfo_dielectric;
}





void FitDielectric::read_all_dielectric_data(const StructureClass& sysVar,
                                             int process,
                                             double Temp_d)
{
    dielectric_freq_overall.clear();
    
    string input;
    vector<vector<double>> dielectric_freq;
    
    /* read in dielectric data source */
    //**************************************************************************
    if (form=="NPIC_ht")
    {
        int spc_id      = get_spc_id().at(0);
        int T_id        = get_T_id().at(0);
        dielectric_freq = dielectric_freq_NPIC.at(spc_id).at(T_id);
        
        if (true)
        {
            cout << "\n";
            cout << "Species Name   "<< storedspecies.at(spc_id)                     <<"\n";
            cout << "Species ID     "<< get_spc_id().at(0)<<"/"<<get_spc_id().at(1)  <<"\n";
            cout << "Temp ID        "<< get_T_id().at(0)+1  <<"/"<<get_T_id().at(1)  <<"\n";
            cout << "Temp(K)        "<< spc_Temps.at(spc_id).at(T_id)                <<"\n";
            cout << "Temp(C)        "<< spc_Temps.at(spc_id).at(T_id)-273            <<"\n";
        }
    }
    else if (form=="NIST"||form=="NPIC_cs")
    {
        if (process==0) {
            input.append(return_AnalysisFolderPath(sysVar));
            input.append("/statistics/dielectric");
            input.append("/cooling_eps_");
            input.append(sysVar.get_usic());
            input.append("_T"+to_string((long long int)Temp_d));
            input.append(".dat");
            
        } else if (process==1) {
            input.append(return_AnalysisFolderPath(sysVar));
            input.append("/statistics/dielectric");
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
                /** avoid contamination by negative signal **/
                //--------------------------------------------------------------
                if (loss>0&&storage>0)
                {
                    dielectric_freq.push_back({afreq,sfreq,loss,storage});
                }
                //--------------------------------------------------------------
            } readFile.close();
        } else {
            cout
            << "FitDielectric::read_all_dielectric_data():\n"
            << "NIST format: dielectric file cannot open.\n";
        }
    }
    else if (form=="manuallySet")
    {
        input.append(return_AnalysisFolderPath(sysVar));
        input.append("/statistics/dielectric");
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
                /** avoid contamination by negative signal **/
                //--------------------------------------------------------------
                if (loss>0&&storage>0)
                {
                    dielectric_freq.push_back({afreq,sfreq,loss,storage});
                }
                //--------------------------------------------------------------
            } readFile.close();
        } else {
            cout
            << "FitDielectric::read_all_dielectric_data():\n"
            << "manuallySet format: dielectric file cannot open.\n";
        }
    }
    /* NOTE:
     * original dielectric data are sorted by angular frequency (low to high) */
    try {
        if (dielectric_freq.size()==0) throw 0;
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
    }
    catch (int i) {
        cout
        << "in FitDielectric::read_all_dielectric_data():\n"
        << "Dielectric data were not able to read!\n";
        exit(EXIT_FAILURE);
    }
    //**************************************************************************
    
    afreq_MIN=dielectric_freq[0][0];
    sfreq_MIN=dielectric_freq[0][1];
    afreq_MAX=dielectric_freq[dielectric_freq.size()-1][0];
    sfreq_MAX=dielectric_freq[dielectric_freq.size()-1][1];
    
    log_afreq_MIN=log10(fabs(afreq_MIN));
    log_afreq_MAX=log10(fabs(afreq_MAX));
    log_sfreq_MIN=log10(fabs(sfreq_MIN));
    log_sfreq_MAX=log10(fabs(sfreq_MAX));
    
    /* establish derivative loss */
    //--------------------------------------------------------------------------
    vector<vector<double>> dielectric_dc_free;
    vector<vector<double>> interpolantdata;
    for (int i=0; i<(int)dielectric_freq.size(); ++i) {
        interpolantdata.push_back
        ({
            log(fabs(dielectric_freq[i][0])), //NOTE: ln(afreq)
            dielectric_freq[i][3]             //NOTE: storage
        });
    } build_penalizedspline_interpolant(interpolantdata,0,1);
    
    double afreqi=0,sfreqi=0,f=0,df=0,d2f=0;
    double storage_dc=0,loss_dc=0;
    for (int i=0; i<(int)dielectric_freq.size(); ++i) {
        afreqi     = dielectric_freq[i][0];
        sfreqi     = dielectric_freq[i][1];
        storage_dc = dielectric_freq[i][3];
        alglib::spline1ddiff(interpolant,log(afreqi),f,df,d2f); // NOTE: ln(afreq)
        loss_dc    = -(pi()/2)*df;
        dielectric_dc_free.push_back({afreqi,sfreqi,loss_dc,storage_dc});
    }
    dielectric_freq_dcfree=dielectric_dc_free;
    /* cubic spline interpolation of derivative loss (log-log) */
    //build_cubicspline_loglog_interpolant(dielectric_freq_dcfree,0,2);
    build_penalizedspline_loglog_interpolant(dielectric_freq_dcfree,0,2);
    interpolantDCLoss=interpolant;
    //--------------------------------------------------------------------------
    
    if (is_use_dc_free) {
        dielectric_freq_overall=dielectric_dc_free;
    } else {
        dielectric_freq_overall=dielectric_freq;
    }
    
    /* cubic spline interpolation of raw loss (log-log) */
    //build_cubicspline_loglog_interpolant(dielectric_freq_overall,0,2);
    build_penalizedspline_loglog_interpolant(dielectric_freq_overall,0,2);
    interpolantRawLoss=interpolant;
    
    /*----------------------------------------------------------------------
     Container Structure
     
     "dielectric_freq_overall"
     --------------------------------
     <I dimension>
     dielectric data
     
     <II dimension>
     dielectric data:
     0: afreq
     1: sfreq
     2: loss
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
    
    for (int i=0; i<(int)dielectric_freq_overall.size(); ++i)
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
    for (int i=cutoff_lo; i<(int)dielectric_freq_overall.size(); ++i)
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
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
    }
    catch (int i)
    {
        cout
        << "in FitDielectric::peak_sampling_byValue1:" << "\n"
        << "dielectric_freq size = " << i << "\n";
        dielectric_freq.push_back({0,0});
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
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
    for (int i=0; i<(int)dielectric_freq_peak.size(); ++i)
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
    for (int i=0; i<(int)dielectric_freq_peak.size(); ++i)
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
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
    }
    catch (int i)
    {
        cout
        << "in FitDielectric::peak_sampling_byValue2:" << "\n"
        << "dielectric_freq size = " << i << "\n";
        dielectric_freq.push_back({0,0});
        std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
    }
    dielectric_freq_peak.clear();
    dielectric_freq_peak=dielectric_freq;
    
    
    dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
}





void FitDielectric::peak_sampling_bySlope(const vector<vector<double>>& data)
{
    /** total # of data points **/
    n_points=(int)data.size();
    
    /** cubic spline interpolantion of the log-log loss curve **/
    build_penalizedspline_loglog_interpolant(data,0,2);
    
    
    /* 1. Distinguish Local Peaks */
    //**************************************************************************
    vector<int> peak_indices;
    
    for (indexii=1; indexii<n_points; ++indexii) { // NOTE: start from 1
        
        afreq_now=log10(fabs(data[indexii][0]));
        afreq_pre=log10(fabs(data[indexii-1][0]));
        alglib::spline1ddiff(interpolant,afreq_now,f_now,df_now,d2f_now);
        alglib::spline1ddiff(interpolant,afreq_pre,f_pre,df_pre,d2f_pre);
        
        if (df_pre*df_now<0) {
            
            if (data[indexii][2]>data[indexii-1][2]) { // NOTE: >
                peak_indexii=indexii;
            } else {
                peak_indexii=indexii-1;
            }
            
            /* NOTE: current peak should NOT repeat or be too close to peaks
             * already sampled */
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
            
            if (n_count==n_total) {
                peak_indices.push_back(peak_indexii);
                peak_indices_original.push_back(peak_indexii);
            }
        }
    }
    //**************************************************************************
    
    
    
    /* 1.1 Reduced Loss Mode */
    //**************************************************************************
    if( n_reduction>0 && peak_indices.size()>0 ) {
        
        /* sort found peaks of redu.loss in decreasing order in loss magnitude */
        double peakLoss=0;
        vector<vector<double>> peakLoss_vec;
        for (whichpeak=0; whichpeak<(int)peak_indices.size(); ++whichpeak) {
            peak_indexii=peak_indices.at(whichpeak);
            peakLoss_vec.push_back({data[peak_indexii][2],(double)peak_indexii});//(loss,index)
        } std::sort(peakLoss_vec.begin(),peakLoss_vec.end(),sortDecreasing0);
        peakLoss=peakLoss_vec[0][0];
        
        /* only use the highest peak to fit */
        for (whichpeak=0; whichpeak<(int)peak_indices.size(); ++whichpeak) {
            peak_indexii=peak_indices.at(whichpeak);
            if (peakLoss==data[peak_indexii][2]) break;
        } peak_indices.clear();
        
        /* the loss peak needs to be above a threshold value
         * NOTE: threshold based on raw loss magnitude */
        if (peakLoss>(lossfrac_cut*dielectric_freq_overall[peak_indexii][2])) {
            peak_indices.push_back(peak_indexii);
        }
    }
    //**************************************************************************
    
    
    
    /* 2. Data Sampling Around Found Local Peaks */
    //**************************************************************************
    if(peak_indices.size()>0)
    {
        for (whichpeak=0; whichpeak<(int)peak_indices.size(); ++whichpeak) {
            
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
            
            /* sample loss data wrt to the current peak (peak_indexii) */
            /* NOTE: backward_sampling() and forward_sampling() edits the content
             * of dielectric_freq that represents the sampled data */
            backward_sampling(data,peak_indexii);
            forward_sampling(data,peak_indexii);
            
            if (dielectric_freq.size()<n_points_min) // # of data too few, not taken the peak
            {
                continue;
            }
            else
            {
                std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
                for (int i=0; i<(int)dielectric_freq.size(); ++i)
                {
                    dielectric_freq_peak.push_back
                    ({
                        log10(fabs(dielectric_freq[i][0])), //log(afreq)
                        log10(fabs(dielectric_freq[i][2])), //log(loss)
                        log10(fabs(dielectric_freq[i][1])), //log(sfreq)
                        dielectric_freq[i][3]               //storage
                    });
                }
                //for (int i=0;i<(int)dielectric_freq_peak.size();++i) {
                //    cout
                //    << dielectric_freq_peak[i][2] << " "
                //    << dielectric_freq_peak[i][1] << "\n";
                //}
                
                /* NOTE: only take the current peak as a final validated peak when
                 * the r2 of HN fit is greater than r2_peak_threshold */
                fitdata_processing_lsfit(dielectric_freq_peak);
                fit_loss_peak();
                
                if (rep.r2>r2_peak_threshold)
                {
                    /** current fit info **/
                    c_losspeaks.push_back(c);
                    rep_losspeaks.push_back(rep);
                    info_losspeaks.push_back(info);
                    cyclecount_losspeaks.push_back(cyclecount);
                    usedt_losspeaks.push_back(usedt);
                    
                    /** domain info **/
                    afreq_max.push_back(data[peak_indexii][0]);
                    sfreq_max.push_back(data[peak_indexii][1]);
                    loss_max.push_back(data[peak_indexii][2]);
                    cutoff_lo_afreq.push_back(dielectric_freq[0][0]);
                    cutoff_lo_sfreq.push_back(dielectric_freq[0][1]);
                    cutoff_hi_afreq.push_back(dielectric_freq[dielectric_freq.size()-1][0]);
                    cutoff_hi_sfreq.push_back(dielectric_freq[dielectric_freq.size()-1][1]);
                    peak_indices_sampled.push_back(peak_indexii);
                    dielectric_freq_sortincreasing.push_back(dielectric_freq_peak);
                    
                    /* NOTE: all vectors shown in this block are containers for
                     * storing domain info; as a new loss domain (excluding the dc
                     * part) is found, it'll be pushed back into these containers
                     * so that we can retrieve data of each domain later */
                }
                else
                {
                    continue;
                }
            }
        }
    }
    //**************************************************************************
}





void FitDielectric::dc_curve_sampling(const vector<vector<double>>& data)
{
    // total # of data points
    n_points=(int)data.size();
    
    int firstLossPeakii=0;
    if (peak_indices_original.size()>0) {
        firstLossPeakii=
        *std::min_element(peak_indices_original.begin(),peak_indices_original.end());
        //firstLossPeakii=peak_indices_original[0];
    }
    
    if (firstLossPeakii>n_check) {
        
        /* 1. Distinguish Valid Base (lowest point in loss w/ dc effect ) */
        //**********************************************************************
        for (indexii=0; indexii<firstLossPeakii; ++indexii) { // NOTE: start from 0
            
            afreq_now=log10(fabs(data.at(indexii)[0]));
            alglib::spline1ddiff(interpolantRawLoss,afreq_now,f_now,df_now,d2f_now);
            
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
        dielectric_freq.clear(); dielectric_freq_peak.clear();
        
        dielectric_freq.push_back
        ({
            data[base_indexii][0],
            data[base_indexii][1],
            data[base_indexii][2],
            data[base_indexii][3]
        });
        /* sample data for dc fitting */
        base_sampling(data,base_indexii);
        
        if (dielectric_freq.size()<n_points_min) {
            return;
        } else {
            std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
            for (int i=0; i<(int)dielectric_freq.size(); ++i) {
                dielectric_freq_peak.push_back
                ({
                    log10(fabs(dielectric_freq[i][0])),       // log(afreq)
                    log10(fabs(fabs(dielectric_freq[i][2]))), // log(loss)
                    log10(fabs(dielectric_freq[i][1])),       // log(sfreq)
                    dielectric_freq[i][3]                     // storage
                });
            } dc_curve_freq_sortincreasing=dielectric_freq_peak; // NOTE!
        }
        //**********************************************************************
    }
    else {
        return;
    }
}





void FitDielectric::combine_domain_data_to_fit()
{
    double logafreq=0,logloss=0;
    vector<double> logafreq_vec,logloss_vec;
    
    /** 1. dc part (dc_curve_freq_sortincreasing) **/
    if (n_dcpoints>1)
    {
        for (int i=0; i<(int)dc_curve_freq_sortincreasing.size(); ++i) {
            logafreq=dc_curve_freq_sortincreasing[i][0];//log(afreq)
            logloss =dc_curve_freq_sortincreasing[i][1];//log(rawLoss)
            logafreq_vec.push_back(logafreq);
            logloss_vec.push_back(logloss);
        }
    }
    
    /** 2. Local relaxations (dielectric_freq_sortincreasing)
     **    NOTE: afreq should not repeat in local domains **/
    if (n_peaks>0)
    {
        for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
            for (int i=0; i<(int)dielectric_freq_sortincreasing.at(whichpeak).size(); ++i) {
                logafreq=dielectric_freq_sortincreasing.at(whichpeak)[i][0];//log(afreq)
                if (find_loss_given_logafreq(logafreq)>0) {
                    logloss=log10(find_loss_given_logafreq(logafreq));//log(rawLoss)
                } else {
                    continue;
                }
                if (!is_afreqRepeat(logafreq,logafreq_vec)) {
                    logafreq_vec.push_back(logafreq);
                    logloss_vec.push_back(logloss);
                }
            }
        }
    }
    
    /** 3. Combine 1.& 2. to give domain range data **/
    if (n_dcpoints>1||n_peaks>0)
    {
        dielectric_freq.clear();
        
        if (logafreq_vec.size()==logloss_vec.size()) {
            for (int i=0; i<(int)logafreq_vec.size(); ++i) {
                dielectric_freq.push_back({logafreq_vec[i],logloss_vec[i]});
            } std::sort(dielectric_freq.begin(),dielectric_freq.end(),sortIncreasing0);
            /** push back combined domain data into container **/
            dielectric_freq_domains.push_back(dielectric_freq);
        } else {
            cout
            << "in FitDielectric::combine_domain_data()" << "\n"
            << "logafreq_vec.size()!=logloss_vec.size()" << "\n";
            exit(EXIT_FAILURE);
        }
    }
}





void FitDielectric::combine_global_data_to_fit()
{
    dielectric_freq.clear();
    vector<vector<double>> raw=dielectric_freq_overall;
    
    double loss=0;
    for (int i=0; i<(int)raw.size(); ++i) {
        afreq=raw[i][0]; loss=raw[i][2];
        if (loss>0) {
            dielectric_freq.push_back({log10(raw[i][0]),log10(loss)});
        } else {
            dielectric_freq.push_back({log10(raw[i][0]),log10(Epsilon)});
        }
    } dielectric_freq_global.push_back(dielectric_freq);
}





void FitDielectric::update_fitinfo(const int whichpeak)
{
    n_points=(int)dielectric_freq_sortincreasing.at(whichpeak).size();
    fitdata_processing_lsfit(dielectric_freq_sortincreasing.at(whichpeak));
    set_fitParams(func);
    
    c=c_losspeaks.at(whichpeak);
    rep=rep_losspeaks.at(whichpeak);
    info=info_losspeaks.at(whichpeak);
    cyclecount=cyclecount_losspeaks.at(whichpeak);
    usedt=usedt_losspeaks.at(whichpeak);
}





double FitDielectric::find_loss_given_logafreq(const double logafreq)
{
    bool is_found=false;
    double rawloss=0;
    for (int i=0;i<(int)dielectric_freq_overall.size();++i) {
        if (logafreq==log10(dielectric_freq_overall[i][0])) {
            rawloss=dielectric_freq_overall[i][2];
            is_found=true;
        }
    }
    if (!is_found) {
        cout
        << "\nin FitDielectric::find_loss_given_logafreq():\n"
        << "can't find loss value at given afreq ("<<afreq<<").\n\n";
        for (int i=0;i<(int)dielectric_freq_overall.size();++i) {
            for (int ii=0;ii<(int)dielectric_freq_overall[i].size();++ii) {
                cout << dielectric_freq_overall[i][ii] << " ";
            } cout << "\n";
        } exit(EXIT_FAILURE);
    } return rawloss;
}





void FitDielectric::fit_dc_loss()
{
    //for (int i=0;i<(int)dc_curve_freq_sortincreasing.size();++i) {
    //    cout
    //    << dc_curve_freq_sortincreasing[i][2] << " "
    //    << dc_curve_freq_sortincreasing[i][1] << "\n";
    //}
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





void FitDielectric::fit_local_domains()
{
    /** sample data in the proximity of found peaks **/
    //----------------------------------------------------------------------
    if (n_reduction==0) {
        /* peak detection of raw loss */
        peak_sampling_bySlope(dielectric_freq_overall);
    } else {
        /* peak detection of reduced loss */
        peak_sampling_bySlope(dielectric_freq_reduced);
    } n_peaks=(int)dielectric_freq_sortincreasing.size();
    
    
    /** sample and fit dc data if non-dc free mode is chosen **/
    //----------------------------------------------------------------------
    if ((!is_use_dc_free)&&(n_reduction==0)) {
        dc_curve_sampling(dielectric_freq_overall);
        n_dcpoints=(int)dc_curve_freq_sortincreasing.size();
        if(n_dcpoints>1) fit_dc_loss();
    }
    
    /*----------------------------------------------------------------------
     Container Structure
     
     "dielectric_freq_sortincreasing (3D double)"
     --------------------------------------------
     <I dimension>
     peaks
     
     <II dimension>
     dielectric data
     
     <III dimension>
     components
     0: log(afreq)
     1: log(loss)
     2: log(sfreq)
     3: storage
     
     "dc_curve_freq_sortincreasing (2D double)"
     --------------------------------------------
     <I dimension>
     dielectric data
     
     <II dimension>
     components
     0: log(afreq)
     1: log(loss)
     2: log(sfreq)
     3: storage
     ---------------------------------------------------------------------*/
}





void FitDielectric::fit_combined_local_domains()
{
    // combine functional forms (cumulative)
    // properly set fit params
    // use starting params from previous fit results
    // gather data in each domain
    // save fitting results for global fitting
    
    combine_domain_data_to_fit();
    fitdata_processing_lsfit(dielectric_freq);
    set_fitParams("combined");
    alglib_lsfit("combined");
    c_combined.push_back(c);
    rep_combined.push_back(rep);
    //cout_r2(); cout_c_parms();
}





void FitDielectric::fit_global_domain()
{
    // use combined functional forms (accumulative)
    // use overall data for fitting
    // use starting params from combined domain fitting
    // save fitting results
    
    combine_global_data_to_fit();
    fitdata_processing_lsfit(dielectric_freq);
    set_fitParams("global");
    alglib_lsfit("global");
    c_global.push_back(c);
    rep_global.push_back(rep);
    //cout_r2(); cout_c_parms();
}





void FitDielectric::calc_taueq_derivatives()
{
    /* establish cubic spline interpolation routine */
    vector<vector<double>> interpolantdata;
    double invT=0;
    double logtau=0;
    for (int i=0; i<(int)taueq_sortdecreasing.size(); ++i)
    {
        invT=(1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0]; // NOTE: inverse T
        logtau=taueq_sortdecreasing[i][1]; // NOTE: log(taueq)
        interpolantdata.push_back({invT,logtau});
    }
    build_penalizedspline_interpolant(interpolantdata,0,1);
    interpolanttaueqCurve=interpolant;
    
    //fitdata_processing_1d(interpolantdata);
    //real_1d_array x_in=fit_xData.c_str();
    //real_1d_array y_in=fit_yData.c_str();
    //alglib::spline1dbuildcubic(x_in,y_in,interpolanttaueqCurve);
}





void FitDielectric::peak_validation(const vector<vector<double>>& data,
                                    const int peak_index)
{
    int n_points=(int)data.size();
    
    /* forward checking: strict monotonic "decreasing" */
    //------------------------------------------------------------------
    n_count=0; n_total=n_check;
    
    for (int count=1; count<=n_check; ++count) {
        if ((peak_index+count)<n_points) {
            bool is_valid =
            data[peak_index+count-1][2]>data[peak_index+count][2];
            if (is_valid) ++n_count;
        }
    }
    
    /* backward checking: strict monotonic "decreasing" */
    //------------------------------------------------------------------
    bool is_at_hifreq_end = abs(peak_index-(n_points-1))<=(int)ceil(npoints_frac*n_points);
    bool is_at_lofreq_end = peak_index<=ceil(npoints_frac*n_points);
    
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
    int n_points=(int)data.size();
    
    /* forward checking */
    //------------------------------------------------------------------
    n_count=0; n_total=n_check;
    
    for (int count=1; count<=n_check; ++count) {
        if ((base_index+count)<n_points && base_index+count-1>=0)
        {
            if (data[base_index+count-1][2]<data[base_index+count][2]) {
                /* NOTE: begins to deflect */
                n_count=n_total; return;
            }
            //bool is_valid =
            //data[base_index+count-1][2]>data[base_index+count][2];
            //if (is_valid) ++n_count;
        }
    }
    /* NOTE:
     for conductivity effect, only forward checking for monotonic decreasing
     should be sufficient */
}





void FitDielectric::base_sampling(const vector<vector<double>>& data,
                                  const int base_index)
{
    int  n_points=(int)data.size();
    
    int  index_lo=0;
    bool is_valid=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        /* backward sampling from base_index-1 */
        index_lo = base_index-j;
        is_valid = index_lo<n_points && index_lo>=0;
        
        if (is_valid)
        {
            afreq_now=log10(fabs(data.at(index_lo)[0]));
            alglib::spline1ddiff(interpolantRawLoss,afreq_now,f_now,df_now,d2f_now);
            
            is_taken = df_now<0 && fabs(df_now+1.0)<dc_threshold;
            
            if (is_taken)
            {
                dielectric_freq.push_back
                ({
                    data[index_lo][0],
                    data[index_lo][1],
                    data[index_lo][2],
                    data[index_lo][3]
                });
            } else {
                break;
            }
        }
    } return;
}





void FitDielectric::backward_sampling(const vector<vector<double>>& data,
                                      const int peak_index)
{
    int  n_points=(int)data.size();
    
    int  index_lo=0;
    bool is_valid=false;
    bool is_at_lofreq_end=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        index_lo = peak_index-j;
        is_valid = index_lo<n_points && index_lo-2>=0;
        
        if (is_valid) {
            
            is_at_lofreq_end =
            index_lo<=ceil(npoints_frac*n_points);
            
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
    int  n_points=(int)data.size();
    
    int  index_hi=0;
    bool is_valid=false;
    bool is_at_hifreq_end=false;
    bool is_taken=false;
    
    for (int j=1; j<=n_each_side; ++j) {
        
        index_hi = peak_index+j;
        is_valid = index_hi+2<n_points && index_hi>=0;
        
        if (is_valid) {
            
            is_at_hifreq_end =
            abs(index_hi-(n_points-1))<=(int)ceil(npoints_frac*n_points);
            
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





bool FitDielectric::is_afreqRepeat(const double logafreq,
                                   std::vector<double>& domain)
{
    bool is_repeat=false;
    //--------------------------------------------------------------
    // "find" returns iterators in the range: [first,last); it would
    // return a pointer to the "last" element if nothing is found.
    //--------------------------------------------------------------
    vector<double>::iterator itr;
    itr=find(domain.begin(),domain.end(),logafreq);
    if (itr!=domain.end()) is_repeat=true;
    else is_repeat=false;
    return is_repeat;
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
    for (int i=0; i<(int)peakindices.size(); ++i) {
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
    
    for (int i=1; i<(int)data.size(); ++i) // NOTE: start from 1
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





void FitDielectric::reduce_raw_loss()
{
    dielectric_freq.clear();
    vector<vector<double>> raw=dielectric_freq_overall;
    
    double loss_reduced=0;
    vector<vector<double>> loss_reduced_vec;
    
    for (int i=0; i<(int)raw.size(); ++i)
    {
        afreq=raw[i][0];
        loss_reduced=raw[i][2]-combined_func_value(c_combined[n_reduction-1],afreq);
        
        /** NOTE: See static_combined_func() for fitting form **/
        
        if (loss_reduced>0) {
            //if (loss_reduced>(lossfrac_cut*raw[i][2])) {
            loss_reduced_vec.push_back({raw[i][0],raw[i][1],loss_reduced,raw[i][3]});
            /** NOTE: if the reduced loss has negative value, it would not be
             ** valid to show on log-log plot, and is not valid for fitting as
             ** well; therefore, only take positive reduced loss points **/
            //}
        } else {
            loss_reduced_vec.push_back({raw[i][0],raw[i][1],Epsilon,raw[i][3]});
            /** if a reduced loss is negative, in order to still keep the size of
             ** reduced loss vector the same as the original loss vector, while not
             ** couning the negative value, assign an epsilon value to it **/
        }
    }
    dielectric_freq_reduced=loss_reduced_vec;
    dielectric_freq_reduced_all.push_back(dielectric_freq_reduced);
    
    //for (int i=0;i<(int)dielectric_freq_reduced.size();++i) {
    //    cout
    //    << dielectric_freq_reduced[i][1] << " "
    //    << dielectric_freq_reduced[i][2] << "\n";
    //}
}





void FitDielectric::lossCurveReduction(const vector<vector<double>>& loss,
                                       const vector<vector<double>>& masterloss)
{
    /* NOTE:
     "dielectric_freq_reduced" should preserve the data structure:
     {afreq,sfreq,lossReduced,storage} */
    double loss_reduced=0; vector<vector<double>> loss_reduced_vec;
    for (int i=0; i<(int)loss.size(); ++i) {
        //loss_reduced=pow(10,masterloss[i][1])-loss[i][2];
        loss_reduced=loss[i][2]-pow(10,masterloss[i][1]);
        loss_reduced_vec.push_back
        ({loss[i][0],loss[i][1],loss_reduced,loss[i][3]});
    } //dielectric_freq_reduced=loss_reduced_vec;
}





void FitDielectric::clear_containers()
{
    //NOTE:
    //containers whose content are having T-dependent information
    //should make clear each time it is used to store new info of a new T's
    
    n_reduction=0;
    
    afreq_max.clear();       afreq_base.clear();
    loss_max.clear();        loss_base.clear();
    sfreq_max.clear();       sfreq_base.clear();
    cutoff_lo_afreq.clear(); cutoff_lo_sfreq.clear();
    cutoff_hi_afreq.clear(); cutoff_hi_sfreq.clear();
    
    c_losspeaks.clear();     rep_losspeaks.clear();
    c_combined.clear();      rep_combined.clear();
    c_global.clear();        rep_global.clear();
    
    info_losspeaks.clear();
    cyclecount_losspeaks.clear();
    usedt_losspeaks.clear();
    
    peak_indices_original.clear();
    peak_indices_sampled.clear();
    dielectric_freq_sortincreasing.clear();
    dielectric_freq_domains.clear();
    dielectric_freq_reduced_all.clear();
    dielectric_freq_global.clear();
    dc_curve_freq_sortincreasing.clear();
    dielectric_freq_overall.clear();
    dielectric_freq_reduced.clear();
    dielectric_freq.clear();
    dielectric_freq_peak.clear();
    masterloss.clear();
}





void FitDielectric::cout_peaksinfo()
{
    if (T_nominal!=0&&form!="NPIC_ht") {
        cout << "\n";
        double TinK=T_nominal*pow(corrFacforT,-1);
        cout << "Temp(K)        "<< TinK     <<"\n";
        cout << "Temp(C)        "<< TinK-273 <<"\n";
    }
    cout << "n_peaks        " << n_peaks << "\n";
    cout << "log(ang.freq)  ";
    for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
        cout << log10(fabs(afreq_max.at(whichpeak))) << " ";
    } cout << "\n";
    cout << "log(lin.freq)  ";
    for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
        cout << log10(fabs(sfreq_max.at(whichpeak))) << " ";
    } cout << "\n";
}





void FitDielectric::write_peaksinfo(const StructureClass& sysVar)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/fit_data");
    o.append("/info_peaks_");
    o.append(sysVar.get_usic());
    o.append(".dat");
    ofstream output(o.c_str(),ofstream::app); // append
    
    output
    << "T_nominal " << T_nominal << "  "
    << "n_peaks "   << n_peaks   << "  ";
    for (whichpeak=0; whichpeak<n_peaks; ++whichpeak) {
        output << log10(fabs(afreq_max.at(whichpeak))) << " ";
    } output << "\n\n";
    output.close();
}





void FitDielectric::write_localFitinfo(const StructureClass& sysVar)
{
    if (dielectric_freq_sortincreasing.size()>0)
    {
        string o;
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/info_localDomains_");
        o.append(sysVar.get_usic());
        o.append(".dat");
        ofstream output(o.c_str(),ofstream::app); // append
        
        output << "T_nominal " << T_nominal << "\n";
        for (whichpeak=0; whichpeak<dielectric_freq_sortincreasing.size(); ++whichpeak) {
            output
            << "peak #"     <<whichpeak+1<< " "
            << "log.afreq " << log10(afreq_max.at(whichpeak)) << "  "
            << dielectric_freq_sortincreasing.at(whichpeak).size() << "  "
            << rep_losspeaks.at(whichpeak).r2 << "  ";
            for (indexii=0; indexii<c_losspeaks.at(whichpeak).length(); ++indexii) {
                output << c_losspeaks.at(whichpeak)[indexii] << " ";
            } output << "\n";
        } output.close();
    }
}





void FitDielectric::write_combinedFitinfo(const StructureClass& sysVar)
{
    if (dielectric_freq_domains.size()>0)
    {
        string o;
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/info_combinedLocals_");
        o.append(sysVar.get_usic());
        o.append(".dat");
        ofstream output(o.c_str(),ofstream::app); // append
        
        output << "T_nominal " << T_nominal << "\n";
        for (int redu=0; redu<=n_reduction; ++redu) { // NOTE: <=
            output
            << "reduction #" << redu << "  "
            << dielectric_freq_domains.at(redu).size() << "  "
            << rep_combined.at(redu).r2 << "  ";
            for (int i=0; i<c_combined.at(redu).length(); ++i) {
                output << c_combined.at(redu)[i] << " ";
            } output << "\n";
        } output << "\n";
        output.close();
    }
}





void FitDielectric::write_globalFitinfo(const StructureClass& sysVar)
{
    if (dielectric_freq_global.size()>0)
    {
        string o;
        o.append(return_AnalysisFolderPath(sysVar));
        o.append("/fit_data");
        o.append("/info_combinedGlobal");
        o.append(sysVar.get_usic());
        o.append(".dat");
        ofstream output(o.c_str(),ofstream::app); // append
        
        output << "T_nominal " << T_nominal << "\n";
        for (int redu=0; redu<=n_reduction; ++redu) { // NOTE: <=
            output
            << "reduction #" << redu << "  "
            << dielectric_freq_global.at(redu).size() << "  "
            << rep_global.at(redu).r2 << "  ";
            for (int i=0; i<c_global.at(redu).length(); ++i) {
                output << c_global.at(redu)[i] << " ";
            } output << "\n";
        } output << "\n";
        output.close();
    }
}





void FitDielectric::fit_dielectric_loss(const StructureClass& sysVar,
                                        int process,
                                        double Temp_d)
{
    time_t tbeg=time(NULL),tend=time(NULL);
    
    is_use_FG=false;
    
    if (Temp_d==0) {
        T_nominal=spc_Temps.at(spc_id[0]).at(T_id[0])*corrFacforT;
        T_actual =spc_Temps.at(spc_id[0]).at(T_id[0]);
    } else {
        T_nominal=Temp_d;
        T_actual =T_nominal*pow(corrFacforT,-1);
    }
    
    string o1,o2,o3,o4,o5,o6;
    
    if (process==0) /** cooling process **/
    {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/cooling_taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/cooling_taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/cooling_fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/cooling_fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)T_nominal));
        o4.append(".dat");
    }
    else if (process==1) /** heating process **/
    {
        o1.clear();
        o1.append(return_AnalysisFolderPath(sysVar));
        o1.append("/fit_data");
        o1.append("/heating_taueq_HN_");
        o1.append(sysVar.get_usic());
        o1.append(".dat");
        o2.clear();
        o2.append(return_AnalysisFolderPath(sysVar));
        o2.append("/fit_data");
        o2.append("/heating_taueq_HN_inverseT_");
        o2.append(sysVar.get_usic());
        o2.append(".dat");
        o3.clear();
        o3.append(return_AnalysisFolderPath(sysVar));
        o3.append("/fit_data");
        o3.append("/heating_fit_HN_params_");
        o3.append(sysVar.get_usic());
        o3.append(".dat");
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/heating_fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append("_T"+to_string((long long int)T_nominal));
        o4.append(".dat");
    }
    else if (process==2) /** manual fit mode **/
    {
        o4.clear();
        o4.append(return_AnalysisFolderPath(sysVar));
        o4.append("/fit_data");
        o4.append("/fit_HN_");
        o4.append(sysVar.get_usic());
        o4.append(".dat");
    }
    if (process!=2) /** store heating and cooling data **/
    {
        o5.clear();
        o5.append(return_AnalysisFolderPath(sysVar));
        o5.append("/fit_data");
        o5.append("/overall_taueq_HN_");
        o5.append(sysVar.get_usic());
        o5.append(".dat");
        o6.clear();
        o6.append(return_AnalysisFolderPath(sysVar));
        o6.append("/fit_data");
        o6.append("/overall_taueq_HN_inverseT_");
        o6.append(sysVar.get_usic());
        o6.append(".dat");
    }
    ofstream write_taueq(o1.c_str(),ofstream::app);      //'append'
    ofstream write_taueqinv(o2.c_str(),ofstream::app);   //'append'
    ofstream write_params(o3.c_str(),ofstream::app);     //'append'
    ofstream write_HN(o4.c_str(),ofstream::app);         //'append'
    ofstream write_alltaueq(o5.c_str(),ofstream::app);   //'append'
    ofstream write_alltaueqinv(o6.c_str(),ofstream::app);//'append'
    
    if (form=="NPIC_ht") {
        /** To add a blank row in-between sample speices **/
        if (spc_id[0]>0) {
            if (T_id[0]==0) {
                write_blankToFile(write_taueq);
                write_blankToFile(write_alltaueq);
                write_blankToFile(write_taueqinv);
                write_blankToFile(write_alltaueqinv);
            }
        }
    }
    
    clear_containers();
    /** containers storing T-dependent infomation need to clear at each use **/
    
    read_all_dielectric_data(sysVar,process,T_nominal);
    /** "read_all_dielectric_data()" gives data container
     ** 'dielectric_freq_overall' that stores complete dielectric reponse of a T;
     ** dielectric_freq_overall has the structure: {afreq,sfreq,loss,storage}
     ** all in its original value (ie. not in log scale) **/
    
    /** avoid contamination by negative dielectric loss response **/
    for (int i=0;i<(int)dielectric_freq_overall.size();++i) {
        if (dielectric_freq_overall[i][2]<0||
            dielectric_freq_overall.size()<n_points_min) {
            cout << "Bad Dielectric Signal. Analysis skipped.\n";
            return;
        }
    }
    
    fit_local_domains(); n_peaks_pre=n_peaks;
    /** "fit_local_domains()" fits sampled data from each individual domain,
     ** including dc part (if found) and apparent relaxation features,
     ** characterized by found peaks from loss response **/
    
    if(n_dcpoints>1||n_peaks>0)
    {
        fit_combined_local_domains(); r2_combined=rep.r2;
        /** "fit_combined_local_domains()" combines the functional forms used in
         ** fit_local_domains() and fits the data in local domains; the initial
         ** parameter values are using those from local fits; this shuold give a
         ** btter fit of the sampled domains **/
        
        fit_global_domain(); r2_global_pre=rep.r2;
        /** fit_global_domain() fits loss data in global frequency domain by the
         ** combined functional form; this would tell how the combined function
         ** decribes the global data; if it describes the data well enough according to
         ** some measures, this would be the final step; if it doesn't, then a looping
         ** structure would be taken, in which the global loss data would first be
         ** reduced by the results from "fit_combined_local_domains()", followed by
         ** routines of "fit_local_domains()", to look for local peaks, and
         ** "fit_combined_domains()", to get better combined local fits, and finally
         ** fit_global_domain(), in order to find a global fit that describes the
         ** overall data nearly perfect **/
        
        counter=0;
        while (r2_global_pre<r2_threshold&&counter<3) { // NOTE: which threshold
            
            ++counter;
            ++n_reduction;
            
            reduce_raw_loss();
            
            fit_local_domains();          n_peaks_now   = n_peaks;
            fit_combined_local_domains(); r2_combined   = rep.r2;
            fit_global_domain();          r2_global_now = rep.r2;
            
            if(r2_global_now<=r2_global_pre) {
                /** if further global fit doesn't increase quality of fit (r2),
                 ** break the loop and use the highest achieved global fit   **/
                n_reduction -= 1;
                n_peaks      = n_peaks_pre;
                break;
            } else {
                r2_global_pre = r2_global_now;
                n_peaks_pre   = n_peaks_now;
            }
        }
    }
    
    /* write DS analysis info to file */
    //--------------------------------------------------------------------------
    if (n_peaks==0)
    {
        write_fit_HN(write_HN);
    }
    else
    {
        for (whichpeak=0; whichpeak<n_peaks; ++whichpeak)
        {
            update_fitinfo(whichpeak);
            write_peak_info(write_HN);
            write_fitModel(write_HN,func);
            write_fitarrays(write_HN);
            write_stopcond(write_HN);
            write_fitinfo(write_HN);
            write_fit_HN(write_HN);
            write_HN_params(write_params);
        }
        write_dc_loss(write_HN);
        write_combined_curve(write_HN);
        write_global_curve(write_HN);
        write_reducedLoss(write_HN);
        
        write_taueqToFile(write_taueq,T_actual);
        write_taueqToFile(write_alltaueq,T_actual);
        write_taueqToFile(write_taueqinv,(1000.0/(corrFacforT/precision))/T_actual);
        write_taueqToFile(write_alltaueqinv,(1000.0/(corrFacforT/precision))/T_actual);
    }
    //--------------------------------------------------------------------------
    cout_peaksinfo();
    write_peaksinfo(sysVar);
    write_localFitinfo(sysVar);
    write_combinedFitinfo(sysVar);
    write_globalFitinfo(sysVar);
    
    write_taueq.close();
    write_taueqinv.close();
    write_params.close();
    write_HN.close();
    write_alltaueq.close();
    write_alltaueqinv.close();
    
    tend=time(NULL);
    double dt=difftime(tend,tbeg);
    if (dt<1) {
        cout << "Time elapsed less than 1 sec.\n";
    } else {
        cout << "Time elapsed  "<<dt<<" sec\n";
    }
}





void FitDielectric::fit_tauFit(const StructureClass& sysVar)
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
        
        get_Theq(sysVar);//cf.get_Theq() of FitData
        set_fitParams(model);
        
        /* cutlowT determined manually */
        if (is_manual_cutlowT) {
            cutlowT=cutlowT_manual;
        }
        /* determine cutlowT by turning point at taueq */
        else {
            cutlowT=0;
            read_all_taueq_data(sysVar);
            calc_taueq_derivatives();
            decide_cutlowT(taueq_sortdecreasing);
        }
        
        read_all_taueq_data(sysVar);
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
        write_tauFit_vs_T(write_fit_tauFit);
        write_tauFit_vs_invT(write_fit_tauFit);
        write_tauFitfit_vs_T(write_fit_tauFit,c_org,Theq,model);
        write_tauFitfit_vs_invT(write_fit_tauFit,c_org,Theq,model);
        
        
        write_fit_tauFit.close();
    }
    else {
        cout
        << "FitDielectric::fit_tauFit: 'fit_tauFit.dat' cannot open." << "\n";
        //exit(EXIT_FAILURE);
    }
}





void FitDielectric::read_all_taueq_data(const StructureClass& sysVar)
{
    taueq_sortdecreasing.clear();
    taueq_sortincreasing.clear();
    vector<vector<double>> taueq;
    
    string i1;
    i1.append(return_AnalysisFolderPath(sysVar));
    i1.append("/fit_data");
    i1.append("/overall_taueq_HN_");
    i1.append(sysVar.get_usic());
    i1.append(".dat");
    
    ifstream readtaueqdata(i1.c_str());
    
    if (readtaueqdata.is_open()) {
        string lineContent;
        double Tdata=0,tauFitdata=0;
        while (getline(readtaueqdata,lineContent))
        {
            if (lineContent=="") continue;
            istringstream iss(lineContent); // (T,taueq)
            iss >> Tdata;
            if (is_manual_cutlowT) {
                if (Tdata>=cutlowT_manual) {// NOTE: >=
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            } else {
                if (Tdata>=cutlowT) { // NOTE: >=
                    iss >> tauFitdata;
                    taueq.push_back({Tdata,tauFitdata});
                }
            }
        } readtaueqdata.close();
    }
    else {
        cout
        << "in FitDielectric::read_all_taueq_data():\n"
        << "taueq datafile cannot open.\n";
    }
    
    try {
        if(taueq.size()==0) throw 0;
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
    }
    catch (int i) {
        cout
        << "in FitDielectric::read_all_taueq_data():\n"
        << "taueq size = " << i << "\n";
        taueq.push_back({0,0});
        std::sort(taueq.begin(),taueq.end(),sortIncreasing0);
        averaging_sorted_data(taueq,taueq_sortincreasing);
        std::sort(taueq.begin(),taueq.end(),sortDecreasing0);
        averaging_sorted_data(taueq,taueq_sortdecreasing);
    }
}





void FitDielectric::write_fit_HN(std::ofstream& outputFile)
{
    if (outputFile.is_open()) {
        
        if (n_peaks>0) {
            
            outputFile
            << "cutoff_lo_afreq " << log10(fabs(cutoff_lo_afreq.at(whichpeak))) << "\n"
            << "cutoff_lo_sfreq " << log10(fabs(cutoff_lo_sfreq.at(whichpeak))) << "\n\n";
            
            outputFile
            << "cutoff_hi_afreq " << log10(fabs(cutoff_hi_afreq.at(whichpeak))) << "\n"
            << "cutoff_hi_sfreq " << log10(fabs(cutoff_hi_sfreq.at(whichpeak))) << "\n\n";
            
            outputFile
            << "afreq_max " << log10(fabs(afreq_max.at(whichpeak))) << "\n"
            << "sfreq_max " << log10(fabs(sfreq_max.at(whichpeak))) << "\n"
            << "loss_max  " << log10(fabs(loss_max.at(whichpeak)))  << "\n\n";
            
            outputFile
            << "tauMax(afreq) "
            <<HN_tauFit_calc(c_losspeaks.at(whichpeak),func)[0] << "\n"
            << "log(tauMax_a) "
            <<log10(fabs(HN_tauFit_calc(c_losspeaks.at(whichpeak),func)[0]))<< "\n"
            << "tauMax(sfreq) "
            <<HN_tauFit_calc(c_losspeaks.at(whichpeak),func)[1] << "\n"
            << "log(tauMax_s) "
            <<log10(fabs(HN_tauFit_calc(c_losspeaks.at(whichpeak),func)[1]))<< "\n\n";
            
            
            double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
            double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
            
            
            if (whichpeak==0)
            {
                // Raw Data
                //--------------------------------------------------------------
                outputFile
                << "Raw Data "
                << "(#.="<< dielectric_freq_overall.size()<<")" << "\n"
                << "log.afreq  log.sfreq  log.loss  storage   " << "\n"
                << "==========================================" << "\n";
                for (size_t i=0; i<(int)dielectric_freq_overall.size(); ++i) {
                    outputFile
                    << log10(fabs(dielectric_freq_overall[i][0])) << " "
                    << log10(fabs(dielectric_freq_overall[i][1])) << " "
                    << log10(fabs(dielectric_freq_overall[i][2])) << " "
                    << dielectric_freq_overall[i][3]        << "\n";
                } outputFile << "\n\n";
                //--------------------------------------------------------------
                
                
                // DC Free Mode
                //--------------------------------------------------------------
                outputFile
                << "DC Free Mode "
                << "(#.="<< dielectric_freq_dcfree.size()<<")" << "\n"
                << "log.afreq  log.sfreq  log.loss  storage   " << "\n"
                << "==========================================" << "\n";
                for (size_t i=0; i<(int)dielectric_freq_dcfree.size(); ++i) {
                    outputFile
                    << log10(fabs(dielectric_freq_dcfree[i][0])) << " "
                    << log10(fabs(dielectric_freq_dcfree[i][1])) << " "
                    << log10(fabs(dielectric_freq_dcfree[i][2])) << " "
                    << dielectric_freq_dcfree[i][3]        << "\n";
                } outputFile << "\n\n";
                //--------------------------------------------------------------
                
                
                // Interpolated Data (Full Range)
                //--------------------------------------------------------------
                outputFile
                << "Cubic Spline Interpolation (Global Range), Raw Loss" << "\n"
                << "log.afreq  log.sfreq  f  df  d2f                   " << "\n"
                << "===================================================" << "\n";
                double f=0,df=0,d2f=0;
                for (int i=0; i<=equal_pieces; ++i) {
                    
                    afreqi=log_afreq_MIN+daf*i;
                    sfreqi=log_sfreq_MIN+dsf*i;
                    
                    alglib::spline1ddiff(interpolantRawLoss,afreqi,f,df,d2f);
                    
                    outputFile
                    << afreqi << " "
                    << sfreqi << " "
                    << f      << " "   // f is log(loss)
                    << df     << " "   // 1st derivative of log(loss) vs log(afreq)
                    << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
                } outputFile << "\n\n";
                
                outputFile
                << "Cubic Spline Interpolation (Global Range), Derivative Loss" << "\n"
                << "log.afreq  log.sfreq  f  df  d2f                          " << "\n"
                << "==========================================================" << "\n";
                for (int i=0; i<=equal_pieces; ++i) {
                    
                    afreqi=log_afreq_MIN+daf*i;
                    sfreqi=log_sfreq_MIN+dsf*i;
                    
                    alglib::spline1ddiff(interpolantDCLoss,afreqi,f,df,d2f);
                    
                    outputFile
                    << afreqi << " "
                    << sfreqi << " "
                    << f      << " "   // f is log(loss)
                    << df     << " "   // 1st derivative of log(loss) vs log(afreq)
                    << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
                } outputFile << "\n\n";
                //--------------------------------------------------------------
            }
            
            
            // Fit Data (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "##"<<whichpeak+1<<" DSFit Global Extrapolation" << "\n"
            << "Fit Data (Full Range)                         " << "\n"
            << "log.afreq  log.sfreq  log.DS_fit              " << "\n"
            << "==============================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << log10(fabs(HN_func_value(c_losspeaks.at(whichpeak),pow(10,afreqi))))
                << "\n";
            } outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            // Fit Data (Fitting Range)
            //------------------------------------------------------------------
            outputFile
            << "##"<<whichpeak+1<<" DSFit Range "
            << "(#.="<< dielectric_freq_sortincreasing.at(whichpeak).size()<<")" << "\n"
            << "log.afreq  log.sfreq  log.loss  log.DS_fit                     " << "\n"
            << "===============================================================" << "\n";
            double afreq_tmp=0,sfreq_tmp=0,loss_tmp=0;
            for (size_t i=0; i<(int)dielectric_freq_sortincreasing.at(whichpeak).size(); ++i) {
                
                afreq_tmp=dielectric_freq_sortincreasing.at(whichpeak)[i][0];
                sfreq_tmp=dielectric_freq_sortincreasing.at(whichpeak)[i][2];
                loss_tmp =dielectric_freq_sortincreasing.at(whichpeak)[i][1];
                
                outputFile
                << afreq_tmp << " "
                << sfreq_tmp << " "
                << loss_tmp  << " "
                << log10(fabs(HN_func_value(c_losspeaks.at(whichpeak),pow(10,afreq_tmp))))
                << "\n";
            } outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            // Distribution of tau (Full Range)
            //------------------------------------------------------------------
            double a=c_losspeaks.at(whichpeak)[1];
            double b=c_losspeaks.at(whichpeak)[2];
            double tau=pow(10,c_losspeaks.at(whichpeak)[3]);
            double t=0,y=0,theta=0,omega=0,f_t=0;
            
            outputFile
            << "##"<<whichpeak+1<<" DSFit                 " << "\n"
            << "Distribution of tau "
            << "(#.="<< dielectric_freq_overall.size()<<")" << "\n"
            << "log.time  log.f(time)                     " << "\n"
            << "==========================================" << "\n";
            for (int i=0; i<(int)dielectric_freq_overall.size(); ++i)
            {
                t=pow(dielectric_freq_overall[i][1],-1);
                y=t/tau;
                theta=atan(sin(pi()*a)/(pow(y,a)+cos(pi()*a)));
                omega=pow(1+2*pow(y,a)*cos(pi()*a)+pow(y,2*a),-b/2);
                f_t=pow(pi(),-1)*pow(y,a*b)*sin(b*theta)*omega;
                
                outputFile
                << log10(t)       << " "
                << log(fabs(f_t)) << "\n";
            } outputFile << "\n\n";
            //------------------------------------------------------------------
            
        }
        else {
            
            // Raw Data
            //------------------------------------------------------------------
            outputFile
            << "Raw Data "
            << "(#.="<< dielectric_freq_overall.size()<<")" << "\n"
            << "log.afreq  log.sfreq  log.loss  storage   " << "\n"
            << "==========================================" << "\n";
            for (size_t i=0; i<dielectric_freq_overall.size(); ++i)
            {
                outputFile
                << log10(fabs(dielectric_freq_overall[i][0])) << " "
                << log10(fabs(dielectric_freq_overall[i][1])) << " "
                << log10(fabs(dielectric_freq_overall[i][2])) << " "
                << dielectric_freq_overall[i][3]              << "\n";
            } outputFile << "\n\n";
            //------------------------------------------------------------------
            
            
            int    equal_pieces=1000;
            double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
            double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
            
            
            // Interpolated Data (Full Range)
            //------------------------------------------------------------------
            outputFile
            << "Cubic Spline Interpolation (Global Range)" << "\n"
            << "log.afreq  log.sfreq  f  df  d2f         " << "\n"
            << "=========================================" << "\n";
            double f=0,df=0,d2f=0;
            for (int i=0; i<=equal_pieces; ++i)
            {
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                
                alglib::spline1ddiff(interpolantRawLoss,afreqi,f,df,d2f);
                
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << f      << " "   // f is log(loss)
                << df     << " "   // 1st derivative of log(loss) vs log(afreq)
                << d2f    << "\n"; // 2nd derivative of log(loss) vs log(afreq)
            } outputFile << "\n\n";
            //------------------------------------------------------------------
        }
    }
    else {
        cout << "in FitDielectric::write_fit_HN(): file cannot open.\n";
    }
}





void FitDielectric::write_dc_loss(std::ofstream& outputFile)
{
    if (n_dcpoints>1)
    {
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        outputFile
        << "DC Loss Fit                            " << "\n"
        << "sig/eps  s                             " << "\n"
        << "=======================================" << "\n";
        for (int i=0; i<c_dc_loss.length(); ++i) outputFile << c_dc_loss[i] << " ";
        outputFile
        << "\n"
        << "r2 " << rep_dc_loss.r2 << "\n\n";
        
        outputFile
        << "DC Loss Fit Range "
        << "(#.="<< n_dcpoints <<")"                 << "\n"
        << "log.afreq  log.loss                    " << "\n"
        << "=======================================" << "\n";
        for (size_t i=0; i<n_dcpoints; ++i) {
            outputFile
            << dc_curve_freq_sortincreasing[i][0] << " "
            << dc_curve_freq_sortincreasing[i][1] << "\n";
        }   outputFile << "\n";
        
        outputFile
        << "Gloabl Range Extrapolation        " << "\n"
        << "log.afreq  log.sfreq  DCLoss_fit  " << "\n"
        << "==================================" << "\n";
        for (int i=0; i<=equal_pieces; ++i) {
            
            afreqi=log_afreq_MIN+daf*i;
            sfreqi=log_sfreq_MIN+dsf*i;
            
            outputFile
            << afreqi << " "
            << sfreqi << " "
            << log10(fabs(dc_loss_value(c_dc_loss,pow(10,afreqi)))) << "\n";
        } outputFile << "\n\n";
    }
    else {
        outputFile << "No DC (conductivity) contribution detected." << "\n";
        outputFile << "\n\n";
    }
}





void FitDielectric::write_combined_curve(std::ofstream& outputFile)
{
    if (dielectric_freq_domains.size()>0)
    {
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        for (int redu=0; redu<=n_reduction; ++redu) { // NOTE: <=
            outputFile
            << "##"<<redu+1<<" Combined Domain(s) Fit     " << "\n"
            << "c[i]                                      " << "\n"
            << "==========================================" << "\n";
            for (int i=0; i<c_combined.at(redu).length(); ++i) {
                outputFile << c_combined.at(redu)[i] << " ";
            }
            outputFile
            << "\n"
            << "r2 " << rep_combined.at(redu).r2 << "\n\n";
            
            outputFile
            << "##"<<redu+1<<" Domain(s) Fit Range "
            << "(#.="<< dielectric_freq_domains.at(redu).size()<<")" << "\n"
            << "log.afreq  log.loss                                " << "\n"
            << "===================================================" << "\n";
            for (size_t i=0; i<(int)dielectric_freq_domains.at(redu).size(); ++i) {
                outputFile
                << dielectric_freq_domains.at(redu)[i][0] << " "
                << dielectric_freq_domains.at(redu)[i][1] << "\n";
            }   outputFile << "\n";
            
            outputFile
            << "##"<<redu+1<<" Domain Fit Global Extrapolation     " << "\n"
            << "log.afreq  log.sfreq  fit_combined                 " << "\n"
            << "===================================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << log10(fabs(combined_func_value(c_combined.at(redu),pow(10,afreqi)))) << " "
                << "\n";
            } outputFile << "\n\n";
        }
    }
}





void FitDielectric::write_global_curve(std::ofstream& outputFile)
{
    if (dielectric_freq_global.size()>0)
    {
        int    equal_pieces=1000;
        double daf=fabs(log_afreq_MAX-log_afreq_MIN)/(double)equal_pieces;
        double dsf=fabs(log_sfreq_MAX-log_sfreq_MIN)/(double)equal_pieces;
        
        for (int redu=0; redu<=n_reduction; ++redu) { // NOTE: <=
            outputFile
            << "##"<<redu+1<<" Attempted Global Fit       " << "\n"
            << "c[i]                                      " << "\n"
            << "==========================================" << "\n";
            for (int i=0; i<c_global.at(redu).length(); ++i) {
                outputFile << c_global.at(redu)[i] << " ";
            }
            outputFile
            << "\n"
            << "r2 " << rep_global.at(redu).r2 << "\n\n";
            
            outputFile
            << "##"<<redu+1<<" Attempted Global Fit Range "
            << "(#.="<< dielectric_freq_global.at(redu).size()<<")" << "\n"
            << "log.afreq  log.loss                               " << "\n"
            << "==================================================" << "\n";
            for (size_t i=0; i<(int)dielectric_freq_global.at(redu).size(); ++i) {
                outputFile
                << dielectric_freq_global.at(redu)[i][0] << " "
                << dielectric_freq_global.at(redu)[i][1] << "\n";
            }   outputFile << "\n";
            
            outputFile
            << "##"<<redu+1<<" Global Fit Interpolation       " << "\n"
            << "log.afreq  log.sfreq  loss_globalFit          " << "\n"
            << "==============================================" << "\n";
            for (int i=0; i<=equal_pieces; ++i) {
                afreqi=log_afreq_MIN+daf*i;
                sfreqi=log_sfreq_MIN+dsf*i;
                outputFile
                << afreqi << " "
                << sfreqi << " "
                << log10(fabs(combined_func_value(c_global.at(redu),pow(10,afreqi))))
                << "\n";
            } outputFile << "\n\n";
        }
    }
}





void FitDielectric::write_reducedLoss(std::ofstream& outputFile)
{
    if (n_reduction>0)
    {
        for (int redu=0; redu<n_reduction; ++redu) {
            outputFile
            << "##"<<redu+1<<" Reduced Loss "
            << "(#.="<< dielectric_freq_reduced_all.at(redu).size()<<")" << "\n"
            << "log.afreq  log.sfreq  log.loss  storage                " << "\n"
            << "=======================================================" << "\n";
            for (int i=0; i<(int)dielectric_freq_reduced_all.at(redu).size(); ++i) {
                outputFile
                << log10(fabs(dielectric_freq_reduced_all.at(redu)[i][0])) << " "   // log(afreq)
                << log10(fabs(dielectric_freq_reduced_all.at(redu)[i][1])) << " "   // log(afreq)
                << log10(fabs(dielectric_freq_reduced_all.at(redu)[i][2])) << " "   // log(loss)
                << dielectric_freq_reduced_all.at(redu)[i][3]              << "\n"; // storage
            } outputFile << "\n\n";
        }
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





void FitDielectric::write_taueqToFile(std::ofstream& writeFile,
                                      const double Temp)
{
    /* write tau_eq data to separate files (cooling or heating) */
    writeFile << Temp << " ";
    double logafreq=0,logsfreq=0;
    for (int i=0; i<(int)c_losspeaks.size(); ++i) {
        logafreq=log10(fabs(HN_tauFit_calc(c_losspeaks[i],func)[0]));
        logsfreq=log10(fabs(HN_tauFit_calc(c_losspeaks[i],func)[1]));
        writeFile << logafreq << " " << logsfreq << " ";
    } writeFile << "(Species Name: "<<storedspecies.at(spc_id[0])<<")";
    writeFile << "\n";
}





void FitDielectric::write_blankToFile(std::ofstream& writeFile)
{
    writeFile << "\n";
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
    double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
    double m_extrp =calc_x_given_y(c,extrp_time,model)[1];
    return {Tg_extrp,m_extrp};
}





double FitDielectric::calc_y_given_x(const real_1d_array& c,
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





vector<double> FitDielectric::calc_x_given_y(const real_1d_array& c,
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
        << "in FitDielectric::calc_x_given_y, "<<model<<" model not found!\n"
        << "Program aborted.\n"; exit(EXIT_FAILURE);
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
        double Tg_extrp = calc_x_given_y(c,extrp_time,model)[0];
        double m_extrp  = calc_x_given_y(c,extrp_time,model)[1];
        //------------------------------------------------------------------
        
        /* computational Tg and fragility */
        //------------------------------------------------------------------
        //double Tg_compu = calc_x_given_y(c,compu_time,model)[0];
        //double m_compu  = calc_x_given_y(c,compu_time,model)[1];
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
                                         const double T_initial,
                                         const std::string& model)
{
    if (outputFile.is_open()) {
        
        //double T_initial=Theq;
        double Ti=0;
        
        /* format: T log10(tauFit_fit) */
        //------------------------------------------------------------------
        double Tg_extrp=calc_x_given_y(c,extrp_time,model)[0];
        
        outputFile << "\n"
        << "T  log10(tauFit_fit)             " << "\n"
        << "=================================" << "\n";
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
            
            outputFile
            << Ti << " "
            << log10(fabs(calc_y_given_x(c,Ti,model)))
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
                                            const double T_initial,
                                            const std::string& model)
{
    if (outputFile.is_open()) {
        
        //double T_initial=Theq;
        double Ti=0;
        
        /* format: 1/T log10(tauFit_fit) */
        //------------------------------------------------------------------
        double Tg_extrp = calc_x_given_y(c,extrp_time,model)[0];
        
        outputFile << "\n"
        << 1000.0/(corrFacforT/precision)
        << "/T  log10(tauFit_fit)            " << "\n"
        << "=================================" << "\n";
        
        for (int i=0; i<=equal_pieces; ++i) {
            
            Ti=T_initial-(fabs(Tg_extrp-T_initial)/(double)equal_pieces)*i;
            
            outputFile
            << (1000.0/(corrFacforT/precision))/Ti << " "
            << log10(fabs(calc_y_given_x(c,Ti,model)))
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
    //cf.get_Theq() of FitData
    read_all_taueq_data(sysVar);
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
    else if (model=="combined")
    {
        alglib::lsfitfit(state, static_combined_func);
    }
    else if (model=="global")
    {
        alglib::lsfitfit(state, static_combined_func);
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
        //HN: shape params are free
        
        //NOTE:
        //tau uses log scale
        
        set_coeffs({       1.0, 0.5, 0.5, -8.0});
        set_coeffs_scale({ 1.0, 0.5, 0.5, -8.0});
        set_coeffs_bndl({  0.0, 0.0, 0.0, -inf});
        set_coeffs_bndu({  inf, 1.0, 1.0,  inf});
    }
    else if (model=="CC_loss") // [A,a,b,tau]
    {
        //Cole-Cole: special case of HN where beta=1
        
        //NOTE:
        //tau uses log scale
        
        set_coeffs({       1.0, 0.5, 1.0, -8.0});
        set_coeffs_scale({ 1.0, 0.5, 1.0, -8.0});
        set_coeffs_bndl({  0.0, 0.0, 0.0, -inf});
        set_coeffs_bndu({  inf, 1.0, 1.0,  inf});
    }
    else if (model=="CD_loss") // [A,a,b,tau]
    {
        //Cole-Davidson: special case of HN where alpha=1
        
        //NOTE:
        //tau uses log scale
        
        set_coeffs({       1.0, 1.0, 0.5, -8.0});
        set_coeffs_scale({ 1.0, 1.0, 0.5, -8.0});
        set_coeffs_bndl({  0.0, 0.0, 0.0, -inf});
        set_coeffs_bndu({  inf, 1.0, 1.0,  inf});;
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
    else if (model=="combined")
    {
        vector<double> coeffs,coeffs_scale,coeffs_bndl,coeffs_bndu;
        
        /* dc part */
        if (n_dcpoints>1) {
            for (int i=0; i<(int)c_dc_loss.length(); ++i)
            {
                coeffs.push_back(c_dc_loss[i]);
                coeffs_scale.push_back(c_dc_loss[i]);
                coeffs_bndl.push_back(-inf);
                coeffs_bndu.push_back(inf);
            }
        }
        /* HN fits */
        if (n_peaks>0) {
            for (int i=0; i<n_peaks; ++i) {
                for (int ii=0; ii<(int)c_losspeaks.at(i).length(); ++ii)
                {
                    coeffs.push_back(c_losspeaks.at(i)[ii]);
                    coeffs_scale.push_back(c_losspeaks.at(i)[ii]);
                    coeffs_bndl.push_back(-inf);
                    coeffs_bndu.push_back(inf);
                }
            }
        }
        set_coeffs(coeffs);
        set_coeffs_scale(coeffs_scale);
        set_coeffs_bndl(coeffs_bndl);
        set_coeffs_bndu(coeffs_bndu);
    }
    else if (model=="global")
    {
        vector<double>
        coeffs,coeffs_scale,coeffs_bndl,coeffs_bndu;
        for (int i=0; i<(int)c_combined.at(n_reduction).length(); ++i)
        {
            coeffs.push_back(c_combined.at(n_reduction)[i]);
            coeffs_scale.push_back(c_combined.at(n_reduction)[i]);
            coeffs_bndl.push_back(-inf);
            coeffs_bndu.push_back(inf);
        }
        set_coeffs(coeffs);
        set_coeffs_scale(coeffs_scale);
        set_coeffs_bndl(coeffs_bndl);
        set_coeffs_bndu(coeffs_bndu);
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
        << "in FitDielectric::set_fitParams, "<<model<<" model not found!\n"
        << "Program aborted.\n"; exit(EXIT_FAILURE);
    }
}



void FitDielectric::fitParams_correction(const std::string& model)
{
    if ((model=="HN_loss")||(model=="CC_loss")||(model=="CD_loss"))
    {
        if (cyclecount>0)
        {
            coeffs_vD[3] += 0.1; // NOTE: log scale
            set_coeffs(coeffs_vD);
            set_coeffs_scale(coeffs_vD);
            c = get_coeffs().c_str();
        }
    }
    if (model=="combined")
    {
        if (cyclecount>0&&coeffs_vD.size()>=4)
        {
            int tauindex=(int)coeffs_vD.size()-1;
            coeffs_vD[tauindex] += 0.1;
            set_coeffs(coeffs_vD);
            set_coeffs_scale(coeffs_vD);
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
    
    //cout << "x \n" << fit_xData.c_str()    << "\n";
    //cout << "y \n" << fit_yData.c_str()    << "\n";
    //cout << "c \n" << get_coeffs().c_str() << "\n";
    
    cyclecount=0;
    usedt.push_back(time(NULL));
    do
    {
        //cout << "c \n" << get_coeffs().c_str() << "\n";
        
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





void FitDielectric::build_function_interpolant(const string& function,
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





void FitDielectric::dc_loss(const alglib::real_1d_array &c,
                            const alglib::real_1d_array &x,
                            double &func,
                            void *ptr)
{
    double A=c[0],s=c[1],afreq=pow(10,x[0]);
    func = log10(fabs(A*pow(afreq,-s))); // log form
}





void FitDielectric::static_combined_func(const alglib::real_1d_array &c,
                                         const alglib::real_1d_array &x,
                                         double &func,
                                         void *ptr)
{
    double afreq=pow(10,x[0]);
    double A_HN=0,a_HN=0,b_HN=0,tau_HN=0,r_w=0,p_w=0;
    double dc=0,peak=0;
    int    n_peaks=0,indexii=0;
    vector<double> peakii;
    
    func = 0;
    
    /* based on the length of c params, there could be following possibilities:
     1. 0 dc part,  0 peak  (length=0)
     2. 0 dc part,  1 peak  (length=4)
     3. 0 dc part, >1 peaks (length=4x)
     4. 1 dc part,  0 peak  (length=2)
     5. 1 dc part,  1 peak  (length=6)
     6. 1 dc part, >1 peaks (length=2+4x) */
    
    if (c.length()==0) // 0 dc part, 0 peak (0)
    {
        func = 0;
    }
    else if (c.length()==2) // 1 dc part, 0 peak (2)
    {
        dc   = log10(fabs(c[0]*pow(afreq,-c[1])));
        func = log10(pow(10,dc));
    }
    else if (c.length()==4) // 0 dc part, 1 peak (4)
    {
        A_HN   = c[0];
        a_HN   = c[1];
        b_HN   = c[2];
        tau_HN = pow(10,c[3]);
        r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
        p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
        peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
        func = log10(pow(10,peak));
    }
    else if (c.length()==6) // 1 dc part, 1 peak (6)
    {
        dc     = log10(fabs(c[0]*pow(afreq,-c[1])));
        A_HN   = c[2];
        a_HN   = c[3];
        b_HN   = c[4];
        tau_HN = pow(10,c[5]);
        r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
        p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
        peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
        func = log10(pow(10,dc)+pow(10,peak));
    }
    else if (c.length()>4 && (c.length()%4)==0) // 0 dc part, >1 peaks (4x)
    {
        n_peaks = (int)c.length()/4;
        for (int i=0; i<n_peaks; ++i)
        {
            indexii= i*4;
            A_HN   = c[indexii];
            a_HN   = c[indexii+1];
            b_HN   = c[indexii+2];
            tau_HN = pow(10,c[indexii+3]);
            r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
            p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
            peak  = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
            peakii.push_back(peak);
        }
        peak=0;
        for (int i=0; i<n_peaks; ++i)
        {
            peak += pow(10,peakii.at(i));
        }
        func = log10(peak);
    }
    else if(c.length()>6) // dc part + at least 2 peaks
    {
        dc = log10(fabs(c[0]*pow(afreq,-c[1])));
        
        n_peaks = (int)(c.length()-2)/4;
        for (int i=0; i<n_peaks; ++i)
        {
            indexii= 2 + i*4;
            A_HN   = c[indexii];
            a_HN   = c[indexii+1];
            b_HN   = c[indexii+2];
            tau_HN = pow(10,c[indexii+3]);
            r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
            p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
            peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
            peakii.push_back(peak);
        }
        peak=pow(10,dc);
        for (int i=0; i<n_peaks; ++i)
        {
            peak += pow(10,peakii.at(i));
        }
        func = log10(peak);
    }
    else
    {
        func = 0;
    }
}





double FitDielectric::HN_func_value(const real_1d_array& c,
                                    const double afreq)
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
        << "\n"; exit(EXIT_FAILURE);
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
        << "\n"; exit(EXIT_FAILURE);
    }
    return {tauFit_afreq,tauFit_sfreq};
}





double FitDielectric::dc_loss_value(const alglib::real_1d_array& c,
                                    const double afreq)
{
    double A=c[0],s=c[1];
    return A*pow(afreq,-s);
}





double FitDielectric::combined_func_value(const alglib::real_1d_array& c,
                                          const double afreq)
{
    double A_HN=0,a_HN=0,b_HN=0,tau_HN=0,r_w=0,p_w=0;
    double dc=0,peak=0;
    int    n_peaks=0,indexii=0;
    vector<double> peakii;
    
    double func = 0;
    
    /* based on the length of c params, there could be following possibilities:
     1. 0 dc part,  0 peak  (length=0)
     2. 0 dc part,  1 peak  (length=4)
     3. 0 dc part, >1 peaks (length=4x)
     4. 1 dc part,  0 peak  (length=2)
     5. 1 dc part,  1 peak  (length=6)
     6. 1 dc part, >1 peaks (length=2+4x) */
    
    if (c.length()==0) // 0 dc part, 0 peak (0)
    {
        func = 0;
    }
    else if (c.length()==2) // 1 dc part, 0 peak (2)
    {
        dc   = log10(fabs(c[0]*pow(afreq,-c[1])));
        func = log10(pow(10,dc));
    }
    else if (c.length()==4) // 0 dc part, 1 peak (4)
    {
        A_HN   = c[0];
        a_HN   = c[1];
        b_HN   = c[2];
        tau_HN = pow(10,c[3]);
        r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
        p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
        peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
        func = log10(pow(10,peak));
    }
    else if (c.length()==6) // 1 dc part, 1 peak (6)
    {
        dc     = log10(fabs(c[0]*pow(afreq,-c[1])));
        A_HN   = c[2];
        a_HN   = c[3];
        b_HN   = c[4];
        tau_HN = pow(10,c[5]);
        r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
        p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
        peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
        func = log10(pow(10,dc)+pow(10,peak));
    }
    else if (c.length()>4 && (c.length()%4)==0) // 0 dc part, >1 peaks (4x)
    {
        n_peaks = (int)c.length()/4;
        for (int i=0; i<n_peaks; ++i)
        {
            indexii= i*4;
            A_HN   = c[indexii];
            a_HN   = c[indexii+1];
            b_HN   = c[indexii+2];
            tau_HN = pow(10,c[indexii+3]);
            r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
            p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
            peak  = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
            peakii.push_back(peak);
        }
        peak=0;
        for (int i=0; i<n_peaks; ++i)
        {
            peak += pow(10,peakii.at(i));
        }
        func = log10(peak);
    }
    else if(c.length()>6) // dc part + at least 2 peaks
    {
        dc = log10(fabs(c[0]*pow(afreq,-c[1])));
        
        n_peaks = (int)(c.length()-2)/4;
        for (int i=0; i<n_peaks; ++i)
        {
            indexii= 2 + i*4;
            A_HN   = c[indexii];
            a_HN   = c[indexii+1];
            b_HN   = c[indexii+2];
            tau_HN = pow(10,c[indexii+3]);
            r_w=pow(1+2*pow(afreq*tau_HN,a_HN)*cos(pi()*a_HN/2)+pow(afreq*tau_HN,2*a_HN),-b_HN/2);
            p_w=atan(sin(pi()*a_HN/2)/(pow(afreq*tau_HN,-a_HN)+cos(pi()*a_HN/2)));
            peak = log10(fabs(A_HN*r_w*sin(b_HN*p_w)));
            peakii.push_back(peak);
        }
        peak=pow(10,dc);
        for (int i=0; i<n_peaks; ++i)
        {
            peak += pow(10,peakii.at(i));
        }
        func = log10(peak);
    }
    else
    {
        func = 0;
    }
    return pow(10,func);
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
    grad[1]= ((((T*u*exponent)/E_inf)+exponent+1)/T)*log10(exp(1));
    grad[2]= ((E_inf*(b-T/E_inf)*exponent)/T)*log10(exp(1));
    grad[3]= ((E_inf*u*exponent)/T)*log10(exp(1));
}




void FitDielectric::initialize()
{
    if(is_manualFitSingleDS)
    {
        form="manual";
        set_is_log_data({true,false});
    }
    if (form=="NPIC_ht")
    {
        set_process("cooling");//"cooling" or "heating"
        
        //if process is set cooling: T will be sorted from high to low
        //if procses is set heating: T will be sorted from low to high
        
        if (process=="cooling") {
            process_id=0;
        } else if (process=="heating") {
            process_id=1;
        } else {
            cout
            <<"in FitDielectric::initialize():\n"
            <<"process ("<<process<<") not found!\n";
            exit(EXIT_FAILURE);
        }
    }
}





void FitDielectric::transfer_raw_data(const StructureClass& sysVar)
{
    /**  Procedure:
     **  1. create folder for dielectric data
     **  2. copy dielectric data from source to destination
     **  3. rename files
     **/
    
    static int countaccess=0;
    ++countaccess;
    
    if (countaccess==1)
    {
        /* create folder for dielectric data (if not exist) */
        string dir,mkdir;
        dir.append(return_AnalysisFolderPath(sysVar));
        dir.append("/statistics/dielectric");
        mkdir=
        "test -e "+dir+";"
        "if [ $? -ne 0 ]; then "
        "   mkdir "+dir+";"
        "else "
        "   cd "+dir+";rm *.dat;"
        "fi"; system(mkdir.c_str());
        
        /* copy dielectric data from source to destination */
        string path,cp;
        path.append(sysVar.get_Path());
        path.append("/dielectric");
        cp="cd "+path+";cp * "+dir;
        system(cp.c_str());
        
        /* rename files (name to ${species}_*) */
        string mv;
        if (form=="NPIC_ht") {
            //
        } else if (form=="NPIC_cs") {
            mv="cd "+dir+";for i in *;do mv $i "+sysVar.get_usic()+".dat;done";
        } else if (form=="NIST") {
            mv="cd "+dir+";for i in *;do mv $i "+sysVar.get_usic()+".dat;done";
        } else if (form=="manuallySet") {
            mv="cd "+dir+";for i in *;do mv $i "+sysVar.get_usic()+".dat;done";
        } system(mv.c_str());
    }
}





/** public setters **/
//------------------------------------------------------------------------------
/* string */
void FitDielectric::set_form(const std::string& str){form=str;}
void FitDielectric::set_func(const std::string& str){func=str;}
void FitDielectric::set_model(const std::string& str){model=str;};
void FitDielectric::set_process(const std::string& str){process=str;}
/* bool */
void FitDielectric::set_is_use_dc_free(const bool b){is_use_dc_free=b;}
void FitDielectric::set_is_manual_cutlowT(const bool b){is_manual_cutlowT=b;}
void FitDielectric::set_is_manualFitSingleDS(const bool b){is_manualFitSingleDS=b;}
void FitDielectric::set_process_id(const int i){process_id=i;}
/* double */
void FitDielectric::set_cutlowT_manual(const double d){cutlowT_manual=d;}
/* STL */
void FitDielectric::set_is_log_data(const std::vector<bool>& vb){is_log_data=vb;}
void FitDielectric::set_spc_id(const vector<int>& vi){spc_id=vi;}
void FitDielectric::set_T_id(const vector<int>& vi){T_id=vi;}
void FitDielectric::set_sfreq_range(const std::vector<double>& vd){sfreq_range=vd;}




