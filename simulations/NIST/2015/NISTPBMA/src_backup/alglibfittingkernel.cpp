//
//  alglibfittingkernel.cpp
//  cppWork
//
//  Created by SJH on 11/11/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#include "alglibfittingkernel.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

AlglibFittingKernel::AlglibFittingKernel()
{
    NewtonGuessTemp=100.0;
}





string AlglibFittingKernel::get_fit_info(const int info_index)
{
    string infostr;
    switch (info_index)
    {
        case (-7):
        {
            infostr="gradient verification failed.";
        }
            break;
        case (1):
        {
            infostr="relative function improvement is no more than EpsF.";
        }
            break;
        case (2):
        {
            infostr="relative step is no more than EpsX.";
        }
            break;
        case (4):
        {
            infostr="gradient norm is no more than EpsG";
        }
            break;
        case (5):
        {
            infostr="MaxIts steps was taken";
        }
            break;
        case (7):
        {
            infostr="stopping conditions are too stringent, further improvement is impossible";
        }
            break;
    }
    return infostr;
}





void AlglibFittingKernel::fitdata_processing_lsfit(const vector<vector<double>>& fitVec)
{
    xdata.clear();
    ydata.clear();
    xraw.clear();
    yraw.clear();
    xorg.clear();
    yorg.clear();
    fit_xData.clear();
    fit_yData.clear();
    
    double x=0,y=0;
    
    xdata.push_back("[[");
    ydata.push_back("[");
    
    for (size_t i=0; i<fitVec.size(); ++i) {
        
        x = fitVec[i][0]; // 1st column data
        y = fitVec[i][1]; // 2nd column data
        
        // x_data to fit
        xdata.push_back(to_string((long double)x));
        xdata.push_back("],[");
        xraw.push_back(to_string((long double)x));
        xorg.push_back(x);
        
        // y_data to fit
        ydata.push_back(to_string((long double)y));
        ydata.push_back(",");
        yraw.push_back(to_string((long double)y));
        yorg.push_back(y);
        
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





void AlglibFittingKernel::fitdata_processing_1d(const vector<vector<double>>& fitVec)
{
    xdata.clear();
    ydata.clear();
    xraw.clear();
    yraw.clear();
    xorg.clear();
    yorg.clear();
    fit_xData.clear();
    fit_yData.clear();
    
    double x=0,y=0;
    
    xdata.push_back("[");
    ydata.push_back("[");
    
    for (size_t i=0; i<fitVec.size(); ++i) {
        
        x = fitVec[i][0]; // 1st column data
        y = fitVec[i][1]; // 2nd column data
        
        // x_data to fit
        xdata.push_back(to_string((long double)x));
        xdata.push_back(",");
        xraw.push_back(to_string((long double)x));
        xorg.push_back(x);
        
        // y_data to fit
        ydata.push_back(to_string((long double)y));
        ydata.push_back(",");
        yraw.push_back(to_string((long double)y));
        yorg.push_back(y);
        
    }
    xdata.pop_back(); // rm appended "],["
    ydata.pop_back(); // rm appended ","
    
    xdata.push_back("]");
    ydata.push_back("]");
    for (size_t i=0; i<xdata.size(); ++i) {
        fit_xData += xdata[i];
        fit_yData += ydata[i];
    }
    
    //cout << "\n" << fit_xData << "\n";
    //cout << "\n" << fit_yData << "\n";
}





void AlglibFittingKernel::write_fitModel(std::ofstream& outputFile,
                                         const std::string& model)
{
    outputFile
    << "Data fit to <" << model << "> functional form"
    << "\n\n";
}





void AlglibFittingKernel::write_fitarrays(std::ofstream& outputFile)
{
    if (outputFile.is_open())
    {
        string xData,yData;
        
        for (size_t i=0; i<xdata.size(); ++i) {
            xData += xdata[i];
            yData += ydata[i];
        }
        
        outputFile
        << "xData (#.=" << (xdata.size()-1)/2 << ")"
        << "\n"
        << xData
        << "\n\n"
        
        << "yData (#.=" << (ydata.size()-1)/2 << ")"
        << "\n"
        << yData
        << "\n\n";
    }
    else {
        cout
        << "write_fitarrays: file cannot open."
        << "\n";
    }
}





void AlglibFittingKernel::write_stopcond(std::ofstream& outputFile)
{
    
    if (outputFile.is_open()) {
        
        outputFile
        << "stop condition:" << "\n"
        << "intial coeffs = ";
        for (int i=0; i<get_coeffs_vD().size(); ++i) {
            outputFile << get_coeffs_vD()[i] << " ";
        }
        outputFile
        << "\n"
        << "scale  = ";
        for (int i=0; i<get_coeffs_scale_vD().size(); ++i) {
            outputFile << get_coeffs_scale_vD()[i] << " ";
        }
        outputFile
        << "\n"
        << "bndl   = ";
        for (int i=0; i<get_coeffs_bndl_vD().size(); ++i) {
            outputFile << get_coeffs_bndl_vD()[i] << " ";
        }
        outputFile
        << "\n"
        << "bndu   = ";
        for (int i=0; i<get_coeffs_bndu_vD().size(); ++i) {
            outputFile << get_coeffs_bndu_vD()[i] << " ";
        }
        outputFile
        << "\n"
        << "epsf   = " << epsf << "\n"
        << "epsx   = " << epsx << "\n"
        << "maxits = " << maxits
        << "\n\n";
        
    }
    else {
        cout
        << "write_stopcond: file cannot open."
        << "\n";
    }
}





void AlglibFittingKernel::write_fitinfo(std::ofstream& outputFile)
{
    if (outputFile.is_open()) {
        
        if(get_is_use_FG()) outputFile << "Operating Mode: FG";
        else outputFile << "Operating Mode: F";
        
        double timediff=difftime(usedt[1],usedt[0]);
        
        outputFile
        << "\n\n"
        
        << "fit info " << int(info)                               << "\n"
        << "(" << get_fit_info(int(info)) << ")"                  << "\n\n"
        
        << "internal iterations = "   << rep.iterationscount      << "\n"
        << "time used = " << timediff << " seconds"               << "\n\n";
        
        outputFile
        << "fit coeffs"                                           << "\n";
        for (size_t i=0; i<c.length(); ++i) {
            outputFile << c[i] << " ";
        }
        outputFile << "\n\n"
        
        << "fit coeffs error"                                     << "\n";
        for (size_t i=0; i<c.length(); ++i) {
            outputFile << rep.errpar[i] << " ";
        }
        outputFile << "\n\n";
        
        outputFile
        << "Statistics"                                           << "\n"
        << "==================================="                  << "\n"
        << "R^2 (coefficient of determination) "                  << "\n"
        << rep.r2 << " ("<<cyclecount<<" cycles)"                 << "\n\n"
        
        << "RMS error "                                           << "\n"
        << rep.rmserror                                           << "\n\n"
        
        << "average error "                                       << "\n"
        << rep.avgerror                                           << "\n\n"
        
        << "average relative error "                              << "\n"
        << rep.avgrelerror                                        << "\n\n"
        
        << "maximum error "                                       << "\n"
        << rep.maxerror                                           << "\n"
        << "==================================="
        << "\n\n";
        
    }
    else {
        cout
        << "write_fitinfo: file cannot open."
        << "\n";
    }
}




void AlglibFittingKernel::write_errcurve(std::ofstream& outputFile,
                                         const std::string& model)
{
    if (outputFile.is_open()) {
        
        outputFile
        << "\n"
        << "T  log(tau_org) log10(tauFit_fit) ErrCurve  " << "\n"
        << "=========================================== " << "\n";
        for (int i=0; i<rep.errcurve.length(); ++i) {
            
            outputFile
            << taueq_sortdecreasing[i][0] << " "
            << taueq_sortdecreasing[i][1] << " "
            << log10(calc_time_given_temp(c,taueq_sortdecreasing[i][0],model)) << " "
            << rep.errcurve[i] << "\n";
        }
        outputFile
        << "\n\n";
        
        outputFile
        << 1000.0/(corrFacforT/precision)
        << "/T  log(tau_org) log10(tauFit_fit) ErrCurve " << "\n"
        << "=========================================== " << "\n";
        for (int i=0; i<rep.errcurve.length(); ++i) {
            
            outputFile
            << (1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0] << " "
            << taueq_sortdecreasing[i][1] << " "
            << log10(calc_time_given_temp(c,taueq_sortdecreasing[i][0],model)) << " "
            << rep.errcurve[i] << "\n";
        }
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "write_errcurve: file cannot open."
        << "\n";
    }
}





void AlglibFittingKernel::write_tauFit_vs_T(std::ofstream& outputFile)
{
    if (outputFile.is_open()) {
        
        /* format: T log10(tauFit_fit) */
        //------------------------------------------------------------------
        outputFile << "\n"
        
        << "T  log10(tauFit)                 " << "\n"
        << "=================================" << "\n";
        for (size_t i=0; i<taueq_sortdecreasing.size(); ++i) {
            outputFile
            << taueq_sortdecreasing[i][0] << " "
            << taueq_sortdecreasing[i][1] << "\n";
        }
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
    }
    else {
        cout
        << "write_tauFit_vs_T: file cannot open."
        << "\n";
    }
}





void AlglibFittingKernel::write_tauFit_vs_invT(std::ofstream& outputFile)
{
    if (outputFile.is_open()) {
        
        /* format: 1/T log10(tauFit) */
        //------------------------------------------------------------------
        outputFile
        << 1000.0/(corrFacforT/precision)
        << "/T  log10(tauFit)                " << "\n"
        << "=================================" << "\n";
        for (size_t i=0; i<taueq_sortdecreasing.size(); ++i) {
            outputFile
            << (1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0] << " "
            << taueq_sortdecreasing[i][1]
            << "\n";
        }
        //------------------------------------------------------------------
        outputFile
        << "\n\n";
        
    }
    else {
        cout
        << "write_tauFit_vs_invT: file cannot open."
        << "\n";
    }
}





void AlglibFittingKernel::write_badR2(std::ofstream& outputFile,
                                      const double r2_threshold)
{
    outputFile
    << "\n\n"
    << "==================================="
    << "\n"
    << "NOTE:"
    << "\n"
    << "r2 = " << rep.r2 << " < " << r2_threshold
    << "\n"
    << "tauFit not written into taueq file!"
    << "\n"
    << "==================================="
    << "\n";
}





void AlglibFittingKernel::write_insufficientDara(std::ofstream& outputFile,
                                                 const int n_fit,
                                                 const int n_threshold)
{
    outputFile
    << "\n\n"
    << "==================================="
    << "\n"
    << "Insufficient data for curve-fitting:"
    << "\n"
    << "Data for fitting = "
    << n_fit << " points" << "\n"
    << "less than " << n_threshold << " points threshold"
    << "\n"
    << "==================================="
    << "\n";
}





void AlglibFittingKernel::set_coeffs(const vector<double>& vD)
{
    /* store as vector<double> type */
    coeffs_vD = vD;
    
    /* vector<double> to string type */
    vector<string> strVec;
    strVec.push_back("[");
    for (size_t i=0; i<vD.size(); ++i) {
        strVec.push_back(to_string((long double)vD[i]));
        strVec.push_back(",");
    }
    strVec.pop_back(); // rm ","
    strVec.push_back("]");
    
    string str;
    for (size_t i=0; i<strVec.size(); ++i) str += strVec[i];
    
    coeffs = str;
}
void AlglibFittingKernel::set_coeffs_scale(const vector<double>& vD)
{
    /* store as vector<double> type */
    coeffs_scale_vD = vD;
    
    /* vector<double> to string type */
    vector<string> strVec;
    strVec.push_back("[");
    for (size_t i=0; i<vD.size(); ++i) {
        strVec.push_back(to_string((long double)vD[i]));
        strVec.push_back(",");
    }
    strVec.pop_back(); // rm ","
    strVec.push_back("]");
    
    string str;
    for (size_t i=0; i<strVec.size(); ++i) str += strVec[i];
    
    coeffs_scale = str;
}
void AlglibFittingKernel::set_coeffs_bndl(const vector<double>& vD)
{
    /* store as vector<double> type */
    coeffs_bndl_vD = vD;
    
    /* vector<double> to string type */
    vector<string> strVec;
    strVec.push_back("[");
    for (size_t i=0; i<vD.size(); ++i) {
        strVec.push_back(to_string((long double)vD[i]));
        strVec.push_back(",");
    }
    strVec.pop_back(); // rm ","
    strVec.push_back("]");
    
    string str;
    for (size_t i=0; i<strVec.size(); ++i) str += strVec[i];
    
    coeffs_bndl=str;
}
void AlglibFittingKernel::set_coeffs_bndu(const vector<double>& vD)
{
    /* store as vector<double> type */
    coeffs_bndu_vD = vD;
    
    /* vector<double> to string type */
    vector<string> strVec;
    strVec.push_back("[");
    for (size_t i=0; i<vD.size(); ++i) {
        strVec.push_back(to_string((long double)vD[i]));
        strVec.push_back(",");
    }
    strVec.pop_back(); // rm ","
    strVec.push_back("]");
    
    string str;
    for (size_t i=0; i<strVec.size(); ++i) str += strVec[i];
    
    coeffs_bndu=str;
}





void AlglibFittingKernel::averaging_sorted_data(const vector<vector<double>>& sorted_org,
                                                vector<vector<double>>& sorted_avg)
{
    //--------------------------------------------------------------------------
    // NOTE:
    // the algorithm is only valid if the 1st column has been sorted in a
    // decreasing or increasing order.
    // i.e. if the data is random, the algorithm will not give a right answer.
    //--------------------------------------------------------------------------
    sorted_avg.clear();
    
    int n_data=(int)sorted_org.size();
    int n_columns=(int)sorted_org[0].size();
    
    vector<double> sum;
    vector<double> finaldata;
    
    if (n_columns>1) { // ex. {T,tau}
        
        // initialization
        for (int j=0; j<(n_columns-1); ++j) sum.push_back(0);
        
        int count=0,i=0,i_beg=0;
        do {
            
            i=i_beg;
            for (int j=0; j<(n_columns-1); ++j) sum[j]=sorted_org[i][j+1];
            // this stores at row i the value of each column except the first column
            
            count=0;
            for (int ii=i+1; ii<n_data; ++ii) {
                if (sorted_org[i][0]==sorted_org[ii][0]) { // if 1st column the same
                    ++count;
                    for (int j=0; j<(n_columns-1); ++j) sum[j] += sorted_org[ii][j+1];
                    i_beg=ii+1;
                }
            }
            if (count==0) ++i_beg;
            for (int j=0; j<(n_columns-1); ++j) sum[j] /= (double)(count+1);
            
            finaldata.clear();
            finaldata.push_back(sorted_org[i][0]);
            for (int j=0; j<(n_columns-1); ++j) finaldata.push_back(sum[j]);
            
            sorted_avg.push_back(finaldata);
            
        } while(i_beg<n_data);
    }
    else { // ex. {T}
        
        int count=0,i=0,i_beg=0;
        do {
            
            i=i_beg;
            
            count=0;
            for (int ii=i+1; ii<n_data; ++ii) {
                if (sorted_org[i][0]==sorted_org[ii][0]) { // if 1st column the same
                    ++count;
                    i_beg=ii+1;
                }
            }
            if (count==0) ++i_beg;
            
            finaldata.clear();
            finaldata.push_back(sorted_org[i][0]);
            
            sorted_avg.push_back(finaldata);
            
        } while(i_beg<n_data);
    }
}





/* setters */
void AlglibFittingKernel::set_is_use_FG(const bool b){is_use_FG=b;}

