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

#include "alglibfittingkernel.h"
#include "sysinfo.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

AlglibFittingKernel::AlglibFittingKernel():

/* bool */
is_use_FG(true), // true: use FG-form; false: use F-form

/* int */
equal_pieces(1000),
cyclecount(0),
indexi(0),
indexii(0),

/* double */
NewtonGuessTemp(100.0),
epsf(0),
epsx(0),
diffstep(0),
corrFacforT(0),
precision(0),
convert(1.0),
extrp_time(0),
compu_time(0),

/* string */
systemUnit("real"),

/** ALGLIB fitting DS **/
//----------------------------------------------------------------------
xdata(),
ydata(),
xraw(),
yraw(),
xorg(),
yorg(),
fit_xData(),
fit_yData(),
usedt(),
x(),
y(),
c(),
s(),
bndl(),
bndu(),
maxits(),
info(),
state(),
rep(),
coeffs(),
coeffs_scale(),
coeffs_bndl(),
coeffs_bndu(),
coeffs_vD(),
coeffs_scale_vD(),
coeffs_bndl_vD(),
coeffs_bndu_vD(),

/** STL containers **/
//----------------------------------------------------------------------
taueq_sortdecreasing(),  // (Teq,taueq)
taueq_sortincreasing(),  // (Teq,taueq)
temps_sortdecreasing(),  // temperatures
temps_sortincreasing(),  // temperatures
taueq2d_sortdecreasing(),
standard2dcontainer()
{
    /* Constructor Assignment of AlglibFittingKernel */
}





const string AlglibFittingKernel::get_fit_info(const int info_index)
{
    string infostr;
    switch (info_index)
    {
        case (-7): {
            infostr="gradient verification failed.";
        } break;
        case (1): {
            infostr="relative function improvement is no more than EpsF.";
        } break;
        case (2): {
            infostr="relative step is no more than EpsX.";
        } break;
        case (4): {
            infostr="gradient norm is no more than EpsG";
        } break;
        case (5): {
            infostr="MaxIts steps was taken";
        } break;
        case (7): {
            infostr="stopping conditions are too stringent, further improvement is impossible";
        } break;
    } return infostr;
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
    for (size_t i=0; i<fitVec.size(); ++i)
    {
        x = fitVec.at(i).at(0); // 1st column data
        y = fitVec.at(i).at(1); // 2nd column data
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
        fit_xData += xdata.at(i);
        fit_yData += ydata.at(i);
    }
    //cout << "\n" << fit_xData << "\n";
    //cout << "\n" << fit_yData << "\n";
}





void AlglibFittingKernel::fitdata_processing_lsfit_2dx(const vector<vector<double>>& fitVec)
{
    /** NOTE:
     ** The final column stores the fitting target (ex. taueq) **/
    
    xdata.clear();
    ydata.clear();
    xraw.clear();
    yraw.clear();
    xorg.clear();
    yorg.clear();
    fit_xData.clear();
    fit_yData.clear();
    
    int    n_column=(int)fitVec.at(0).size();
    double y=0;
    string tmp;
    
    xdata.push_back("[[");
    ydata.push_back("[");
    for (size_t i=0; i<fitVec.size(); ++i)
    {
        // x_data to fit
        tmp.clear();
        for (indexi=0; indexi<(n_column-1); ++indexi) {
            if (indexi==0) {
                tmp += to_string((long double)fitVec.at(i).at(indexi));
            } else {
                tmp += ","+to_string((long double)fitVec.at(i).at(indexi));
            }
        }
        xdata.push_back(tmp);
        xdata.push_back("],[");
        // y_data to fit
        y = fitVec.at(i).at(n_column-1);
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
    for (size_t i=0; i<fitVec.size(); ++i)
    {
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





void AlglibFittingKernel::write_relaxation_target(std::ofstream& outputFile,
                                                  const std::string& reltarget)
{
    outputFile << "Relaxation target is <"<<reltarget<<">\n";
}





void AlglibFittingKernel::write_message(std::ofstream& outputFile,
                                        const std::string& msg)
{
    outputFile << msg << "\n\n";
}





void AlglibFittingKernel::write_fitModel(std::ofstream& outputFile,
                                         const std::string& model)
{
    outputFile << "Data fit to <"<<model<<"> functional form\n\n";
}





void AlglibFittingKernel::write_fitarrays(std::ofstream& outputFile)
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
    << "\n\n";
    
    outputFile
    << "yData (#.=" << (ydata.size()-1)/2 << ")"
    << "\n"
    << yData
    << "\n\n";
}





void AlglibFittingKernel::write_stopcond(std::ofstream& outputFile)
{
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
    outputFile << "\n"
    << "epsf   = " << epsf << "\n"
    << "epsx   = " << epsx << "\n"
    << "maxits = " << maxits << "\n\n";
}





void AlglibFittingKernel::write_fitinfo(std::ofstream& outputFile)
{
    if(get_is_use_FG()) outputFile << "Operating Mode: FG";
    else outputFile << "Operating Mode: F";
    
    outputFile << "\n\n"
    << "fit info " << int(info) << "\n"
    << "(" << get_fit_info(int(info)) << ") \n\n";
    
    outputFile
    << "internal iterations = " << rep.iterationscount << "\n"
    << "time used = "
    << difftime(usedt[1],usedt[0]) << " seconds \n\n";
    
    outputFile << "fit coeffs \n";
    for (size_t i=0; i<c.length(); ++i) {
        outputFile << c[i] << " ";
    } outputFile << "\n\n";
    
    outputFile << "fit coeffs error \n";
    for (size_t i=0; i<c.length(); ++i) {
        outputFile << rep.errpar[i] << " ";
    } outputFile << "\n\n";
    
    outputFile
    << "Statistics \n"
    << "=================================== \n"
    << "R^2 (coefficient of determination)  \n"
    << rep.r2 << " ("<<cyclecount<<" cycles) \n\n"
    << "RMS error \n"
    << rep.rmserror << "\n\n"
    << "average error \n"
    << rep.avgerror << "\n\n"
    << "average relative error \n"
    << rep.avgrelerror << "\n\n"
    << "maximum error \n"
    << rep.maxerror << "\n"
    << "=================================== \n\n";
}





void AlglibFittingKernel::write_yvsxby2dvecdub(std::ofstream& outputFile,
                                               const std::vector<std::vector<double>>& xy)
{
    outputFile << "\n"
    << "x  y                                         " << "\n"
    << "=============================================" << "\n";
    for (int i=0;i<(int)xy.size();++i) {
        outputFile
        << xy[i][0] << " "
        << xy[i][1] << "\n";
    } outputFile    << "\n";
}





void AlglibFittingKernel::write_errcurve(std::ofstream& outputFile,
                                         const std::string& model)
{
    outputFile << "\n"
    << "T  log10(tau_org) log10(tauFit_fit) ErrCurve " << "\n"
    << "=============================================" << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << taueq_sortdecreasing[i][0] << " "
        << taueq_sortdecreasing[i][1] << " "
        << log10(calc_y_given_x(c,taueq_sortdecreasing[i][0],model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n";
    
    outputFile << "\n"
    << 1000.0/(corrFacforT/precision)
    << "/T  log10(tau_org) log10(tauFit_fit) ErrCurve " << "\n"
    << "==============================================" << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << (1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0] << " "
        << taueq_sortdecreasing[i][1] << " "
        << log10(calc_y_given_x(c,taueq_sortdecreasing[i][0],model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n\n";
}





void AlglibFittingKernel::write_errcurve(std::ofstream& outputFile,
                                         const std::string& model,
                                         const std::vector<std::vector<double>>& vvd)
{
    outputFile << "\n"
    << "T  log10(tau_org) log10(tauFit_fit) ErrCurve " << "\n"
    << "=============================================" << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << vvd[i][0] << " "
        << vvd[i][1] << " "
        << log10(calc_y_given_x(c,vvd[i][0],model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n";
    
    outputFile << "\n"
    << 1000.0/(corrFacforT/precision)
    << "/T  log10(tau_org) log10(tauFit_fit) ErrCurve " << "\n"
    << "==============================================" << "\n";
    for (int i=0; i<rep.errcurve.length(); ++i) {
        outputFile
        << (1000.0/(corrFacforT/precision))/vvd[i][0] << " "
        << vvd[i][1] << " "
        << log10(calc_y_given_x(c,vvd[i][0],model)) << " "
        << rep.errcurve[i] << "\n";
    } outputFile << "\n\n";
}





void AlglibFittingKernel::write_tauFit_vs_T(std::ofstream& outputFile)
{
    /* format: T log10(tauFit_fit) */
    //------------------------------------------------------------------
    outputFile << "\n"
    << "T  log10(tauFit)                 " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<taueq_sortdecreasing.size(); ++i) {
        outputFile
        << taueq_sortdecreasing[i][0] << " "
        << taueq_sortdecreasing[i][1] << " "
        << sorted_avg[i][1]           << " "
        << sorted_stdev[i][1]         << "\n";
    } outputFile << "\n\n";
}





void AlglibFittingKernel::write_tauFit_vs_invT(std::ofstream& outputFile)
{
    /* format: 1/T log10(tauFit) */
    //------------------------------------------------------------------
    outputFile << "\n"
    << 1000.0/(corrFacforT/precision)
    << "/T  log10(tauFit)                " << "\n"
    << "=================================" << "\n";
    for (size_t i=0; i<taueq_sortdecreasing.size(); ++i) {
        outputFile
        << (1000.0/(corrFacforT/precision))/taueq_sortdecreasing[i][0] << " "
        << taueq_sortdecreasing[i][1] << " "
        << sorted_avg[i][1]           << " "
        << sorted_stdev[i][1]         << "\n";
    } outputFile << "\n\n";
}





void AlglibFittingKernel::write_badFit(std::ofstream& outputFile,
                                       const int code)
{
    if (code==0) {
        outputFile << "\n\n"
        << "=========================================================\n"
        << "NOTE:\n"
        << "r2 < threshold, so tauFit is not written into taueq file!\n"
        << "=========================================================\n";
    }
    else if (code==1) {
        outputFile << "\n\n"
        << "=================================================\n"
        << "NOTE:\n"
        << "Current relaxation function is not fully relaxed.\n"
        << "=================================================\n";
    } else {
        cout
        << "in AlglibFittingKernel::write_badFit(): code value not valid.\n";
        exit(EXIT_FAILURE);
    }
}





void AlglibFittingKernel::write_insufficientData(std::ofstream& outputFile,
                                                 const int n_fit,
                                                 const int n_threshold)
{
    outputFile << "\n\n"
    << "===================================\n"
    << "Insufficient data for curve-fitting:\n"
    << "Data for fitting = " << n_fit << " points\n"
    << "less than "<<n_threshold<<" points threshold \n"
    << "===================================\n";
}





void AlglibFittingKernel::build_cubicspline_interpolant(const vector<vector<double>>& raw,
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
    fitdata_processing_1d(data);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    spline1dbuildcubic(x_in,y_in,interpolant);
}





void AlglibFittingKernel::build_cubicspline_loglog_interpolant(const vector<vector<double>>& raw,
                                                               const int x_index,
                                                               const int y_index)
{
    vector<vector<double>> data;
    for (int i=0; i<(int)raw.size(); ++i) {
        data.push_back
        ({
            log10(fabs(raw.at(i).at(x_index))),
            log10(fabs(raw.at(i).at(y_index)))
        });
    }
    fitdata_processing_1d(data);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    spline1dbuildcubic(x_in,y_in,interpolant);
}





void AlglibFittingKernel::build_polynomial_interpolant(const vector<vector<double>>& raw,
                                                       const int x_index,
                                                       const int y_index,
                                                       const ae_int_t m)
{
    /** fit to mth order polynomial form **/
    if (m<1) {
        cout
        << "in AlglibFittingKernel::build_polynomial_interpolant():\n"
        << "please specify polynomail order > 0.\n"; exit(EXIT_FAILURE);
    }
    vector<vector<double>> data;
    for (int i=0; i<(int)raw.size(); ++i) {
        data.push_back
        ({
            raw.at(i).at(x_index),
            raw.at(i).at(y_index)
        });
    }
    fitdata_processing_1d(data);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    ae_int_t info;
    polynomialfitreport rep;
    //barycentricinterpolant p;
    polynomialfit(x_in,y_in,m+1,info,interpolantpolynom,rep);
}





void AlglibFittingKernel::build_penalizedspline_interpolant(const vector<vector<double>>& raw,
                                                            const int x_index,
                                                            const int y_index,
                                                            const ae_int_t m,
                                                            const double rho)
{
    vector<vector<double>> data;
    for (int i=0; i<(int)raw.size(); ++i) {
        data.push_back
        ({
            raw.at(i).at(x_index),
            raw.at(i).at(y_index)
        });
    }
    fitdata_processing_1d(data);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    ae_int_t info;
    spline1dfitreport rep;
    //ae_int_t n_basis_functions=50;
    //double rho=3.0;
    spline1dfitpenalized(x_in,y_in,m,rho,info,interpolant,rep);
}





void AlglibFittingKernel::build_penalizedspline_loglog_interpolant(const vector<vector<double>>& raw,
                                                                   const int x_index,
                                                                   const int y_index,
                                                                   const ae_int_t m,
                                                                   const double rho)
{
    vector<vector<double>> data;
    for (int i=0; i<(int)raw.size(); ++i) {
        data.push_back
        ({
            log10(fabs(raw.at(i).at(x_index))),
            log10(fabs(raw.at(i).at(y_index)))
        });
    }
    fitdata_processing_1d(data);
    real_1d_array x_in=fit_xData.c_str();
    real_1d_array y_in=fit_yData.c_str();
    ae_int_t info;
    spline1dfitreport rep;
    //ae_int_t n_basis_functions=50;
    //double rho=3.0;
    spline1dfitpenalized(x_in,y_in,m,rho,info,interpolant,rep);
}





double AlglibFittingKernel::tfunction_at_given_time(const vector<vector<double>>& tfunction,
                                                    const double time,
                                                    bool& flag)
{
    /** function returns value at a specified time
     ** NOTE: for <r2>, need to use log-log data **/
    
    double tfunc_val=0;
    size_t max=tfunction.size();
    double mintime=tfunction.at(0).at(0);
    double maxtime=tfunction.at(max-1).at(0);
    if (time<=maxtime&&time>=mintime) {
        build_penalizedspline_interpolant(tfunction,0,1,50,2.0);
        //build_cubicspline_interpolant(tfunction,0,1);
        tfunc_val=spline1dcalc(interpolant,time);
    } else {
        //cout
        //<< "in AlglibFittingKernel::tfunction_at_given_time():\n"
        //<< "Specfied time ("<<time<<") exceeds time boundary "
        //<< "("<<mintime<<","<<maxtime<<"). Please check.\n";
        flag=true;
        //exit(EXIT_FAILURE);
    } return tfunc_val;
}





double AlglibFittingKernel::tfunction_at_given_valu(const vector<vector<double>>& tfunction,
                                                    const double valu,
                                                    bool& flag)
{
    /** function returns time at a specified value
     ** NOTE: for <r2>, need to use log-log data **/
    
    //double tfunc_val=0;
    double time=0;
    size_t max=tfunction.size();
    double mintime=tfunction.at(0).at(0);
    double maxtime=tfunction.at(max-1).at(0);
    
    build_penalizedspline_interpolant(tfunction,1,0,50,2.0);//reverse domains
    time=spline1dcalc(interpolant,valu);
    
    if (time<=maxtime&&time>=mintime) {
        flag=false;
    } else {
        flag=true;
    } return time;
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
                                                vector<vector<double>>& sorted_avg_local)
{
    if (sorted_org.size()==0) {
        cout
        << "in AlglibFittingKernel::averaging_sorted_data():\n"
        << "input vector size = 0, please check."; exit(EXIT_FAILURE);
    }
    
    //--------------------------------------------------------------------------
    // NOTE:
    // this algorithm is only valid when the 1st column has been sorted, either
    // increasing or decreasing; if 1st column data is random, the algorithm
    // will not be able to give the right answer.
    //--------------------------------------------------------------------------
    /** Modified by SJH (20170609) **/
    sorted_avg.clear();
    sorted_stdev.clear();
    SysInfo stats;
    deque<double> avg,stdev;
    vector<vector<double>> tmpvvd;
    
    sorted_avg_local.clear();
    vector<double> sum_columns;//stores data of 2nd to final columns of a row
    vector<double> finaldata;
    
    
    int n_data=(int)sorted_org.size();
    int n_columns=(int)sorted_org.at(0).size();
    for (int i=1;i<sorted_org.size();++i) {
        if (sorted_org.at(i).size()>n_columns) {
            n_columns=(int)sorted_org.at(i).size();
        }
    }
    
    if (n_columns>1) //ex. {T,tau,...}
    {
        //initialize container
        for (int j=0; j<(n_columns-1); ++j)
        {
            sum_columns.push_back(0);
            tmpvvd.push_back({0});
        }
        int count=0,i=0,i_beg=0;
        do
        {
            i=i_beg;
            //put current row values (2nd to final column) to container
            for (int j=0; j<(n_columns-1); ++j)
            {
                sum_columns.at(j)=sorted_org.at(i).at(j+1);//NOTE: j+1 since from 2nd column
                tmpvvd.at(j)={sorted_org.at(i).at(j+1)};
                avg.push_back(0);
                stdev.push_back(0);
            }
            count=0;
            for (int ii=i+1; ii<n_data; ++ii)
            {
                //if 1st column has the same value (ex. same T), add following column values
                //to current container at current row; count number of repeats.
                if (sorted_org.at(i).at(0)==sorted_org.at(ii).at(0)) {
                    ++count;
                    for (int j=0; j<(n_columns-1); ++j) {
                        sum_columns.at(j)+=sorted_org.at(ii).at(j+1);
                        tmpvvd.at(j).push_back(sorted_org.at(ii).at(j+1));
                    } i_beg=ii+1;
                }
            }
            if (count==0) ++i_beg;
            
            //if first colimn value has repeated, do arithmetic average in
            //corresponding columns by repeat times
            for (int j=0; j<(n_columns-1); ++j) {
                sum_columns.at(j)/=(double)(count+1);
                stats.set_calcvector(tmpvvd.at(j));
                avg.at(j)=stats.calc_mean();
                stdev.at(j)=stats.get_sample_stdev();
            }
            finaldata.clear();
            finaldata.push_back(sorted_org.at(i).at(0));//1st column value
            
            //push back (averaged) values from 2nd to the final column
            for (int j=0; j<(n_columns-1); ++j) {
                finaldata.push_back(sum_columns.at(j));
            } sorted_avg_local.push_back(finaldata);
            
            /** Modified by SJH (20170609) **/
            avg.push_front(sorted_org.at(i).at(0));
            stdev.push_front(sorted_org.at(i).at(0));
            sorted_avg.push_back(avg);
            sorted_stdev.push_back(stdev);
            avg.clear();
            stdev.clear();
            
        } while(i_beg<n_data);
    }
    else //ex. {T}
    {
        int count=0,i=0,i_beg=0;
        do {
            i=i_beg;
            count=0;
            for (int ii=i+1; ii<n_data; ++ii) {
                if (sorted_org.at(i).at(0)==sorted_org.at(ii).at(0)) {
                    ++count;
                    i_beg=ii+1;
                }
            }
            if (count==0) ++i_beg;
            finaldata.clear();
            finaldata.push_back(sorted_org.at(i).at(0));
            sorted_avg_local.push_back(finaldata);
        } while(i_beg<n_data);
    }
}





void AlglibFittingKernel::cout_c_parms()
{
    for (indexi=0; indexi<c.length(); ++indexi) {
        cout << c[indexi] << " ";
    } cout << "\n";
}





void AlglibFittingKernel::cout_r2()
{
    cout << rep.r2 << "\n";
}





/** setters **/
void AlglibFittingKernel::set_is_use_FG(const bool b){is_use_FG=b;}




