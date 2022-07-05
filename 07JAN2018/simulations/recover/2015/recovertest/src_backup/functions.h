//
//  Created by Sean Jh H. on 7/9/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "sysinfo.h"
#include "structureclass.h"
#include "lmpscripts.h"
#include "amdatanalysis.h"

namespace autoWork {
    
    //=============================================================
    // The interface function interacts with external master code
    // and handles automation;
    // This function returns path to the final report file when
    // all automation is completed.
    //=============================================================
    const std::string submit_jobs();
    
    //=============================================================
    // This is a subroutine dwelling in submit_jobs() that handles
    // automated simulations and analysis,
    // including jobs submission (simulation and analysis),
    // and interal calculation (ALGLIB curve-fitting)
    //
    // NOTE:
    // This function contains "regime" information which will be
    // changed according to previously set up in the submit_job()
    // function, ex. # of T's in each regime, ts in each regime;
    // also, temperatures for different trials and systems in
    // different regimes will be automatically updated as regime
    // develops. All these are handled in this function.
    //=============================================================
    void automated_batch(StructureClass&,const LmpScripts&,bool&);
    
    //=============================================================
    // Watch list functionality
    //=============================================================
    void watch_hold(const StructureClass&,
                    const std::vector<std::string>&,
                    const std::string&,
                    const int);
    
    //=============================================================
    // Write out progress info to file
    //=============================================================
    void write_progress(StructureClass&);
    
    //=============================================================
    // Write out simulation results to file
    //=============================================================
    void write_report(StructureClass&);
    
    //=============================================================
    // System command to interact with bash shell
    //=============================================================
    void call_system_bash(std::string&);
    
    //=============================================================
    // Get linearly spaced temperatures defined by user
    //=============================================================
    std::vector<std::vector<double>> linearTemperatures(const SysInfo&);
    
    
    //=============================================================
    // Convenience functions used throughout the program
    //=============================================================
    const std::string return_year();
    const std::string return_date(); // format: M/D/Y
    const std::string return_datetime();
    const std::string make_SimulationsFolders(const StructureClass&);
    const std::string test_AnalysisFolder(const StructureClass&,
                                          const AmdatAnalysis&);
    const std::string move_InnerLoopFiles(const StructureClass&,
                                          const int,const int);
    const std::string move_OuterLoopFiles(const StructureClass&,
                                          const int);
    const std::string make_InitialFolders(const SysInfo&);
    const std::string archive_InitialFiles(const StructureClass&);
    const std::string return_SimulationFolderPath(const SysInfo&);
    const std::string return_AnalysisFolderPath(const SysInfo&);
    const std::string echo_NodesUsage(const std::vector<int>&,
                                      const std::vector<int>&);
    void system_initialization(StructureClass&);
    void error_checking(const StructureClass&);
    void clean_tmpfiles(const StructureClass&);
    void delete_existing_file(const std::string& fileName);
    double error_propagation(const std::vector<double>& dFi,
                             const std::vector<double>& errori);
    
    std::vector<double> thermocalc_equ(const StructureClass&,
                                       const int n_sys,
                                       const double Temp_d,
                                       const std::string& keyWord);
    
    std::vector<double> thermocalc_equ_avgT(const StructureClass&,
                                            const int n_sys,
                                            const double Temp_d);
    
    std::vector<double> thermocalc_equ_avgL(const StructureClass&,
                                            const int n_sys,
                                            const double Temp_d);
    
    void thermodataprocess_prd(const StructureClass&,
                               const int n_trl,
                               const int n_sys,
                               const double Temp_d);
    
    void thermocalc_prd(const StructureClass&,
                        const int n_trl,
                        const int n_sys,
                        const double Temp_d);
    
    void test_RepeatednessByTwoBases(const int,const int);
}
#endif
