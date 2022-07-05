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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "userinterface.h"
#include "sysinfo.h"
#include "structureclass.h"
#include "lmpscripts.h"
#include "workscripts.h"
#include "amdatanalysis.h"
#include "fitdata.h"

#include "moltemplatelmpdata.h"
#include "fitdielectric.h"
#include "theorytest.h"

namespace autoWork
{
    const std::string submit_jobs(const UserInterface * const inpVar=nullptr);
    //--------------------------------------------------------------------------
    /** submit_jobs() serves as an interfacee that interacts with external code
     ** and handles automation, including but not limited to file structure
     ** genreration, simulation automation, amdat calculation automation, and
     ** automated data analysis; this function returns the path to a final
     ** report file when automation is completed **/
    //--------------------------------------------------------------------------
    
    void automated_batch(StructureClass&,const LmpScripts&,bool&,const UserInterface * const inpVar=nullptr);
    //--------------------------------------------------------------------------
    /** automated_batch() is a subroutine that dwells in submit_jobs() that
     ** handles batch automation of simulation and analysis; since the
     ** automation design is based on a scheme of iterative relaxation regimes,
     ** all regime-specific information will be automatically and properly
     ** handled in this subroutine, such as temperatures, simulation time,
     ** number of temperatures, and so forth; automated data analysis
     ** is also handled at the end of this subroutine to continue the
     ** bootstrapping process until the final regime has finished **/
    //--------------------------------------------------------------------------
    
    /** Watch for targets and hold process **/
    void watch_hold(const StructureClass&,const std::vector<std::string>&,const std::string&,const int);
    
    /** Initialization of system variables **/
    void system_initialization(StructureClass&,const UserInterface * const inpVar=nullptr);
    void moltemplate_initialization(StructureClass&,MoltemplateLmpData&,const UserInterface * const inpVar=nullptr);
    void fitDielectric_initialization(StructureClass&,FitDielectric&,const UserInterface * const inpVar=nullptr);
    void simulation_initialization(StructureClass&,WorkScripts&,const UserInterface * const inpVar=nullptr);
    void analysis_initialization(StructureClass&,AmdatAnalysis&,const UserInterface * const inpVar=nullptr);
    void fd_initialization(StructureClass&,FitData&,const UserInterface * const inpVar=nullptr);
    void tt_initialization(StructureClass&,TheoryTest&,const UserInterface * const inpVar=nullptr);
    
    /** Write out progress info to file **/
    void write_progress(StructureClass&);
    
    /** Write out simulation results to file **/
    void write_report(StructureClass&);
    
    /** System command to interact with Linux Bash **/
    void call_system_bash(const std::string&);
    
    /** Get linearly spaced temperatures **/
    std::vector<std::vector<double>> linearTemperatures(const SysInfo&);
    
    /** Invoke independent ALGLIB fit routine **/
    //---------------------------------------------
    void alglibfit_routine(const StructureClass&);
    //---------------------------------------------
    
    /** Convenience functions  **/
    void use_prepScripts(const StructureClass&);
    void error_checking(const StructureClass&);
    void clean_tmpfiles(const StructureClass&);
    void delete_existing_file(const std::string&);
    void system_wait(const int waitsec=3);
    void system_countdown(const int waitsec=3);
    void makeNewAnanlysisFolder(const StructureClass&);
    void backupfitdata(const StructureClass&);
    void backupstatistics(const StructureClass&);
    void copy_from_to(const StructureClass&);
    void coreAllocation(const StructureClass&,const int,int&);
    void makeltfilefrompdb();
    void makeTacltfile();
    void check_finished_blocks(const StructureClass&,const WorkScripts&,const std::vector<std::vector<double>>&);
    void cutoff_trajectory(const StructureClass&,const WorkScripts&,const std::vector<std::vector<double>>&,const int);
    void make_SimulationsFolders(const StructureClass&);
    void move_InnerLoopFiles(const StructureClass&,const int,const int);
    void move_OuterLoopFiles(const StructureClass&,const int);
    void backup_simulation_traj(const StructureClass&,const std::string&);
    void backup_analysis_stats(const StructureClass&);
    bool is_path_exist(const std::string&);
    bool randbool(const double p=0.5);
    int  randint(const int, const int);
    double error_propagation(const std::vector<double>&,const std::vector<double>&);
    const std::string return_year();
    const std::string return_date();
    const std::string return_datetime();
    const std::string return_SimulationFolderPath(const SysInfo&);
    const std::string return_AnalysisFolderPath(const SysInfo&);
    const std::string extract_amdat_svn(const std::string&);
    
    template <typename T>
    void transpose(std::vector<std::vector<T>>& vec) {
        std::vector<std::vector<T>>
        ret(vec[0].size(),std::vector<T>(vec.size()));
        for (int i=0;i<vec[0].size();++i) {
            for (int j=0;j<vec.size();++j) {
                ret[i][j]=vec[j][i];
            }
        } vec=ret;
    }
}
#endif /* FUNCTIONS_H */




