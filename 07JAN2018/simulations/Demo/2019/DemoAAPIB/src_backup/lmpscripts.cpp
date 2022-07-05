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

#include "lmpscripts.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

LmpScripts::LmpScripts(const StructureClass& sysVar):

/* bool */
is_tampering(false),     // true: use parallel tampering after equil.
is_respa(false),         // true: RESPA; false: Verlet
is_resize(false),        // true: resize and fix box size after equil. run
is_fixMomentum(true),    // true: fixp in gen. thru equil. run
is_fixPrdMomentum(true), // true: fixp in production run
is_fixShake(false),
is_bondswap(false),

/* int */
fixpsteps(0),
fixpsteps_prd(0),
delay(0),
shakeIter(20),
print_shake(5000),

/* double */
shakeTol(0.0001),
Tdamp(1000.0),
Pdamp(1000.0),
pressure(1.0),

/* string */
GPU_mode("force/neigh")
{
    if (sysVar.get_systemUnit()=="real") {
        Tdamp         = 100.0;
        Pdamp         = 1000.0;
        pressure      = 1.0;
        fixpsteps     = 1000; // fix steps in equilibration
        fixpsteps_prd = 1000; // fix steps in production
    } else if (sysVar.get_systemUnit()=="lj") {
        Tdamp         = 2.0;
        Pdamp         = 2.0;
        pressure      = 0.0;
        fixpsteps     = 200;
        fixpsteps_prd = 2000;
    } else if (sysVar.get_systemUnit()=="metal") {
        Tdamp         = 0.1;
        Pdamp         = 1.0;
        pressure      = 1.0;
        fixpsteps     = 1000;
        fixpsteps_prd = 1000;
    }
    /* get the is_fixShake information from the StructureClass (VERY OLD) */
    int counter=0;
    for (size_t i=0; i<sysVar.get_is_shakeBonds().size(); ++i) {
        if (sysVar.get_is_shakeBonds()[i]) ++counter;
    }
    for (size_t i=0; i<sysVar.get_is_shakeAngles().size(); ++i) {
        if (sysVar.get_is_shakeAngles()[i]) ++counter;
    }
    if (counter>0) {
        set_is_fixShake(true);
    } else {
        set_is_fixShake(false);
    }
}





void LmpScripts::make_GenerationInputFile(const StructureClass& sysVar) const
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/generation");
    o.append("/generation.inp");
    //cout << o << "\n";
    
    ofstream gen(o.c_str());
    if (gen.is_open()) {
        
        //======================================================================
        // Header
        //======================================================================
        gen
        << "# LAMMPS Script for Generation"
        << "\n\n";
        
        if (sysVar.get_is_GPU()) {
            gen
            << "package gpu "
            << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
        }
        
        gen << "\n";
        
        if (sysVar.get_systemUnit()=="lj") {
            gen << "units lj" << "\n";
        }
        else if (sysVar.get_systemUnit()=="real") {
            gen << "units real" << "\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            gen << "units metal" << "\n";
        }
        
        if (sysVar.get_angleType()==0) {
            gen << "atom_style bond" << "\n";
        }
        else {
            gen << "atom_style angle" << "\n";
        }
        
        gen << "\n";
        
        
        //======================================================================
        // READ_DATA
        //======================================================================
        gen << "# ----------------------------------"  << "\n";
        gen << "## READ_DATA"                          << "\n";
        gen << "# ----------------------------------"  << "\n";
        
        if (1) {
            gen
            << "read_data ./lammps_inputs/start_data/start"
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << ".data"
            << "\n\n";
        }
        else {
            gen
            << "read_data ./lammps_inputs/start_data/input"
            << "_${trial}"
            << ".data"
            << "\n\n";
        }
        
        
        //======================================================================
        // Bond_style and Bond Coeffs
        //======================================================================
        gen << "# ----------------------------------"  << "\n";
        gen << "## BOND"                               << "\n";
        gen << "# ----------------------------------"  << "\n";
        ifstream bP("bonds_Params.txt");
        if (bP.is_open()) {
            string lineContent;
            string stringData;
            //double doubleData;
            //int intData;
            
            getline(bP,lineContent);
            
            gen << lineContent << "\n"; // bond_style harmonic
            
            while (getline(bP,lineContent)) {
                
                gen << lineContent << "\n";
                
                /* < old method >
                 
                 // bond_coeff 1 956.02 2.50
                 istringstream lineStream(lineContent);
                 lineStream >> stringData;
                 gen << stringData;
                 lineStream >> intData;
                 gen << " " << intData;
                 lineStream >> doubleData;
                 gen << " " << fixed << setprecision(2) << doubleData;
                 lineStream >> doubleData;
                 gen << " " << doubleData << "\n";
                 
                 */
            }
            bP.close();
        }
        
        gen << "\n";
        
        
        //======================================================================
        // Angle_style and Angle Coeffs
        //======================================================================
        if (sysVar.get_trueAngles()!=0) {
            
            gen << "# ----------------------------------" << "\n";
            gen << "## ANGLE" << "\n";
            gen << "# ----------------------------------" << "\n";
            ifstream aP("angles_Params.txt");
            if (aP.is_open()) {
                string lineContent;
                string stringData;
                //double doubleData;
                //int intData;
                
                getline(aP,lineContent);
                
                gen << lineContent << "\n"; // ex. angle_style cosine/squared
                
                while (getline(aP,lineContent)) {
                    
                    
                    gen << lineContent << "\n";
                    
                    /* < old method >
                     
                     istringstream lineStream(lineContent);
                     // angle_coeff 1 2.99 170.0
                     lineStream >> stringData;
                     gen << stringData;
                     lineStream >> intData;
                     gen << " " << intData;
                     lineStream >> doubleData;
                     gen << " " << fixed << setprecision(2) << doubleData;
                     lineStream >> doubleData;
                     gen << " " << fixed << setprecision(1) << doubleData << "\n";
                     
                     */
                }
                aP.close();
            }
            gen << "\n";
        }
        
        
        //======================================================================
        // Pair_style and Pair Coeffs
        //======================================================================
        gen << "# ----------------------------------" << "\n";
        gen << "## PAIR INTERACTION" << "\n";
        gen << "# ----------------------------------" << "\n";
        ifstream lj("pairs_Params.txt");
        if (lj.is_open()) {
            
            string lineContent;
            string stringData;
            //double doubleData;
            //int intData;
            
            getline(lj,lineContent);
            
            gen << lineContent << "\n"; // ex. pair_style lj/gromacs 9.0 12.0
            
            /* < old method >
             
             istringstream liness(lineContent);
             // pair_style lj/gromacs 9.0_Ang. 12.0_Ang.
             liness >> stringData;
             gen << stringData;
             liness >> stringData;
             gen << " " << stringData;
             liness >> doubleData;
             gen << " " << fixed << setprecision(1) << doubleData;
             liness >> doubleData;
             gen << " " << doubleData << "\n";
             
             */
            
            while (getline(lj,lineContent)) {
                
                gen << lineContent << "\n";
                
                /* < old method >
                 
                 istringstream lineStream(lineContent);
                 // pair_coeff 1 1 0.6274 4.30 9.0 12.0
                 lineStream >> stringData;
                 gen << stringData;
                 lineStream >> intData;
                 gen << " " << intData;
                 lineStream >> intData;
                 gen << " " << intData;
                 lineStream >> doubleData;
                 gen << " " << fixed << setprecision(4) << doubleData;
                 lineStream >> doubleData;
                 gen << " " << fixed << setprecision(2) << doubleData;
                 lineStream >> doubleData;
                 gen << " " << fixed << setprecision(1) << doubleData;
                 lineStream >> doubleData;
                 gen << " " << doubleData << "\n";
                 
                 */
            }
            lj.close();
        }
        
        gen << "\n";
        
        
        //======================================================================
        // SPECIAL_BONDS
        //======================================================================
        gen
        << "# ----------------------------------" << "\n"
        << "## PAIR"                              << "\n"
        << "# ----------------------------------" << "\n";
        
        if (sysVar.get_trueAngles()!=0) {
            
            // With angle forcefield
            gen
            << "special_bonds lj 0.0 0.0 1.0"
            << "\n";
        }
        else {
            
            // No angle forcefield
            gen
            << "special_bonds lj 0.0 1.0 1.0"
            << "\n";
        }
        
        
        //======================================================================
        // PAIR_MODIFY
        //======================================================================
        gen
        << "pair_modify shift yes"
        << "\n";
        
        
        //======================================================================
        // NEIGH_MODIFY
        //======================================================================
        gen
        << "neigh_modify delay 1 every 1 check yes"
        << "\n\n";
        
        //======================================================================
        // NOTE: Group ID and Types of Atoms
        //======================================================================
        
        /*
         For current version, this is dependent on methods
         in the StructureClass, i.e. if they were not executed,
         there would be no information shown in the Group IDs */
        
        gen << "# ----------------------------------" << "\n";
        gen << "## Group IDs (Note the Atom types)"   << "\n";
        gen << "# ----------------------------------" << "\n";
        
        int overalltypes=0,backbonetype=0,sidegrouptype=0;
        for (int i=0; i<sysVar.get_trueAtoms(); ++i)
            if(sysVar.get_is_backboneType()[i]) {++overalltypes; ++backbonetype;}
        for (int i=0; i<sysVar.get_trueAtoms();++i)
            if(sysVar.get_is_sidegroupType()[i]) {++overalltypes; ++sidegrouptype;}
        
        if(overalltypes>0) {
            
            gen << "group polymer type";
            for (int i=0; i<sysVar.get_trueAtoms(); ++i) gen << " " << i+1;
            gen << "\n";
            
            if (backbonetype>0) {
                gen << "group backbone type";
                for (int i=0; i<sysVar.get_trueAtoms(); ++i)
                    if(sysVar.get_is_backboneType()[i]) gen << " " << i+1;
                gen << "\n";
            }
            
            if (sidegrouptype>0) {
                gen << "group sidegroup type";
                for (int i=0; i<sysVar.get_trueAtoms(); ++i)
                    if(sysVar.get_is_sidegroupType()[i]) gen << " " << i+1;
                gen << "\n";
            }
        }
        
        gen << "\n";
        
        
        //======================================================================
        // DUMP
        //======================================================================
        gen
        << "# ----------------------------------"       << "\n"
        << "## DUMP"                                    << "\n"
        << "# ----------------------------------"       << "\n"
        << "dump 1 all xyz 1000 "
        << "./generation/trajectory/trajectory"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".gen.xyz"
        << "\n"
        << "dump_modify 1 append yes"
        << "\n\n";
        
        
        //======================================================================
        // FIX
        //======================================================================
        int fix=1;
        gen << "# ----------------------------------"   << "\n";
        gen << "## FIX"                                 << "\n";
        gen << "# ----------------------------------"   << "\n";
        if (get_is_fixMomentum()) {
            gen
            << "fix " << fix << " "
            << "all momentum " << get_fixpsteps() << " "
            << "linear 1 1 1 angular"
            << "\n";
            ++fix;
        }
        if (get_is_bondswap()) {
            gen
            << "fix " << fix << " "
            << "backbone bond/swap 50 0.5 4.5 ${vseed}"
            << "\n";
            ++fix;
        }
        //======================================================================
        if (get_is_fixShake()) {
            
            int overallshake=0,bondshake=0,angleshake=0;
            for (int i=0; i<sysVar.get_trueBonds(); ++i)
                if(sysVar.get_is_shakeBonds()[i]) {++overallshake; ++bondshake;}
            for (int i=0; i<sysVar.get_trueAngles();++i)
                if(sysVar.get_is_shakeAngles()[i]) {++overallshake; ++angleshake;}
            
            if(overallshake>0) {
                gen
                << "fix " << fix << " "
                << "polymer shake " << fixed << setprecision(4)
                << shakeTol  << " "
                << shakeIter << " "
                << "${ts_save}";
                
                if (bondshake>0) {
                    gen << " b";
                    for (int i=0; i<sysVar.get_trueBonds(); ++i)
                        if(sysVar.get_is_shakeBonds()[i])
                            gen << " " << i+1;
                }
                if (angleshake>0) {
                    gen << " a";
                    for (int i=0; i<sysVar.get_trueAngles(); ++i)
                        if(sysVar.get_is_shakeAngles()[i])
                            gen << " " << i+1;
                }
                gen << "\n";
                ++fix;
                
                gen.unsetf(ios_base::floatfield);
            }
        }
        //======================================================================
        gen << "\n";
        
        
        //======================================================================
        // START WITH HIGH TEMPERATURE
        //======================================================================
        gen
        << "# ----------------------------------"     << "\n"
        << "## START WITH HIGH TEMPERATURE"           << "\n"
        << "# ----------------------------------"     << "\n";
        gen
        << fixed << setprecision(1)
        << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
        
        gen.unsetf(ios_base::floatfield);
        
        gen
        << "fix " << fix << " "
        << "all nvt "
        << "temp "
        << "${T} ${T} "
        << get_Tdamp()
        << "\n\n";
        
        
        //======================================================================
        // THERMO_STYLE
        //======================================================================
        gen << "# ----------------------------------" << "\n";
        gen << "## THERMO_STYLE"                      << "\n";
        gen << "# ----------------------------------" << "\n";
        
        if (sysVar.get_is_GPU()) {
            gen
            << "thermo_style custom "
            << "step temp press pe ke etotal ebond eangle epair "
            << "lx ly lz vol density "
            << "dt time cpu tpcpu spcpu cpuremain"
            << "\n\n";
        }
        else {
            gen
            << "thermo_style custom "
            << "step temp press pe ke etotal lx ly lz vol"
            << "\n\n";
        }
        
        
        //======================================================================
        // RUN_STYLE
        //======================================================================
        gen << "# ----------------------------------" << "\n";
        gen << "## RUN_STYLE"                         << "\n";
        gen << "# ----------------------------------" << "\n";
        if (get_is_respa()) {
            
            if (sysVar.get_sideGroupInWords()=="stdFENE") {
                gen
                << "run_style respa 2 3 bond 1 pair 2"
                << "\n";
            }
            else {
                gen
                << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                << "\n";
            }
        }
        else {
            gen
            << "run_style verlet"
            << "\n";
        }
        
        gen << "\n";
        
        
        //======================================================================
        // VELOCITY (different seeds for different trials)
        //======================================================================
        if (1) {
            gen << "# ----------------------------------" << "\n";
            gen << "## VELOCITY"                          << "\n";
            gen << "# ----------------------------------" << "\n";
            
            gen
            << "velocity all create ${T} "
            << "${vseed} mom yes rot yes dist gaussian"
            << "\n";
            
            gen << "\n";
        }
        
        
        //======================================================================
        // RUN
        //======================================================================
        if (sysVar.get_systemUnit()=="real") {
            
            gen
            << "# ----------------------------------"      << "\n"
            << "## PRE-RUN"                                << "\n"
            << "# ----------------------------------"      << "\n"
            << "timestep 0.5"                              << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n"
            
            << "timestep 1.0"                              << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
            
            gen
            << "unfix " << fix << " "
            << "\n";
            
            gen
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            << get_pressure() << " " << get_pressure() << " "
            << get_Pdamp()
            << "\n\n";
            
            gen
            << "timestep ${ts}"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            
            gen
            << "# ----------------------------------"      << "\n"
            << "## PRE-RUN"                                << "\n"
            << "# ----------------------------------"      << "\n"
            << "timestep 0.0005"                           << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n"
            
            << "timestep 0.001"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
            
            gen
            << "unfix " << fix << " "
            << "\n";
            
            gen
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            << get_pressure() << " " << get_pressure() << " "
            << get_Pdamp()
            << "\n\n";
            
            gen
            << "timestep ${ts}"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
        }
        else if (sysVar.get_systemUnit()=="lj") {
            
            gen
            << "# ----------------------------------"      << "\n"
            << "## PRE-RUN"                                << "\n"
            << "# ----------------------------------"      << "\n"
            << "timestep 0.0005"                           << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n"
            
            << "timestep 0.001"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
            
            gen
            << "unfix " << fix << " "
            << "\n";
            
            // use high pressure of 1 here
            gen
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            <<"1.0 1.0 "
            << get_Pdamp()
            << "\n\n";
            
            gen
            << "timestep ${ts}"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
            
            // ramp pressure from 1 to final
            gen
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            << "1.0 " << get_pressure() << " "
            << get_Pdamp()
            << "\n\n";
            
            gen
            << "timestep ${ts}"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
            
            // run fianl pressure for a short period
            gen
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            << get_pressure() << " " << get_pressure() << " "
            << get_Pdamp()
            << "\n\n";
            
            gen
            << "timestep ${ts}"                            << "\n"
            << "thermo 1000"                               << "\n"
            << "run 10000"                                 << "\n\n";
        }
        
        
        //======================================================================
        // Saving Frequency
        //======================================================================
        gen
        << "# ----------------------------------"            << "\n"
        << "## Saving Frequency"                             << "\n"
        << "# ----------------------------------"            << "\n"
        << "variable ts_save equal floor(${steps_gen}*0.02)" << "\n"
        << "dump_modify 1 every ${ts_save}"
        << "\n\n";
        
        gen
        << "reset_timestep 0"
        << "\n\n";
        
        gen
        << "# ----------------------------------"      << "\n"
        << "## GENERATION RUN"                         << "\n"
        << "# ----------------------------------"      << "\n"
        << "timestep ${ts}"                            << "\n"
        << "thermo ${ts_save}"                         << "\n"
        << "run ${steps_gen}"                          << "\n"
        << "# ----------------------------------"      << "\n";
        //======================================================================
        
        gen << "\n";
        
        //======================================================================
        // GENERATION RESTART
        //======================================================================
        gen
        << "# ----------------------------------" << "\n"
        << "## GENERATION RESTART"                << "\n"
        << "# ----------------------------------" << "\n"
        << "write_restart "
        << "./generation/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".gen.restart"
        << "\n";
        //======================================================================
        gen << "\n";
        gen.close();
        cout << "lmpscripts(gen) made!\n\n";system_wait(1);
    }
    else cout << "Generation: 'lmp_gen.inp' cannot open." << "\n";
}





void LmpScripts::make_QuenchInputFile(const StructureClass& sysVar) const
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/quench");
    o.append("/quench.inp");
    //cout << o << "\n";
    
    ofstream qch(o.c_str());
    if (qch.is_open()) {
        
        //======================================================================
        // Header
        //======================================================================
        qch
        << "# LAMMPS Script for Quench"
        << "\n\n";
        
        if (sysVar.get_is_GPU()) {
            qch
            << "package gpu "
            << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
        }
        
        qch << "\n";
        
        //======================================================================
        // READ_RESTART
        //======================================================================
        qch
        << "# ----------------------------------" << "\n"
        << "## READ_RESTART"                      << "\n"
        << "# ----------------------------------" << "\n"
        << "read_restart ${readRestartPath}"
        << "\n";
        //======================================================================
        qch << "\n";
        
        qch << "reset_timestep 0" << "\n";
        
        qch << "\n";
        //======================================================================
        // NEIGH_MODIFY
        //======================================================================
        qch
        << "# ----------------------------------" << "\n"
        << "## PAIR"                              << "\n"
        << "# ----------------------------------" << "\n"
        << "neigh_modify "
        << "delay "<< get_delay() << " "
        << "every 1 check yes"
        << "\n\n";
        
        
        //======================================================================
        // FIX
        //======================================================================
        int fix=1;
        qch << "# ----------------------------------" << "\n";
        qch << "## FIX"                               << "\n";
        qch << "# ----------------------------------" << "\n";
        if (is_fixMomentum) {
            qch
            << "fix " << fix << " "
            << "all momentum " << get_fixpsteps() << " "
            << "linear 1 1 1 angular"
            << "\n";
            ++fix;
        }
        
        qch << "\n";
        
        //======================================================================
        // THERMO_STYLE
        //======================================================================
        qch << "# ----------------------------------" << "\n";
        qch << "## THERMO_STYLE"                      << "\n";
        qch << "# ----------------------------------" << "\n";
        
        if (sysVar.get_is_GPU()) {
            qch
            << "thermo_style custom "
            << "step temp press pe ke etotal ebond eangle epair "
            << "lx ly lz vol density "
            << "dt time cpu tpcpu spcpu cpuremain"
            << "\n\n";
        }
        else {
            qch
            << "thermo_style custom "
            << "step temp press pe ke etotal lx ly lz vol"
            << "\n\n";
        }
        
        
        //======================================================================
        // RUN_STYLE
        //======================================================================
        qch
        << "# ----------------------------------"       << "\n"
        << "## RUN_STYLE"                               << "\n"
        << "# ----------------------------------"       << "\n";
        if (get_is_respa()) {
            
            if (sysVar.get_sideGroupInWords()=="stdFENE") {
                qch
                << "run_style respa 2 3 bond 1 pair 2"
                << "\n";
            }
            else {
                qch
                << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                << "\n";
            }
        }
        else {
            
            qch
            << "run_style verlet" << "\n";
        }
        
        qch << "\n";
        
        //======================================================================
        // quench_steps & ts_save
        //======================================================================
        qch
        << "variable quench_steps equal "
        << "floor(${steps_qch}/${n_intervals})"
        << "\n"
        << "variable ts_save equal floor(${quench_steps}*0.2)"
        << "\n\n";
        
        
        //======================================================================
        // DUMP
        //======================================================================
        qch
        << "# ----------------------------------" << "\n"
        << "## DUMP"                              << "\n"
        << "# ----------------------------------" << "\n"
        
        << "dump 1 "
        << "all xyz ${ts_save} "
        << "./quench/trajectory/trajectory"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_Regime${Regime}"
        << ".qch.xyz"
        << "\n"
        
        << "dump_modify 1 append yes"
        << "\n\n";
        
        
        //======================================================================
        // save highest T of quench
        //======================================================================
        qch
        << "# ----------------------------------" << "\n"
        << "## save highest T of quench"          << "\n"
        << "# ----------------------------------" << "\n"
        << "write_restart "
        << "./quench/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${TempA}"
        << ".qch.restart"
        << "\n\n";
        
        
        //======================================================================
        // Looping Over Quenching Temperatures
        //======================================================================
        qch
        << "label loopOverT"
        << "\n\n";
        
        qch
        << fixed << setprecision(1)
        << "variable TA equal ${TempA}/" << sysVar.get_corrFacforT() << "\n"
        << "variable TB equal ${TempB}/" << sysVar.get_corrFacforT() << "\n\n";
        
        qch.unsetf(ios_base::floatfield); // clear decision setting
        
        qch
        << "fix " << fix << " "
        << "all npt "
        << "temp "
        << "${TA} ${TB} "
        << get_Tdamp() << " "
        << "iso "
        << get_pressure() << " " << get_pressure() << " "
        << get_Pdamp()
        << "\n\n";
        
        qch
        << "# ----------------------------------" << "\n"
        << "## QUENCH RUN"                        << "\n"
        << "# ----------------------------------" << "\n"
        << "timestep ${ts}"                       << "\n"
        << "thermo ${ts_save}"                    << "\n"
        << "run ${quench_steps}"                  << "\n"
        << "# ----------------------------------" << "\n";
        
        qch << "\n";
        
        qch
        << "# ----------------------------------" << "\n"
        << "## QUENCH RESTARTS"                   << "\n"
        << "# ----------------------------------" << "\n"
        << "write_restart "
        << "./quench/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${TempB}"
        << ".qch.restart"
        << "\n";
        
        qch << "\n";
        
        qch
        << "next TempA" << "\n"
        << "next TempB" << "\n"
        << "\n"
        << "jump SELF loopOverT";
        //======================================================================
        qch << "\n";
        qch.close();
        cout << "lmpscripts(qch) made!\n\n";system_wait(1);
    }
    else cout << "Quench: 'lmp_qch.inp' cannot open." << "\n";
}





void LmpScripts::make_EquilibrationInputFile(const StructureClass& sysVar) const
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/equilibration");
    o.append("/equilibration.inp");
    //cout << o << "\n";
    
    ofstream equ(o.c_str());
    if (equ.is_open()) {
        
        //======================================================================
        // Header
        //======================================================================
        equ
        << "# LAMMPS Script for Equilibration"
        << "\n\n";
        
        if (sysVar.get_is_GPU()) {
            equ
            << "package gpu "
            << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
        }
        
        equ << "\n";
        //======================================================================
        // READ_RESTART
        // NOTE: Read "Which" restart file
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## READ_RESTART"                      << "\n"
        << "# ----------------------------------" << "\n";
        /*
         equ
         << "read_restart "
         << "./quench/restart/restart"
         << "_${usic}"
         << "_${trial}"
         << "_${namestring}"
         << "_T${Temp}"
         << ".qch.restart"
         << "\n";
         */
        equ << "read_restart ${read_res}" << "\n";
        //======================================================================
        equ << "\n";
        
        equ
        << "if \"${run_phase} == 1\" then "
        << "\"reset_timestep 0\"" << "\n";
        
        equ << "\n";
        
        //======================================================================
        // NEIGH_MODIFY
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## PAIR"                              << "\n"
        << "# ----------------------------------" << "\n"
        << "neigh_modify delay " << get_delay() << " "
        << "every 1 check yes"
        << "\n";
        
        equ << "\n";
        
        //======================================================================
        // Saving Frequency
        //======================================================================
        if (sysVar.get_systemUnit()=="lj") {
            if (sysVar.get_n_equ_blocks()>=5) {
                equ
                << "# ----------------------------------"            << "\n"
                << "## Saving Frequency (50 frames)"                 << "\n"
                << "# ----------------------------------"            << "\n"
                << "variable ts_save equal floor(${steps_equ}*0.02)" << "\n"
                << "\n";
            }
            else
            {
                equ
                << "# ----------------------------------"            << "\n"
                << "## Saving Frequency (every frames)"              << "\n"
                << "# ----------------------------------"            << "\n"
                << "variable ts_save equal floor(${steps_equ})"      << "\n"
                << "\n";
            }
        }
        else if ((sysVar.get_systemUnit()=="real")||(sysVar.get_systemUnit()=="metal")) {
            equ
            << "# ----------------------------------"            << "\n"
            << "## Saving Frequency (50 frames)"                 << "\n"
            << "# ----------------------------------"            << "\n"
            << "variable ts_save equal floor(${steps_equ}*0.02)" << "\n"
            << "\n";
        }
        
        //======================================================================
        // DUMP
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## DUMP"                              << "\n"
        << "# ----------------------------------" << "\n"
        << "dump 1 all xyz ${ts_save} "
        << "./equilibration/trajectory/trajectory"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".equ.xyz"
        << "\n"
        
        << "dump_modify 1 append yes" << "\n";
        //======================================================================
        equ << "\n";
        
        
        if(1) {
            //==================================================================
            // RESTART CHECKPOINT
            //==================================================================
            equ
            << "# ----------------------------------" << "\n"
            << "## CHECKPOINT"                        << "\n"
            << "# ----------------------------------" << "\n"
            << "if \"${set_CheckPoint} == 1\" then "
            << "\"restart ${restartf} "
            << "./equilibration/restart/restart"
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".equ.restart\""
            << "\n";
            //==================================================================
            equ << "\n";
        }
        
        //======================================================================
        // FIX
        //======================================================================
        int fix=1;
        equ
        << "# ----------------------------------" << "\n"
        << "## FIX"                               << "\n"
        << "# ----------------------------------" << "\n";
        if (is_fixMomentum) {
            equ
            << "fix " << fix << " "
            << "all momentum " << get_fixpsteps() << " "
            << "linear 1 1 1 angular"
            << "\n";
            ++fix;
        }
        //======================================================================
        if (get_is_fixShake()) {
            
            int overallshake=0,bondshake=0,angleshake=0;
            for (int i=0; i<sysVar.get_trueBonds(); ++i)
                if(sysVar.get_is_shakeBonds()[i]) {++overallshake; ++bondshake;}
            for (int i=0; i<sysVar.get_trueAngles();++i)
                if(sysVar.get_is_shakeAngles()[i]) {++overallshake; ++angleshake;}
            
            if(overallshake>0) {
                equ
                << "fix " << fix << " "
                << "polymer shake " << fixed << setprecision(4)
                << shakeTol  << " "
                << shakeIter << " "
                << "${ts_save}";
                
                if (bondshake>0) {
                    equ << " b";
                    for (int i=0; i<sysVar.get_trueBonds(); ++i)
                        if(sysVar.get_is_shakeBonds()[i])
                            equ << " " << i+1;
                }
                if (angleshake>0) {
                    equ << " a";
                    for (int i=0; i<sysVar.get_trueAngles(); ++i)
                        if(sysVar.get_is_shakeAngles()[i])
                            equ << " " << i+1;
                }
                equ << "\n";
                ++fix;
                
                equ.unsetf(ios_base::floatfield);
            }
        }
        //======================================================================
        
        equ
        << "\n"
        << fixed << setprecision(1)
        << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
        
        equ.unsetf(ios_base::floatfield); // reset float precision
        
        equ
        << "fix " << fix << " "
        << "all npt "
        << "temp "
        << "${T} ${T} "
        << get_Tdamp() << " "
        << "iso "
        << get_pressure() << " "
        << get_pressure() << " "
        << get_Pdamp()
        << "\n\n";
        
        
        //======================================================================
        // THERMO_STYLE
        //======================================================================
        equ << "# ----------------------------------" << "\n";
        equ << "## THERMO_STYLE"                      << "\n";
        equ << "# ----------------------------------" << "\n";
        
        if (sysVar.get_is_GPU()) {
            equ
            << "thermo_style custom "
            << "step temp press pe ke etotal ebond eangle epair "
            << "lx ly lz vol density "
            << "dt time cpu tpcpu spcpu cpuremain"
            << "\n\n";
        }
        else {
            equ
            << "thermo_style custom "
            << "step temp press pe ke etotal lx ly lz vol"
            << "\n\n";
        }
        
        
        //======================================================================
        // RUN_STYLE
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## RUN_STYLE"                         << "\n"
        << "# ----------------------------------" << "\n";
        if (get_is_respa()) {
            
            if (sysVar.get_sideGroupInWords()=="stdFENE") {
                equ
                << "run_style respa 2 3 bond 1 pair 2"
                << "\n";
            }
            else {
                equ
                << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                << "\n";
            }
        }
        else {
            
            equ
            << "run_style verlet" << "\n";
        }
        equ << "\n";
        
        
        //======================================================================
        // RUN
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## EQUILIBRATION RUN"                 << "\n"
        << "# ----------------------------------" << "\n"
        << "timestep ${ts}"                       << "\n"
        << "thermo ${ts_save}"                    << "\n"
        << "run ${steps_equ}"                     << "\n"
        << "# ----------------------------------" << "\n";
        //======================================================================
        equ << "\n";
        
        
        //======================================================================
        // WRITE_RESTART
        //======================================================================
        equ
        << "# ----------------------------------" << "\n"
        << "## WRITE_RESTART"                     << "\n"
        << "# ----------------------------------" << "\n";
        /*
         equ
         << "write_restart "
         << "./equilibration/restart/restart"
         << "_${usic}"
         << "_${trial}"
         << "_${namestring}"
         << "_T${Temp}"
         << ".equ.restart"
         << "\n";
         */
        equ << "write_restart ${write_res}";
        //======================================================================
        equ.close();
        cout << "lmpscripts(equ) made!\n\n";system_wait(1);
    }
    else cout << "Equilibration: 'lmp_equ.inp' cannot open." << "\n";
}





void LmpScripts::make_TamperingInputFile(const StructureClass& sysVar) const
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/tampering");
    o.append("/tampering.inp");
    //cout << o << "\n";
    
    ofstream res(o.c_str());
    if (res.is_open()) {}
    else cout << "tampering: 'lmp_tampering.inp' cannot open." << "\n";
}





void LmpScripts::make_ResizeInputFile(const StructureClass& sysVar) const
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/resize");
    o.append("/resize.inp");
    //cout << o << "\n";
    
    ofstream res(o.c_str());
    if (res.is_open()) {
        
        //======================================================================
        // Header
        //======================================================================
        res
        << "# LAMMPS Script for Resize"
        << "\n\n";
        
        if (sysVar.get_is_GPU()) {
            res
            << "package gpu "
            << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
        }
        
        res << "\n";
        //======================================================================
        // READ_RESTART
        // NOTE: Read "Which" restart file
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## READ_RESTART"                      << "\n";
        res << "# ----------------------------------" << "\n";
        
        res
        << "read_restart "
        << "./equilibration/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".equ.restart"
        << "\n";
        //======================================================================
        res << "\n";
        
        res << "reset_timestep 0" << "\n";
        
        res << "\n";
        //======================================================================
        // NEIGH_MODIFY
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## PAIR"                              << "\n";
        res << "# ----------------------------------" << "\n";
        res
        << "neigh_modify delay " << get_delay() << " "
        << "every 1 check yes"
        << "\n\n";
        
        
        //======================================================================
        // Saving Frequency
        //======================================================================
        res << "# ----------------------------------"            << "\n";
        res << "## Saving Frequency (50 frames)"                 << "\n";
        res << "# ----------------------------------"            << "\n";
        res << "variable ts_save equal floor(${steps_res}*0.02)" << "\n";
        res << "\n";
        
        
        //======================================================================
        // DUMP
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## DUMP"                              << "\n";
        res << "# ----------------------------------" << "\n";
        
        res
        << "dump 1 all xyz ${ts_save} "
        << "./resize/trajectory/trajectory"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".res.xyz" << "\n";
        
        res << "dump_modify 1 append yes" << "\n";
        //======================================================================
        res << "\n";
        
        
        
        //======================================================================
        // FIX
        //======================================================================
        int fix=1;
        res << "# ----------------------------------" << "\n";
        res << "## FIX"                               << "\n";
        res << "# ----------------------------------" << "\n";
        
        if (false) {
            res
            << "fix " << fix << " "
            << "all momentum " << get_fixpsteps() << " "
            << "linear 1 1 1 angular"
            << "\n";
            ++fix;
        }
        //======================================================================
        if (get_is_fixShake()) {
            
            int overallshake=0,bondshake=0,angleshake=0;
            for (int i=0; i<sysVar.get_trueBonds(); ++i)
                if(sysVar.get_is_shakeBonds()[i]) {++overallshake; ++bondshake;}
            for (int i=0; i<sysVar.get_trueAngles();++i)
                if(sysVar.get_is_shakeAngles()[i]) {++overallshake; ++angleshake;}
            
            if(overallshake>0) {
                res
                << "fix " << fix << " "
                << "polymer shake " << fixed << setprecision(4)
                << shakeTol  << " "
                << shakeIter << " "
                << "${ts_save}";
                
                if (bondshake>0) {
                    res << " b";
                    for (int i=0; i<sysVar.get_trueBonds(); ++i)
                        if(sysVar.get_is_shakeBonds()[i])
                            res << " " << i+1;
                }
                if (angleshake>0) {
                    res << " a";
                    for (int i=0; i<sysVar.get_trueAngles(); ++i)
                        if(sysVar.get_is_shakeAngles()[i])
                            res << " " << i+1;
                }
                res << "\n";
                ++fix;
                
                res.unsetf(ios_base::floatfield);
            }
        }
        //======================================================================
        
        
        res
        << "fix " << fix << " "
        << "all deform 100 "
        << "x final -${halfX} ${halfX} "
        << "y final -${halfY} ${halfY} "
        << "z final -${halfZ} ${halfZ} "
        << "remap x "
        << "units box "
        << "\n\n";
        
        
        //======================================================================
        // THERMO_STYLE
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## THERMO_STYLE"                      << "\n";
        res << "# ----------------------------------" << "\n";
        
        if (sysVar.get_is_GPU()) {
            res
            << "thermo_style custom "
            << "step temp press pe ke etotal ebond eangle epair "
            << "lx ly lz vol density "
            << "dt time cpu tpcpu spcpu cpuremain"
            << "\n\n";
        }
        else {
            res
            << "thermo_style custom "
            << "step temp press pe ke etotal lx ly lz vol"
            << "\n\n";
        }
        
        
        //======================================================================
        // RUN_STYLE
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## RUN_STYLE"                         << "\n";
        res << "# ----------------------------------" << "\n";
        if (get_is_respa()) {
            
            if (sysVar.get_sideGroupInWords()=="stdFENE") {
                res
                << "run_style respa 2 3 bond 1 pair 2"
                << "\n";
            }
            else {
                res
                << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                << "\n";
            }
        }
        else {
            
            res
            << "run_style verlet" << "\n";
        }
        
        res << "\n";
        
        
        //======================================================================
        // RESIZE RUN
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## RESIZE RUN"                        << "\n";
        res << "# ----------------------------------" << "\n";
        res << "timestep ${ts}"                       << "\n";
        res << "thermo ${ts_save}"                    << "\n";
        res << "run ${steps_res}"                     << "\n";
        res << "# ----------------------------------" << "\n";
        //======================================================================
        res << "\n";
        
        
        
        //======================================================================
        // POST-RESIZE RUN
        //======================================================================
        res
        << "unfix " << fix << "\n\n";
        
        res
        << fixed << setprecision(1)
        << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
        
        res.unsetf(ios_base::floatfield);
        
        res
        << "fix " << fix << " "
        << "all nvt "
        << "temp "
        << "${T} ${T} "
        << get_Tdamp()
        << "\n\n";
        
        res << "# ----------------------------------" << "\n";
        res << "## POST-RESIZE RUN"                   << "\n";
        res << "# ----------------------------------" << "\n";
        res << "timestep ${ts}"                       << "\n";
        res << "thermo ${ts_save}"                    << "\n";
        res << "run ${steps_res}"                     << "\n";
        res << "# ----------------------------------" << "\n";
        
        res << "\n";
        //======================================================================
        // WRITE_RESTART
        //======================================================================
        res << "# ----------------------------------" << "\n";
        res << "## WRITE_RESTART"                     << "\n";
        res << "# ----------------------------------" << "\n";
        
        res
        << "write_restart "
        << "./resize/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".res.restart"
        << "\n";
        //======================================================================
        res.close();
        cout << "lmpscripts(res) made!\n\n";system_wait(1);
    }
    else cout << "Resize: 'lmp_res.inp' cannot open." << "\n";
}





void LmpScripts::make_ProductionInputFile(const StructureClass& sysVar,
                                          const int n_sys) const
{
    
    /* Exponential Timestepping */
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/lammps_inputs/production");
    o.append("/production.inp");
    //cout << o << "\n";
    
    //const int n_poly_d  = sysVar.get_n_polyVec()[n_sys];
    
    ofstream prd(o.c_str());
    if (prd.is_open()) {
        
        //======================================================================
        // Header
        //======================================================================
        prd
        << "# LAMMPS Script for Production (Exponential Timestepping)"
        << "\n\n";
        
        if (sysVar.get_is_GPU()) {
            prd
            << "package gpu "
            << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
        }
        
        prd << "\n";
        
        //======================================================================
        // READ_RESTART
        // NOTE: Read "Which" restart file
        //======================================================================
        prd
        << "# ----------------------------------"  << "\n"
        << "## READ_RESTART"                       << "\n"
        << "# ----------------------------------"  << "\n";
        
        if (sysVar.get_is_singleTempTest()) {
            
            prd
            << "read_restart "
            << "./generation/restart/restart"
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".gen.restart"
            << "\n";
        }
        else {
            
            if (get_is_resize()) {
                
                prd
                << "read_restart "
                << "./resize/restart/restart"
                << "_${usic}"
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << ".res.restart"
                << "\n";
            }
            else {
                
                prd
                << "read_restart "
                << "./equilibration/restart/restart"
                << "_${usic}"
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << ".equ.restart"
                << "\n";
            }
        }
        
        //======================================================================
        prd << "\n";
        
        prd << "reset_timestep 0" << "\n";
        
        prd << "\n";
        //======================================================================
        // NEIGH_MODIFY
        //======================================================================
        prd
        << "# ----------------------------------"  << "\n"
        << "## PAIR"                               << "\n"
        << "# ----------------------------------"  << "\n"
        << "neigh_modify delay " << get_delay() << " "
        << "every 1 check yes"
        << "\n\n";
        
        
        //======================================================================
        // DUMP
        // ** The file path affects the AMDAT inputfile
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## DUMP"                              << "\n"
        << "# ----------------------------------" << "\n"
        << "dump 1 all custom 1 "
        << "./production/trajectory/trajectory"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".prd.custom "
        << "type x y z ix iy iz"
        << "\n"
        << "dump_modify 1 first no sort id"
        << "\n";
        //======================================================================
        prd << "\n";
        
        
        
        //======================================================================
        // COMPUTE
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## COMPUTE"                           << "\n"
        << "# ----------------------------------" << "\n"
        << "compute 1 polymer gyration/molecule"  << "\n";
        //======================================================================
        prd << "\n";
        
        
        //======================================================================
        // FIX
        //======================================================================
        int fix=1;
        prd
        << "# ----------------------------------" << "\n"
        << "## FIX"                               << "\n"
        << "# ----------------------------------" << "\n";
        if (get_is_fixPrdMomentum()) {
            prd
            << "fix " << fix << " "
            << "all momentum " << get_fixpsteps_prd() << " "
            << "linear 1 1 1 angular"
            << "\n";
            ++fix;
        }
        //======================================================================
        if (get_is_fixShake()) {
            
            int overallshake=0,bondshake=0,angleshake=0;
            for (int i=0; i<sysVar.get_trueBonds(); ++i)
                if(sysVar.get_is_shakeBonds()[i]) {++overallshake; ++bondshake;}
            for (int i=0; i<sysVar.get_trueAngles();++i)
                if(sysVar.get_is_shakeAngles()[i]) {++overallshake; ++angleshake;}
            
            if(overallshake>0) {
                prd
                << "fix " << fix << " "
                << "polymer shake " << fixed << setprecision(4)
                << shakeTol  << " "
                << shakeIter << " "
                << "0";
                
                if (bondshake>0) {
                    prd << " b";
                    for (int i=0; i<sysVar.get_trueBonds(); ++i)
                        if(sysVar.get_is_shakeBonds()[i])
                            prd << " " << i+1;
                }
                if (angleshake>0) {
                    prd << " a";
                    for (int i=0; i<sysVar.get_trueAngles(); ++i)
                        if(sysVar.get_is_shakeAngles()[i])
                            prd << " " << i+1;
                }
                prd << "\n";
                ++fix;
                
                prd.unsetf(ios_base::floatfield);
            }
        }
        //======================================================================
        
        if (get_is_resize()) {
            
            prd
            << "\n"
            << fixed << setprecision(1)
            << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
            
            prd.unsetf(ios_base::floatfield);
            
            prd
            << "fix " << fix << " "
            << "all nvt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp()
            << "\n";
        }
        else {
            
            prd
            << "\n"
            << fixed << setprecision(1)
            << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
            
            prd.unsetf(ios_base::floatfield);
            
            prd
            << "fix " << fix << " "
            << "all npt "
            << "temp "
            << "${T} ${T} "
            << get_Tdamp() << " "
            << "iso "
            << get_pressure() << " "
            << get_pressure() << " "
            << get_Pdamp()
            << "\n\n";
        }
        
        
        //======================================================================
        // THERMO_STYLE
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## THERMO_STYLE"                      << "\n"
        << "# ----------------------------------" << "\n";
        //<< "thermo_style custom "
        //<< "step temp press pe ke etotal ebond eangle epair "
        //<< "lx ly lz vol density "
        //<< "dt time cpu tpcpu spcpu ";
        //for (int i=1; i<=n_poly_d; ++i) prd << "c_1[" << i << "]" << " ";
        if (sysVar.get_is_GPU()) {
            prd
            << "thermo_style custom "
            << "step temp press pe ke etotal ebond eangle epair "
            << "lx ly lz vol density "
            << "dt time cpu tpcpu spcpu cpuremain"
            << "\n\n";
        }
        else {
            prd
            << "thermo_style custom "
            << "step temp press pe ke etotal lx ly lz vol"
            << "\n\n";
        }
        
        
        //======================================================================
        // RUN_STYLE
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## RUN_STYLE"                         << "\n"
        << "# ----------------------------------" << "\n";
        if (get_is_respa()) {
            
            if (sysVar.get_sideGroupInWords()=="stdFENE") {
                prd
                << "run_style respa 2 3 bond 1 pair 2"
                << "\n";
            }
            else {
                prd
                << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                << "\n";
            }
        }
        else {
            
            prd
            << "run_style verlet" << "\n";
        }
        
        prd << "\n";
        
        //======================================================================
        // TIMESTEP
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## TIMESTEP"                          << "\n"
        << "# ----------------------------------" << "\n"
        << "timestep ${ts}"                       << "\n";
        
        prd << "\n";
        
        //======================================================================
        // RUN
        //======================================================================
        prd << "# ----------------------------------------------------" << "\n";
        prd << "## Exponential Timestepping"                            << "\n";
        prd << "# ----------------------------------------------------" << "\n";
        prd << "variable EXP equal ${blocksize}-17"                     << "\n";
        prd << "variable s equal 0"                                     << "\n";
        prd << "variable l loop ${n_prd_blocks}  ## NUMBER OF BLOCKS"   << "\n";
        prd << "# ----------------------------------------------------" << "\n";
        prd << "label outerloop"                                        << "\n";
        
        prd << "\n";
        
        prd
        << "if \"${divide_prd} > 1\" then "
        << "\"reset_timestep 0\""
        << "\n"
        << "if \"${divide_prd} > 1\" then "
        << "\"variable s equal 0\""
        << "\n";
        
        prd << "\n";
        
        prd << "  # ----------------------------"                       << "\n";
        prd << "  ## 1. \"PSEUDO-LINEAR\""                              << "\n";
        prd << "  # ----------------------------"                       << "\n";
        prd << "  variable a loop 16"                                   << "\n";
        prd << "  label linearloop"                                     << "\n";
        prd << "    variable s equal $s+1"                              << "\n";
        prd << "    dump_modify 1 every $s append yes sort id"          << "\n";
        prd << "    thermo $s"                                          << "\n";
        prd << "    run 1"                                              << "\n";
        prd << "    print \"$s\""                                       << "\n";
        prd << "  next a"                                               << "\n";
        prd << "  jump SELF linearloop"                                 << "\n";
        prd << "  # ----------------------------"                       << "\n";
        
        prd << "\n";
        
        prd << "  # ----------------------------"                       << "\n";
        prd << "  ## 2. \"BRIDGING\""                                   << "\n";
        prd << "  # ----------------------------"                       << "\n";
        prd << "  variable m equal floor(${exp_base}^(16))-16"          << "\n";
        prd << "  variable s equal $s+$m"                               << "\n";
        prd << "  dump_modify 1 append yes every $s sort id"            << "\n";
        prd << "  thermo $s"                                            << "\n";
        prd << "  run $m"                                               << "\n";
        prd << "  print \"$s\""                                         << "\n";
        prd << "  # ----------------------------"                       << "\n";
        
        prd << "\n";
        
        prd << "  # ----------------------------"                       << "\n";
        prd << "  ## 3. \"EXPONENTIAL TIMESTEPS\""                      << "\n";
        prd << "  # ----------------------------"                       << "\n";
        prd << "  variable c loop ${EXP}"                               << "\n";
        prd << "  label innerloop"                                      << "\n";
        prd << "    variable m equal floor(${exp_base}^($c+16))-"
        << "floor(${exp_base}^($c+16-1))"                               << "\n";
        prd << "    variable s equal $s+$m"                             << "\n";
        prd << "    dump_modify 1 append yes every $s sort id"          << "\n";
        prd << "    thermo $s"                                          << "\n";
        prd << "    run $m"                                             << "\n";
        prd << "    print \"$s\""                                       << "\n";
        prd << "  next c"                                               << "\n";
        prd << "  jump SELF innerloop"                                  << "\n";
        prd << "  # ----------------------------"                       << "\n";
        
        prd << "\n";
        
        prd
        << "if \"${set_CheckPoint} == 1\" then "
        << "\"write_restart "
        << "./production/restart/restart"
        << "_${l}_${n_prd_blocks}"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".prd.restart\""
        << "\n";
        
        prd << "\n";
        
        prd << "next l"                                                 << "\n";
        prd << "jump SELF outerloop"                                    << "\n";
        prd << "# ----------------------------------------------------" << "\n";
        prd << "\n";
        //======================================================================
        
        
        //======================================================================
        // WRITE_RESTART
        // The file path will be watched for
        //======================================================================
        prd
        << "# ----------------------------------" << "\n"
        << "## WRITE_RESTART"                     << "\n"
        << "# ----------------------------------" << "\n"
        << "write_restart "
        << "./production/restart/restart"
        << "_${usic}"
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".prd.restart"
        << "\n";
        //======================================================================
        prd.close();
        cout << "lmpscripts(prd) made!\n\n";system_wait(1);
    }
    else cout << "Production: 'lmp_prd.inp' cannot open." << "\n";
    
    
    
    
    /* Linear Timestepping */
    // For # of timesteps > 2^31 in a single run
    
    if (0) {
        
        o.clear();
        o.append(return_SimulationFolderPath(sysVar));
        o.append("/lammps_inputs/production");
        o.append("/production_linear.inp");
        //cout << o << "\n";
        
        ofstream lin(o.c_str());
        if (lin.is_open()) {
            
            //==================================================================
            // Header
            //==================================================================
            lin
            << "# LAMMPS Script for Production (Linear Timestepping)"
            << "\n\n";
            
            if (sysVar.get_is_GPU()) {
                lin
                << "package gpu "
                << get_GPU_mode() << " ${GPU} ${GPU} -1.0" << "\n";
            }
            
            lin << "\n";
            
            //==================================================================
            // READ_RESTART
            // NOTE: Read "Which" restart file
            //==================================================================
            lin
            << "# ----------------------------------"  << "\n"
            << "## READ_RESTART"                       << "\n"
            << "# ----------------------------------"  << "\n"
            << "read_restart ${read_res}" << "\n";
            //==================================================================
            lin << "\n";
            
            lin
            << "if \"${run_phase} == 1\" then "
            << "\"reset_timestep 0\""
            << "\n";
            
            lin << "\n";
            
            //==================================================================
            // NEIGH_MODIFY
            //==================================================================
            lin
            << "# ----------------------------------"  << "\n"
            << "## PAIR"                               << "\n"
            << "# ----------------------------------"  << "\n"
            << "neigh_modify delay " << get_delay() << " "
            << "every 1 check yes"
            << "\n";
            
            lin << "\n";
            
            //==================================================================
            // Saving Frequency
            //==================================================================
            lin
            << "# ----------------------------------"            << "\n"
            << "## Saving Frequency (50 frames)"                 << "\n"
            << "# ----------------------------------"            << "\n"
            << "variable ts_save equal floor(${steps_prd}*0.02)" << "\n";
            
            lin << "\n";
            
            //==================================================================
            // DUMP
            // ** The file path affects the AMDAT inputfile
            //==================================================================
            lin
            << "# ----------------------------------" << "\n"
            << "## DUMP"                              << "\n"
            << "# ----------------------------------" << "\n"
            << "dump 1 all custom ${ts_save} "
            << "./production/trajectory/trajectory"
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".prd.custom "
            << "type x y z ix iy iz"
            << "\n"
            << "dump_modify 1 append yes first no sort id"
            << "\n";
            //==================================================================
            
            lin << "\n";
            
            //==================================================================
            // COMPUTE
            //==================================================================
            lin
            << "# ----------------------------------" << "\n"
            << "## COMPUTE"                           << "\n"
            << "# ----------------------------------" << "\n"
            << "compute 1 polymer gyration/molecule"
            << "\n";
            //==================================================================
            
            lin << "\n";
            
            //==================================================================
            // FIX
            //==================================================================
            int fix=1;
            lin
            << "# ----------------------------------" << "\n"
            << "## FIX"                               << "\n"
            << "# ----------------------------------" << "\n";
            if (get_is_fixPrdMomentum()) {
                lin
                << "fix " << fix << " "
                << "all momentum " << get_fixpsteps_prd() << " "
                << "linear 1 1 1 angular"
                << "\n";
                ++fix;
            }
            //==================================================================
            if (get_is_fixShake()) {
                
                int overallshake=0,bondshake=0,angleshake=0;
                for (int i=0; i<sysVar.get_trueBonds(); ++i)
                    if(sysVar.get_is_shakeBonds()[i]) {++overallshake; ++bondshake;}
                for (int i=0; i<sysVar.get_trueAngles();++i)
                    if(sysVar.get_is_shakeAngles()[i]) {++overallshake; ++angleshake;}
                
                if(overallshake>0) {
                    lin
                    << "fix " << fix << " "
                    << "polymer shake " << fixed << setprecision(4)
                    << shakeTol  << " "
                    << shakeIter << " "
                    << "0";
                    
                    if (bondshake>0) {
                        lin << " b";
                        for (int i=0; i<sysVar.get_trueBonds(); ++i)
                            if(sysVar.get_is_shakeBonds()[i])
                                lin << " " << i+1;
                    }
                    if (angleshake>0) {
                        lin << " a";
                        for (int i=0; i<sysVar.get_trueAngles(); ++i)
                            if(sysVar.get_is_shakeAngles()[i])
                                lin << " " << i+1;
                    }
                    lin << "\n";
                    ++fix;
                    
                    lin.unsetf(ios_base::floatfield);
                }
            }
            //==================================================================
            
            if (get_is_resize()) {
                
                lin
                << "\n"
                << fixed << setprecision(1)
                << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
                
                lin.unsetf(ios_base::floatfield);
                
                lin
                << "fix " << fix << " "
                << "all nvt "
                << "temp "
                << "${T} ${T} "
                << get_Tdamp()
                << "\n";
            }
            else {
                
                lin
                << "\n"
                << fixed << setprecision(1)
                << "variable T equal ${Temp}/" << sysVar.get_corrFacforT() << "\n\n";
                
                lin.unsetf(ios_base::floatfield);
                
                lin
                << "fix " << fix << " "
                << "all npt "
                << "temp "
                << "${T} ${T} "
                << get_Tdamp() << " "
                << "iso "
                << get_pressure() << " "
                << get_pressure() << " "
                << get_Pdamp()
                << "\n\n";
            }
            
            lin << "\n";
            
            //==================================================================
            // THERMO_STYLE
            //==================================================================
            //lin
            //<< "# ----------------------------------" << "\n"
            //<< "## THERMO_STYLE"                      << "\n"
            //<< "# ----------------------------------" << "\n"
            //<< "thermo_style custom "
            //<< "step temp press pe ke etotal ebond eangle epair "
            //<< "lx ly lz vol density "
            //<< "dt time cpu tpcpu spcpu ";
            //for (int i=1; i<=n_poly_d; ++i) lin << "c_1[" << i << "]" << " ";
            
            if (sysVar.get_is_GPU()) {
                lin
                << "thermo_style custom "
                << "step temp press pe ke etotal ebond eangle epair "
                << "lx ly lz vol density "
                << "dt time cpu tpcpu spcpu cpuremain"
                << "\n\n";
            }
            else {
                lin
                << "thermo_style custom "
                << "step temp press pe ke etotal lx ly lz vol"
                << "\n\n";
            }
            
            
            //==================================================================
            // RUN_STYLE
            //==================================================================
            lin
            << "# ----------------------------------" << "\n"
            << "## RUN_STYLE"                         << "\n"
            << "# ----------------------------------" << "\n";
            if (get_is_respa()) {
                
                if (sysVar.get_sideGroupInWords()=="stdFENE") {
                    lin
                    << "run_style respa 2 3 bond 1 pair 2"
                    << "\n";
                }
                else {
                    lin
                    << "run_style respa 3 2 2 bond 1 angle 2 pair 3"
                    << "\n";
                }
            }
            else {
                
                lin
                << "run_style verlet" << "\n";
            }
            
            lin << "\n";
            
            //==================================================================
            // RUN
            //==================================================================
            lin
            << "# ----------------------------------" << "\n"
            << "## Linear Timestepping"               << "\n"
            << "# ----------------------------------" << "\n"
            << "timestep ${ts}"                       << "\n"
            << "thermo ${ts_save}"                    << "\n"
            << "run ${steps_equ}"                     << "\n"
            << "\n";
            //==================================================================
            
            
            //==================================================================
            // WRITE_RESTART
            // The file path will be watched for
            //==================================================================
            lin
            << "# ----------------------------------" << "\n"
            << "## WRITE_RESTART"                     << "\n"
            << "# ----------------------------------" << "\n"
            << "write_restart ${write_res}";
            //==================================================================
            
            lin.close();
        }
        else cout << "Production_linear: 'lmp_prd.inp' cannot open." << "\n";
    }
    
}





const string LmpScripts::work_target(const StructureClass& sysVar,
                                     const string& state,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d) const
{
    string target;
    
    if (state=="equ") {
        target.append(return_SimulationFolderPath(sysVar));
        target.append("/equilibration/restart");
        target.append("/restart_");
        target.append(sysVar.get_usic());
        target.append("_00"+to_string((long long int)n_trl));
        target.append("_"+sysVar.get_nameString(n_sys));
        target.append("_T"+to_string((long long int)Temp_d));
        target.append(".equ.restart");
        //cout << target << "\n";
    }
    else if (state=="prd") {
        target.append(return_SimulationFolderPath(sysVar));
        target.append("/production/restart");
        target.append("/restart_");
        target.append(sysVar.get_usic());
        target.append("_00"+to_string((long long int)n_trl));
        target.append("_"+sysVar.get_nameString(n_sys));
        target.append("_T"+to_string((long long int)Temp_d));
        target.append(".prd.restart");
        //cout << target << "\n";
    }
    
    
    return target;
}





/** public setters **/
//----------------------------------------------------------------------
/* bool */
void LmpScripts::set_is_tampering(const bool b){is_tampering=b;}
void LmpScripts::set_is_resize(const bool b){is_resize=b;}
void LmpScripts::set_is_fixMomentum(const bool b){is_fixMomentum=b;}
void LmpScripts::set_is_fixPrdMomentum(const bool b){is_fixPrdMomentum=b;}
void LmpScripts::set_is_fixShake(const bool b){is_fixShake=b;}
void LmpScripts::set_is_respa(const bool b){is_respa=b;}
void LmpScripts::set_is_bondswap(const bool b){is_bondswap=b;}
/* int */
void LmpScripts::set_fixpsteps(const int i){fixpsteps=i;}
void LmpScripts::set_fixpsteps_prd(const int i){fixpsteps_prd=i;}
void LmpScripts::set_delay(const int i){delay=i;}
void LmpScripts::set_shakeIter(const int i){shakeIter=i;}
void LmpScripts::set_print_shake(const int i){print_shake=i;}
/* double */
void LmpScripts::set_shakeTol(const double f){shakeTol=f;}
void LmpScripts::set_Tdamp(const double f){Tdamp=f;}
void LmpScripts::set_Pdamp(const double f){Pdamp=f;}
void LmpScripts::set_pressure(const double d){pressure=d;}
/* string */
void LmpScripts::set_GPU_mode(const string& s){GPU_mode=s;}




