//
//  Created by Sean Jh H. on 7/4/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "amdatanalysis.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

AmdatAnalysis::AmdatAnalysis(const StructureClass& sysVar,
                             const WorkScripts& ws):

/* bool */
is_changeAtomType(false),
is_changeFrameSteps(false),
is_NPT(true),
is_fixed_q(false),
is_fixed_maxL(false),

is_strFac(false),
is_msd(false),
is_ngp(false),
is_isfs(false),

is_composition(false),
is_u2dist(false),
is_stiffness_dist(false),
is_isf(false),
is_mbodies(false),

is_strings(false),
is_peak_frame(true),

/* int */
amdat_cores(1),
amdat_priority(0),
fullblock(0),
vhs_nbins(75),
u2dist_nbins(100),
stiffness_nbins(100),

/* double */
DWF_time(0),
vhs_rangebinned(1.5),
strings_threshold(0.6),
max_u2(1.5),
max_stiffness(1.0),

threshold_percentile(93.5),
threshold_percentile_lo(45.0),
threshold_percentile_hi(55.0),

/* string */
amdat_exe(),
symmetry("symmetric"),
geometry("xyz"),
analysispart(),
threshold_keyword("greater"),
centertypename("centroid")

{
    
    //amdat_exe="/home/jh148/AMDAT/AMDAT";   // amdat version 4.2.0
    //amdat_exe="/home/jh148/AMDAT81/AMDAT"; // amdat subversion 81
    //amdat_exe="/home/jh148/AMDAT84/AMDAT"; // amdat subversion 84
    amdat_exe="/home/jh148/AMDAT89/AMDAT"; // amdat subversion 89
    
    if (sysVar.get_systemUnit()=="real")
    {
        DWF_time = 1e+3; // fs
    }
    else if (sysVar.get_systemUnit()=="lj")
    {
        DWF_time = 1.0;  // tau
    }
    else if (sysVar.get_systemUnit()=="metal")
    {
        DWF_time = 1.0;  // ps
    }
    
    prd_blocksize = ws.get_prd_blocksize();
    prd_exp_base  = ws.get_prd_exp_base();
    n_prd_blocks  = ws.get_n_prd_blocks();
}





void AmdatAnalysis::systemBlock_definition(ofstream& ai,
                                           const string& path)
{
    
    /*------------------------
     species list:
     
     1. KGPolymer
     2. binaryLJ
     3. Cu4Ag6
     4. SiO2
     5. MartiniPS30
     ------------------------*/
    string species="KGPolymer";
    
    
    if (species=="KGPolymer") { // KGPolymer
        
        strings_threshold=0.6;
        
        is_mbodies=false;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << path                                           << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base
        << " 0 0 ${ts} "                                  << "\n"
        << "polymer 80"                                   << "\n"
        << "1"                                            << "\n"
        << "20"                                           << "\n"
        
        << "create_list all"        << "\n"
        << "all"                    << "\n\n";
        
        ai << "\n\n";
    }
    
    else if (species=="binaryLJ") { // binaryLJ
        
        strings_threshold *= 0.89;
        
        is_mbodies=false;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << path                                           << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base
        << " 0 0 ${ts} "                                  << "\n"
        << "type1 800 type2 200"                          << "\n"
        << "1 2"                                          << "\n"
        << "1 0"                                          << "\n"
        << "0 1"                                          << "\n\n"
        
        << "create_list all"        << "\n"
        << "all"                    << "\n\n"
        
        << "create_list type1_list" << "\n"
        << "type_species type1 1"   << "\n"
        
        << "create_list type2_list" << "\n"
        << "type_species type2 2"   << "\n";
        
        ai << "\n\n";
    }
    
    else if (species=="Cu4Ag6") { // Cu4Ag6
        
        strings_threshold *= 0.89;
        
        is_mbodies=false;
        
        ai
        << "system_nv"                                    << "\n" // NOTE: system_nv
        << "custom"                                       << "\n"
        << path                                           << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base
        << " 0 0 ${ts} "                                  << "\n"
        << "system 1"                                     << "\n"
        << "1 2"                                          << "\n"
        << "400 600"                                      << "\n\n"
        
        << "create_list all"        << "\n"
        << "all"                    << "\n\n"
        
        << "create_list type1_list" << "\n"
        << "type_species system 1"  << "\n"
        
        << "create_list type2_list" << "\n"
        << "type_species system 2"  << "\n";
        
        ai << "\n\n";
    }
    
    else if (species=="SiO2") { // SiO2
        
        strings_threshold *= 0.89;
        
        is_mbodies=false;
        
        ai
        << "system_nv"                                    << "\n" // NOTE: system_nv
        << "custom"                                       << "\n"
        << path                                           << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base
        << " 0 0 ${ts} "                                  << "\n"
        << "system 1"                                     << "\n"
        << "1 2"                                          << "\n"
        << "512 1024"                                     << "\n\n"
        
        << "create_list all"        << "\n"
        << "all"                    << "\n\n"
        
        << "create_list type1_list" << "\n"
        << "type_species system 1"  << "\n"
        
        << "create_list type2_list" << "\n"
        << "type_species system 2"  << "\n";
        
        ai << "\n\n";
    }
    
    else if (species=="MartiniPS30") { // MartiniPS30
        
        strings_threshold *= 0.89;
        
        is_mbodies=true;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << path                                           << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base
        << " 0 0 ${ts} "                                  << "\n"
        << "monomer 960"                                  << "\n"
        << "1 2"                                          << "\n"
        << "1 3"                                          << "\n\n"
        
        << "create_list all"        << "\n"
        << "all"                    << "\n\n"
        
        << "create_multibodies mbodies all_molecule"
        << "\n"
        << "#create_multibodies mbodies species_molecule monomer"
        << "\n"
        << "#create_multibodies mbodies species_type monomer 1"
        << "\n\n"
        
        << "trajectories_from_multibodies mbodies_list "
        << "mbodies "
        << centertypename
        << "\n\n";
        
        ai << "\n\n";
    }
}





void AmdatAnalysis::make_amdatInputFile(const StructureClass& sysVar,
                                        const WorkScripts& ws)
{
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_inputs");
    
    if (is_strings) {
        o.append("/amdat_strings.inp");
    }
    else {
        o.append("/amdat.inp");
    }
    ofstream ai(o.c_str());
    
    
    string i;
    if (get_is_changeAtomType()||get_is_deleteZerothFrames()) {
        
        i.append(return_SimulationFolderPath(sysVar));
        i.append("/production/trajectory");
        
        i.append("/new_trajectory_"); // NOTE: "new_"
        i.append(sysVar.get_usic());
        i.append("_${trial}");
        i.append("_${namestring}");
        i.append("_T${Temp}");
        i.append(".prd.custom");
        //cout << i << "\n";
    }
    else {
        
        i.append(return_SimulationFolderPath(sysVar));
        i.append("/production/trajectory");
        
        i.append("/trajectory_");
        i.append(sysVar.get_usic());
        i.append("_${trial}");
        i.append("_${namestring}");
        i.append("_T${Temp}");
        i.append(".prd.custom");
        //cout << i << "\n";
    }
    
    
    if (true) { // manually specify
        
        systemBlock_definition(ai,i);
        
    }
    
    else { // use programmed structure
        
        
        if(get_is_NPT())
            ai << "system_np" << "\n";
        else
            ai << "system_nv" << "\n";
        
        ai << "custom" << "\n";
        
        ai << i << "\n"; // trajectory file from production
        
        ai
        << "exponential ${n_prd_blocks} ${blocksize} "
        << prd_exp_base << " 0 0 ${ts} "
        << "\n";
        
        ai
        << "polymer ${n_poly}"
        << "\n";
        
        for (int i=1; i<=sysVar.get_trueAtoms(); ++i) {
            ai << i << " ";
        }
        
        ai << "\n";
        ai << "${backLen} ";
        if (sysVar.get_trueAtoms()>1) {
            
            if (sysVar.get_typeBackbone()!=sysVar.linear) {
                
                if (sysVar.get_trueAtoms()!=2) {
                    
                    for (int i=1; i<sysVar.get_trueAtoms()-1; ++i) {
                        ai << "${backLen} ";
                    }
                    ai << "${numSide}" << "\n";
                }
                else {
                    ai << "${backLen}" << "\n";
                }
                
            }
            else {
                ai << "${numSide}" << "\n";
            }
        }
        /*
         if (sysVar.get_typeBackbone()==sysVar.linear) {
         ai << "1 2" << "\n";
         ai << "${backLen} ${numSide}" << "\n";
         }
         else {
         ai << "1 2 3" << "\n";
         ai << "${backLen} ${backLen} ${numSide}" << "\n";
         }
         */
        //======================================================================
        ai << "\n";
        
        
        ai << "create_list all" << "\n";
        ai << "all" << "\n";
        ai << "\n";
        
        
        //======================================================================
        // Define List of Backbone
        //======================================================================
        ai << "create_list backbone" << "\n";
        ai << "type_species polymer ";
        if (sysVar.get_typeBackbone()==sysVar.linear)
            ai << "1";
        else
            ai << "1";
        //======================================================================
        ai << "\n\n";
        
        
        //======================================================================
        // Define List of Side Group
        //======================================================================
        if ((sysVar.get_typeBackbone()==sysVar.linear) &&
            (sysVar.get_trueAtoms()>1))
        {
            ai << "create_list sideGroup" << "\n";
            ai << "type_species polymer 2";
        }
        //======================================================================
        
        ai << "\n\n";
    }
    
    
    
    //======================================================================
    // NOTE:
    // Which AMDAT analysis to be taken as the major target
    // will affect the returned path strings
    //======================================================================
    
    /* Mean Square Displacement */
    if (is_msd) {
        
        if (is_mbodies) {
            
            ai
            << "msd "
            << "./statistics/msd/all/msd_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n";
            
            ai << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "msd "
            << "./statistics/msd/all/msd_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n";
            
            ai << "list all" << "\n";
            
            ai
            << "msd "
            << "./statistics/msd/backbone/msd_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << "\n";
            
            ai << "list backbone" << "\n";
            
            ai << "\n";
        }
    }
    
    /* Calculates non-Gaussian parameter of the mean squared displacement */
    if (is_ngp) {
        
        if (is_mbodies) {
            
            ai
            << "ngp "
            << "./statistics/ngp/all/ngp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n";
            
            ai << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "ngp "
            << "./statistics/ngp/all/ngp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n";
            
            ai << "list all" << "\n";
            
            ai
            << "ngp "
            << "./statistics/ngp/backbone/ngp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << "\n";
            
            ai << "list backbone" << "\n";
            
            ai << "\n";
        }
    }
    
    /* Structure Factor */
    if (is_strFac) {
        
        if (is_mbodies) {
            
            ai
            << "structure_factor "
            << "./statistics/strFac/all/strFac_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << " "
            
            << get_symmetry()   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n";
            
            ai << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "structure_factor "
            << "./statistics/strFac/all/strFac_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << " "
            
            << get_symmetry()   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n";
            
            ai << "list all" << "\n";
            
            ai
            << "structure_factor "
            << "./statistics/strFac/backbone/strFac_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << " "
            << get_symmetry()   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n";
            
            ai << "list backbone" << "\n";
            
            ai << "\n";
        }
    }
    
    /* Calculates the composition and number density */
    if (is_composition) {
        
        if (is_mbodies) {
            
            ai
            << "composition "
            << "./statistics/composition/all/comp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n"
            
            << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "composition "
            << "./statistics/composition/all/comp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << "\n"
            << "list all"
            << "\n";
            
            ai
            << "composition "
            << "./statistics/composition/backbone/comp_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << "\n"
            << "list backbone"
            << "\n";
            
            ai << "\n";
        }
    }
    
    /* Calculates distribution of square displacements at a specified time */
    if (is_u2dist) {
        
        if (is_mbodies) {
            
            ai
            << "u2dist "
            << "./statistics/u2dist/all/u2dist_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            << u2dist_nbins
            << " "
            << max_u2
            << " "
            << "${dwf_frame}"
            << "\n"
            
            << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "u2dist "
            << "./statistics/u2dist/all/u2dist_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            << u2dist_nbins
            << " "
            << max_u2
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list all"
            << "\n";
            
            ai
            << "u2dist "
            << "./statistics/u2dist/backbone/u2dist_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << " "
            << u2dist_nbins
            << " "
            << max_u2
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list backbone"
            << "\n";
            
            ai << "\n";
        }
    }
    
    /* Calculates distribution of inverse Debye-Waller factor values 1/u2 */
    if (is_stiffness_dist) {
        
        if (is_mbodies) {
            
            ai
            << "stiffness_dist "
            << "./statistics/stiffness/all/stiffness_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            << stiffness_nbins
            << " "
            << max_stiffness
            << " "
            << "${dwf_frame}"
            << "\n"
            
            << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "stiffness_dist "
            << "./statistics/stiffness/all/stiffness_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            << stiffness_nbins
            << " "
            << max_stiffness
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list all"
            << "\n";
            
            ai
            << "stiffness_dist "
            << "./statistics/stiffness/backbone/stiffness_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << " "
            << stiffness_nbins
            << " "
            << max_stiffness
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list backbone"
            << "\n";
            
            ai << "\n";
        }
    }
    
    /* Self Part of Intermediate Scattering Function */
    if (is_isfs) {
        
        if (is_mbodies) {
            
            ai
            << "isfs "
            << "./statistics/isfs/all/isfs_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << " "
            
            << "${waveindex}"   << " "
            << "${waveindex}"   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n"
            
            << "list mbodies_list" << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "isfs "
            << "./statistics/isfs/all/isfs_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << " "
            
            << "${waveindex}"   << " "
            << "${waveindex}"   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n"
            << "list all"
            << "\n";
            
            ai
            << "isfs "
            << "./statistics/isfs/backbone/isfs_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"  << " "
            
            << "${waveindex}"   << " "
            << "${waveindex}"   << " "
            << get_geometry()   << " "
            << "${maxLenScale}" << " "
            << get_fullblock()
            << "\n"
            << "list backbone"
            << "\n";
            
            ai << "\n";
        }
    }
    /* Full Intermediate Scattering Function */
    if (is_isf) {
        
        if (is_mbodies) {
            
            ai
            << "isf_list "
            << "./statistics/isf/all/isf_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            
            << "mbodies_list"
            
            << " "
            << get_geometry()
            << " "
            << "${waveindex}"
            << " "
            << "${waveindex}"
            << " "
            << "\n";
            
            ai << "\n";
        }
        else {
            
            ai
            << "isf_list "
            << "./statistics/isf/all/isf_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            
            << "all"
            << " "
            << get_geometry()
            << " "
            << "${waveindex}"
            << " "
            << "${waveindex}"
            << " "
            << "\n";
            
            ai
            << "isf_list "
            << "./statistics/isf/backbone/isf_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".backbone.dat"
            << " "
            
            << "backbone"
            << " "
            << get_geometry()
            << " "
            << "${waveindex}"
            << " "
            << "${waveindex}"
            << " "
            << "\n";
            
            ai << "\n";
        }
    }
    
    /* string analysis */
    if (is_strings) {
        
        if (true) { // use "value_list" to find mobile particles
            
            ai
            /* create list analysis file (WITH .dat appended) */
            << "displacement_list "
            << "./statistics/strings/all/displacement_stats_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            /* create list obj. (WITHOUT .dat appended) */
            << "displacement_list_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << " "
            << "${frame}" // timegap_index
            << "\n";
            
            /* specify list obj. */
            if (is_mbodies) {
                ai
                << "list mbodies_list"
                << "\n\n";
            }
            else {
                ai
                << "list all"
                << "\n\n";
            }
            
            ai
            << "value_list "
            << "write_pdb"
            << " "
            /* use list object */
            << "displacement_list_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << " "
            /* create pdb filename stem */
            << "./statistics/strings/all/displacement_pdb_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << " "
            << "0 0"
            << "\n\n";
            
            if (threshold_keyword=="between") { // "between"
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi << "_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile_lo
                << " "
                << threshold_percentile_hi
                << "\n\n";
            }
            else { // "greater" or "less"
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile
                << "\n\n";
            }
        }
        
        
        if (false) { // use "compare_gaussian" & "find_fast" for mobile particles
            
            ai
            << "vhs "
            << "./statistics/strings/all/vhs_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat"
            << " "
            << vhs_rangebinned
            << " "
            << vhs_nbins
            << "\n"
            
            << "all" << "\n";
            
            ai << "\n";
            
            ai
            << "compare_gaussian ./statistics/strings/all/gaussian_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << " "
            << "${frame}" << "\n"
            
            << "all" << "\n";
            
            ai << "\n";
            
            ai
            << "find_fast "
            << "fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << " "
            << "./statistics/strings/all/fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << ".all.dat" << "\n"
            
            << "all" << "\n";
            
        }
        
        ai << "\n";
        
        ai
        << "strings "
        << "./statistics/strings/all/strings_"
        << sysVar.get_usic()
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << ".all.dat"
        << " "
        << "${frame}"
        << " "
        << strings_threshold
        << " "
        << "./statistics/strings/sigma_matrix.txt"
        << " "
        << "stringlist_"
        << sysVar.get_usic()
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        
        << "\n";
        
        if (true) { // based on "value_list"
            
            if (threshold_keyword=="between") { // "between"
                
                ai
                << "list "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi << "_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << "\n\n";
            }
            else { // "greater" or "less"
                
                ai
                << "list "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << sysVar.get_usic()
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << "\n\n";
            }
            
        }
        else { // based on "find_fast"
            
            ai
            << "list "
            << "fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "\n\n";
        }
        
    }
    
    ai.close();
    
}





void AmdatAnalysis::make_amdatInputFile_strings_loop(const StructureClass& sysVar,
                                                     const WorkScripts& ws,
                                                     const int blockframe)
{
    string i;
    i.append(return_SimulationFolderPath(sysVar));
    i.append("/production/trajectory");
    i.append("/trajectory_");
    i.append(sysVar.get_usic());
    i.append("_${trial}");
    i.append("_${namestring}");
    i.append("_T${Temp}");
    i.append(".prd.custom");
    //cout << i << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_inputs");
    o.append("/amdat_strings_"+to_string((long long int)blockframe));
    o.append(".inp");
    ofstream ai(o.c_str());
    
    /* System block of AMDAT here needs to be manually hard coded */
    //==========================================================================
    
    ai.close();
}





std::vector<int> AmdatAnalysis::make_amdatInputFile_strings_loop(const StructureClass& sysVar,
                                                                 const WorkScripts& ws,
                                                                 const int n_trl,
                                                                 const int n_sys,
                                                                 const double Temp_d,
                                                                 const int peakframe,
                                                                 const int blockframe)
{
    string i;
    i.append(return_SimulationFolderPath(sysVar));
    i.append("/production/trajectory");
    i.append("/trajectory_");
    i.append(sysVar.get_usic());
    i.append("_${trial}");
    i.append("_${namestring}");
    i.append("_T${Temp}");
    i.append(".prd.custom");
    //cout << i << "\n";
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_inputs");
    o.append("/amdat_strings_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".inp");
    ofstream ai(o.c_str());
    
    int pre_frames=10; // # of frames before NGP peak
    int pos_frames=10; // # of frames after "
    
    int beg_frame=peakframe-pre_frames;
    int end_frame=peakframe+(pos_frames-1);
    
    if (beg_frame<1) { // should at least start from 1st frame
        beg_frame=1;
    }
    if (end_frame>=(blockframe-2)) { // avoid last frames
        end_frame=blockframe-2;
    }
    
    /* System Block Definition */
    systemBlock_definition(ai,i);
    
    ai
    << "## peakframe  " << peakframe  << "\n"
    << "## blockframe " << blockframe << "\n"
    << "## beg_frame  " << beg_frame  << "\n"
    << "## end_frame  " << end_frame  << "\n";
    
    ai << "\n\n";
    
    
    bool is_use_findfast=false;
    
    
    for (int frame_index=beg_frame; frame_index<=end_frame; ++frame_index) {
        
        if(is_use_findfast){
            
            // vhs
            //------------------------------------------------------------------
            ai
            << "vhs "
            << "./statistics/strings/all/vhs_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_frame"
            << frame_index
            << ".all.dat"
            << " "
            << vhs_rangebinned
            << " "
            << vhs_nbins
            << "\n"
            
            << "list all" << "\n";
            //------------------------------------------------------------------
            
            ai << "\n";
            
            // compare_gaussian (need frame)
            //------------------------------------------------------------------
            ai
            << "compare_gaussian ./statistics/strings/all/gaussian_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_frame"
            << frame_index
            << ".all.dat"
            << " "
            << frame_index
            << "\n"
            
            << "list all" << "\n";
            //------------------------------------------------------------------
            
            ai << "\n";
            
            // find_fast
            //------------------------------------------------------------------
            ai
            << "find_fast "
            << "fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_frame"
            << frame_index
            << " "
            << "./statistics/strings/all/fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_frame"
            << frame_index
            << ".all.dat"
            << "\n"
            
            << "list all" << "\n";
            //------------------------------------------------------------------
            
            ai << "\n";
            
        }
        
        // strings (need frame)
        //----------------------------------------------------------------------
        ai
        << "strings "
        << "./statistics/strings/all/strings_"
        << sysVar.get_usic()
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << "_frame"
        << frame_index
        << ".all.dat"
        << " "
        << frame_index
        << " "
        << strings_threshold
        << " "
        << "./statistics/strings/sigma_matrix.txt"
        << " "
        << "./statistics/strings/all/stringlist/stringlist_"
        << sysVar.get_usic()
        << "_${trial}"
        << "_${namestring}"
        << "_T${Temp}"
        << "_frame"
        << frame_index
        << ".all.dat"
        << "\n";
        
        if(1){
            ai
            << "list "
            << "fast_"
            << sysVar.get_usic()
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_frame"
            << frame_index;
        }
        else {
            ai
            << "list all";
        }
        ai << "\n";
        //----------------------------------------------------------------------
        
        ai << "\n\n";
    }
    
    ai.close();
    
    return {beg_frame,end_frame};
}





void AmdatAnalysis::make_amdatSubFile(const StructureClass& sysVar,
                                      const WorkScripts& ws,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d,
                                      const int frame)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_submission_files");
    o.append("/AMDAT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".qsub");
    //cout << o << "\n";
    
    const int n_poly_d  = sysVar.get_n_polyVec()[n_sys];
    const int backLen_d = sysVar.get_backLenVec()[n_sys];
    const int sideLen_d = sysVar.get_sideLenVec()[n_sys];
    
    ofstream as0(o.c_str());
    if(as0.is_open())	{
        
        as0
        << "#!/bin/bash"                                                << "\n"
        << "#$ -V"                                                      << "\n"
        << "#$ -cwd"                                                    << "\n"
        << "#$ -j y"                                                    << "\n"
        << "#$ -pe orte " << get_amdat_cores()                          << "\n"
        << "#$ -p " << get_amdat_priority()                             << "\n";
        
        streamsize ss=as0.precision();
        as0 << fixed << setprecision(0);
        
        as0
        << "#$ -hold_jid prd_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        as0
        << "#$ -N amdat_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << "\n";
        
        as0
        << "#$ -o ./AMDAT_submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".amdat.o"
        << "\n";
        
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        ////////////////////////////////////////////////////////////////////////
        //======================================================================
        // Nodes Allocation
        //======================================================================
        /*
         For compute nodes alloccation: within a node, jobs are filled on a half
         node with the same GPU, and then to the next half with the other GPU;
         the same procedure continues in the next requested node. */
        
        int numNode0 = 0; // number of nodes on node0
        int numNode1 = 0; // number of nodes on node1
        
        /*
         using both node0 and node1 */
        if (ws.get_is_node0() && ws.get_is_node1()) {
            numNode0 = (int)(ws.get_node0().size());
            numNode1 = (int)(ws.get_node1().size());
            
        }
        /*
         only using node0 */
        else if (ws.get_is_node0()) {
            numNode0 = (int)(ws.get_node0().size());
            numNode1 = 0;
        }
        /*
         only using node1 */
        else if (ws.get_is_node1()) {
            numNode0 = 0;
            numNode1 = (int)(ws.get_node1().size());
        }
        
        int sumNode = numNode0 + numNode1;
        
        static int counter=0;
        
        /*
         nonimal number of job allocation on a node;
         when this nominal number increases, it means change to the next node
         because the current node is filled completely */
        int noml=counter/(ws.get_CoresPerNode()/get_amdat_cores());
        /*
         index of the node */
        int tmpNode=noml%sumNode;
        int tmpNode1=0;
        /*
         using both node0 and node1 */
        if (ws.get_is_node0() && ws.get_is_node1()) {
            if (tmpNode<numNode0) {
                as0
                << "#$ -q all.q@compute-0-"
                << ws.get_node0()[tmpNode] << ".local"
                << "\n";
            }
            else {
                tmpNode1=(tmpNode-numNode0);
                as0
                << "#$ -q all.q@compute-1-"
                << ws.get_node1()[tmpNode1] << ".local"
                << "\n";
            }
        }
        /*
         only using node0 */
        else if (ws.get_is_node0()) {
            as0
            << "#$ -q all.q@compute-0-"
            << ws.get_node0()[tmpNode] << ".local"
            << "\n";
        }
        /*
         only using node1 */
        else if (ws.get_is_node1()) {
            as0
            << "#$ -q all.q@compute-1-"
            << ws.get_node1()[tmpNode] << ".local"
            << "\n";
        }
        //======================================================================
        // GPU Allocation
        //
        // NOTE: this section is copied from workscripts.cpp
        // AMDAT is usually not utilizing GPUs
        //
        //======================================================================
        /*
         gpuNum is switched between 0 and 1 */
        int gpuNum=0;
        /*
         number of half node cores */
        int hNode=ws.get_CoresPerNode()/2;
        /*
         individual jobs on a node */
        int indv=counter%(ws.get_CoresPerNode()/get_amdat_cores());
        /*
         assign GPU to individual jobs */
        if(sysVar.get_is_GPU()) {
            if (indv<(hNode/get_amdat_cores())) gpuNum=0; else gpuNum=1;
        }
        /*
         Mycroft compute-0-16's GPU1 is broken */
        if((ws.get_is_node0())||
           (ws.get_is_node0() && ws.get_is_node1()))
            if(ws.get_node0()[tmpNode]==16) gpuNum=0;
        
        ++counter;
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        as0 << "\n";
        
        
        /* AMDAT Executable */
        //----------------------------------------------------------------------
        as0
        << "mpirun -np "
        << get_amdat_cores() << " "
        << get_amdat_exe();
        //----------------------------------------------------------------------
        
        as0 << fixed << setprecision(0);
        
        /* AMDAT Inputfile */
        //----------------------------------------------------------------------
        if (is_strings) {
            
            if (is_peak_frame) {
                as0
                << " -i ./AMDAT_inputs/amdat_strings.inp";
            }
            else {
                as0
                << " -i ./AMDAT_inputs/amdat_strings_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".inp";
            }
        }
        else {
            as0
            << " -i ./AMDAT_inputs/amdat.inp";
        }
        //----------------------------------------------------------------------
        
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /* AMDAT temporary file */
        //----------------------------------------------------------------------
        as0
        << fixed << setprecision(0)
        << " -t ./AMDAT_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".tmp";
        //----------------------------------------------------------------------
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        const int numSide=sideLen_d*backLen_d;
        double maxLenScale=0;
        if (get_is_fixed_maxL()){maxLenScale=sysVar.get_maxLenScale()[n_sys];}
        
        
        /* input varaibles */
        //======================================================================
        as0
        << " -c trial 00"       << n_trl
        << " -c namestring "    << sysVar.get_nameString(n_sys)
        << " -c Regime "        << sysVar.get_current_regime()
        << " -c n_prd_blocks "  << ws.get_n_prd_blocks()
        << " -c blocksize "     << ws.get_prd_blocksize()
        << " -c n_poly "        << n_poly_d
        << " -c backLen "       << backLen_d
        << " -c numSide "       << numSide
        << " -c maxLenScale "   << maxLenScale
        << " -c ts "            << ws.get_timestep_size()
        << fixed << setprecision(0)
        << " -c Temp "          << Temp_d;
        
        if (is_strings) {
            if (is_peak_frame) {
                as0
                << " -c frame "     << frame;
            }
        }
        //======================================================================
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /* to find the first waveindex from strFac file */
        if (get_is_isfs() || get_is_isf()) {
            as0
            << " -c waveindex ";
            if (get_is_fixed_q()) {
                as0 << sysVar.get_waveindex()[n_sys];}
            else {
                as0 << find_waveindex(sysVar,n_trl,n_sys,Temp_d);}
        }
        
        
        /* screen file */
        as0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        as0
        << " > ./screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".amdat.screen"
        << "\n";
        //----------------------------------------------------------------------
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        as0.close();
    }
    else cout << "amdatSubFile: 'amdat.sh' cannot open." << "\n";
    
}





int AmdatAnalysis::find_waveindex(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d)
{
    /*
     find wavenumber used overall structure factor information */
    string analysispart="all";
    
    //==============================================
    // READ: The path to the structure factor file
    //==============================================
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/strFac/");
    i.append(analysispart);
    i.append("/strFac_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    //cout << i << "\n";
    
    
    double q_threshold=0; // x-direction cutoff
    
    if (sysVar.get_systemUnit()=="real") {
        q_threshold=1.0;
    }
    else if (sysVar.get_systemUnit()=="lj") {
        q_threshold=5.0;
    }
    
    int waveindex=0;
    
    ifstream readstrFacFile(i.c_str());
    if (readstrFacFile.is_open()) {
        
        string lineContent;
        vector<double> wavevector;
        vector<double> strFac;
        double wavetmp;
        double strFactmp;
        
        /* This is the first trash line in structure factor file */
        getline(readstrFacFile,lineContent);
        
        /*
         1st column: wavevector, k;
         2nd column: structure factor, S(k)( */
        while (getline(readstrFacFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> wavetmp;   wavevector.push_back(wavetmp);
            iss >> strFactmp; strFac.push_back(strFactmp);
            
            if (strFac.size()>=3) {
                
                int midindex=(int)strFac.size()-2; // 0-based index
                double q_mid=wavevector[midindex];
                double forwdiff=strFac[midindex]-strFac[midindex+1];
                double backdiff=strFac[midindex]-strFac[midindex-1];
                
                if ((q_mid>q_threshold) &&
                    (forwdiff>0) && (backdiff>0))
                {
                    waveindex = midindex; // 0-based index
                    //cout << strFac[midindex] << "\n";
                    break;
                }
            }
        }
        readstrFacFile.close();
    }
    else cout << "find_firstwaveindex: 'strFac' file cannot open." << "\n";
    
    return waveindex;
}





int AmdatAnalysis::find_ngp_frame(const StructureClass& sysVar,
                                  const int n_trl,
                                  const int n_sys,
                                  const double Temp_d)
{
    string analysispart="all";
    
    //==============================================
    // READ: The path to the structure factor file
    //==============================================
    string i;
    i.append(return_AnalysisFolderPath(sysVar));
    i.append("/statistics/ngp/");
    i.append(analysispart);
    i.append("/ngp_");
    i.append(sysVar.get_usic());
    i.append("_00"+to_string((long long int)n_trl));
    i.append("_"+sysVar.get_nameString(n_sys));
    i.append("_T"+to_string((long long int)Temp_d));
    i.append("."+analysispart+".dat");
    //cout << i << "\n";
    
    double time_threshold=0; // in-block time cutoff
    
    const int current_regime=sysVar.get_current_regime();
    double previous_teq,current_teq,teq;
    if (current_regime>0) {
        previous_teq = sysVar.get_equilibration_times()[current_regime-1];
        current_teq  = sysVar.get_equilibration_times()[current_regime];
        teq=max(previous_teq,current_teq);
    }
    else {
        teq=sysVar.get_equilibration_times()[current_regime];
    }
    time_threshold=teq/sysVar.get_n_equ_blocks();
    
    
    vector<double> data_decreasing;
    
    int frame=0;
    
    ifstream readFile(i.c_str());
    if (readFile.is_open()) {
        
        string lineContent;
        vector<double> time;
        vector<double> ngp;
        double timetmp;
        double ngptmp;
        
        /* This is the first trash line in file */
        getline(readFile,lineContent);
        
        /*
         1st column: time
         2nd column: ngp
         */
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> timetmp; time.push_back(timetmp);
            iss >> ngptmp;  ngp.push_back(ngptmp);
            
            if (ngp.size()>=3) {
                
                int midindex=(int)ngp.size()-2; // 0-based index
                double mid=ngp[midindex];
                double forwdiff=ngp[midindex]-ngp[midindex+1];
                double backdiff=ngp[midindex]-ngp[midindex-1];
                
                if ((mid>time_threshold) &&
                    (forwdiff>0) && (backdiff>0))
                {
                    frame = midindex; // 0-based index
                    //cout << strFac[midindex] << "\n";
                    break;
                }
            }
        }
        readFile.close();
    }
    else cout << "find_ngp_frame: 'ngp' file cannot open." << "\n";
    
    return frame;
}





void AmdatAnalysis::make_amdatSubScript(const StructureClass& sysVar,
                                        const std::vector<std::vector<double>> tinfo) const
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/amdatSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    
    const int n_trial  = sysVar.get_n_trial();
    const int n_system = sysVar.get_n_system();
    
    
    ofstream as1(o.c_str());
    if (as1.is_open()) {
        
        as1 << "#!/bin/bash" << "\n";
        as1 << "\n";
        
        /* Trial */
        //----------------------------------------------------------------------
        as1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            as1 << i;
            if (i!=(n_trial-1)) as1 << " ";
        }
        as1 << ")" << "\n";
        
        
        /* namestring */
        //----------------------------------------------------------------------
        as1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            as1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) as1 << " ";
        }
        as1 << ")" << "\n";
        
        
        /* Temperature */
        streamsize ss=as1.precision();
        as1 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        as1 << "Temp=(";
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=0; n_sys<n_system; ++n_sys) {
                
                int index=n_trl*n_system+n_sys;
                
                int n_Temp=(int)tinfo[index].size();
                
                as1 << "\"";
                for (int T=0; T<n_Temp; ++T) {
                    as1
                    << tinfo[index][T];
                    if (T!=(n_Temp-1)) as1 << " ";
                }
                if (n_trl!=(n_trial-1)) as1 << "\"" << " ";
            }
        }
        as1
        << "\""
        << ")"
        << "\n\n";
        //----------------------------------------------------------------------
        as1.precision(ss);
        as1 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        as1
        << "cd "
        << return_AnalysisFolderPath(sysVar)                            << "\n";
        
        as1 << "\n";
        
        as1
        << "for i in ${trial[@]}; do"                                   << "\n"
        << "  for ii in ${namestring}; do"                              << "\n"
        << "    for iii in ${Temp[$i]}; do"                             << "\n"
        << "      qsub_File=./AMDAT_submission_files/AMDAT_"
        << sysVar.get_usic()
        << "_00${i}"
        << "_${ii}"
        << "_T${iii}"
        << ".qsub"                                                      << "\n"
        
        << "      test -e ${qsub_File}"                                 << "\n"
        << "      if [ \"$?\" -eq \"0\" ]; then"                        << "\n"
        << "        qsub ${qsub_File}"                                  << "\n"
        << "      else"                                                 << "\n"
        << "        continue"                                           << "\n"
        << "      fi"                                                   << "\n"
        << "    done"                                                   << "\n"
        << "  done"                                                     << "\n"
        << "done"                                                       << "\n";
        
        as1.close();
    }
    else cout << "amdatSubScript: 'amdat.sh' cannot open." << "\n";
    
    
    //==========================================================================
    // Direct submission of all AMDAT jobs
    //==========================================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_AnalysisFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/amdatSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::change_atom_order(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    
    /* get the path to the production trajectory file */
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/trajectory/");
    
    file.append("trajectory_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.custom");
    //cout << i << "\n";
    
    /* input file */
    input=folder+file;
    ifstream read_input(input.c_str());
    
    /* output file */
    output=folder+"newOrder_"+file;
    ofstream write_output(output.c_str());
    
    int n_poly=32;
    int backlen=30;
    int back_per_monomer=1;
    int side_per_monomer=3;
    int monomersize=back_per_monomer+side_per_monomer;
    
    int n_columns_of_data=7; // # of columns of data to be stored
    
    if (read_input.is_open()) {
        
        string lineContent;
        string strData0,strData1,strData2,strData3;
        
        bool is_lineForAtomNum=false;
        //int  intData0=0;
        int  countStepkeyword=0;
        int  n_atoms=0;
        int  line_index=0;
        
        vector<double> vecData; // used for storing line content
        
        vector<vector<double>> trajectoryinfo;
        
        while (getline(read_input,lineContent)) {
            
            ++line_index;
            
            istringstream iss(lineContent);
            
            // write out "type x y z ix iy iz" in custom file
            //--------------------------------------------------------------
            if ((countStepkeyword!=0)&&(countStepkeyword<=n_atoms)) {
                
                for (int i=0; i<n_columns_of_data; ++i) {
                    vecData.push_back(0);
                    iss >> vecData[i];
                }
                trajectoryinfo.push_back(vecData);
                
                /* NOTE: clean container content after each use */
                vecData.clear();
                
                
                // when all trajectories are read into storing vector
                // write out the intended format to file
                if (countStepkeyword==n_atoms) {
                    
                    for (int i=0; i<n_poly; ++i) {
                        
                        int offset_index=i*(backlen*monomersize);
                        
                        for (int ii=offset_index; ii<(offset_index+backlen); ++ii) {
                            
                            /* backbone bead */
                            int back_index=ii;
                            
                            for (int j=0; j<back_per_monomer; ++j) {
                                for (int iii=0; iii<trajectoryinfo[back_index+j].size(); ++iii) {
                                    write_output
                                    << trajectoryinfo[back_index+j][iii] << " ";
                                }
                                write_output << "\n";
                            }
                            
                            /* side group beads */
                            int sum_sides=side_per_monomer*(back_index-offset_index);
                            int side_index=(offset_index+backlen)+sum_sides;
                            
                            for (int j=0; j<side_per_monomer; ++j) {
                                for (int iii=0; iii<trajectoryinfo[side_index+j].size(); ++iii) {
                                    write_output
                                    << trajectoryinfo[side_index+j][iii] << " ";
                                }
                                write_output << "\n";
                            }
                        }
                    }
                    /* NOTE: clean container content after each use */
                    trajectoryinfo.clear();
                }
                
            }
            else {
                
                write_output << lineContent << "\n";
            }
            //--------------------------------------------------------------
            
            iss >> strData0;
            
            if (is_lineForAtomNum) {
                
                n_atoms=stoi(strData0); // use # of atoms defined in file
                //cout << n_atoms << "\n";
                
                is_lineForAtomNum=false;
            }
            
            if (strData0=="ITEM:") {
                
                iss >> strData1;
                
                if (strData1=="NUMBER") {
                    
                    is_lineForAtomNum=true;
                    //cout << line_index << "\n";
                }
                else if (strData1=="ATOMS") {
                    
                    // NOTE:
                    // if the read-in data type (ex. double) is not the same as
                    // the receiving variable type (ex. string), then the
                    // receiving varibale will not store the read-in data;
                    // therefore, the memory of the receiving variable would
                    // remain as the previously stored value where the read-in
                    // data type is the same as the receiving varaible.
                    
                    ++countStepkeyword;
                    //cout << line_index << "\n";
                }
                else if (strData1=="TIMESTEP") {
                    
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
            
        }
        read_input.close();
    }
    else {
        cout
        << "change_atom_order: input file cannot open." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    write_output.close();
}





void AmdatAnalysis::change_productionAtomTypes(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d)
{
    
    /* get the path to the production trajectory file */
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/trajectory/");
    
    file.append("trajectory_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.custom");
    //cout << i << "\n";
    
    /* input file */
    input=folder+file;
    ifstream read_input(input.c_str());
    
    /* output file */
    output=folder+"new_"+file;
    ofstream write_output(output.c_str());
    
    vector<int> backbone=sysVar.get_backboneTypes();
    vector<int> sideGroup=sysVar.get_sideGroupTypes();
    vector<int>::iterator itr0,itr1;
    
    int n_columns_of_data=7; // # of columns of data to be stored
    
    if (read_input.is_open()) {
        
        string lineContent;
        string strData0,strData1,strData2,strData3;
        
        bool is_lineForAtomNum=false;
        int  intData0=0;
        int  countStepkeyword=0;
        int  n_atoms=0;
        int  line_index=0;
        
        vector<double> vecData; // used for storing line content
        
        while (getline(read_input,lineContent)) {
            
            ++line_index;
            
            istringstream iss(lineContent);
            
            // write out "type x y z ix iy iz" in custom file
            //--------------------------------------------------------------
            if ((countStepkeyword!=0)&&(countStepkeyword<=n_atoms)) {
                
                for (int i=0; i<n_columns_of_data; ++i) {
                    vecData.push_back(0);
                    iss >> vecData[i];
                }
                
                intData0=(int)vecData[0]; // this is the atom type
                
                // if found matching backbone type,
                // change to assigned backbone type
                itr0=find(backbone.begin(),backbone.end(),intData0);
                if (itr0!=backbone.end()) {
                    vecData[0]=1;
                }
                
                // if found matching sideGroup type,
                // change to assigned sideGroup type
                itr1=find(sideGroup.begin(),sideGroup.end(),intData0);
                if (itr1!=sideGroup.end()) {
                    vecData[0]=2;
                }
                
                // write rest of content in the same row
                for (int i=0; i<n_columns_of_data; ++i) {
                    write_output << vecData[i] << " ";
                }
                write_output << "\n";
                
                
                // clean container content after each use
                vecData.clear();
            }
            else {
                
                write_output << lineContent << "\n";
            }
            //--------------------------------------------------------------
            
            
            iss >> strData0;
            
            if (is_lineForAtomNum) {
                
                n_atoms=stoi(strData0); // use # of atoms defined in file
                //cout << n_atoms << "\n";
                
                is_lineForAtomNum=false;
            }
            
            if (strData0=="ITEM:") {
                
                iss >> strData1;
                
                if (strData1=="NUMBER") {
                    
                    is_lineForAtomNum=true;
                    //cout << line_index << "\n";
                }
                else if (strData1=="ATOMS") {
                    
                    ++countStepkeyword;
                    //cout << line_index << "\n";
                }
                else if (strData1=="TIMESTEP") {
                    
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
            
            
        }
        read_input.close();
    }
    else {
        cout
        << "change_productionAtomTypes: input file cannot open." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    write_output.close();
}





void AmdatAnalysis::change_prdcustomfilesteps(const StructureClass& sysVar,
                                              const int n_trl,
                                              const int n_sys,
                                              const double Temp_d)
{
    /* get the path to the production trajectory file */
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/trajectory/");
    
    file.append("trajectory_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.custom");
    //cout << i << "\n";
    
    /* input file */
    input=folder+file;
    ifstream read_input(input.c_str());
    
    /* output file */
    output=folder+"new_"+file;
    ofstream write_output(output.c_str());
    
    if (write_output.is_open()) {
        
        if (read_input.is_open()) {
            
            string lineContent;
            string strData0,strData1,strData2,strData3;
            
            int steps_zero=0;
            int steps_corr=0;
            int countStepkeyword=0;
            int line_index=0;
            int block_index=0;
            int frame_index=0;
            int countaccess=0;
            
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                
                // NOTE:
                // the very first frame is at timestep=0;
                // in the using exponential timestepping scheme,
                // 0th block has N+1 frames, 1st-final blocks
                // each has N frames.
                if ((countStepkeyword==0)&&(countaccess>1)) {
                    
                    if (block_index>0) {
                        iss >> strData0;
                        steps_zero=stoi(strData0);
                        steps_corr=steps_zero+block_index*(int)pow(prd_exp_base,prd_blocksize-1);
                        write_output << steps_corr << "\n";
                    }
                    else {
                        write_output << lineContent << "\n";
                    }
                    
                    ++frame_index;
                    if (frame_index==prd_blocksize) {
                        ++block_index;
                        frame_index=0;
                    }
                }
                else {
                    write_output << lineContent << "\n";
                }
                
                ++countStepkeyword;
                
                iss >> strData0;
                
                if (strData0=="ITEM:") {
                    
                    iss >> strData1;
                    
                    if (strData1=="TIMESTEP") {
                        ++countaccess;
                        countStepkeyword=0;
                        //cout << line_index << "\n";
                    }
                }
            }
            read_input.close();
        }
        else {
            cout
            << "change_prdcustomfilesteps: input file cannot open." << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
        write_output.close();
    }
    else {
        cout
        << "change_prdcustomfilesteps: output file cannot open." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::delete_zerothFrames(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys,
                                        const double Temp_d)
{
    
    /* get the path to the production trajectory file */
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/trajectory/");
    
    file.append("trajectory_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.custom");
    //cout << i << "\n";
    
    /* input file */
    input=folder+file;
    ifstream read_input(input.c_str());
    
    /* output file */
    output=folder+"new_"+file;
    ofstream write_output(output.c_str());
    
    if (write_output.is_open()) {
        
        if (read_input.is_open()) {
            
            string lineContent,lineContent1;
            string strData0,strData1,strData2,strData3;
            
            bool is_readinframe=true;
            int  line_index=0;
            int  countaccess=0;
            
            vector<double> vecData; // used for storing line content
            
            // Read every line one at a time of the file
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                
                iss >> strData0;
                
                if (strData0=="ITEM:") {
                    
                    iss >> strData1;
                    
                    if (strData1=="TIMESTEP") {
                        
                        ++countaccess;
                        //cout << line_index << "\n";
                        getline(read_input,lineContent1);
                        if ((stoi(lineContent1)==0)&&(countaccess>1)) is_readinframe=false;
                        else is_readinframe=true;
                        
                        if (is_readinframe) {
                            write_output
                            << lineContent  << "\n"
                            << lineContent1 << "\n";
                        }
                    }
                    else {
                        if (is_readinframe) write_output << lineContent << "\n";
                    }
                }
                else {
                    if (is_readinframe) write_output << lineContent << "\n";
                }
                
                
            }
            read_input.close();
        }
        else {
            cout
            << "delete_zerothFrames: input file cannot open." << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
        write_output.close();
    }
    else {
        cout
        << "delete_zerothFrames: output file cannot open." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    // rm and rename file
    //string rm="rm "+input;
    //system(rm.c_str());
    //string mv="mv "+output+" "+input;
    //system(mv.c_str());
}





const string AmdatAnalysis::amdat_target(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys,
                                         const double Temp_d,
                                         const string& phase) const
{
    string target;
    target.append(return_AnalysisFolderPath(sysVar));
    target.append("/statistics");
    target.append("/"+phase);
    target.append("/"+get_analysispart());
    target.append("/"+phase+"_");
    target.append(sysVar.get_usic());
    target.append("_00"+to_string((long long int)n_trl));
    target.append("_"+sysVar.get_nameString(n_sys));
    target.append("_T"+to_string((long long int)Temp_d));
    target.append("."+get_analysispart()+".dat");
    //cout << target << "\n";
    
    return target;
}

const string AmdatAnalysis::amdat_target(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys,
                                         const double Temp_d,
                                         const string& phase,
                                         const int frame) const
{
    string target;
    target.append(return_AnalysisFolderPath(sysVar));
    target.append("/statistics");
    target.append("/"+phase);
    target.append("/"+get_analysispart());
    target.append("/"+phase+"_");
    target.append(sysVar.get_usic());
    target.append("_00"+to_string((long long int)n_trl));
    target.append("_"+sysVar.get_nameString(n_sys));
    target.append("_T"+to_string((long long int)Temp_d));
    target.append("_frame"+to_string((long long int)frame));
    target.append("."+get_analysispart()+".dat");
    //cout << target << "\n";
    
    return target;
}




/*==( public setters )==*/

/* bool */
void AmdatAnalysis::set_is_changeAtomOrder(const bool b){is_changeAtomOrder=b;}
void AmdatAnalysis::set_is_changeAtomType(const bool b){is_changeAtomType=b;}
void AmdatAnalysis::set_is_changeFrameSteps(const bool b){is_changeFrameSteps=b;}
void AmdatAnalysis::set_is_deleteZerothFrames(const bool b){is_deleteZerothFrames=b;}
void AmdatAnalysis::set_is_NPT(const bool b){is_NPT=b;}
void AmdatAnalysis::set_is_fixed_q(const bool b){is_fixed_q=b;}
void AmdatAnalysis::set_is_fixed_maxL(const bool b){is_fixed_maxL=b;}
void AmdatAnalysis::set_is_strFac(const bool b){is_strFac=b;}
void AmdatAnalysis::set_is_msd(const bool b){is_msd=b;}
void AmdatAnalysis::set_is_ngp(const bool b){is_ngp=b;}
void AmdatAnalysis::set_is_isfs(const bool b){is_isfs=b;}
void AmdatAnalysis::set_is_composition(const bool b){is_composition=b;}
void AmdatAnalysis::set_is_u2dist(const bool b){is_u2dist=b;}
void AmdatAnalysis::set_is_stiffness_dist(const bool b){is_stiffness_dist=b;}
void AmdatAnalysis::set_is_isf(const bool b){is_isf=b;}
void AmdatAnalysis::set_is_mbodies(const bool b){is_mbodies=b;}
void AmdatAnalysis::set_is_strings(const bool b){is_strings=b;}
void AmdatAnalysis::set_is_peak_frame(const bool b){is_peak_frame=b;}

/* int */
void AmdatAnalysis::set_amdat_cores(const int i){amdat_cores=i;}
void AmdatAnalysis::set_amdat_priority(const int i){amdat_priority=i;}
void AmdatAnalysis::set_fullblock(const int i){fullblock=i;}

/* string */
void AmdatAnalysis::set_amdat_exe(const string& s){amdat_exe=s;}
void AmdatAnalysis::set_symmetry(const std::string& str){symmetry=str;}
void AmdatAnalysis::set_geometry(const std::string& str){geometry=str;}
void AmdatAnalysis::set_analysispart(const std::string& str){analysispart=str;}
void AmdatAnalysis::set_centertypename(const std::string& str){centertypename=str;}

