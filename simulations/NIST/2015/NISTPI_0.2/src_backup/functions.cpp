//
//  Created by Sean Jh H. on 7/9/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "functions.h"
#include "structureclass.h"
#include "lmpscripts.h"
#include "workscripts.h"
#include "amdatanalysis.h"
#include "fitdata.h"
#include "fitdielectric.h"
#include <algorithm>

#include "theorytest.h"

using namespace std;

const string autoWork::submit_jobs()
{
    time_t begin_overall=time(NULL);
    string startdatetime=return_datetime();
    
    /*
     Supported Backbones
     *****************************************************
     linear  = 0
     typeB1  = 1  Methyl Acrylate, (C4H6O2)n
     typeB2  = 2  Methyl MethAcrylate, (C5H8O2)n
     typeB3  = 3  Methyl Acrylamide, (C4H7NO)n
     typeB4  = 4  1-ethyl-4-vinylbenzene, (C10H12)n
     LJliqB  = 5  LJ liquid
     *****************************************************
     
     Supported Side Groups
     *****************************************************
     phenyl  = 0  C6H5-
     alkyl   = 1  CnH2n+1-
     tButyl  = 2  tert-butyl, (CH3)3C–
     clProp  = 3  1-chloropropane, C3H7Cl-
     C3H8O   = 4  methoxyethane, C3H8O-
     stdFENE = 5  linear stnadard FENE chain
     LJliqS  = 6  LJ liquid
     *****************************************************
     
     Example Combination
     *****************************************************
     (5,6): Lennard-Jones particles
     (0,5): linear bead-spring FENE chains
     (0,0): linear+phenyl = polystyrene
     (2,1): typeB2+alkyl  = PxMA (x=length of side alkyls)
     *****************************************************/
    
    /*====================================================
     Declare an object of 'StructureClass' class
     Constructor arguments: (typeB,typeS);
     This object will be used as the simulation information
     carrier throughout the automation process
     ====================================================*/
    StructureClass species(0,5);
    
    
    //==========================================================================
    // Master Directory Path
    //==========================================================================
    //species.set_Path("/home/jh148/AUTOWORK");
    species.set_Path("/Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork");
    
    //==========================================================================
    // USIC = simType + (year) + usicID
    //==========================================================================
    species.set_simType("NIST");
    species.set_year("2015");
    species.set_usicID("PI_0.2");
    
    //==========================================================================
    // Number of Systems and Trials
    //==========================================================================
    species.set_n_trial(4);
    species.set_n_sys_beg(0); // inclusive
    species.set_n_sys_end(0); // inclusive
    
    //==========================================================================
    // System Unit and Starting Temperatures
    //==========================================================================
    species.set_systemUnit("real");  // "real" or "metal" or "lj"
    species.set_startTemp(500.0);    // highest starting T
    species.set_crossTemp(300.0);    // crossing T
    species.set_hTres(50.0);         // resolution above crossTemp
    species.set_lTres(25.0);         // resolution below crossTemp
    
    // NOTE:
    // The starting temperature needs to be "guaranteed high temperature"
    // i.e. well above TA and close to T0
    
    
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    // Equilibration Cutoff Times
    //==========================================================================
    const double n_equ_blocks=10;
    species.set_n_equ_blocks(n_equ_blocks);
    
    const double n_prd_blocks=10;
    species.set_n_prd_blocks(n_prd_blocks);
    
    const double r=pow(10.0,0.5);
    species.set_equilibration_times
    ({
        n_equ_blocks*  1e+3,
        n_equ_blocks*r*1e+3,
        n_equ_blocks*  1e+4,
        n_equ_blocks*  1e+5,
        n_equ_blocks*  1e+6,
        n_equ_blocks*  1e+7,
        n_equ_blocks*  1e+8,
        n_equ_blocks*  1e+9
    });
    // NOTE:
    // - First 3 regimes are suggested to cover just 1 decade in high T regime
    // - Overall number of regimes should > 3 regimes to have good quality fit
    // - in real units: time unit is femto second
    // - in  lj  units: time unit is tau (~pico second in real unit)
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    /* set dynamic range covered by simulation */
    //--------------------------------------------------------------------------
    const bool   is_useDynamicRange=false;
    const double dynamicRange_on_tau0=4; // NOTE: in log scale
    //--------------------------------------------------------------------------
    
    
    /* set number of T's used in each regime */
    //--------------------------------------------------------------------------
    // NOTE:
    // number of T's should >= 2 for extrapolation to be valid
    // also, number of element should = number of regimes
    //--------------------------------------------------------------------------
    species.set_n_regime_temps
    ({8,8,8,8,8,8,4,4});
    
    
    /* set timestep size used in each regime */
    //--------------------------------------------------------------------------
    // NOTE:
    // 1. number of element should = number of regimes
    // 2. should use correct time unit
    //    real  fs
    //    metal ps
    //    lj    tau
    //--------------------------------------------------------------------------
    species.set_ts_regime
    ({1,1,1,1,1,1,1,1});
    
    
    /* set [start/end] regime */
    //--------------------------------------------------------------------------
    // NOTE:
    // starting and ending regime indices are both inclusive
    //--------------------------------------------------------------------------
    species.set_regime_beg(0); // inclusive
    species.set_regime_end(7); // inclusive
    
    
    //--------------------------------------------------------------------------
    // Module Switches
    //--------------------------------------------------------------------------
    species.set_is_makeFileSystem(1);  // File System Generation
    species.set_is_DataPhase(0);       // Data File Generation
    species.set_is_LammpsPhase(0);     // Simulation Scripts Generation
    species.set_is_Simulations(0);     // Bash Scripts Generation
    species.set_is_AMDAT(0);           // AMDAT Scripts Generation
    species.set_is_fitData(0);         // Curve-fitting toolkit with ALGLIB
    species.set_is_fit_dielectric(1);  // Fit Experimetal Dielectric Data
    species.set_is_theorytest(0);      // Test the Established Theories of G.F.
    
    //--------------------------------------------------------------------------
    // Direct Submission Switch
    //--------------------------------------------------------------------------
    species.set_is_directSub(0);       // Automatic Job Submissions
    
    //--------------------------------------------------------------------------
    // Watch Functionality Switch
    //--------------------------------------------------------------------------
    species.set_is_watch_hold(0);      // Watch for Target Files
    
    //--------------------------------------------------------------------------
    // other useful switches
    //--------------------------------------------------------------------------
    //species.set_is_GPU(false);  // Enable GPU Features
    //species.set_is_fullquench(true); // Straight or Stepwise Quench
    //species.set_is_singleTempTest(true); // singleT run: gen --> prd --> amdat
    
    
    
    
    
    /* System Initialiation */
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    int n_digits=3;
    // NOTE: precision
    // real units: temperature precise to "n_digits" below decimal point
    // (ex. n_digits=3 --> 1234.567; n_digits=0 --> 1234)
    // lj units: temperature precise to "n_digits+3" below decimal point
    // (ex. n_digits=3 --> 1.234567; n_digits=0 --> 1.234)
    double precision=pow(10,n_digits);
    species.set_n_digits(n_digits);
    species.set_precision(precision);
    
    /* define system unit related proprties */
    if (species.get_systemUnit()=="real")
    {
        /* time unit: fs */
        species.set_corrFacforT(precision);
        species.set_extrp_time(1e+17);       // 10^17 fs = 100 seconds
        species.set_is_use_viscosity(false); // true: using viscosity data
    }
    else if (species.get_systemUnit()=="metal")
    {
        /* time unit: ps */
        species.set_corrFacforT(precision);
        species.set_extrp_time(1e+14);       // 10^14 ps = 100 seconds
        species.set_is_use_viscosity(false); // true: using viscosity data
    }
    else if (species.get_systemUnit()=="lj")
    {
        /* time unit: tau ~ ps */
        species.set_corrFacforT(precision*1000.0);
        species.set_extrp_time(1e+14);       // 10^14 tau ~ 100 seconds
        species.set_is_use_viscosity(false);
    }
    
    
    /* store a cpoy of original equilibration time */
    species.set_equilibration_times_org(species.get_equilibration_times());
    
    
    /* dynamics range in simulation */
    species.set_is_useDynamicRange(is_useDynamicRange);
    if (is_useDynamicRange)
    {species.set_dynamicRange_tau0(dynamicRange_on_tau0);}
    
    
    /* the starting time */
    species.set_startdatetime(startdatetime);
    
    
    /* setup force fields */
    species.set_structureParam();
    
    
    /* n_trial */
    const int n_trial=species.get_n_trial();
    
    
    /* [start/end] regimes */
    const int regime_beg=species.get_regime_beg();
    const int regime_end=species.get_regime_end();
    
    
    /* n_regime */
    const int n_regime=(int)species.get_equilibration_times().size();
    species.set_n_regime(n_regime);
    const int n_sys_beg=species.get_n_sys_beg();
    const int n_sys_end=species.get_n_sys_end();
    
    
    /* Temperature containers */
    //--------------------------------------------------
    // linearTemperatures function returns:
    // the 1st element [0] is the overall T's;
    // the 2nd element [1] is the maually set hiT's
    // the 3rd element [2] is the maually set loT's
    //--------------------------------------------------
    species.set_n_Temp(species.get_n_regime_temps()[0]);
    vector<double> tinfo=linearTemperatures(species)[0];
    for (int i=0; i<tinfo.size(); ++i) {tinfo[i] *= precision;}
    vector<vector<double>> initialTs;
    vector<vector<double>> tinfo2D={tinfo,tinfo};
    
    /* set the cutoff T for Arrhenius fit */
    //==================================================
    double cutT=0;
    if (true)
    {
        /* lowest T in first regime */
        size_t s=tinfo.size();
        cutT=floor(tinfo[s-1]); // already expanded
    }
    else
    {
        /* arbitrary T based on experience */
        cutT=1000.0;
        cutT=floor(cutT*precision); // use expanded form
    }
    species.set_cutTforArrhenius(cutT/species.get_corrFacforT());
    //==================================================
    
    
    //--------------------------------------------------
    // containers below have structure:
    // row:    trials and systems
    // column: elemnts per trial ans system
    //--------------------------------------------------
    for (int n_trl=0; n_trl<n_trial; ++n_trl)
    {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
        {
            initialTs.push_back(tinfo);
            species.get_temperaturesInfo().push_back(tinfo); // 2D double
            species.get_equilibratedTs().push_back(tinfo);   // 2D double
            species.get_quenchingTs().push_back(tinfo);      // 2D double
            species.get_qchTdiff().push_back(0);             // 1D double
            species.get_tauFit().push_back(tinfo2D);         // 3D double
            species.get_tauEqu().push_back(tinfo2D);         // 3D double
            
        }
    }
    species.set_initialtemps(initialTs);
    
    
    //--------------------------------------------------
    // containers below have structure:
    // row:    regime
    // column: trials and systems
    //--------------------------------------------------
    vector<double> tmpvd;
    for (int i=0; i<n_trial; ++i)
    {
        for (int ii=n_sys_beg; ii<=n_sys_end; ++ii)
        {
            tmpvd.push_back(0.0);
        }
    }
    for (int i=0; i<n_regime; ++i)
    {
        /* container has current regime data of all trials of systems */
        species.get_Tg_extrp().push_back(tmpvd);
        species.get_Tg_compu().push_back(tmpvd);
        species.get_m_extrp().push_back(tmpvd);
        species.get_m_compu().push_back(tmpvd);
        
        /* container has elapsed time of current regime */
        species.get_elapsedtime().push_back(0);
    }
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                       I. Data Files Generation                         //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    string str;
    
    if (species.get_is_makeFileSystem())
    {
        /* Construct File System */
        //-------------------------------------------
        str=make_SimulationsFolders(species);
        call_system_bash(str);
        //-------------------------------------------
        
    }
    if (species.get_is_DataPhase())
    {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) // species
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl) // trial
            {
                /* Construct Single Chain */
                //----------------------------------
                species.make_SingleChainXYZ(n_sys);
                species.make_BondFile(n_sys);
                species.make_AngleFile(n_sys);
                species.write_ForcefieldToFile();
                //----------------------------------
                
                /* Invoke PACKMOL */
                species.invoke_PACKMOL(n_sys);
                
                /* Make LAMMPS Data Files */
                species.make_LammpsDataFile(n_sys);
                
                /* Change File Names (per trial) */
                str=move_InnerLoopFiles(species,n_sys,n_trl);
                call_system_bash(str);
                
                /* Change File Names (per system) */
                str=move_OuterLoopFiles(species,n_sys);
                call_system_bash(str);
            }
        }
        cout
        << "\n\n"
        << "Finished Data File Generation!"
        << "\n\n";
        
        /* print single chain info */
        //species.echoInformation();
    }
    else
    {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) // species
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl) // trial
            {
                //--------------------------------------------
                // if data file generation is skipped,
                // single chain info is still neccesary;
                // ex. it will be used in AMDAT input script.
                //--------------------------------------------
                species.make_SingleChainXYZ(n_sys);
                species.make_BondFile(n_sys);
                species.make_AngleFile(n_sys);
                species.write_ForcefieldToFile();
                //--------------------------------------------
            }
        }
    }
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                    II. LAMMPS InputFiles Generation                    //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    LmpScripts lmp(species);
    
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    // 'LmpScripts' class setters
    //==========================================================================
    lmp.set_is_tampering(false); // true: use parallel tampering after equil.
    lmp.set_is_respa(false);  // true: RESPA; false: Verlet
    lmp.set_is_resize(false); // resize after equil. and fix box size afterwards
    lmp.set_is_fixMomentum(true);    // fixp in generation thru equilibration
    lmp.set_is_fixPrdMomentum(true); // fixp in production
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    if (species.get_is_LammpsPhase())
    {
        /* Generation */
        lmp.make_GenerationInputFile(species);
        
        /* Quench */
        lmp.make_QuenchInputFile(species);
        
        /* Equilibration */
        lmp.make_EquilibrationInputFile(species);
        
        /* Parallel Tempering */
        if (lmp.get_is_tampering()) {
            lmp.make_TamperingInputFile(species);}
        
        /* Resizing */
        if (lmp.get_is_resize()) {
            lmp.make_ResizeInputFile(species);}
        
        /* Production */
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
        {
            lmp.make_ProductionInputFile(species,n_sys);
        }
    }
    /* Clean up residual temporary files if found */
    clean_tmpfiles(species);
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                   Fit Experimental Dielectric Data                     //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    FitDielectric ds(species);
    
    ds.set_is_use_dc_free(false);
    ds.set_is_manualFitSingleDS(true);
    
    /*--------------------------------------------------------------------------
     form -- Format of Dielectric Data:
     1. NPIC
     2. NIST
     
     func -- Fitting Function:
     1. HN_loss (Loss Part of Havriliak-Negami)-- shape params are free
     2. CC_loss (Loss Part of Cole-Cole)       -- beta=1  (symmetric broadening)
     3. CD_loss (Loss Part of Cole-Davidson)   -- alpha=1 (asymmetric '')
     
     model -- Model describing T Dependence of Relaxation Time:
     1. VFT
     2. COOP
     -------------------------------------------------------------------------*/
    ds.set_form("NIST");
    ds.set_func("HN_loss");
    ds.set_model("COOP");
    
    if (species.get_is_fit_dielectric())
    {
        time_t startds=time(NULL);
        
        vector<vector<double>> tinfo_dielectric;
        
        /* make folder for and transfer dielectric data */
        ds.make_dielectric_folder(species);
        
        /* intialized internal values */
        ds.initialize();
        
        /* data pre-processing */
        tinfo_dielectric=ds.dielectricData_preprocess(species);
        
        if (!ds.get_is_manualFitSingleDS()) {
            
            /* fit relaxation profile to the chosen func */
            for (int process=0; process<tinfo_dielectric.size(); ++process) {
                // NOTE:
                // if there are both cooling and heating processes:
                // process=0: container has cooling temperatures
                // process=1: container has heating temperatures
                for (int T=0; T<tinfo_dielectric[process].size(); ++T) {
                    ds.fit_dielectric_loss
                    (species,process,tinfo_dielectric[process][T]);
                }
            }
            /* fit temperature denpendence of relaxation to the chosen model */
            ds.fit_tauFit(species);
            
        } else {
            
            /* fit independent DS data manually */
            // special case: process=2
            ds.fit_dielectric_loss(species);
        }
        
        time_t endds=time(NULL);
        double timediff=difftime(endds,startds);
        cout << "\nTime elapsed: " << timediff << " sec\n";
        return "\nfinished!\n";
    }
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                 Automate Over All Simulation Regimes                   //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    for (int i=regime_beg; i<=regime_end; ++i)
    {
        time_t begin_regime=time(NULL);
        
        species.set_current_regime(i);
        
        double t_eq=species.get_equilibration_times()[i];
        double tauFit_cut=t_eq/n_equ_blocks;
        species.set_compu_time(tauFit_cut);
        
        /* Simulations & Analyses & Curve-Fitting */
        //-------------------------------------------------
        species.max_times_retry=2;
        species.times_retry=0;
        species.times_write_tempsinfo=0;
        int times_retry=0;
        bool is_retry=false;
        do
        {
            species.times_retry=times_retry;
            autoWork::automated_batch(species,lmp,is_retry);
            
            if (is_retry) {
                ++times_retry;
                cout
                << "\n"
                << "retry Regime_" << species.get_current_regime()
                << "\n"
                << "#" << times_retry << " retry"
                << "\n";
            }
            
        } while (is_retry);
        species.get_simulation_retry().push_back({i,times_retry});
        //-------------------------------------------------
        
        time_t end_regime=time(NULL);
        
        double timediff=difftime(end_regime,begin_regime);
        
        /* store elapsed time of each regime */
        species.get_elapsedtime()[i-regime_beg]=timediff;
        
        cout
        << "\n"
        << "Regime_" << i << " Completed."
        << "\n\n";
        system("date");
        cout << "\n";
        
        cout
        << "Time Elapsed: "
        << timediff << " seconds."
        << "\n\n";
    }
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    // write progress info of final regime to file
    //--------------------------------------------------------
    autoWork::write_progress(species);
    //--------------------------------------------------------
    
    time_t end_overall=time(NULL);
    double timediff_overall=difftime(end_overall,begin_overall);
    species.get_elapsedtime().push_back(timediff_overall);
    species.set_timediff_overall(timediff_overall);
    
    // write final report to file
    //--------------------------------------------------------
    autoWork::write_report(species);
    //--------------------------------------------------------
    
    
    // return path to the final report signifying completion
    //--------------------------------------------------------
    return species.get_sys_targets();
}










////////////////////////////////////////////////////////////////////////////////
//==============================================================================
void autoWork::automated_batch(StructureClass& species,
                               const LmpScripts& lmp,
                               bool& is_retry)
{
    const int n_trial   = species.get_n_trial();
    const int n_system  = species.get_n_system();
    const int n_sys_beg = species.get_n_sys_beg();
    const int n_sys_end = species.get_n_sys_end();
    const int n_regime  = species.get_n_regime();
    const int current_regime = species.get_current_regime();
    
    /************************************************************
     Very Important:
     The temperature containers.
     All are 2D vectors: vector<vector<double>>
     
     in class SysInfo,
     
     -	temperaturesInfo:
     storing temperatures with a size of number of T's
     specified for the current T regime.
     
     -	equilibratedTs:
     storing the in-equilibrium T's; size of container may be smaller
     than the numnber of T's of current regime.
     
     -	quenchingTs:
     storing T's used for equilibration and prodution;
     these T's should be smaller than the lowest in-equilibrium T
     in the previous regime.
     *************************************************************/
    
    vector<vector<double>> tinfo;
    
    if (species.get_is_fullquench()) tinfo = species.get_temperaturesInfo();
    else tinfo = species.get_quenchingTs();
    
    species.set_is_updatetinfo(true);
    
    //#define MANUAL
#ifdef MANUAL
    if (current_regime==0)
    {
        tinfo=
        {
            {},
            {},
            {},
            {}
        };
    }
    else if (current_regime==1)
    {
        tinfo=
        {
            {},
            {},
            {},
            {}
        };
    }
    else if (current_regime==2)
    {
        tinfo=
        {
            {},
            {},
            {},
            {}
        };
    }
    else if (current_regime==3)
    {
        tinfo=
        {
            {},
            {},
            {},
            {}
        };
    }
    else if (current_regime==4)
    {
        tinfo=
        {
            {},
            {},
            {},
            {}
        };
    }
    species.set_is_updatetinfo(false);
    species.get_temperaturesInfo()=tinfo;
    species.get_quenchingTs()=tinfo;
#endif
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                          III. Do Simulations                           //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    WorkScripts ws(species);
    
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    // 'WorkScripts' class setters
    //==========================================================================
    
    /* Generation (@ highest specified T) */
    //--------------------------------------------------------------------------
    const double steps_generation=1e+6;
    ws.set_steps_gen(steps_generation);
    //--------------------------------------------------------------------------
    
    
    /* Quench */
    //--------------------------------------------------------------------------
    // NOTE:
    // unit of quench rate:
    // "real"  'kelvins per nano second' (default: 10)
    // "metal" 'kelvins per nano second' (default: 10)
    // "lj"    'T per tau' (default: 1e-4)
    // ws.set_quenchRate(1e-4);
    //--------------------------------------------------------------------------
    species.set_quenchRate(ws.get_quenchRate());
    //--------------------------------------------------------------------------
    
    
    //--------------------------------------------------------------------------
    // NOTE:
    // Generation,Equilibration,Resize,Production were set in class constructor.
    //--------------------------------------------------------------------------
    
    
    /*==( Setup Compute Nodes )==*/
    
    
    // compute-0-nodes
    //--------------------------------------------------------------------------
    ws.set_is_node0(0);
    ws.set_node0({14,15});
    //--------------------------------------------------------------------------
    
    // compute-1-nodes
    //--------------------------------------------------------------------------
    ws.set_is_node1(1);
    ws.set_node1({8,9});
    //--------------------------------------------------------------------------
    
    species.set_is_node0_info(ws.get_is_node0());
    species.set_is_node1_info(ws.get_is_node1());
    species.set_node0_info(ws.get_node0());
    species.set_node1_info(ws.get_node1());
    
    
    // NOTE:
    // c is how many cores are occupied by one job
    // ex. c=8 means 1 job occupies 8 cores so it uses an entire GPU;
    // this is based on that one compute node has 16 cores and 2 GPUs
    int c=(8/species.get_n_regime_temps()[current_regime]);
    if (c<1) c=1; // c cannot be less than one
    //--------------------------------------------------------------------------
    // order of elements corresponds to simulation stages of:
    // 0: generation    (1cpu_1gpu)
    // 1: quenching     (1cpu_1gpu)
    // 2: equilibration (1cpu_1/8gpu)
    // 3: resizing      (1cpu_1/8gpu)
    // 4: production    (1cpu_1/8gpu)
    //--------------------------------------------------------------------------
    ws.set_numCores({1,1,1,1,1});  // # of total requsted cpu cores
    ws.set_run_cores({1,1,1,1,1}); // # of cpu cores to run jobs on
    ws.set_cores({8,8,c,c,c});     // # used for jobs allocation per half-node
    ws.set_priority({0,0,0,0,0});  // # priority of simulation phases
    //--------------------------------------------------------------------------
    
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    /* write out progress info to file */
    //--------------------------------------------------------------------------
    autoWork::write_progress(species);
    //--------------------------------------------------------------------------
    
    
    if((species.get_is_Simulations())&&
       (species.times_retry<=species.max_times_retry)) // NOTE: "<="
    {
        // target files of simulation
        vector<string> work_targets;
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl)
        {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
            {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                species.set_n_Temp((int)(tinfo[index].size()));
                
                /* Generation Subfile (only in first regime) */
                if ((current_regime==0)&&(is_retry==false)) // NOTE: is_retry
                {
                    ws.make_GenerationSubFiles
                    (species,lmp,n_trl,n_sys,tinfo[index][0]);
                }
                
                if (!species.get_is_singleTempTest())
                {
                    /* Quench Subfile */
                    if (species.get_is_fullquench())
                    {
                        if (current_regime==0) {
                            ws.make_QuenchSubFiles(species,n_trl,n_sys);}
                    }
                    else
                    {
                        ws.make_QuenchSubFiles(species,n_trl,n_sys);
                    }
                    
                    /* Equilibration Subfile */
                    for (int T=0; T<species.get_n_Temp(); ++T)
                    {
                        ws.make_EquilibrationSubfiles
                        (species,lmp,n_trl,n_sys,tinfo[index][T]);
                        
                        /* Specify target files to be watched for */
                        //------------------------------------------------------
                        work_targets.push_back
                        (lmp.work_target(species,"equ",n_trl,n_sys,tinfo[index][T]));
                        //------------------------------------------------------
                    }
                }
            }
        }
        
        /* Generation Jobs Submissions */
        if ((current_regime==0)&&(is_retry==false)) // NOTE: is_retry
        {
            ws.make_GenerationSubScripts(species);
        }
        
        
        if (!species.get_is_singleTempTest())
        {
            /* Quench Jobs Submissions */
            if (species.get_is_fullquench())
            {
                if (current_regime==0) {ws.make_QuenchSubScripts(species);}
            }
            else
            {
                ws.make_QuenchSubScripts(species);
            }
            
            /* Equilibration Jobs Submissions */
            ws.make_EquilibraitionSubScripts(species,tinfo);
            
            /* Watch and Hold */
            //------------------------------------------------------------------
            if (species.get_is_watch_hold())
            {
                autoWork::watch_hold
                (species,work_targets,"simulation",current_regime);
            }
            //------------------------------------------------------------------
            work_targets.clear();
        }
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl)
        {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
            {
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                species.set_n_Temp((int)(tinfo[index].size()));
                
                for (int T=0; T<species.get_n_Temp(); ++T)
                {
                    /* Resizing Subfile */
                    if (!species.get_is_singleTempTest())
                    {
                        if (lmp.get_is_resize())
                        {
                            ws.make_ResizeSubfiles
                            (species,lmp,n_trl,n_sys,tinfo[index][T]);
                        }
                    }
                    
                    /* Production Subfile */
                    ws.make_ProductionSubFiles
                    (species,lmp,n_trl,n_sys,tinfo[index][T]);
                    
                    /* Specify target files to be watched for */
                    //----------------------------------------------------------
                    work_targets.push_back
                    (lmp.work_target(species,"prd",n_trl,n_sys,tinfo[index][T]));
                    //----------------------------------------------------------
                }
            }
        }
        
        if (!species.get_is_singleTempTest())
        {
            /* Resize Jobs Submissions */
            if (lmp.get_is_resize()) {ws.make_ResizeSubScripts(species,tinfo);}
        }
        
        /* Production Jobs Submissions */
        ws.make_ProductionSubScripts(species,tinfo);
        
        /* Watch and Hold */
        //----------------------------------------------------------------------
        if (species.get_is_watch_hold())
        {
            autoWork::watch_hold
            (species,work_targets,"simulation",current_regime);
        }
        //----------------------------------------------------------------------
        work_targets.clear();
    }
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                           IV. AMDAT Analyses                           //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    AmdatAnalysis aa(species,ws);
    
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    // 'AmdatAnalysis' class setters
    //==========================================================================
    aa.set_is_deleteZerothFrames(false); // delete the 0th frame from 2nd to final block
    aa.set_is_changeFrameSteps(false);   // change frame timesteps in custom file
    aa.set_is_changeAtomType(false);     // change atom types in custom file
    aa.set_is_changeAtomOrder(false);    // change orders of atoms in custom file
    aa.set_analysispart("all");          // ex. "all" or "backbone" for analysis
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    /* setup AMDAT jobs allocation configuration */
    aa.set_amdat_cores(1);
    aa.set_amdat_priority(0);
    
    
    if((species.get_is_AMDAT())&&
       (species.times_retry<=species.max_times_retry)) // NOTE: "<="
    {
        
        if (lmp.get_is_resize())
        {
            aa.set_is_NPT(false);
        }
        
        // target files of amdat
        vector<string> amdat_targets;
        
        aa.set_is_strFac(true);
        aa.set_is_msd(true);
        aa.set_is_ngp(true);
        aa.set_is_isfs(true);
        //aa.set_is_composition(true);
        //aa.set_is_u2dist(true);
        //aa.set_is_stiffness_dist(true);
        
        aa.set_is_isf(false);
        aa.set_is_strings(false);
        
        string o=test_AnalysisFolder(species,aa);
        call_system_bash(o);
        
        /* make AMDAT inputfile */
        aa.make_amdatInputFile(species,ws);
        
        for (int n_trl=0; n_trl<n_trial; ++n_trl)
        {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
            {
                // locate position of current system
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                // set # of current running temperatures
                species.set_n_Temp((int)(tinfo[index].size()));
                
                /* set for the box size */
                //----------------------------------------------------------
                if (species.get_maxLenScale()[n_sys]!=0) {
                    aa.set_is_fixed_maxL(true);}
                else {
                    aa.set_is_fixed_maxL(false);}
                //----------------------------------------------------------
                
                /* set for the q index */
                //----------------------------------------------------------
                if (species.get_waveindex()[n_sys]!=0) {
                    aa.set_is_fixed_q(true);}
                else {
                    aa.set_is_fixed_q(false);}
                //----------------------------------------------------------
                
                for (int T=0; T<species.get_n_Temp(); ++T) {
                    
                    /* delete the 0th frame from 2nd to final block */
                    if (species.get_is_set_CheckPoint())
                    {
                        if (aa.get_is_deleteZerothFrames())
                        {
                            aa.delete_zerothFrames
                            (species,n_trl,n_sys,tinfo[index][T]);
                        }
                    }
                    
                    /* change frame timesteps in custom file */
                    if (aa.get_is_changeFrameSteps())
                    {
                        aa.change_prdcustomfilesteps
                        (species,n_trl,n_sys,tinfo[index][T]);
                    }
                    
                    /* change atom type in custom file if it's needed */
                    if (aa.get_is_changeAtomType())
                    {
                        species.set_backboneTypes({1,2,3});
                        species.set_sideGroupTypes({5,7});
                        aa.change_productionAtomTypes
                        (species,n_trl,n_sys,tinfo[index][T]);
                    }
                    
                    /* change the order of atoms in custom file */
                    if (aa.get_is_changeAtomOrder())
                    {
                        aa.change_atom_order
                        (species,n_trl,n_sys,tinfo[index][T]);
                    }
                    
                    /* AMDAT submission files */
                    //------------------------------------------------------
                    aa.make_amdatSubFile(species,ws,n_trl,n_sys,tinfo[index][T]);
                    
                    /* Specify target files to be watched for */
                    //------------------------------------------------------
                    amdat_targets.push_back
                    (aa.amdat_target(species,n_trl,n_sys,tinfo[index][T],"isfs"));
                    //NOTE: "strFac" || "isfs"
                    //------------------------------------------------------
                }
            }
        }
        
        /* AMDAT Jobs Submissions */
        //------------------------------------------------------------------
        aa.make_amdatSubScript(species,tinfo);
        
        /* Watch and Hold */
        //------------------------------------------------------------------
        if (species.get_is_watch_hold())
        {
            autoWork::watch_hold
            (species,amdat_targets,"analysis",current_regime);
        }
        //------------------------------------------------------------------
        amdat_targets.clear();
    }
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    //                         V. Fit Data (by ALGLIB)                        //
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    FitData fd(species,ws,aa);
    
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    // 'FitData' class setters
    //==========================================================================
    
    /* Curve-Fitting Switches */
    fd.set_is_use_FG(true);
    fd.set_is_applytauFitcut(true);
    fd.set_is_imposeStrictTeq(false);
    fd.set_is_fit_full_alpha(true);
    fd.set_is_use_gammafunc(false);
    
    /* Fit Different Relaxation Mechanisms */
    fd.set_is_fit_sExp(true);         // true: fit KWW function
    
    /* Fit Temperature Dependence of Relaxation Time */
    fd.set_is_fit_Arrhenius(true);    // true: fit Arrhemius function
    fd.set_is_fit_tauFit(true);       // true: fit the chosen relaxatiom model
    
    /* Set Customized Analyses */
    fd.set_is_find_DWF(true);         // true: find the Debye-Waller factor
    fd.set_is_calc_thermoData(true);  // true: calculate average thermo data
    fd.set_is_fitByEveryPoint(true);  // true: fit and report Tg by every point
    fd.set_is_fit_by_Tc(false);       // true: fit only T<Tc in final regime
    fd.set_is_fit_by_TA(true);        // true: fit only T<TA in final regime
    fd.set_is_avgFitParams(false);    // true: shoot by avg; false: by trial
    
    
    /* Fitting Models */
    /*--------------------------------------------------------------------------
     Models for Fitting Different Relaxation Mechanisms
     ---------------------------------------------------------------------------
     Density Correlation Function: (fit by Stretched Exponential Function)
     KWW            F = A*exp(-(t/tau)^beta)
     KWW_pwr        F = A*(t^-a)*exp(-(t/tau)^beta)
     mKWW           F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
     mKWW_pwr       F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
     
     ---------------------------------------------------------------------------
     Models for Fitting Temperature Dependence of Relaxation Time:
     ---------------------------------------------------------------------------
     3-PARAMETER:
     ------------
     AM             tau = tau0*exp((p/T)^alpha)
     DG             tau = tau0*exp((A/T)*exp(-B*T))
     VFT            tau = tau0*exp((D*T0)/(T-T0))
     Mauro          tau = tau0*exp((K/T)exp(C/T))
     
     4-PARAMETER:
     ------------
     CG             tau = tau0*exp(B/{(T-T0)+[(T-T0)^2+C*T]^0.5})
     COOP           tau = tau0*exp(E_inf*(1+exp(-u*((T/E_inf)-b)))/T)
     DEAG           tau = tau0*exp(((A-B*T)/T)*exp(C/T))
     FourParamVFT   tau = tau0*exp([(D*T0)/(T-T0)]^alpha)
     
     OTHER:
     ------------
     MCT            tau = A*((T-Tc)/Tc)^(-r)
     SOU            tau = A*((T-Tc)/T)^(-r)
     ArrheniusII    tau = tau0*exp((A/T)+(B/T^2))
     ArrheniusIII   tau = tau0*exp((A/T)+(B/T^2)+(C/T^3))
     
     -------------------------------------------------------------------------*/
    
    
    /* Model for Fitting Relxation Profile */
    //==========================================================================
    string tcalc,format;
    if (fd.get_is_fit_sExp()){
        tcalc="KWW";
    }
    
    
    /* Model for Fitting the T Dependence of Relaxation Time */
    //==========================================================================
    string model="COOP";
    
    
    /* Model for Predicting New Temperatures */
    //==========================================================================
    string extrp;
    if (current_regime==0) {
        extrp="VFT";
    } else {
        extrp="COOP";
    }
    
    
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    species.set_relaxationModel(model);
    species.set_extrpolateModel(extrp);
    
    vector<string> relmodel;
    relmodel.push_back(model);
    relmodel.push_back(extrp);
    
    if (species.get_is_fitData())
    {
        /* clear content of relavent containers used in class FitData */
        //**********************************************************************
        species.get_Tg_extrp()[current_regime].clear();
        species.get_Tg_compu()[current_regime].clear();
        species.get_m_extrp()[current_regime].clear();
        species.get_m_compu()[current_regime].clear();
        //**********************************************************************
        
        if (species.times_retry<=species.max_times_retry) // NOTE: "<="
        {
            if (fd.get_is_fit_sExp())
            {
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                        
                        /* the positioning index */
                        //------------------------------------------------------
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        species.set_n_Temp((int)(tinfo[index].size()));
                        //------------------------------------------------------
                        
                        /* clear content of contianers for new values */
                        //******************************************************
                        species.get_equilibratedTs()[index].clear();
                        species.get_tauFit()[index].clear();
                        species.get_tauEqu()[index].clear();
                        //******************************************************
                        for (int T=0; T<species.get_n_Temp(); ++T) {
                            
                            /* fit relaxation profile to KWW funtion */
                            fd.fit_sExp(species,n_trl,n_sys,tinfo[index][T],tcalc);
                        }
                        /* write temperatures to file */
                        //------------------------------------------------------
                        fd.write_qchTs(species,n_trl,n_sys); // simulated T's
                        fd.write_equTs(species,n_trl,n_sys); // in-equilibrium T's
                        //------------------------------------------------------
                    }
                }
                if (true)
                {
                    //**************************************************************
                    // check the number of in-equilibrium T's of current regime;
                    // if further simulation of current regime is needed, it will be
                    // based on a more aggressive interpolation of temperatures to
                    // ensure the number of in-equilibrium temperatures in current
                    // regime is larger than a threshold number.
                    //**************************************************************
                    is_retry=fd.check_is_retry(species);
                    if (is_retry) {
                        
                        if (species.times_retry<species.max_times_retry) return; // NOTE: "<"
                        else is_retry=false;
                    }
                }
            }
        }
        
        if (fd.get_is_find_DWF())
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl)
            {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                {
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    species.set_n_Temp((int)(tinfo[index].size()));
                    
                    for (int T=0; T<species.get_n_Temp(); ++T) {
                        
                        fd.find_individual_msd_data
                        (species,n_trl,n_sys,tinfo[index][T]);
                    }
                }
            }
            if (current_regime==(n_regime-1))
            {
                fd.write_DWF_equ(species);
                fd.write_DWF_equ_avg(species);
            }
        }
        
        if (fd.get_is_fit_Arrhenius())
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl)
            {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                {
                    /* fit relaxation data to Arrhenius */
                    fd.fit_Arrhenius(species,n_trl,n_sys);
                }
            }
        }
        
        if (fd.get_is_fit_tauFit())
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl)
            {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                {
                    /* fit relaxation to the chosen model */
                    fd.fit_tauFit(species,n_trl,n_sys,relmodel);
                }
            }
        }
        
        if (fd.get_is_calc_thermoData())
        {
            for (int n_trl=0; n_trl<n_trial; ++n_trl)
            {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                {
                    /* the positioning index */
                    //----------------------------------------------------------
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    species.set_n_Temp((int)(tinfo[index].size()));
                    //----------------------------------------------------------
                    for (int T=0; T<species.get_n_Temp(); ++T)
                    {
                        thermodataprocess_prd(species,n_trl,n_sys,tinfo[index][T]);
                        thermocalc_prd(species,n_trl,n_sys,tinfo[index][T]);
                    }
                }
            }
        }
    }
    
    
    /* Theory Test (Not Included in Standard Automation) */
    if(species.get_is_theorytest())
    {
        TheoryTest tt(species,ws,aa,tinfo,model);
        
        tt.set_is_use_ngp_peak_frame(true);
        tt.set_is_data_smoothing(false);
        
        tt.amdatstringanalysis(species,ws,aa);
        
        if (current_regime==(n_regime-1))
        {
            tt.AGtest(species);
            tt.GLMtest(species);
        }
    }
    
}
//==============================================================================
////////////////////////////////////////////////////////////////////////////////










void autoWork::watch_hold(const StructureClass& sysVar,
                          const vector<std::string>& targets,
                          const string& current_state,
                          const int current_regime)
{
    int counter=0;
    int display=0;
    int n_total=(int)targets.size();
    
    /* default wait time = 100s */
    int waitTime=100;
    string system_wait="read -t "+to_string((long long int)waitTime);
    
    do {
        
        counter=0;
        
        /* test if target files exist */
        for (size_t i=0; i<targets.size(); ++i) {
            
            ifstream file(targets[i].c_str());
            //cout << targets[i].c_str() << "\n";
            if (file.good()) ++counter;
        }
        
        string i,test;
        /*
         if there're core.* files in simulations folder,
         stop automation, and delete the core.* files */
        if (current_state=="simulation") {
            
            i.append(return_SimulationFolderPath(sysVar));
            i.append("/core.*");
            test.append("for file in "+i+"; ");
            test.append("do test -e $file; done");
            //cout << test.c_str() << "\n";
        }
        
        /*
         if there're core.* files in analysis folder,
         stop automation, and delete the core.* files */
        else if (current_state=="analysis") {
            
            i.append(return_AnalysisFolderPath(sysVar));
            i.append("/core.*");
            test.append("for file in "+i+"; ");
            test.append("do test -e $file; done");
            //cout << test.c_str() << "\n";
        }
        
        int ret=system(test.c_str());
        if (ret==0) {
            cout
            << "\n"
            << "core.* files detected and now deleted. Automation stopped."
            << "\n";
            system("date");
            string rm="rm -rf "+i;
            system(rm.c_str());
            
            /* qdel the jobs with the same USIC */
            string qdel;
            qdel.append("qdel ");
            qdel.append("*_"+sysVar.get_usic()+"_*");
            system(qdel.c_str());
            system("read -t5");exit(EXIT_FAILURE);
        }
        
        if (counter<n_total)
        {
            /*
             Every 1000 waitTime's, show USIC and automation section */
            if ((display%1000)==0) {
                cout
                << "\n"
                << "Regime ["
                << (current_regime+1) << "/" << sysVar.get_n_regime() << "]: "
                << "in " << current_state
                << "\n";
                cout
                << "\n"
                << "USIC = " << sysVar.get_usic()
                << "\n";
            }
            
            /*
             Every 10 waitTime's, show how many files found */
            if ((display%10)==0) {
                cout
                << "(" << counter << "/" << n_total << ") files found; "
                << "\n";
            }
            system(system_wait.c_str());
        }
        
        ++display;
        
    } while(counter<n_total);
}





void autoWork::write_progress(StructureClass& sysVar)
{
    
    /* Write progress to file stored at the USIC folder */
    string o;
    o.append(return_AnalysisFolderPath(sysVar)+"/..");
    o.append("/progress_");
    o.append(sysVar.get_usic());
    o.append(".txt");
    //cout << o << "\n";
    
    ofstream progress(o.c_str());
    if (progress.is_open()) {
        
        progress
        << "Starting_Time    " << sysVar.get_startdatetime()
        << "\n"
        << "Final_Write_Time " << return_datetime()
        << "\n\n";
        
        progress
        << "tau_alpha Cutoffs ";
        if (sysVar.get_systemUnit()=="real") {
            progress << "(timeUnit: fs)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            progress << "(timeUnit: ps)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="lj") {
            progress << "(timeUnit: tau)" << "\n";
        }
        progress
        << "------------------------------------------"
        << "\n";
        for (int i=0; i<sysVar.get_n_regime(); ++i) {
            
            double previous_teq,current_teq,teq;
            if (i>0) {
                previous_teq = sysVar.get_equilibration_times()[i-1];
                current_teq  = sysVar.get_equilibration_times()[i];
                teq = max(previous_teq,current_teq);
            }
            else {
                teq = sysVar.get_equilibration_times()[i];
            }
            double n_equ_blocks=sysVar.get_n_equ_blocks();
            progress
            << "Regime_" << i << " " << teq/n_equ_blocks
            << "\n";
        }
        progress
        << "------------------------------------------"
        << "\n"
        << "Current_Regime "
        << "\n"
        << (sysVar.get_current_regime()+1) << " of "
        << sysVar.get_n_regime()
        << "\n\n";
        
        
        progress
        << "<Compute Nodes Info>"
        << "\n";
        if (sysVar.get_is_node0_info()) {
            progress
            << "Compute_0: ";
            for (size_t i=0; i<sysVar.get_node0_info().size(); ++i) {
                progress
                << sysVar.get_node0_info()[i] << " ";
            }
            progress << "\n";
        }
        if (sysVar.get_is_node1_info()) {
            progress
            << "Compute_1: ";
            for (size_t i=0; i<sysVar.get_node1_info().size(); ++i) {
                progress
                << sysVar.get_node1_info()[i] << " ";
            }
            progress << "\n";
        }
        progress << "\n";
        
        
        vector<double> teq=sysVar.get_equilibration_times();
        int current_regime=sysVar.get_current_regime();
        double n_prd_blocks=sysVar.get_n_prd_blocks();
        double n_equ_blocks=sysVar.get_n_equ_blocks();
        
        
        progress
        << "System Unit:    " << sysVar.get_systemUnit() << "\n\n"
        
        << "Quenching Rate: " << sysVar.get_quenchRate() << " units\n"
        << "Equilibration:  " << n_equ_blocks << " tau_alpha\n"
        << "Production:     " << n_prd_blocks << " tau_alpha\n\n";
        
        
        progress
        << "Time Elapsed up to Current Regime         " << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  "
            << sysVar.get_elapsedtime()[i] << " "
            << "seconds"
            << "\n";
        }
        
        progress
        << "\n"
        << "Tg_extrp(trials) avg(Tg) stdev(Tg) stdev/avg " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress
                << sysVar.get_Tg_extrp()[i][ii] << " ";
            }
            sysVar.set_calcvector(sysVar.get_Tg_extrp()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            progress
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        
        progress
        << "\n"
        << "m_extrp(trials) avg(m) stdev(m) stdev/avg    " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress
                << sysVar.get_m_extrp()[i][ii] << " ";
            }
            sysVar.set_calcvector(sysVar.get_m_extrp()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            progress
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        
        progress
        << "\n"
        << "Tg_compu(trials) avg(Tg) stdev(Tg) stdev/avg " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            
            progress
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_Tg_compu()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_Tg_compu()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            progress
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        
        progress
        << "\n"
        << "m_compu(trials) avg(m) stdev(m) stdev/avg    " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            
            progress
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_m_compu()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_m_compu()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            progress
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        
    }
    else cout << "progress: progress file cannot open." << "\n";
}





void autoWork::write_report(StructureClass& sysVar)
{
    
    const int n_sys_beg = sysVar.get_n_sys_beg();
    const int n_regime  = sysVar.get_n_regime();
    
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar)+"/..");
    o.append("/report_");
    o.append(sysVar.get_usic());
    o.append(".txt");
    //cout << o << "\n";
    
    sysVar.set_sys_targets(o);
    
    cout
    << "Total Time Elapsed: "
    << sysVar.get_timediff_overall() << " seconds."
    << "\n\n";
    
    ofstream report(o.c_str());
    if (report.is_open()) {
        
        report
        << "Started  " << sysVar.get_startdatetime()
        << "\n"
        << "Finished " << return_datetime()
        << "\n\n"
        
        
        /* report general information */
        //----------------------------------------------------------------------
        << "General Info" << "\n"
        << "------------------------------------------" << "\n"
        << "System_Unit               " << sysVar.get_systemUnit()
        << "\n\n"
        
        << "<Structure>"
        << "\n"
        << "Backbone_type             "
        << sysVar.get_typeB()
        << " (" << sysVar.get_backboneInWords() << ")"
        << "\n"
        << "SideGroup_type            "
        << sysVar.get_typeS()
        << " (" << sysVar.get_sideGroupInWords() << ")"
        << "\n"
        << "Monomers_Per_Chain        " << sysVar.get_backLenVec()[n_sys_beg]
        << "\n"
        << "Number_of_Chains          " << sysVar.get_n_polyVec()[n_sys_beg]
        << "\n"
        << "Number_of_Beads           " << sysVar.get_polyBeads()
        << "\n\n"
        
        << "<Simulation>"
        << "\n"
        << "QuenchRate                " << sysVar.get_quenchRate()
        << "\n"
        << "Equilibration_Blocks      " << sysVar.get_n_equ_blocks()
        << "\n"
        << "Production_Blocks         " << sysVar.get_n_prd_blocks()
        << "\n\n"
        
        << "<Analysis>"
        << "\n"
        << "Relaxation_Model          " << sysVar.get_relaxationModel()
        << "\n"
        << "Extrapolation_Model       " << sysVar.get_extrpolateModel()
        << "\n"
        << "Wave_Number(q)            " << sysVar.get_wavenumber()
        << "\n"
        << "------------------------------------------"
        << "\n\n\n";
        //----------------------------------------------------------------------
        
        
        vector<double> teq=sysVar.get_equilibration_times();
        double n_equ_blocks=sysVar.get_n_equ_blocks();
        
        
        /* report Extrapolated Tg */
        //----------------------------------------------------------------------
        /*
         get_Tg_extrp() is a 2D vector;
         rows store per regime data, columns store per trial data */
        report
        << "Extrapolated Tg [@("<<sysVar.get_extrp_time()<<")timeUnit]"
        << "\n"
        << "------------------------------------------"
        << "\n";
        for (size_t i=0; i<sysVar.get_Tg_extrp()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_extrp()[n_regime-1][i] << "\n";}
        
        report << "\n";
        
        sysVar.set_calcvector(sysVar.get_Tg_extrp()[n_regime-1]);
        double Tg_mean=sysVar.calc_mean();
        double Tg_stdev=sysVar.get_sample_stdev();
        report
        << "mean(Tg)          " << Tg_mean              << "\n"
        << "stdev(Tg)         " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------" << "\n\n";
        
        report
        << "Tg(trials) avg(Tg) stdev(Tg) stdev/avg"     << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_Tg_extrp()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_Tg_extrp()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            report
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        /* report fragility @ Tg_extrp */
        //----------------------------------------------------------------------
        /*
         get_m_extrp() is a 2D vector;
         rows store per regime data, columns store per trial data */
        report
        << "Fragility @ Tg_extrp                      " << "\n"
        << "------------------------------------------" << "\n";
        for (size_t i=0; i<sysVar.get_m_extrp()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_extrp()[n_regime-1][i] << "\n";}
        
        report << "\n";
        
        sysVar.set_calcvector(sysVar.get_m_extrp()[n_regime-1]);
        double m_mean=sysVar.calc_mean();
        double m_stdev=sysVar.get_sample_stdev();
        report
        << "mean(m_extrp)   " << m_mean                 << "\n"
        << "stdev(m_extrp)  " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------" << "\n\n";
        
        report
        << "m(trials) avg(m) stdev(m) stdev/avg"        << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_m_extrp()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_m_extrp()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            report
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        report << "\n\n";
        
        
        /* report Computational Tg */
        //----------------------------------------------------------------------
        /*
         get_Tg_compu() is a 2D vector;
         rows store per regime data, columns store per trial data */
        report
        << "Computational Tg [@("<<sysVar.get_compu_time()<<")timeUnit]" << "\n"
        << "------------------------------------------" << "\n";
        for (size_t i=0; i<sysVar.get_Tg_compu()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_compu()[n_regime-1][i] << "\n";}
        
        report << "\n";
        
        sysVar.set_calcvector(sysVar.get_Tg_compu()[n_regime-1]);
        Tg_mean=sysVar.calc_mean();
        Tg_stdev=sysVar.get_sample_stdev();
        report
        << "mean(Tg_compu)    " << Tg_mean              << "\n"
        << "stdev(Tg_compu)   " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------" << "\n\n";
        
        report
        << "Tg(trials) avg(Tg) stdev(Tg) stdev/avg"     << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_Tg_compu()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_Tg_compu()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            report
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        /* report fragility @ Tg_compu */
        //----------------------------------------------------------------------
        /*
         get_m_compu() is a 2D vector;
         rows store per regime data, columns store per trial data */
        report
        << "Fragility @ Tg_compu                      " << "\n"
        << "------------------------------------------" << "\n";
        for (size_t i=0; i<sysVar.get_m_compu()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_compu()[n_regime-1][i] << "\n";}
        
        report << "\n";
        
        sysVar.set_calcvector(sysVar.get_m_compu()[n_regime-1]);
        m_mean=sysVar.calc_mean();
        m_stdev=sysVar.get_sample_stdev();
        report
        << "mean(m_compu)   " << m_mean                 << "\n"
        << "stdev(m_compu)  " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------" << "\n\n";
        
        report
        << "m(trials) avg(m) stdev(m) stdev/avg"        << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  ";
            
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_m_compu()[i][ii] << " ";}
            
            sysVar.set_calcvector(sysVar.get_m_compu()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            report
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        report << "\n\n";
        
        
        /* report number of T's simulated in each regime */
        //----------------------------------------------------------------------
        report
        << "Number of T's Simulated in Each Regime    " << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  "
            << sysVar.get_n_regime_temps()[i]
            << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        /* report timestep size used in each regime */
        //----------------------------------------------------------------------
        report
        << "Timestep size in each regime ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        }
        report
        << "------------------------------------------"
        << "\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "ts[" << i << "]=  "
            << sysVar.get_ts_regime()[i]
            << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        /* report tau_alpha cutoffs */
        //----------------------------------------------------------------------
        report
        << "tau_alpha Cutoffs ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        }
        report
        << "------------------------------------------"
        << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            double previous_teq,current_teq,teq;
            if (i>0) {
                previous_teq = sysVar.get_equilibration_times()[i-1];
                current_teq  = sysVar.get_equilibration_times()[i];
                teq = max(previous_teq,current_teq);
            }
            else {
                teq = sysVar.get_equilibration_times()[i];
            }
            double n_equ_blocks=sysVar.get_n_equ_blocks();
            report
            << "Regime_" << i << " " << teq/n_equ_blocks
            << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        
        /* report equilibration times */
        //----------------------------------------------------------------------
        report
        << "Equilibration Times ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        }
        else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        }
        report
        << "------------------------------------------"
        << "\n";
        for (int i=0; i<n_regime; ++i) {
            
            double previous_teq,current_teq,teq;
            if (i>0) {
                previous_teq = sysVar.get_equilibration_times()[i-1];
                current_teq  = sysVar.get_equilibration_times()[i];
                teq = max(previous_teq,current_teq);
            }
            else {
                teq = sysVar.get_equilibration_times()[i];
            }
            
            report
            << "Regime_" << i << " " << teq
            << "\n";
        }
        report
        << "------------------------------------------"
        << "\n\n";
        //----------------------------------------------------------------------
        
        
        /* report if simulations were retried */
        //----------------------------------------------------------------------
        bool is_retry=false;
        for (int i=0; i<sysVar.get_simulation_retry().size(); ++i) {
            if (sysVar.get_simulation_retry()[i][1]!=0) {
                is_retry=true;
                break;
            }
        }
        if (is_retry) {
            report
            << "Regime(tau_alpha)  Number of Times Retried" << "\n"
            << "------------------------------------------" << "\n";
            for (int i=0; i<n_regime; ++i) {
                report
                << "Regime_" << sysVar.get_simulation_retry()[i][0] << " "
                << teq[i]/n_equ_blocks << "  "
                << sysVar.get_simulation_retry()[i][1]
                << "\n";
            }
            report
            << "------------------------------------------"
            << "\n\n";
        }
        //----------------------------------------------------------------------
        
        
        /* report elapased time */
        //----------------------------------------------------------------------
        report
        << "Time Elapsed for Each Temperature Regime  " << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "Regime_" << i << " "
            << teq[i]/n_equ_blocks << "  "
            << sysVar.get_elapsedtime()[i]
            << " seconds"
            << "\n";
        }
        report
        << "------------------------------------------"
        << "\n";
        
        report
        << "\n"
        << "Total Elapsed Time: "
        << sysVar.get_elapsedtime()[n_regime]
        << " seconds"
        << "\n";
        //----------------------------------------------------------------------
        
        
        report.close();
    }
    else cout << "report: final report file cannot open." << "\n";
}






//=============================================================
// Function that uses system command to invoke bash
//=============================================================
void autoWork::call_system_bash(string& strVar)
{
    string str0,str1;
    str0="bash ./";
    str1=str0+strVar;
    /* invoke bash */
    int ret=system(str1.c_str());
    
    vector<string> strVec;
    /* This is the list of files that should not be deleted when success */
    //strVec.push_back("pkm.inp");
    
    int count=0;
    for (size_t i=0; i<strVec.size(); ++i) {
        if(strVar==strVec[i]) ++count;
    }
    
    /*
     If not in the above list,
     the temporary files are deleted after use */
    if (count==0) {
        str0="rm ./";
        str1=str0+strVar;
        if (ret!=0) exit(EXIT_FAILURE);
        else system(str1.c_str());
    }
}





//=============================================================
// Get overall temperatures, highT and loT regime temperatures
//=============================================================
vector<vector<double>> autoWork::linearTemperatures(const SysInfo& sysVar)
{
    vector<double> hiT,loT,tinfo;
    vector<vector<double>> sortedTemps;
    
    //==============================================
    // use the upper and lower temperature bounds
    //==============================================
    if (sysVar.get_n_Temp()==0) {
        
        /* manually set above crossTemp */
        int n_hiT=(sysVar.get_startTemp()-sysVar.get_crossTemp())/sysVar.get_hTres();
        for (int i=0; i<n_hiT; ++i)	{
            tinfo.push_back(sysVar.get_startTemp()-sysVar.get_hTres()*i);
            hiT.push_back(sysVar.get_startTemp()-sysVar.get_hTres()*i);
        }
        /* manually set below crossTemp */
        int n_loT=((sysVar.get_crossTemp()-sysVar.get_finalTemp())/sysVar.get_lTres())+1;
        for (int i=0; i<n_loT; ++i){
            tinfo.push_back(sysVar.get_crossTemp()-sysVar.get_lTres()*i);
            loT.push_back(sysVar.get_crossTemp()-sysVar.get_lTres()*i);
        }
        
        sortedTemps.push_back(tinfo); // index 0: tinfo
        sortedTemps.push_back(hiT);   // index 1: hiT
        sortedTemps.push_back(loT);   // index 2: loT
    }
    
    //==============================================
    // use N starting temperatures
    //==============================================
    else
    {
        double T=sysVar.get_startTemp();
        
        tinfo.push_back(T);
        if (T>sysVar.get_crossTemp()) {hiT.push_back(T);}
        else {loT.push_back(T);}
        
        for (int i=1; i<sysVar.get_n_Temp(); ++i){
            
            
            /* manually set above crossTemp */
            if (T>sysVar.get_crossTemp()) {
                T -= sysVar.get_hTres();
                tinfo.push_back(T);
                hiT.push_back(T);
            }
            /* manually set below crossTemp */
            else if (T<=sysVar.get_crossTemp()) {
                T -= sysVar.get_lTres();
                tinfo.push_back(T);
                loT.push_back(T);
            }
        }
        sortedTemps.push_back(tinfo); // index 0: tinfo
        sortedTemps.push_back(hiT);   // index 1: hiT
        sortedTemps.push_back(loT);   // index 2: loT
    }
    
    return sortedTemps;
}





//=============================================================
// Convenience functions used throughout the program
//=============================================================
const string autoWork::return_year()
{
    time_t t=time(NULL);tm* timePtr=localtime(&t);
    stringstream ss; ss << timePtr->tm_year+1900;
    string year = ss.str();
    
    return year;
}





/*==( format: M/D/Y )==*/
const string autoWork::return_date()
{
    //cout << "day of month = " << timePtr->tm_mday << endl;
    //cout << "month of year = " << timePtr->tm_mon << endl;
    //cout << "year = " << timePtr->tm_year+1900 << endl;
    //cout << "weekday = " << timePtr->tm_wday << endl;
    //cout << "day of year = " << timePtr->tm_yday << endl;
    //cout << "daylight savings = " << timePtr->tm_isdst << endl;
    
    time_t t=time(NULL);
    tm* timePtr=localtime(&t);
    
    stringstream year,month,day;
    month<< timePtr->tm_mon+1;
    day  << timePtr->tm_mday;
    year << timePtr->tm_year+1900;
    
    string date;
    date.append(month.str());
    date.append("/");
    date.append(day.str());
    date.append("/");
    date.append(year.str());
    
    return date;
}





const string autoWork::return_datetime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y/%m/%d %X", &tstruct);
    
    return buf;
}





const string autoWork::make_SimulationsFolders(const StructureClass& sysVar)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    
    ofstream sf("simFolder.sh");
    if ( sf.is_open() )
    {
        sf
        << "#!/bin/bash"                                                << "\n"
        << "Automation=" << path_d                                      << "\n"
        << "cd ${Automation}"                                           << "\n";
        
        sf
        /* make usic directory */
        << "usic_dir="
        << "${Automation}/simulations/"
        << simType_d << "/" << year_d << "/" << usic_d                  << "\n";
        
        
        /* whether to check if 'simulations' folder exists */
        if (1) {
            sf
            << "test -e \"${usic_dir}\""                                << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "	mkdir -p ${usic_dir}"                               << "\n"
            << "else"                                                   << "\n"
            << "	read -t10 -p \"" << usic_d
            << " folder exists, delete and make new? (y/n)\" VAR"       << "\n"
            << "	if [ \"$VAR\" == \"y\" ]; then"                     << "\n"
            << "		rm -r ${usic_dir}"                              << "\n"
            << "		wait"                                           << "\n"
            << "		mkdir -p ${usic_dir}"                           << "\n"
            << "	else"                                               << "\n"
            << "	    rm simFolder.sh"                                << "\n"
            << "		echo -e \"\\nEXIT!\\n\"; exit 3"                << "\n"
            << "	fi"                                                 << "\n"
            << "fi"                                                     << "\n";
        }
        else sf	<< " mkdir -p ${usic_dir}"                              << "\n";
        
        
        /* make the "simulations" directry in the usic folder */
        sf
        << "simulations_dir=${usic_dir}/simulations"                    << "\n"
        << "mkdir -p ${simulations_dir}"                                << "\n"
        
        /* make packmol folder */
        << "packmol=${simulations_dir}/packmol"                         << "\n"
        << "mkdir -p ${packmol}"                                        << "\n"
        
        /* make lammps_inputs folders */
        << "lammps_inputs=${simulations_dir}/lammps_inputs"             << "\n"
        << "mkdir -p ${lammps_inputs}/start_data"                       << "\n"
        << "mkdir -p ${lammps_inputs}/generation"                       << "\n"
        << "mkdir -p ${lammps_inputs}/quench"                           << "\n"
        << "mkdir -p ${lammps_inputs}/equilibration"                    << "\n"
        << "mkdir -p ${lammps_inputs}/tampering"                        << "\n"
        << "mkdir -p ${lammps_inputs}/resize"                           << "\n"
        << "mkdir -p ${lammps_inputs}/production"                       << "\n"
        
        /* theses are the folders to be put in each simulation phase */
        << "mkdir -p ${usic_dir}/packet/log"                            << "\n"
        << "mkdir -p ${usic_dir}/packet/screen"                         << "\n"
        << "mkdir -p ${usic_dir}/packet/trajectory"                     << "\n"
        << "mkdir -p ${usic_dir}/packet/restart"                        << "\n"
        << "mkdir -p ${usic_dir}/packet/submission_files"               << "\n"
        << "mkdir -p ${usic_dir}/packet/submission_files/cluster_out"   << "\n"
        
        /* make simualtion phase directories */
        << "mkdir -p ${simulations_dir}/generation"                     << "\n"
        << "mkdir -p ${simulations_dir}/quench"                         << "\n"
        << "mkdir -p ${simulations_dir}/equilibration"                  << "\n"
        << "mkdir -p ${simulations_dir}/tampering"                      << "\n"
        << "mkdir -p ${simulations_dir}/resize"                         << "\n"
        << "mkdir -p ${simulations_dir}/production"                     << "\n"
        << "mkdir -p ${simulations_dir}/submission_scripts"             << "\n"
        
        /* copy content of packet folder to each phase */
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/generation/"   << "\n"
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/quench/"       << "\n"
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/equilibration/"<< "\n"
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/tampering/"    << "\n"
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/resize/"       << "\n"
        << "cp -r ${usic_dir}/packet/* ${simulations_dir}/production/"   << "\n"
        
        /* make the "analysis" directry in the usic folder */
        << "analysis_dir=${usic_dir}/analysis"                           << "\n"
        << "mkdir -p ${analysis_dir}"                                    << "\n"
        
        /* make corresponding directories in the Analysis folder */
        << "mkdir -p ${analysis_dir}/AMDAT_inputs"                       << "\n"
        << "mkdir -p ${analysis_dir}/AMDAT_submission_files"             << "\n"
        << "mkdir -p ${analysis_dir}/AMDAT_submission_files/cluster_out" << "\n"
        << "mkdir -p ${analysis_dir}/submission_scripts"                 << "\n"
        << "mkdir -p ${analysis_dir}/screen"                             << "\n"
        << "mkdir -p ${analysis_dir}/statistics"                         << "\n"
        << "mkdir -p ${analysis_dir}/fit_data"                           << "\n"
        
        /* source files backup */
        << "src_backup=${usic_dir}/src_backup"                           << "\n"
        << "mkdir -p ${src_backup}"                                      << "\n"
        << "cp *.cpp *.h ${src_backup}/"                                 << "\n"
        << "cp ./packmol ${src_backup}/"                                << "\n"
        
        /* remove packet folder */
        << "rm -rf ${usic_dir}/packet" << "\n"
        
        /* remove archive folder */
        << "rm -rf ${Automation}/archive" << "\n";
        
        sf.close();
    }
    else cout << "simFolders: 'simFolder.sh' cannot open." << "\n";
    
    return ("simFolder.sh");
}





const string autoWork::test_AnalysisFolder(const StructureClass& sysVar,
                                           const AmdatAnalysis& aa)
{
    const string path_d=return_AnalysisFolderPath(sysVar);
    
    ofstream t("amdatfolder.sh");
    if(t.is_open()) {
        
        t
        << "#!/bin/bash"                                                << "\n"
        
        << "\n"
        
        /* test the "analysis" folder */
        << "test -e " << path_d                                         << "\n"
        << "if [ \"$?\" -ne \"0\" ]; then"                              << "\n"
        << "  echo -e \"\\n" << path_d << " does NOT exist!\\n\""       << "\n"
        << "  exit 3"                                                   << "\n"
        << "fi"                                                         << "\n"
        << "\n"
        
        /* create the "screen" folder */
        << "test -e " << path_d  << "/screen"                           << "\n"
        << "if [ \"$?\" -ne \"0\" ]; then"                              << "\n"
        << "  mkdir " << path_d << "/screen"                            << "\n"
        << "fi"                                                         << "\n"
        << "\n"
        
        /* create the "submission_script" folder */
        << "test -e " << path_d  << "/submission_scripts"               << "\n"
        << "if [ \"$?\" -ne \"0\" ]; then"                              << "\n"
        << "  mkdir " << path_d << "/submission_scripts"                << "\n"
        << "fi"                                                         << "\n"
        << "\n"
        
        /* create the "fit_data" folder */
        << "test -e " << path_d  << "/fit_data"                         << "\n"
        << "if [ \"$?\" -ne \"0\" ]; then"                              << "\n"
        << "  mkdir " << path_d << "/fit_data"                          << "\n"
        << "fi"                                                         << "\n"
        
        << "\n";
        
        if (aa.get_is_strFac()) {
            t
            << "strFac_all="
            << path_d << "/statistics/strFac/all"                       << "\n"
            << "strFac_backbone="
            << path_d << "/statistics/strFac/backbone"                  << "\n"
            << "test -e ${strFac_all}"                                  << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${strFac_all}"                               << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${strFac_backbone}"                             << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${strFac_backbone}"                          << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_msd()) {
            t
            << "msd_all="
            << path_d << "/statistics/msd/all"                          << "\n"
            << "msd_backbone="
            << path_d << "/statistics/msd/backbone"                     << "\n"
            << "test -e ${msd_all}"                                     << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${msd_all}"                                  << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${msd_backbone}"                                << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${msd_backbone}"                             << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_ngp()) {
            t
            << "ngp_all="
            << path_d << "/statistics/ngp/all"                          << "\n"
            << "ngp_backbone="
            << path_d << "/statistics/ngp/backbone"                     << "\n"
            << "test -e ${ngp_all}"                                     << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${ngp_all}"                                  << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${ngp_backbone}"                                << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${ngp_backbone}"                             << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_isfs()) {
            t
            << "isfs_all="
            << path_d << "/statistics/isfs/all"                         << "\n"
            << "isfs_backbone="
            << path_d << "/statistics/isfs/backbone"                    << "\n"
            << "test -e ${isfs_all}"                                    << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${isfs_all}"                                 << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${isfs_backbone}"                               << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${isfs_backbone}"                            << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_composition()) {
            t
            << "comp_all="
            << path_d << "/statistics/composition/all"                  << "\n"
            << "comp_backbone="
            << path_d << "/statistics/composition/backbone"             << "\n"
            << "test -e ${comp_all}"                                    << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${comp_all}"                                 << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${comp_backbone}"                               << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${comp_backbone}"                            << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_u2dist()) {
            t
            << "u2dist_all="
            << path_d << "/statistics/u2dist/all"                       << "\n"
            << "u2dist_backbone="
            << path_d << "/statistics/u2dist/backbone"                  << "\n"
            << "test -e ${u2dist_all}"                                  << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${u2dist_all}"                               << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${u2dist_backbone}"                             << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${u2dist_backbone}"                          << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_stiffness_dist()) {
            t
            << "stiffness_dist_all="
            << path_d << "/statistics/stiffness_dist/all"               << "\n"
            << "stiffness_dist_backbone="
            << path_d << "/statistics/stiffness_dist/backbone"          << "\n"
            << "test -e ${stiffness_dist_all}"                          << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${stiffness_dist_all}"                       << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${stiffness_dist_backbone}"                     << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${stiffness_dist_backbone}"                  << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_isf()) {
            t
            << "isf_all="
            << path_d << "/statistics/isf/all"                          << "\n"
            << "isf_backbone="
            << path_d << "/statistics/isf/backbone"                     << "\n"
            << "test -e ${isf_all}"                                     << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${isf_all}"                                  << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${isf_backbone}"                                << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${isf_backbone}"                             << "\n"
            << "fi"                                                     << "\n";
        }
        if (aa.get_is_strings()) {
            t
            << "strings_all="
            << path_d << "/statistics/strings/all"                      << "\n"
            << "strings_backbone="
            << path_d << "/statistics/strings/backbone"                 << "\n"
            
            << "test -e ${strings_all}"                                 << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${strings_all}"                              << "\n"
            << "fi"                                                     << "\n"
            << "test -e ${strings_backbone}"                            << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "  mkdir -p ${strings_backbone}"                         << "\n"
            << "fi"                                                     << "\n";
        }
        
        t.close();
    }
    else cout << "amdatfolder: 'amdatfolder.sh' cannot open." << "\n";
    
    return ("amdatfolder.sh");
}





const string autoWork::move_InnerLoopFiles(const StructureClass& sysVar,
                                           const int n_sys,
                                           const int n_trl)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    
    ofstream mv("mvInner.sh");
    if ( mv.is_open() )
    {
        //======================================================================
        // Level of loops:
        //======================================================================
        // "trial is the inner most loop"
        // Files that need to change names/move are:
        // pkm_output.xyz, start.data, pkm.inp
        //======================================================================
        mv
        << "#!/bin/bash"                                                << "\n"
        
        /* master folder */
        << "Automation=" << path_d                                      << "\n"
        
        /* usic folder */
        << "simulations_dir=${Automation}/simulations/"
        << simType_d << "/" << year_d << "/" << usic_d << "/simulations"<< "\n"
        
        /* packmol folder */
        << "packmolDir=${simulations_dir}/packmol"                      << "\n"
        
        /* lammps_inputs folder */
        << "lammps_inputs=${simulations_dir}/lammps_inputs"             << "\n";
        
        /* change names */
        mv
        << "mv "
        << "./pkm_output.xyz "
        << "./pkm_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << ".xyz"                                                       << "\n";
        
        if (1) {
            mv
            << "mv "
            << "./start.data "
            << "./start_" << usic_d
            << "_00" << n_trl << "_"
            << sysVar.get_nameString(n_sys)
            << ".data"                                                  << "\n";
        }
        else {
            mv
            << "mv "
            << "./start.data "
            << "./start_"
            << "00" << n_trl
            << ".data"                                                  << "\n";
        }
        
        mv
        << "mv "
        << "./pkm.inp "
        << "./pkmInput_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << ".inp"                                                       << "\n";
        
        /* move to folders */
        mv
        << "mv ./pkm_*      ${packmolDir}/"                             << "\n"
        << "mv ./start_*    ${lammps_inputs}/start_data/"               << "\n"
        << "mv ./pkmInput_* ${packmolDir}/"                             << "\n";
        
        mv.close();
    }
    else cout << "mvInner: 'mv.sh' cannot open." << "\n";
    
    return ("mvInner.sh");
}





const string autoWork::move_OuterLoopFiles(const StructureClass& sysVar,
                                           const int n_sys)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    
    ofstream mv("mvOuter.sh");
    if ( mv.is_open() )
    {
        //======================================================================
        // In the outer loop,
        // files that need to change names/move are:
        // singleChain.xyz, Bonds.txt, Angles.txt, sys_Params.txt
        //======================================================================
        mv
        << "#!/bin/bash"                                                << "\n"
        
        /* master folder */
        << "Automation=" << path_d                                      << "\n"
        
        /* usic folder */
        << "usic_dir=${Automation}/simulations/"
        << simType_d << "/" << year_d << "/" << usic_d                  << "\n"
        
        /* singleChain folder */
        << "singleChain=${usic_dir}/singleChain"                        << "\n"
        << "mkdir -p ${singleChain}"                                    << "\n";
        
        /* change names */
        mv
        << "mv "
        << "./singleChain.xyz "
        << "./single_"
        << sysVar.get_nameString(n_sys)
        << ".xyz"                                                       << "\n"
        
        << "mv "
        << "./Bonds.txt "
        << "./Bonds_"
        << sysVar.get_nameString(n_sys)
        << ".txt"                                                       << "\n"
        
        << "mv "
        << "./Angles.txt "
        << "./Angles_"
        << sysVar.get_nameString(n_sys)
        << ".txt"                                                       << "\n"
        
        << "mv "
        << "./sys_Params.txt "
        << "./sys_Params_"
        << sysVar.get_nameString(n_sys)
        << ".txt"                                                       << "\n";
        
        /* move to folders */
        mv
        << "mv ./single_*     ${singleChain}/"                          << "\n"
        << "mv ./Bonds_*      ${singleChain}/"                          << "\n"
        << "mv ./Angles_*     ${singleChain}/"                          << "\n"
        << "mv ./sys_Params_* ${singleChain}/"                          << "\n";
        
        mv.close();
    }
    else cout << "mvOuter: 'mv.sh' cannot open." << "\n";
    
    return ("mvOuter.sh");
}





const string autoWork::make_InitialFolders(const SysInfo& sysVar)
{
    const string path_d = sysVar.get_Path();
    const string usic_d = sysVar.get_usic();
    
    ofstream init("mkdir_init.sh");
    if (init.is_open())
    {
        init
        << "#!/bin/bash"                                                << "\n"
        << "Automation=" << path_d                                      << "\n"
        << "packmol=${Automation}/packmol"                              << "\n"
        << "inputFolder=${Automation}/input_" << usic_d                 << "\n";
        
        /* whether to check if 'inputFolder' exists */
        if (1) {
            init
            << "test -e \"${inputFolder}\""                             << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "	mkdir ${inputFolder}"                               << "\n"
            << "else"                                                   << "\n"
            << "	read -t10 -p \"input_" << usic_d
            << " exists, delete and make new? (y/n)\" VAR"              << "\n"
            << "	if [ \"$VAR\" == \"y\" ]; then"                     << "\n"
            << "		rm -r ${inputFolder}"                           << "\n"
            << "		wait"                                           << "\n"
            << "		mkdir ${inputFolder}"                           << "\n"
            << "	else"                                               << "\n"
            << "		echo -e \"\\nEXIT!\\n\" ; exit 3"               << "\n"
            << "	fi"                                                 << "\n"
            << "fi"                                                     << "\n";
        }
        /* always delete the existing and replace with new one */
        else {
            init
            << "test -e \"${inputFolder}\""                             << "\n"
            << "if [ \"$?\" -ne \"0\" ]; then"                          << "\n"
            << "	mkdir ${inputFolder}"                               << "\n"
            << "else"                                                   << "\n"
            << "    rm -r ${inputFolder}"                               << "\n"
            << "	wait"                                               << "\n"
            << "	mkdir ${inputFolder}"                               << "\n"
            << "fi"                                                     << "\n";
        }
        
        init
        << "archive=${Automation}/archive"                              << "\n"
        << "test -e \"${archive}\""                                     << "\n"
        << "if [ \"$?\" -ne \"0\" ]; then"                              << "\n"
        << "	mkdir ${archive}"                                       << "\n"
        << "fi"                                                         << "\n"
        
        /* specify lammps_inputs folders */
        << "lammps_inputs=${inputFolder}/lammps_inputs"                 << "\n"
        << "start_data=${lammps_inputs}/start_data"                     << "\n"
        << "generation=${lammps_inputs}/generation"                     << "\n"
        << "quench=${lammps_inputs}/quench"                             << "\n"
        << "equilibration=${lammps_inputs}/equilibration"               << "\n"
        << "production=${lammps_inputs}/production"                     << "\n"
        /* make lammps_inputs folders */
        << "mkdir ${lammps_inputs}"                                     << "\n"
        << "mkdir ${start_data}"                                        << "\n"
        << "mkdir ${generation}"                                        << "\n"
        << "mkdir ${quench}"                                            << "\n"
        << "mkdir ${equilibration}"                                     << "\n"
        << "mkdir ${production}"                                        << "\n"
        
        /* specify other relavent folders */
        << "other_files=${inputFolder}/other_files"                     << "\n"
        << "src_backup=${other_files}/src_backup"                       << "\n"
        << "packmolDir=${inputFolder}/packmol"                          << "\n"
        /* make other relavent folders */
        << "mkdir ${other_files}"                                       << "\n"
        << "mkdir ${src_backup}"                                        << "\n"
        << "mkdir ${packmolDir}"                                        << "\n"
        /* copy source codes to the "src_backup" folder */
        << "## copy source codes to ./other_files"                      << "\n"
        << "cp *.cpp *.h ${src_backup}/"                                << "\n"
        << "cp ${packmol} ${src_backup}/"                               << "\n";
        
        init.close();
    }
    else cout << "mkdir_init: 'mkdir_init.sh' cannot open." << "\n";
    
    return ("mkdir_init.sh");
}





const string autoWork::archive_InitialFiles(const StructureClass& sysVar)
{
    
    /*
     The purpose of archiving the input folder is preserving the initial
     system configuration, which is stored in the start.data file.
     Every time the program is invoked, a new system configuration will be
     created. If for instance, one needs the original data file information of
     a particular usic simulation, he then can refer to the start.data file
     in the archive with the corresponsing usic. */
    
    const string path_d = sysVar.get_Path();
    const string usic_d = sysVar.get_usic();
    
    ofstream arc("arc.sh");
    if ( arc.is_open() )
    {
        arc
        << "#!/bin/bash"                                                << "\n"
        << "Automation=" << path_d                                      << "\n"
        << "inputFolder=${Automation}/input_" << usic_d                 << "\n"
        << "archive=${Automation}/archive"                              << "\n"
        << "arcinput=${archive}/input_" << usic_d                       << "\n"
        << "other_files=${inputFolder}/other_files"                     << "\n"
        << "singleChain=${other_files}/singleChain"                     << "\n"
        
        << "cd ${Automation}/archive"                                   << "\n";
        
        /* whether to check if any existing input in the archive */
        if (0) {
            arc
            << "test -e \"${arcinput}\""                                << "\n"
            << "if [ \"$?\" -eq \"0\" ]; then"                          << "\n"
            << "  read -t10 -p \"input_" << usic_d
            << " exists in the Archive, replace?(y/n)\" VAR"            << "\n"
            << "	if [ \"$VAR\" == \"y\" ]; then"                     << "\n"
            << "		rm -r ${arcinput}"                              << "\n"
            << "		wait"                                           << "\n"
            << "	else"                                               << "\n"
            << "		echo -e \"\\nEXIT!\\n\" ; exit 3"               << "\n"
            << "	fi"                                                 << "\n"
            << "fi"                                                     << "\n"
            << "mv ${inputFolder} ${archive}"                           << "\n";
        }
        /* always delete the existing input folder in the archive */
        else {
            arc
            << "test -e \"${arcinput}\""                                << "\n"
            << "if [ \"$?\" -eq \"0\" ]; then"                          << "\n"
            << "  rm -r ${arcinput}"                                    << "\n"
            << "  wait"                                                 << "\n"
            << "fi"                                                     << "\n"
            << "mv ${inputFolder} ${archive}"                           << "\n";
        }
        
        arc
        << "echo -e \"\\nFinished Data Files Generation!\\n\""          << "\n";
        
        arc.close();
    }
    else cout << "arc: 'arc.sh' cannot open." << "\n";
    
    return ("arc.sh");
}





/*==( get the path to the Simulation folder )==*/
const string autoWork::return_SimulationFolderPath(const SysInfo& sysVar)
{
    string o;
    o.append(sysVar.get_Path());
    o.append("/simulations/");
    o.append(sysVar.get_simType());
    o.append("/");
    o.append(sysVar.get_year());
    o.append("/");
    o.append(sysVar.get_usic());
    o.append("/simulations");
    
    return o;
}





/*==( get the path to the Analysis folder )==*/
const string autoWork::return_AnalysisFolderPath(const SysInfo& sysVar)
{
    string o;
    o.append(sysVar.get_Path());
    o.append("/simulations/");
    o.append(sysVar.get_simType());
    o.append("/");
    o.append(sysVar.get_year());
    o.append("/");
    o.append(sysVar.get_usic());
    o.append("/analysis");
    
    return o;
}





const string autoWork::echo_NodesUsage(const vector<int> &node0_d,
                                       const vector<int> &node1_d)
{
    ofstream sh("num.sh");
    if(sh.is_open())
    {
        sh
        << "#!/bin/bash"                                                << "\n"
        
        << "node0=(";
        for (size_t i=0; i<node0_d.size(); ++i) {
            sh << node0_d[i];
            if(i!=(node0_d.size()-1)) sh << " ";
        }
        sh << ")"                                                       << "\n";
        
        sh << "node1=(";
        for (size_t i=0; i<node1_d.size(); ++i) {
            sh << node1_d[i];
            if(i!=(node1_d.size()-1)) sh << " ";
        }
        sh << ")"                                                       << "\n";
        
        sh
        << "numNode0=${#node0[@]}"                                      << "\n"
        << "numNode1=${#node1[@]}"                                      << "\n"
        << "sumNode=$[${numNode0}+${numNode1}]"                         << "\n"
        << "hlfNode=$[${sumNode}/2]"                                    << "\n"
        << "numGPU=$[${sumNode}*2]"                                     << "\n"
        << "numCPU=$[${sumNode}*16]"                                    << "\n"
        
        << "\n"
        
        << "echo 'numNode0 = '\"${#node0[@]}\""                         << "\n"
        << "echo 'node0    = ('\"${node0[@]}\"')'"                      << "\n"
        << "echo 'numNode1 = '\"${#node1[@]}\""                         << "\n"
        << "echo 'node1    = ('\"${node1[@]}\"')'"                      << "\n"
        << "echo 'sumNode  = '\"${sumNode}\""                           << "\n"
        << "echo 'hlfNode  = '\"${hlfNode}\""                           << "\n"
        << "echo 'numGPU   = '\"${numGPU}\""                            << "\n"
        << "echo 'numCPU   = '\"${numCPU}\""                            << "\n"
        
        << "\n";
        
        //sh << "read -p \"press [Enter] to continue\"" << "\n";
        
        sh.close();
    }
    else cout << "num: 'num.sh' cannot open." << "\n";
    
    return ("num.sh");
}





void autoWork::clean_tmpfiles(const StructureClass& sysVar)
{
    vector<string> filesToClean;
    filesToClean.push_back("singleChain.xyz");
    filesToClean.push_back("Bonds.txt");
    filesToClean.push_back("Angles.txt");
    filesToClean.push_back("sys_Params.txt");
    filesToClean.push_back("lmp_header.txt");
    filesToClean.push_back("tmp_mass.txt");
    filesToClean.push_back("bonds_Params.txt");
    filesToClean.push_back("angles_Params.txt");
    filesToClean.push_back("pairs_Params.txt");
    
    string rm="rm ./";
    // rm all specifically defined tmp files
    for (size_t i=0; i<filesToClean.size(); ++i) {
        string tmpStr="test -e "+filesToClean[i];
        int ret=system(tmpStr.c_str());
        if (ret==0) {
            tmpStr=rm+filesToClean[i];
            system(tmpStr.c_str());
        }
    }
}





void autoWork::delete_existing_file(const std::string& fileName)
{
    string rm="rm ";
    string tmpStr="test -e "+fileName;
    int ret=system(tmpStr.c_str());
    if (ret==0) {
        tmpStr=rm+fileName;
        system(tmpStr.c_str());
    }
}





double autoWork::error_propagation(const std::vector<double>& dFi,
                                   const std::vector<double>& errori)
{
    if (dFi.size()==errori.size()) {
        double errorF=0;
        for (size_t i=0; i<dFi.size(); ++i) {
            errorF += pow(dFi[i],2)*pow(errori[i],2);
        }
        errorF=sqrt(errorF);
        return errorF;
    }
    else {
        // number of elements doesn't match
        return (-1.0);
    }
}





vector<double> autoWork::thermocalc_equ(const StructureClass& sysVar,
                                        const int n_sys,
                                        const double Temp_d,
                                        const std::string& keyWord)
{
    vector<double> stats;
    vector<double> calc_tmp;
    double value=0;
    double avgValue=0;
    SysInfo calc_stats;
    
    //======================================================
    // get the path to the equilibration log file
    //======================================================
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/equilibration/log/");
    
    int countFiles=0; // count of files used for averaging
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        file.append(sysVar.get_usic());
        file.append("_00"+to_string((long long int)n_trl));
        file.append("_"+sysVar.get_nameString(n_sys));
        file.append("_T"+to_string((long long int)Temp_d));
        file.append(".equ.log");
        //cout << i << "\n";
        
        /* input file */
        input=folder+file;
        ifstream read_input(input.c_str());
        
        if (read_input.is_open()) {
            
            ++countFiles; // +1 if the file exists
            
            string lineContent;
            string strData0;
            bool inside_thermoSection=false;
            int line_index=0;
            
            /* # of columns of data to be stored */
            int n_columns_of_data=14;
            
            // Format of thermo section
            //--------------------------
            // 0:  Step
            // 1:  Temp
            // 2:  Press
            // 3:  PotEng
            // 4:  KinEng
            // 5:  TotEng
            // 6:  E_bond
            // 7:  E_angle
            // 8:  E_pair
            // 9:  Lx
            // 10: Ly
            // 11: Lz
            // 12: Volume
            // 13: Density
            // 14: Dt
            // 15: Time
            // 16: CPU
            // 17: T/CPU
            // 18: S/CPU
            // 19: CPULeft
            //--------------------------
            int keyWordIndex=0;
            if (keyWord=="Step"){keyWordIndex=0;}
            else if (keyWord=="Temp"){keyWordIndex=1;}
            else if (keyWord=="Press"){keyWordIndex=2;}
            else if (keyWord=="PotEng"){keyWordIndex=3;}
            else if (keyWord=="KinEng"){keyWordIndex=4;}
            else if (keyWord=="TotEng"){keyWordIndex=5;}
            else if (keyWord=="E_bond"){keyWordIndex=6;}
            else if (keyWord=="E_angle"){keyWordIndex=7;}
            else if (keyWord=="E_pair"){keyWordIndex=8;}
            else if (keyWord=="Lx"){keyWordIndex=9;}
            else if (keyWord=="Ly"){keyWordIndex=10;}
            else if (keyWord=="Lz"){keyWordIndex=11;}
            else if (keyWord=="Volume"){keyWordIndex=12;}
            else if (keyWord=="Density"){keyWordIndex=13;}
            else if (keyWord=="Dt"){keyWordIndex=14;}
            else if (keyWord=="Time"){keyWordIndex=15;}
            else if (keyWord=="CPU"){keyWordIndex=16;}
            else if (keyWord=="T/CPU"){keyWordIndex=17;}
            else if (keyWord=="S/CPU"){keyWordIndex=18;}
            else if (keyWord=="CPULeft"){keyWordIndex=19;}
            
            vector<double> vecData; // used for storing line content
            
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                iss >> strData0; // the 1st deliminated string of line
                
                if (inside_thermoSection) {
                    
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0);
                        iss >> vecData[i];
                    }
                    // NOTE:
                    // offset by 1 because the first column is for strData0
                    // so basically the first item is wasted, i.e. "Step"
                    value=vecData[keyWordIndex-1];
                    
                    /* clean container content after each use */
                    vecData.clear();
                }
                
                /* This is the next line to the thermo section */
                if (strData0=="Loop") inside_thermoSection=false;
                /* This is the previous line to the thermo section */
                else if (strData0=="Step") inside_thermoSection=true;
            }
            
            /* add the "final" value of thermo over different trials of same T */
            calc_tmp.push_back(value);
            avgValue += value;
            
            read_input.close();
        }
        else {
            //cout
            //<< "thermocalc_integrated: input file cannot open." << "\n";
            //system("read -t5");exit(EXIT_FAILURE);
        }
        
        /* clear content of filename container */
        file.clear();
        
    } /* end block for looping over all trials */
    
    if (countFiles==0) return {0};
    
    /* average by the real count number which may not necessarily = n_trials */
    avgValue /= (double)(countFiles);
    calc_stats.set_calcvector(calc_tmp);
    double avg=calc_stats.calc_mean();
    if (avgValue!=avg) {
        cout
        << "thermocalc_equ: avgValue != calc_avg." << "\n"
        << "Something is wrong!"
        << "\n";
    }
    double stdev=calc_stats.get_sample_stdev();
    
    stats.push_back(Temp_d);
    stats.push_back(avg);
    stats.push_back(stdev);
    
    return stats;
}





vector<double> autoWork::thermocalc_equ_avgT(const StructureClass& sysVar,
                                             const int n_sys,
                                             const double Temp_d)
{
    vector<double> stats;
    vector<double> calc_tmp;
    double value=0;
    double avgValue=0;
    SysInfo calc_stats;
    
    //======================================================
    // get the path to the equilibration log file
    //======================================================
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/equilibration/log/");
    
    int countFiles=0; // count of files used for averaging
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        file.append(sysVar.get_usic());
        file.append("_00"+to_string((long long int)n_trl));
        file.append("_"+sysVar.get_nameString(n_sys));
        file.append("_T"+to_string((long long int)Temp_d));
        file.append(".equ.log");
        //cout << i << "\n";
        
        /* input file */
        input=folder+file;
        ifstream read_input(input.c_str());
        
        if (read_input.is_open()) {
            
            ++countFiles; // +1 if the file exists
            
            string lineContent;
            string strData0;
            bool inside_thermoSection=false;
            int line_index=0;
            
            /* # of columns of data to be stored */
            int n_columns_of_data=14;
            
            // Format of thermo section
            //--------------------------
            // 0:  Step
            // 1:  Temp
            // 2:  Press
            // 3:  PotEng
            // 4:  KinEng
            // 5:  TotEng
            // 6:  E_bond
            // 7:  E_angle
            // 8:  E_pair
            // 9:  Lx
            // 10: Ly
            // 11: Lz
            // 12: Volume
            // 13: Density
            // 14: Dt
            // 15: Time
            // 16: CPU
            // 17: T/CPU
            // 18: S/CPU
            // 19: CPULeft
            //--------------------------
            int keyWordIndex=1; // Temp
            
            vector<double> vecData; // used for storing line content
            
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                iss >> strData0; // the 1st deliminated string of line
                
                if (inside_thermoSection) {
                    
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0);
                        iss >> vecData[i];
                    }
                    // NOTE:
                    // offset by 1 because the first column is for strData0
                    value=vecData[keyWordIndex-1];
                    calc_tmp.push_back(value); // store "every" value
                    
                    /* clean container content after each use */
                    vecData.clear();
                }
                
                /* This is the next line to the thermo section */
                if (strData0=="Loop") inside_thermoSection=false;
                /* This is the previous line to the thermo section */
                else if (strData0=="Step") inside_thermoSection=true;
                // NOTE:
                // this should be put behind because it's the "next" line
                // that it's "inside" the thermo section
            }
            
            /* add the "last line" of thermo over different trials of same T */
            //avgValue += value;
            read_input.close();
        }
        
        /* clear content of filename container */
        file.clear();
        
    } /* end block for looping over all trials */
    
    if (countFiles==0) return {0};
    
    /* average by the real count number which may not necessarily = n_trials */
    avgValue /= (double)(countFiles);
    calc_stats.set_calcvector(calc_tmp);
    double avg=calc_stats.calc_mean();
    double stdev=calc_stats.get_sample_stdev();
    
    stats.push_back(Temp_d);
    stats.push_back(avg);
    stats.push_back(stdev);
    
    return stats;
}





vector<double> autoWork::thermocalc_equ_avgL(const StructureClass& sysVar,
                                             const int n_sys,
                                             const double Temp_d)
{
    
    vector<double> avgboxL;
    double Lx=0,Ly=0,Lz=0;
    double avgLx=0,avgLy=0,avgLz=0;
    
    
    //======================================================
    // get the path to the equilibration log file
    //======================================================
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/equilibration/log/");
    
    int countFiles=0; // count of files used for averaging
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        
        file.append(sysVar.get_usic());
        file.append("_00"+to_string((long long int)n_trl));
        file.append("_"+sysVar.get_nameString(n_sys));
        file.append("_T"+to_string((long long int)Temp_d));
        file.append(".equ.log");
        //cout << i << "\n";
        
        /* input file */
        input=folder+file;
        ifstream read_input(input.c_str());
        
        if (read_input.is_open()) {
            
            ++countFiles; // +1 if the file exists
            
            string lineContent;
            string strData0;
            bool inside_thermoSection=false;
            int line_index=0;
            
            /* # of columns of data to be stored */
            int n_columns_of_data=14;
            
            // Format of thermo section
            //--------------------------
            // 0:  Step
            // 1:  Temp
            // 2:  Press
            // 3:  PotEng
            // 4:  KinEng
            // 5:  TotEng
            // 6:  E_bond
            // 7:  E_angle
            // 8:  E_pair
            
            // 9:  Lx
            // 10: Ly
            // 11: Lz
            
            // 12: Volume
            // 13: Density
            // 14: Dt
            // 15: Time
            // 16: CPU
            // 17: T/CPU
            // 18: S/CPU
            // 19: CPULeft
            //--------------------------
            
            vector<double> vecData; // used for storing line content
            
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                iss >> strData0; // the 1st deliminated string of line
                
                if (inside_thermoSection) {
                    
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0);
                        iss >> vecData[i];
                    }
                    // NOTE:
                    // offset by 1 because the first column is for strData0
                    Lx=vecData[9-1];
                    Ly=vecData[10-1];
                    Lz=vecData[11-1];
                    
                    /* clean container content after each use */
                    vecData.clear();
                }
                
                /* This is the next line to the thermo section */
                if (strData0=="Loop") inside_thermoSection=false;
                /* This is the previous line to the thermo section */
                else if (strData0=="Step") inside_thermoSection=true;
                // NOTE:
                // this should be put behind because it's the "next" line
                // that it's "inside" the thermo section
            }
            
            /* add the "last line" of thermo over different trials of same T */
            avgLx += Lx;
            avgLy += Ly;
            avgLz += Lz;
            
            read_input.close();
        }
        
        /* clear content of filename container */
        file.clear();
        
    } /* end block for looping over all trials */
    
    /* average by the real count number which may not necessarily = n_trials */
    avgLx /= (double)(countFiles);
    avgLy /= (double)(countFiles);
    avgLz /= (double)(countFiles);
    
    avgboxL.push_back(avgLx);
    avgboxL.push_back(avgLy);
    avgboxL.push_back(avgLz);
    
    return avgboxL;
}





void autoWork::thermodataprocess_prd(const StructureClass& sysVar,
                                     const int n_trl,
                                     const int n_sys,
                                     const double Temp_d)
{
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/log/");
    
    /* input file */
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.log");
    input=folder+file;
    ifstream read_input(input.c_str());
    
    /* output file */
    file.clear();
    file.append("new_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.log");
    output=folder+file;
    ofstream write_output(output.c_str());
    
    if (read_input.is_open()) {
        if (write_output.is_open()) {
            
            string lineContent;
            string strData0,strData1,strData2;
            int    countStepkeyword=0;
            int    linesforuse=0;
            int    step0=0,step1=0;
            int    countaccess=0;
            
            const int n_linesforuse=2; // NOTE!
            
            // Read every line one at a time of the file
            while (getline(read_input,lineContent)) {
                
                ++linesforuse;
                
                istringstream iss(lineContent);
                iss >> strData0; // the 1st deliminated string of line
                iss >> strData1; // ""  2nd ""
                iss >> strData2; // ""  3rd ""
                
                // write out thermo values
                if ((countStepkeyword!=0) &&
                    (linesforuse<=n_linesforuse)) {
                    
                    if ((linesforuse%2)==0)      step0=stoi(strData0);
                    else if ((linesforuse%2)==1) step1=stoi(strData0);
                    
                    if ((countaccess==1)||((step0!=step1)&&(stoi(strData0)!=0))) {
                        write_output << lineContent << "\n";
                    }
                    
                }
                
                // This enters the thermo output region
                if (strData0=="Step") {
                    
                    linesforuse=0;
                    ++countaccess;
                    
                    // write out the thermo_style items only one time
                    if (countStepkeyword==0) {
                        write_output << lineContent << "\n";
                    }
                    
                    ++countStepkeyword;
                }
                
                // put a space between the boundary of 2 production blocks
                if ((countStepkeyword!=0) &&
                    (strData0=="##") &&
                    (strData1=="1.") &&
                    (strData2=="\"PSEUDO-LINEAR\"")) {
                    
                    //write_output << "\n";
                }
                
            } // ends while getline
            
            write_output.close();
        }
        
        read_input.close();
    }
    else {
        //cout
        //<< "thermoprocess_prd: input file cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
    }
    
    /* clear content of filename container */
    file.clear();
}





void autoWork::thermocalc_prd(const StructureClass& sysVar,
                              const int n_trl,
                              const int n_sys,
                              const double Temp_d)
{
    vector<double> stats;
    vector<double> calc_tmp;
    double value=0;
    SysInfo calc_stats;
    
    string keyWord;
    vector<string> allkeyWord;
    allkeyWord.push_back("Density");
    allkeyWord.push_back("Volume");
    allkeyWord.push_back("TotEng");
    allkeyWord.push_back("Lx");
    allkeyWord.push_back("Ly");
    allkeyWord.push_back("Lz");
    allkeyWord.push_back("Temp");
    
    /* input file */
    string folder,file,input,output;
    folder.append(return_SimulationFolderPath(sysVar));
    folder.append("/production/log/");
    
    file.append("new_");
    file.append(sysVar.get_usic());
    file.append("_00"+to_string((long long int)n_trl));
    file.append("_"+sysVar.get_nameString(n_sys));
    file.append("_T"+to_string((long long int)Temp_d));
    file.append(".prd.log");
    input=folder+file;
    ifstream read_input(input.c_str());
    
    auto start_position=read_input.tellg();
    
    if (read_input.is_open()) {
        
        for (int i=0; i<allkeyWord.size(); ++i) {
            
            keyWord=allkeyWord[i];
            
            /* output file */
            folder.clear();file.clear();
            
            folder.append(return_AnalysisFolderPath(sysVar));
            folder.append("/fit_data/");
            
            file.append("thermo_");
            file.append(keyWord+"_");
            file.append(sysVar.get_usic());
            file.append("_00"+to_string((long long int)n_trl));
            file.append("_"+sysVar.get_nameString(n_sys));
            file.append(".dat");
            output=folder+file;
            ofstream write_output(output.c_str(),ofstream::app); // 'append'
            
            string lineContent;
            string strData0;
            bool inside_thermoSection=false;
            int line_index=0;
            
            /* # of columns of data to be stored */
            int n_columns_of_data=14;
            
            // Format of thermo section
            //--------------------------
            // 0:  Step
            // 1:  Temp
            // 2:  Press
            // 3:  PotEng
            // 4:  KinEng
            // 5:  TotEng
            // 6:  E_bond
            // 7:  E_angle
            // 8:  E_pair
            // 9:  Lx
            // 10: Ly
            // 11: Lz
            // 12: Volume
            // 13: Density
            // 14: Dt
            // 15: Time
            // 16: CPU
            // 17: T/CPU
            // 18: S/CPU
            //--------------------------
            int keyWordIndex=0;
            if (keyWord=="Step"){keyWordIndex=0;}
            else if (keyWord=="Temp"){keyWordIndex=1;}
            else if (keyWord=="Press"){keyWordIndex=2;}
            else if (keyWord=="PotEng"){keyWordIndex=3;}
            else if (keyWord=="KinEng"){keyWordIndex=4;}
            else if (keyWord=="TotEng"){keyWordIndex=5;}
            else if (keyWord=="E_bond"){keyWordIndex=6;}
            else if (keyWord=="E_angle"){keyWordIndex=7;}
            else if (keyWord=="E_pair"){keyWordIndex=8;}
            else if (keyWord=="Lx"){keyWordIndex=9;}
            else if (keyWord=="Ly"){keyWordIndex=10;}
            else if (keyWord=="Lz"){keyWordIndex=11;}
            else if (keyWord=="Volume"){keyWordIndex=12;}
            else if (keyWord=="Density"){keyWordIndex=13;}
            else if (keyWord=="Dt"){keyWordIndex=14;}
            else if (keyWord=="Time"){keyWordIndex=15;}
            else if (keyWord=="CPU"){keyWordIndex=16;}
            else if (keyWord=="T/CPU"){keyWordIndex=17;}
            else if (keyWord=="S/CPU"){keyWordIndex=18;}
            
            vector<double> vecData; // used for storing line content
            
            read_input.seekg(start_position);
            
            while (getline(read_input,lineContent)) {
                
                ++line_index;
                
                istringstream iss(lineContent);
                iss >> strData0; // 1st deliminated string of input
                
                if (inside_thermoSection) {
                    
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0);
                        iss >> vecData[i];
                    }
                    // NOTE:
                    // offset by 1 because the first column is for strData0
                    // so basically the first item is wasted, i.e. "Step"
                    value=vecData[keyWordIndex-1];
                    calc_tmp.push_back(value);
                    
                    /* clean container content after each use */
                    vecData.clear();
                }
                
                /* This is the previous line to the thermo section */
                if (strData0=="Step") inside_thermoSection=true;
                // NOTE:
                // this should be put behind because it's the "next" line
                // that it's "inside" the thermo section
            }
            
            calc_stats.set_calcvector(calc_tmp);
            double avg=calc_stats.calc_mean();
            double stdev=calc_stats.get_sample_stdev();
            
            write_output
            << Temp_d*pow(sysVar.get_corrFacforT(),-1) << " "
            << avg << " "
            << stdev << "\n";
            
            calc_tmp.clear();
            read_input.clear();
            
        }
        
        read_input.close();
    }
    else {
        //cout
        //<< "thermocalc_integrated: input file cannot open." << "\n";
        //system("read -t5");exit(EXIT_FAILURE);
        return;
    }
    
    
    bool is_keep_thermo_prd=true;
    
    if(is_keep_thermo_prd){
        string rm="rm "+input;
        system(rm.c_str());
    }
}





/*==( Function that tests repetition by two base numbers )==*/
void autoWork::test_RepeatednessByTwoBases(const int b1, const int b2)
{
    double	base1 = (double)b1;
    double	base2 = (double)b2;
    int		index = 0;
    vector<double> series;
    ofstream File("Series.txt");
    if ( File.is_open() )
    {
        for (int i=0; i<17; ++i)
        {
            for (int j=0; j<5; ++j)
            {
                double tmp = pow(base1,(double)i) + pow(base2,(double)j);
                series.push_back(tmp);
                File << tmp << "\n";
                ++index;
            }
        }
        File.close();
        cout << "Series size = " << index << "\n";
    }
    
    // check repeat
    int tmp=0;
    for (int i=0; i<index-1; ++i)
    {
        for (int j=i+1; j<index; ++j)
            if ( series[i] == series[j] )
            {
                ++tmp;
                cout << series[i] << " repeats" << "\n";
            }
    }
    if ( tmp == 0 )
        cout << "No Repeats." << "\n";
    else
        cout << "Repeats." << "\n";
}
