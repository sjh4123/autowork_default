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

#include "functions.h"

using namespace std;
using namespace autoWork;

#define TESTAMDAT 0

const string autoWork::submit_jobs(const UserInterface * const inpVar)
{
    time_t begin_overall=time(NULL);
    
    int typeB=0;
    int typeS=0;
    //**************************************************************************
    // typeB and typeS should be one of the below listed combinations
    // Please pick one that is most similar to the system to be simulated
    //--------------------------------------------------------------------------
    // (typeB,typeS)  Targeted System
    // (5,6)          Lennard-Jones Particles
    // (0,5)          Bead-Spring FENE Chains
    // (0,0)          Polystyrene Chains
    //**************************************************************************
    if (inpVar!=nullptr) {
        typeB=inpVar->get_sysVar().get_typeB();
        typeS=inpVar->get_sysVar().get_typeS();
    }
    
    StructureClass species(typeB,typeS);
    
    /** Master Directory Path **/
    //**************************************************************************
    species.set_userName("jh148");
    species.set_masterFolder("AUTOWORK");
    //species.set_Path("/home/"+species.get_userName()+"/"+species.get_masterFolder());
    species.set_Path("/Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork");
    //**************************************************************************
    
    /** USIC = simType + (year) + usicID **/
    //**************************************************************************
    species.set_simType("Auto");
    species.set_year("2017");
    species.set_usicID("test");
    
    /** Namestring **/
    //**************************************************************************
    // naming convention: *_${usic}_${trl}_${namestring}_T${Temp}.*
    // if not specified, namestring is defualted based on typeB & typeS
    //**************************************************************************
    //species.set_nameString("namestr");
    
    /** Number of Duplicates for the Created Species **/
    //**************************************************************************
    species.set_n_trial(2);
    
    /** System Unit and Starting Temperatures **/
    //**************************************************************************
    species.set_systemUnit("real");//"real" or "metal" or "lj"
    species.set_startTemp(1500.0); //highest starting T
    species.set_crossTemp(1300.0); //crossing T
    species.set_hTres(50.0);       //resolution above crossTemp
    species.set_lTres(25.0);       //resolution below crossTemp
    
    /** Equilibration & Data Collection Time (Mapping to Each Regime) **/
    ////////////////////////////////////////////////////////////////////////////
    //**************************************************************************
    const int    n_khunsegs=1;
    const double n_equ_blocks=10;
    const double n_prd_blocks=10;
    const double n_relaxation=n_equ_blocks*pow(n_khunsegs,2);
    const double r=pow(10.0,0.5);
    //**************************************************************************
    species.set_n_kuhnsegs(n_khunsegs);
    species.set_n_equ_blocks(n_equ_blocks);
    species.set_n_prd_blocks(n_prd_blocks);
    species.set_n_relaxation(n_relaxation);
    species.set_equilibration_times
    ({
        n_relaxation*  1e+3,
        n_relaxation*r*1e+3,
        n_relaxation*  1e+4,
        n_relaxation*  1e+5,
        n_relaxation*  1e+6,
        n_relaxation*  1e+7,
        n_relaxation*  1e+8,
        n_relaxation*  1e+9
    });
    /** First 3 regimes are suggested to cover l.e.(<=) 1 decade in time
     ** Becuase T domain is much wider in highT that that in lowT regime;
     ** Overall number of regimes should be >= 3 to have good data quality
     ** real : time unit is fs
     ** lj   : time unit is tauLJ (~ps in real unit)
     ** metal: time unit is ps **/
    //**************************************************************************
    ////////////////////////////////////////////////////////////////////////////
    
    /** Set number of T's used in each regime (should = number of regimes) **/
    /** NOTE: each regime should have # of Ts >= 2 **/
    //**************************************************************************
    species.set_n_regime_temps
    ({8,8,8,8,8,8,4,4});
    
    /** Set timestep size used in each regime (should = number of regimes) **/
    //**************************************************************************
    species.set_ts_regime
    ({1,1,1,1,1,1,1,1});
    
    /** Set [start,end] regime **/
    //**************************************************************************
    species.set_regime_beg(0);//inclusive
    species.set_regime_end(0);//inclusive
    
    /** Program Module Switch **/
    ////////////////////////////////////////////////////////////////////////////
    //**************************************************************************
    species.set_is_makeFileSystem(0);  // Generate Folder Structure
    species.set_is_defaultLmpData(0);  // Generate Data File by StructureClass
    species.set_is_moltempLmpData(0);  // Generate Data File by Moltemplate
    species.set_is_LammpsPhase(0);     // Generate Default LAMMPS Scripts (OLD)
    species.set_is_use_prepScripts(0); // Use User-Prepared LAMMPS/AMDAT Scripts
    species.set_is_Simulations(0);     // Generate Simulation Bash Scripts
    species.set_is_AMDAT(0);           // Generate AMDAT Bash Scripts
    species.set_is_fitData(0);         // Perform Data Analyses
    species.set_is_fit_dielectric(0);  // Analyze Dielectric Data
    species.set_is_theorytest(0);      // Test Theories of Glass Formation
    species.set_is_cutoff_traject(0);  // Cutoff Production Trajectory
    species.set_is_alglibfit(0);       // Execute ALGLIB routine manually
    //**************************************************************************
    ////////////////////////////////////////////////////////////////////////////
    
    /** Job Submission Switch **/
    //**************************************************************************
    species.set_is_directSub(0);       // Automatic Submission
    species.set_is_watch_hold(0);      // Watch for Targets & Hold Process
    
    /** Other Switches **/
    //**************************************************************************
    //species.set_is_GPU(false);           // false: CPU only simulation
    //species.set_is_fullquench(true);     // true:  traditional linear quench
    //species.set_is_singleTempTest(true); // true:  isothermal (gen/prd/amdat)
    //species.set_which_GPU(1);            // specify GPU for the only trial
    //species.set_is_makeNewAnanlysisFolder(true); // only for test purpose
    //species.set_is_backupfitdata(true);  // only for test purpose
    
    /** System Initialiation **/
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    species.set_n_digits(3);
    // NOTE: precision
    // real units: temperature precise to "n_digits" below decimal point
    // (ex. n_digits=3 --> 1234.567; n_digits=0 --> 1234)
    // lj units: temperature precise to "n_digits+3" below decimal point
    // (ex. n_digits=3 --> 1.234567; n_digits=0 --> 1.234)
    
    /** assign user-input system data **/
    system_initialization(species,inpVar);
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////
    
    static bool is_show_usic=true;
    if (is_show_usic) {
        cout << "\nUSIC = "<<species.get_usic()<<"\n"; system_wait(1);
        is_show_usic=false;
    }
    
    /** Make Folder Structure **/
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
    if (species.get_is_makeFileSystem())
    {
        make_SimulationsFolders(species);
    }
    if (species.get_is_makeNewAnanlysisFolder())
    {
        makeNewAnanlysisFolder(species);
    }
    if (species.get_is_backupfitdata())
    {
        backupfitdata(species);
    }
    if (species.get_is_copy_from_to())
    {
        copy_from_to(species);
    }
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                       I. Data File Generation                        **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    
    /** LAMMPS datafile by StructureClass **/
    //--------------------------------------------------------------------------
    int n_sys=0,n_trl=0;
    if (species.get_is_defaultLmpData())
    {
        for (n_sys=species.get_n_sys_beg(); n_sys<=species.get_n_sys_end();
             ++n_sys)//species
        {
            for (n_trl=0; n_trl<species.get_n_trial(); ++n_trl)//trial
            {
                /** Construct Single Chain **/
                species.make_SingleChainXYZ(n_sys);
                species.make_BondFile(n_sys);
                species.make_AngleFile(n_sys);
                species.write_ForcefieldToFile();
                /** PACKMOL **/
                species.invoke_PACKMOL(n_sys);
                /** Make LAMMPS Data Files **/
                species.make_LammpsDataFile(n_sys);
                /** Change File Names (trial) **/
                move_InnerLoopFiles(species,n_sys,n_trl);
                /** Change File Names (system) **/
                move_OuterLoopFiles(species,n_sys);
            }
        } cout << "\ndefault lmpdata made!\n";system_wait(1);
    }
    else
    {
        for (n_sys=species.get_n_sys_beg(); n_sys<=species.get_n_sys_end();
             ++n_sys)//species
        {
            for (n_trl=0; n_trl<species.get_n_trial(); ++n_trl)//trial
            {
                //--------------------------------------------
                // if data file generation is skipped,
                // single chain info is still neccesary;
                // ie. it will be used in AMDAT input
                //--------------------------------------------
                species.make_SingleChainXYZ(n_sys);
                species.make_BondFile(n_sys);
                species.make_AngleFile(n_sys);
                species.write_ForcefieldToFile();
            }
        }
    } clean_tmpfiles(species);
    
    /** LAMMPS datafile by Moltemplate **/
    //--------------------------------------------------------------------------
    MoltemplateLmpData moldata(species);
    
    bool is_restrictbynAtoms=true;
    bool is_optimizeCharge=true;
    bool is_normalize_mass=true;
    bool is_correctAlkylDihedral=true;
    
    /** Degree of Polymerization **/
    int dop=1;
    
    /** Monomer Set & Monomer Composition **/
    vector<string> merSet={"S010"};
    vector<int> merComp={1};
    
    /** (Co)Polymer Type **/
    string copolymerType="random";//random,block
    
    /** Tacticity Type **/
    string tacticity="atactic";//atactic,syndiotactic,isotactic
    
    /** Program Setup Section **/
    //--------------------------------------------------------------------------
    moldata.set_dop(dop);
    moldata.set_merSet(merSet);
    moldata.set_merComp(merComp);
    moldata.set_copolymerType(copolymerType);
    moldata.set_tacticity(tacticity);
    moldata.set_is_restrictbynAtoms(is_restrictbynAtoms);
    moldata.set_is_optimizeCharge(is_optimizeCharge);
    moldata.set_is_mass_normalization(is_normalize_mass);
    moldata.set_is_correctAlkylDihedral(is_correctAlkylDihedral);
    moldata.set_n_total_atoms(15000);
    moldata.set_sequenceNum(24);
    
    species.set_chainLen(dop);//NOTE
    
    /** NOTE:
     ** if is_applyTacticity is true, sequence set and tacticity are both
     ** handled here; otherwise, final sequence set needs to be specified
     ** in the below "set_sequenceSet()" which takes vector<vector<string>> **/
    
    /** assign user-input MoltemplateLmpData data **/
    moltemplate_initialization(species,moldata,inpVar);
    
    if (false)
    {
        makeltfilefrompdb();
        //makeTacltfile();
        exit(EXIT_SUCCESS);
    }
    if (species.get_is_moltempLmpData())
    {
        if (true) { /** Create auto-gen sequenceSet with tacticity **/
            
            if (moldata.get_is_restrictbynAtoms()) {
                /** System size controlled by sequenceLen & n_total_atoms **/
                //moldata.set_n_total_atoms(15000);
                moldata.set_n_total_atoms(moldata.get_n_total_atoms());
            } else {
                /** System size controlled by sequenceLen & sequenceNum **/
                //moldata.set_sequenceNum(24);
                moldata.set_sequenceNum(moldata.get_sequenceNum());
            }
            moldata.set_sequenceLen(moldata.get_dop());
            moldata.set_chem_mono(moldata.get_merSet());
            moldata.set_mono_beads(moldata.get_merComp());
            
            //moldata.set_sequenceSet(moldata.make_atacticSequenceSet());//byVenkat
            moldata.set_sequenceSet(moldata.make_SequenceSet_adhoc());
            
        } else { /** replace by established sequenceSet **/
            moldata.set_sequenceSet(moldata.read_sequenceSetfromfile(species));
        }
        /** Use sequenceSet to make LAMMPS datafile **/
        moldata.make_LmpDataFilebyMoltemplate(species);
    }
    //--------------------------------------------------------------------------
    moldata.set_atomTypes(species);
    moldata.set_atomNumbers(species);
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                   II. LAMMPS InputFiles Generation                   **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    LmpScripts lmp(species);
    
    //lmp.set_is_resize(true);
    //lmp.set_is_fixMomentum(false);
    //lmp.set_is_fixPrdMomentum(false);
    
    /** Generate LAMMPS scripts by LmpScripts class **/
    if (species.get_is_LammpsPhase())
    {
        lmp.make_GenerationInputFile(species);
        lmp.make_QuenchInputFile(species);
        lmp.make_EquilibrationInputFile(species);
        if(lmp.get_is_tampering())lmp.make_TamperingInputFile(species);
        if(lmp.get_is_resize())lmp.make_ResizeInputFile(species);
        for(n_sys=species.get_n_sys_beg();n_sys<=species.get_n_sys_end();
            ++n_sys)lmp.make_ProductionInputFile(species,n_sys);
    } clean_tmpfiles(species);
    
    /** Use prepared scripts in 'scripts' folder **/
    if (species.get_is_use_prepScripts())
    {
        use_prepScripts(species);
    }
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                  Fit (Experimental) Dielectric Data                  **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    FitDielectric ds(species);
    
    //--------------------------------------------------------------------------
    /** Operating Procedure:
     **
     ** At master path, create a folder named dielectric, and put raw dielectric
     ** files into it.
     **
     ** Naming convention for raw dielectric files:
     ** <NIST>        form: arbitrary name ending with .dat, ex. sample.dat
     ** <NPIC>        form: species-TinC.ext, ex. 1-85C.csv
     ** <manuallySet> form: arbitrary name ending with .dat, ex. sample.dat
     **
     ** NOTE:
     ** for <NIST> and <manuallySet> forms, it only allows one raw file to be
     ** in the designated folder at a time;
     ** <NPIC> form does not have this restriction, ie. there can be multiple
     ** files in the folder for same batch of analysis
     **
     ** <form> -- format for dielectric files:
     ** "NPIC" or "NIST" or "manuallySet"
     **
     ** <func> -- fitting function:
     ** "HN_loss" (Loss Part of Havriliak-Negami)-- free shape parameters
     ** "CC_loss" (Loss Part of Cole-Cole)       -- beta=1  (symmetric  broadening)
     ** "CD_loss" (Loss Part of Cole-Davidson)   -- alpha=1 (asymmetric broadening)
     **
     ** <model> -- model for T-dependence of relaxation time:
     ** "VFT" or "COOP"
     **/
    //--------------------------------------------------------------------------
    ds.set_form("NPIC");
    ds.set_func("HN_loss");
    ds.set_model("COOP");
    ds.set_is_use_dc_free(false);
    ds.set_is_manualFitSingleDS(false);
    ds.set_sfreq_range({1.0,1.0e+5});//{min,max}; only applies to NPIC form
    
    /** assign user-input FitDielectric data **/
    fitDielectric_initialization(species,ds,inpVar);
    
    if (species.get_is_fit_dielectric())
    {
        time_t startds=time(NULL);
        vector<vector<double>> tinfo_dielectric;
        
        /** transfer source data to designated folder and rename files **/
        ds.transfer_raw_data(species);
        
        /** intialized internal values **/
        ds.initialize();
        
        /** process raw data **/
        tinfo_dielectric=ds.dielectricData_preprocess(species);
        
        if (!ds.get_is_manualFitSingleDS())
        {
            /** fit dielectric response **/
            if (ds.get_form()=="NPIC")
            {
                for (int spc=0; spc<tinfo_dielectric.size(); ++spc) {//species
                    for (int T=0; T<tinfo_dielectric.at(spc).size(); ++T) {//T
                        ds.set_spc_id({spc,(int)tinfo_dielectric.size()});
                        ds.set_T_id({T,(int)tinfo_dielectric.at(spc).size()});
                        ds.fit_dielectric_loss(species,ds.get_process_id());
                    }
                } ds.fit_tauFit(species);
            }
            else if (ds.get_form()=="NIST")
            {
                for (int process=0; process<tinfo_dielectric.size(); ++process) {
                    //if there are cooling and heating processes:
                    //process=0: 'tinfo_dielectric' has cooling temperatures
                    //process=1: 'tinfo_dielectric' has heating temperatures
                    double Temp=0;
                    for (int T=0; T<tinfo_dielectric.at(process).size(); ++T) {
                        Temp=tinfo_dielectric.at(process).at(T);
                        ds.fit_dielectric_loss(species,process,Temp);
                    }
                } ds.fit_tauFit(species);
            }
        } else {
            /** fit user-specified DS data **/
            // special case: process=2
            ds.fit_dielectric_loss(species);
        }
        time_t endds=time(NULL);
        cout << "\nTotal time elapsed: " << difftime(endds,startds) << " sec\n";
        return "\n\nDielectric Anaylsis Finished!\n\n";
    }
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    /** For Standalone ALGLIB Routines **/
    if (species.get_is_alglibfit())
    {
        alglibfit_routine(species);
    }
    
    
    
    /** END of SYS_CONFIG **/
    //--------------------------------------------------------------------------
    if (inpVar!=nullptr) {
        if (inpVar->get_sysVar().get_is_return()) return "END_SYS_CONFIG";
    }
    //--------------------------------------------------------------------------
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                 Automate Over All Simulation Regimes                 **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    int regime=0;
    for (regime=species.get_regime_beg();regime<=species.get_regime_end();++regime)
    {
        time_t begin_regime=time(NULL);
        species.set_current_regime(regime);
        
        double t_eq=species.get_equilibration_times().at(regime);
        double compu_time=0;//compu_time is used to model segmental relax time
        
        if (inpVar!=nullptr)
        {
            compu_time=t_eq/species.get_n_relaxation();
            species.set_compu_time(compu_time);
            
            if (inpVar->get_sysVar().get_n_equ_blocks_plus().size()>0) {
                species.set_n_equ_blocks_plus(inpVar->get_sysVar().get_n_equ_blocks_plus());
                species.set_n_relaxation_plus(inpVar->get_sysVar().get_n_relaxation_plus());
                species.set_n_equ_blocks(species.get_n_equ_blocks_plus().at(regime));
                compu_time=t_eq/species.get_n_relaxation_plus().at(regime);
                species.set_compu_time(compu_time);
            }
            if (inpVar->get_sysVar().get_n_prd_blocks_plus().size()>0) {
                species.set_n_prd_blocks_plus(inpVar->get_sysVar().get_n_prd_blocks_plus());
                species.set_n_prd_blocks(species.get_n_prd_blocks_plus().at(regime));
            }
        } else {
            compu_time=t_eq/species.get_n_relaxation();
            species.set_compu_time(compu_time);
        }
        
        //----------------------------------------------------------------------
        /** Automated_batch() with a retry algorithm
         ** Algorithm implemented in check_is_retry() in class FitData **/
        //----------------------------------------------------------------------
        int  times_retry=0;
        bool is_retry=false;
        species.set_n_cum_taueq(0);
        species.set_max_times_retry(2);
        do {
            species.set_times_retry(times_retry);
            autoWork::automated_batch(species,lmp,is_retry,inpVar);
            if (is_retry) {
                ++times_retry;
                if (times_retry>species.get_max_times_retry()) {
                    cout << "\n\n"
                    << "Max time of retry at Regime_"<< species.get_current_regime()<<" "
                    << "has been reached. Automation terminated.\n";
                    system("date");cout<<"\n";exit(EXIT_FAILURE);
                } else {
                    cout << "\n\n"
                    << "retry Regime_" << species.get_current_regime() << "\n"
                    << "#" << times_retry << " retry\n";
                    system("date");cout<<"\n";
                }
            }
        } while (is_retry);
        species.get_simulation_retry().push_back({regime,times_retry});
        //----------------------------------------------------------------------
        time_t end_regime=time(NULL);
        double timediff=difftime(end_regime,begin_regime);
        species.get_elapsedtime().at(regime-species.get_regime_beg())=timediff;
        cout << "\nRegime_" << regime << " Completed.\n\n";
        system("date");
        cout << "\nTime Elapsed: " << timediff << " seconds.\n\n";
    }
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    
    
    /** write progress info of final regime to file **/
    //--------------------------------------------------------------------------
    if (species.get_regime_beg()==0) autoWork::write_progress(species);
    
    time_t end_overall=time(NULL);
    double timediff_overall=difftime(end_overall,begin_overall);
    species.get_elapsedtime().push_back(timediff_overall);
    species.set_timediff_overall(timediff_overall);
    
    /** write final report to file **/
    //--------------------------------------------------------------------------
    autoWork::write_report(species);
    
    /** return the path of the report file **/
    //--------------------------------------------------------------------------
    return species.get_sys_targets();
}










////////////////////////////////////////////////////////////////////////////////
//==============================================================================
void autoWork::automated_batch(StructureClass& species,
                               const LmpScripts& lmp,
                               bool& is_retry,
                               const UserInterface * const inpVar)
{
    const int n_trial    = species.get_n_trial();
    const int n_system   = species.get_n_system();
    const int n_sys_beg  = species.get_n_sys_beg();
    const int n_sys_end  = species.get_n_sys_end();
    const int regime_end = species.get_regime_end();
    const int current_regime = species.get_current_regime();
    
    vector<vector<double>> tinfo;
    if (species.get_is_fullquench()) {
        tinfo = species.get_temperaturesInfo();
    } else {
        tinfo = species.get_quenchingTs();
    }
    /** true: update "temperaturesInfo" & "quenchingTs" containers **/
    species.set_is_updatetinfo(true);
    species.set_max_initialized_regime(species.get_n_regime()-1);//0-based
    
    if (inpVar!=nullptr)
    {
        if (inpVar->get_sysVar().get_inptinfo().size()!=0)
        {
            int n_regimes=species.get_n_regime();
            int n_regimes_initialized=(int)inpVar->get_sysVar().get_inptinfo().size();
            species.set_max_initialized_regime(n_regimes_initialized-1);//0-based
            if (current_regime<=species.get_max_initialized_regime()) {
                tinfo=inpVar->get_sysVar().get_inptinfo().at(current_regime);
                species.get_temperaturesInfo()=tinfo;
                species.get_quenchingTs()=tinfo;
            }
            if (n_regimes_initialized==n_regimes) {
                species.set_is_updatetinfo(false);
            } else {
                if (current_regime<species.get_max_initialized_regime()) {
                    //if number of initialized regimes of tinfo is smaller than
                    //actual n_regime, that means tinfo for the regime next to
                    //the max initialized regime has to be determined by PreSQ;
                    //therefore, all regimes prior to the max initialized regime
                    //shuold NOT update tinfo, but only at and after max regime.
                    species.set_is_updatetinfo(false);
                } else {
                    species.set_is_updatetinfo(true);
                }
            }
        }
    }
    
    //#define MANUALSET
#ifdef MANUALSET
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
    /**                         III. Do Simulations                          **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    WorkScripts ws(species);
    
    /** Generation (@highestT) **/
    //--------------------------------------------------------------------------
    ws.set_steps_gen(1e+6);
    
    /** Quench **/
    //--------------------------------------------------------------------------
    // unit of quench rate:
    // "real"  'kelvins per nano second' (default: 1000)
    // "metal" 'kelvins per nano second' (default: 1000)
    // "lj"    'T per tau' (default: 1e-4)
    //ws.set_quenchRate(1e-4);
    
    /**==( Setup Compute Nodes )==**/
    
    /** Compute Node-0 **/
    //--------------------------------------------------------------------------
    ws.set_is_node0(0);
    ws.set_node0({14,15});
    
    /** Compute Node-1 **/
    //--------------------------------------------------------------------------
    ws.set_is_node1(0);
    ws.set_node1({8,9});
    
    /** Compute Node-2 **/
    //--------------------------------------------------------------------------
    ws.set_is_node2(0);
    ws.set_node2({0,1});
    
    /** Job Allocation **/
    //--------------------------------------------------------------------------
    int n=8;//# of CPUs assigned to a batch (trial) of jobs
    int c=0;//# of CPUs occupied by a job (therefore, c controls GPU/CPU ratio)
    if (ws.get_is_node2()) n=6;
    if (inpVar!=nullptr) if (inpVar->get_ws().get_is_node2()) n=6;
    coreAllocation(species,n,c);
    //--------------------------------------------------------------------------
    // order of elements corresponds to simulation stages of:
    // 0: generation    (1cpu_1gpu)
    // 1: quench        (1cpu_1gpu)
    // 2: equilibration (1cpu_(1/n)gpu)
    // 3: resize        (1cpu_(1/n)gpu)
    // 4: production    (1cpu_(1/n)gpu)
    //--------------------------------------------------------------------------
    ws.set_numCores({1,1,1,1,1}); //number of requsted CPUs per job
    ws.set_run_cores({1,1,1,1,1});//number of CPUs to run a job on
    ws.set_cores({n,n,c,c,c});    //number of CPUs occupied by a job
    ws.set_priority({0,0,0,0,0}); //priority for each simulation phase
    //==========================================================================
    ////////////////////////////////////////////////////////////////////////////
    
    if (species.get_regime_beg()==0) autoWork::write_progress(species);
    
    int sim_restart=0;
    int gen_restart=0,qch_restart=0,equ_restart=0,res_restart=0,prd_restart=0;
    if (inpVar!=nullptr) {
        sim_restart=inpVar->get_sysVar().get_sim_restart();
        gen_restart=inpVar->get_sysVar().get_gen_restart();
        qch_restart=inpVar->get_sysVar().get_qch_restart();
        equ_restart=inpVar->get_sysVar().get_equ_restart();
        res_restart=inpVar->get_sysVar().get_res_restart();
        prd_restart=inpVar->get_sysVar().get_prd_restart();
    }
    
    /** assign user-input WorkScripts data **/
    simulation_initialization(species,ws,inpVar);
    
    if (current_regime>=sim_restart)
    {
        if(species.get_is_Simulations())
        {
            vector<string> work_targets;
            
            /** Simulation Phases: [generation, quench, eequilibration] **/
            //------------------------------------------------------------------
            for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    species.set_n_Temp((int)(tinfo.at(index).size()));
                    
                    /** Generation subfile (only in first regime) **/
                    if ((current_regime==0)&&(is_retry==false)) {//NOTE:is_retry
                        if (current_regime>=gen_restart) {
                            ws.make_GenerationSubFiles
                            (species,lmp,n_trl,n_sys,tinfo.at(index).at(0));
                        }
                    }
                    if (!species.get_is_singleTempTest()) {
                        /** Quench subfile **/
                        if (species.get_is_fullquench()) {
                            if (current_regime==0) {
                                ws.make_QuenchSubFiles(species,n_trl,n_sys);}
                        } else {
                            if (current_regime>=qch_restart) {
                                ws.make_QuenchSubFiles(species,n_trl,n_sys);
                            }
                        }
                        /** Equilibration subfile **/
                        if (current_regime>=equ_restart) {
                            for (int T=0; T<species.get_n_Temp(); ++T) {
                                ws.make_EquilibrationSubfiles
                                (species,lmp,n_trl,n_sys,tinfo.at(index).at(T));
                                /** target files: equilibration restart files **/
                                work_targets.push_back
                                (lmp.work_target(species,"equ",n_trl,n_sys,tinfo.at(index).at(T)));
                            }
                        }
                    }
                }
            }
            /** Jobs Submission: [generation, quench, eequilibration] **/
            //------------------------------------------------------------------
            /** Generation Jobs Submissions **/
            if ((current_regime==0)&&(is_retry==false)) { // NOTE: is_retry
                if (current_regime>=gen_restart) {
                    ws.make_GenerationSubScripts(species);
                }
            }
            if (!species.get_is_singleTempTest()) {
                /** Quench Jobs Submissions **/
                if (species.get_is_fullquench()) {
                    if (current_regime==0) ws.make_QuenchSubScripts(species);
                } else {
                    if (current_regime>=qch_restart) {
                        ws.make_QuenchSubScripts(species);
                    }
                }
                /** Equilibration Jobs Submissions **/
                if (current_regime>=equ_restart) {
                    ws.make_EquilibraitionSubScripts(species,tinfo);
                    /** Watch and Hold **/
                    ////////////////////////////////////////////////////////////
                    //----------------------------------------------------------
                    if (species.get_is_watch_hold()) {
                        autoWork::watch_hold
                        (species,work_targets,"simulation",current_regime);
                    } work_targets.clear();
                    //----------------------------------------------------------
                    ////////////////////////////////////////////////////////////
                }
            }
            /** Simulation Phases: [resize, production] **/
            //------------------------------------------------------------------
            for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys){
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    species.set_n_Temp((int)(tinfo.at(index).size()));
                    for (int T=0; T<species.get_n_Temp(); ++T) {
                        /** resize subfile **/
                        if (!species.get_is_singleTempTest()) {
                            if (lmp.get_is_resize()) {
                                if (current_regime>=res_restart) {
                                    ws.make_ResizeSubfiles
                                    (species,lmp,n_trl,n_sys,tinfo.at(index).at(T));
                                }
                            }
                        }
                        /** Production subfile **/
                        if (current_regime>=prd_restart) {
                            ws.make_ProductionSubFiles
                            (species,lmp,n_trl,n_sys,tinfo.at(index).at(T));
                            /** target files: production restart files **/
                            work_targets.push_back
                            (lmp.work_target(species,"prd",n_trl,n_sys,tinfo.at(index).at(T)));
                        }
                    }
                }
            }
            /** Jobs Submission: [resize, production] **/
            //------------------------------------------------------------------
            if (!species.get_is_singleTempTest()) {
                if (current_regime>=res_restart) {
                    if (lmp.get_is_resize()) ws.make_ResizeSubScripts(species,tinfo);
                }
            }
            if (current_regime>=prd_restart) {
                ws.make_ProductionSubScripts(species,tinfo);
                /** Watch and Hold **/
                ////////////////////////////////////////////////////////////////
                //--------------------------------------------------------------
                if (species.get_is_watch_hold()) {
                    autoWork::watch_hold
                    (species,work_targets,"simulation",current_regime);
                } work_targets.clear();
                //--------------------------------------------------------------
                ////////////////////////////////////////////////////////////////
            }
            
            if (species.get_is_aging())
            {
                backup_simulation_traj(species,"production");
            }
        }
    } /* end current_regime>=sim_restart */
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                         IV. AMDAT Analyses                           **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    AmdatAnalysis aa(species,ws);
    
    /** Class Setters **/
    //--------------------------------------------------------------------------
    aa.set_is_mbodies(false);          // multibody mapping
    aa.set_is_monoStruct(true);        // whether use monomer bank structure
    aa.set_is_changeAtomType(false);   // change atom types in custom file
    aa.set_is_recoverAtomType(false);  // recover modified atom types
    aa.set_is_changeFrameSteps(false); // change frame timesteps in prd xyz
    aa.set_is_changeAtomOrder(false);  // change orders of atoms in prd xyz
    aa.set_analysispart("all");        // analysis target (ex. "all")
    aa.set_relaxation_target("isfs");  // which relaxation: "isfs" or "baf"
    
    /** setup AMDAT jobs allocation configuration **/
    //--------------------------------------------------------------------------
    int c_amdat=1;
    aa.set_amdat_numCores(c_amdat);
    aa.set_amdat_run_cores(1);
    aa.set_amdat_priority(0);
    
    int ana_restart=0;
    if (inpVar!=nullptr) {
        ana_restart=inpVar->get_sysVar().get_ana_restart();
    }
    
    /** assign user-input AmdatAnalysis data **/
    analysis_initialization(species,aa,inpVar);
    
    if (current_regime>=ana_restart)
    {
        if(species.get_is_AMDAT())
        {
            /** use system_np or system_nv **/
            //------------------------------------------------------------------
            if (lmp.get_is_resize()) aa.set_is_NPT(false);
            
            /** typical analyses **/
            //------------------------------------------------------------------
            aa.set_is_strFac(true);
            aa.set_is_msd(true);
            aa.set_is_ngp(true);
            aa.set_is_isfs(true);
            //aa.set_is_baf(false);
            //aa.set_is_rdf(true);
            
            /** other analyses **/
            //------------------------------------------------------------------
            //aa.set_is_baf(false);
            //aa.set_is_composition(true);
            //aa.set_is_u2dist(true);
            //aa.set_is_stiffness_dist(true);
            //aa.set_is_isf(false);
            aa.set_is_strings(false);
            
            /** check and make corresponding statistics folders **/
            //------------------------------------------------------------------
            test_AnalysisFolder(species,aa);
            
            /* test section */
            if (TESTAMDAT)
            {
                aa.set_is_mbodies(0);
                aa.set_is_monoStruct(0);
                aa.set_is_strings(1);
                aa.make_sigmamatrix(species);
                aa.make_amdatInputFile(species,ws);
                aa.set_analysispart("all");
                aa.set_is_use_voroNeighbors(1);
                aa.make_amdatInputFile_strings_loop
                (species,ws,0,0,tinfo.at(0).at(0),80,100);
                exit(EXIT_SUCCESS);
            }
            
            /** make AMDAT inputfile **/
            //------------------------------------------------------------------
            if (species.get_is_amdatinp()) aa.make_amdatInputFile(species,ws);
            
            vector<string> amdat_targets;
            for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                    int index=n_trl*n_system+(n_sys-n_sys_beg);
                    species.set_n_Temp((int)(tinfo.at(index).size()));
                    
                    /** to get wavenumber (q) of the first peak of strFac
                     ** true:  automatic determination by internal function
                     ** false: q value then needs to be manually specified **/
                    aa.set_is_auto_qvalue(false);
                    //double TstrFac=species.get_cutTforArrhenius();
                    //aa.find_waveindex(species,TstrFac);
                    
                    for (int T=0; T<species.get_n_Temp(); ++T) {
                        
                        /** delete 0th frames from 2nd to final block (if any) **/
                        if (ws.get_divide_prd()>1) {
                            aa.delete_zerothFrames
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        /** change frame timesteps in custom file (less useful) **/
                        if (aa.get_is_changeFrameSteps()) {
                            aa.change_prdcustomfilesteps
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        /** change atom type in production custom file **
                         ** NEED to specify the types to be changed below **/
                        //species.set_backboneTypes({1,3});
                        //species.set_sideGroupTypes({2});
                        species.set_backboneTypes(species.get_types_heavy());
                        species.set_sideGroupTypes(species.get_types_light());
                        
                        if (aa.get_is_changeAtomType()) {
                            aa.set_is_keep_new_trajectory(false);
                            aa.change_productionAtomTypes
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        /** recover atom types in custom file to their original **/
                        if(aa.get_is_recoverAtomType()) {
                            aa.recover_productionAtomTypes
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        /** change the order of atoms in custom file **/
                        if (aa.get_is_changeAtomOrder()) {
                            aa.change_atom_order
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                        
                        /** make AMDAT submission files **/
                        //------------------------------------------------------
                        aa.make_amdatSubFile
                        (species,ws,n_trl,n_sys,tinfo.at(index).at(T));
                        
                        /** target files of AMDAT **/
                        amdat_targets.push_back
                        (aa.amdat_target(species,n_trl,n_sys,tinfo.at(index).at(T),
                                         aa.get_relaxation_target()));
                    }
                }
            }
            /** AMDAT Jobs Submissions **/
            //------------------------------------------------------------------
            aa.make_amdatSubScript(species,tinfo);
            
            /** Watch and Hold **/
            ////////////////////////////////////////////////////////////////////
            //------------------------------------------------------------------
            if (species.get_is_watch_hold()) {
                autoWork::watch_hold
                (species,amdat_targets,"analysis",current_regime);
            } amdat_targets.clear();
            //------------------------------------------------------------------
            ////////////////////////////////////////////////////////////////////
            
            if (species.get_is_aging())
            {
                backup_analysis_stats(species);
            }
        }
    } /* end current_regime>=ana_restart */
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                      V. Fit Data (with ALGLIB)                       **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    FitData fd(species,ws,aa);
    
    /** Curve-Fitting Switches **/
    //--------------------------------------------------------------------------
    fd.set_is_use_FG(true);
    fd.set_is_fit_Fs_by_spline(false);
    fd.set_is_fit_full_alpha(true);
    fd.set_is_use_gammafunc(false);
    
    /** Fit Relaxation Profile by Stretched Exponential Function **/
    //--------------------------------------------------------------------------
    fd.set_is_fit_sExp(true);         // true: fit KWW function
    
    /** Fit Temperature Dependence of Relaxation Time **/
    //--------------------------------------------------------------------------
    fd.set_is_fit_Arrhenius(true);    // true: fit Arrhemius function
    fd.set_is_fit_tauFit(true);       // true: fit the chosen relaxatiom model
    
    /** Customized Analysis Switches **/
    //--------------------------------------------------------------------------
    fd.set_is_applytauFitcut(true);   // true: apply current taucut for data
    fd.set_is_use_KWWassist(true);    // true: use KWW coeffs from previous fit
    fd.set_is_find_DWF(true);         // true: find the Debye-Waller factor
    fd.set_is_find_NGP(true);         // true: find the Non-Gaussian Parameter
    fd.set_is_calc_thermoData(false); // true: calculate average thermo data
    
    /** How the onset temperature of G.F. (TA) is defined **/
    /*--------------------------------------------------------------------------
     1. "ArrFit" -- X% deviation from highT Arrhenius fit (def. X=10)
     2. "ArrNor" -- X% deviation from normalized Arrhenius relation (def. X=15)
     3. "MsdCag" -- onset temperature where caging emerges from msd
     4. "FsqCag" -- onset temperature where caging emerges from isfs
     -------------------------------------------------------------------------*/
    fd.set_definitionofTA("ArrNor");
    
    /** Fitting Models **/
    /*--------------------------------------------------------------------------
     Models for Fitting Different Relaxation Functions
     ---------------------------------------------------------------------------
     KWW            F = A*exp(-(t/tau)^beta)
     KWW_full       F = exp(lnA*(1-exp(-(t/tauf)^(2p)))^(1/p)-(t/tau)^beta)
     KWW_pwr        F = A*(t^-a)*exp(-(t/tau)^beta)
     mKWW           F = (1-A)*exp(-t/taulib)+A*exp(-(t/tau)^beta)
     mKWW_pwr       F = (1-A)*exp(-t/taulib)+A*(t^-a)*exp(-(t/tau)^beta)
     
     ---------------------------------------------------------------------------
     Models for Fitting Temperature Dependence of Relaxation Time:
     ---------------------------------------------------------------------------
     <GLMVFT>       tau = tau0*exp(A*[(2*T0)/(T-T0)]^(alpha/2))
     GLMVFT1        param: alpha
     GLMVFT3        param: tau0, T0, alpha
     GLMVFT4        param: tau0, A, T0, alpha
     
     3-PARAMETER
     -------------
     AM             tau = tau0*exp((p/T)^alpha)
     DG             tau = tau0*exp((A/T)*exp(-B*T))
     VFT            tau = tau0*exp((D*T0)/(T-T0))
     Mauro          tau = tau0*exp((K/T)exp(C/T))
     
     4-PARAMETER
     -------------
     CG             tau = tau0*exp(B/{(T-T0)+[(T-T0)^2+C*T]^0.5})
     COOP           tau = tau0*exp(E_inf*(1+exp(-u*((T/E_inf)-b)))/T)
     DEAG           tau = tau0*exp(((A-B*T)/T)*exp(C/T))
     FourParamVFT   tau = tau0*exp([(D*T0)/(T-T0)]^alpha)
     
     OTHERS:
     ------------
     MCT            tau = A*((T-Tc)/Tc)^(-r)
     SOU            tau = A*((T-Tc)/T)^(-r)
     ArrheniusII    tau = tau0*exp((A/T)+(B/T^2))
     ArrheniusIII   tau = tau0*exp((A/T)+(B/T^2)+(C/T^3))
     -------------------------------------------------------------------------*/
    
    /** Model for Fitting Correlation Function, ex. Fs(q,t) **/
    //--------------------------------------------------------------------------
    fd.set_fcorr_model("KWW");
    
    /** Model for Fitting T Dependence of Relaxation Time **/
    //--------------------------------------------------------------------------
    fd.set_extrp_model("COOP");
    
    /** Model for Predicting New Temperatures **/
    //--------------------------------------------------------------------------
    if (current_regime==0) {
        fd.set_presq_model("VFT");
    } else {
        fd.set_presq_model("VFT");
    }
    
    int fd_res=0;
    if (inpVar!=nullptr) {
        fd_res=inpVar->get_sysVar().get_fd_res();
    }
    
    /** assign user-input FitData data **/
    fd_initialization(species,fd,inpVar);
    
    if (current_regime>=fd_res)
    {
        if (species.get_is_fitData())
        {
            /** clear content of containers used in FitData class **/
            //******************************************************************
            species.get_Tg_extrp().at(current_regime).clear();
            species.get_Tg_compu().at(current_regime).clear();
            species.get_m_extrp().at(current_regime).clear();
            species.get_m_compu().at(current_regime).clear();
            species.get_Tg_VFT().at(current_regime).clear();
            species.get_Tg_COOP().at(current_regime).clear();
            species.get_m_VFT().at(current_regime).clear();
            species.get_m_COOP().at(current_regime).clear();
            //******************************************************************
            if (fd.get_is_fit_sExp())
            {
                string fcorr=fd.get_fcorr_model();
                if (fd.get_is_fit_lnFs()) fcorr+="_lnFs";
                
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                    {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        species.set_n_Temp((int)(tinfo.at(index).size()));
                        
                        /** clear content of contianers for new values **/
                        //******************************************************
                        species.get_equilibratedTs().at(index).clear();
                        species.get_tauFit().at(index).clear();
                        species.get_tauEqu().at(index).clear();
                        //******************************************************
                        for (int T=0; T<species.get_n_Temp(); ++T)
                        {
                            /** fit relaxation profile to KWW function **/
                            fd.fit_sExp(species,n_trl,n_sys,tinfo.at(index).at(T),fcorr);
                        }
                        /** write temperatures to file **/
                        //------------------------------------------------------
                        fd.write_qchTs(species,n_trl,n_sys); // simulated T's
                        fd.write_equTs(species,n_trl,n_sys); // in-equilibrium T's
                        
                        /** write averaged sExp params to file **/
                        //------------------------------------------------------
                        if (current_regime==regime_end)
                        {
                            fd.write_sExp_params_avg(species,n_sys,fcorr);
                        }
                    }
                }
                /** check if current regime needs a retry **/
                //--------------------------------------------------------------
                if (fd.check_is_retry(species)) {is_retry=true; return;}
                else is_retry=false;
                //--------------------------------------------------------------
            }
            if (fd.get_is_fit_Arrhenius())
            {
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                    {
                        fd.fit_Arrhenius(species,n_trl,n_sys);
                    }
                }
            }
            if (fd.get_is_find_DWF())
            {
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        species.set_n_Temp((int)(tinfo.at(index).size()));
                        for (int T=0; T<species.get_n_Temp(); ++T)
                        {
                            fd.write_DWF
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                            fd.write_diffusion_coeff
                            (species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                    }
                }
                if (current_regime==regime_end)
                {
                    fd.write_MSDt_equ(species);
                    fd.write_MSDt_equ_avg(species);
                    fd.write_DWF_equ(species);
                    fd.write_DWF_equ_avg(species);
                }
            }
            if (fd.get_is_find_NGP())
            {
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                        int index=n_trl*n_system+(n_sys-n_sys_beg);
                        species.set_n_Temp((int)(tinfo.at(index).size()));
                        for (int T=0; T<species.get_n_Temp(); ++T)
                        {
                            fd.write_NGP(species,n_trl,n_sys,tinfo.at(index).at(T));
                        }
                    }
                }
                if (current_regime==regime_end)
                {
                    fd.write_avg_NGP(species);
                }
            }
            if (fd.get_is_normalModes())
            {
                if (aa.get_relaxation_target()=="isfs")
                {
                    if (fd.get_is_qspectrum())
                    {
                        if (aa.get_analysispart()!="mode12")
                        {
                            if (current_regime==regime_end) {
                                fd.write_avg_qmatch(species);
                                fd.write_avg_qmatch_KWW(species);
                                fd.write_smo_qmatch_KWW(species);
                            }
                        }
                        /** local min/max of bKWW vs T **/
                        if (current_regime==regime_end) {
                            fd.write_bKWW_vs_T(species);
                        }
                    } else {
                        if (current_regime==regime_end) {
                            fd.write_avg_normal_modes_decoupling(species);
                        }
                    }
                }
            }
            if (fd.get_is_fit_tauFit())
            {
                /** fit relaxation to chosen relaxation model **/
                string extrp=fd.get_extrp_model();
                string presq=fd.get_presq_model();
                for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                    for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys)
                    {
                        fd.fit_tauFit(species,n_trl,n_sys,extrp,presq);
                    }
                }
            }
            if (species.get_is_GPU())
            {
                if (fd.get_is_calc_thermoData())
                {
                    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
                        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                            int index=n_trl*n_system+(n_sys-n_sys_beg);
                            species.set_n_Temp((int)(tinfo.at(index).size()));
                            for (int T=0; T<species.get_n_Temp(); ++T)
                            {
                                ws.thermodataprocess_prd
                                (species,n_trl,n_sys,tinfo.at(index).at(T));
                                ws.thermocalc_prd
                                (species,n_trl,n_sys,tinfo.at(index).at(T));
                            }
                        }
                    }
                }
            }
        }
    }/* if (current_regime>=fd_res) */
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**                          VI. Theory Test                             **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    TheoryTest tt(species,ws,aa,fd);
    
    int tt_res=0;
    if (inpVar!=nullptr) {
        tt_res=inpVar->get_sysVar().get_tt_res();
    }
    
    /** assign user-input TheoryTest data **/
    tt_initialization(species,tt,inpVar);
    
    if (current_regime>=tt_res)
    {
        if(species.get_is_theorytest())
        {
            if (tt.get_is_AGtest()) {
                tt.amdatstringanalysis(species,ws,aa,tinfo);
            }
            if (current_regime==regime_end) {
                if (tt.get_is_AGtest()) tt.AGtest(species);
                if (tt.get_is_GLMtest()) tt.GLMtest(species);
                if (tt.get_is_HWtest()) tt.HWtest(species);
                if (tt.get_is_Leporinitest()) tt.Leporinitest(species);
            }
        }
    }/* if (current_regime>=tt_res) */
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //                                                                        //
    /**      Cutoff Broken Production (Exponential Timestep) Trajectory      **/
    //                                                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //==========================================================================
    if (species.get_is_cutoff_traject())
    {
        check_finished_blocks(species,ws,tinfo);
        exit(EXIT_SUCCESS);
        int cutoffBlock=5;//keep 1st to (cutoffBlock-1)th block
        cutoff_trajectory(species,ws,tinfo,cutoffBlock);
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
    
    /** default wait time = 100s **/
    int waitsec=100;
    
    string file,test;
    do {
        counter=0;
        /** see if target files exist **/
        for (size_t i=0; i<targets.size(); ++i) {
            ifstream tarfile(targets.at(i).c_str());
            //cout << targets.at(i).c_str() << "\n";
            if (tarfile.good()) ++counter;
        } file.clear(); test.clear();
        /** check whether core files are generated **/
        if (current_state=="simulation") {
            file.append(return_SimulationFolderPath(sysVar));
            file.append("/core.*");
            test.append("for file in "+file+"; ");
            test.append("do test -e $file; done");
        } else if (current_state=="analysis") {
            file.append(return_AnalysisFolderPath(sysVar));
            file.append("/core.*");
            test.append("for file in "+file+"; ");
            test.append("do test -e $file; done");
        }
        /** if there're core files, then delete the jobs and quit automation **/
        if (system(test.c_str())==0) {
            cout
            << "\ncore.* files detected and now deleted. Automation stopped.\n";
            system("date");
            string rm="rm -rf "+file;
            system(rm.c_str());
            string qdel;
            qdel.append("qdel ");
            qdel.append("*_"+sysVar.get_usic()+"_*");
            system(qdel.c_str());
            system("sleep 5");exit(EXIT_FAILURE);
        }
        if (counter<n_total) {
            /** Every 500 waitTime's, show USIC and automation section **/
            if ((display%500)==0) {
                cout << "\nRegime ["
                << (current_regime+1) << "/" << sysVar.get_n_regime() << "]: "
                << "in " << current_state
                << "\n";
                cout << "\nUSIC = " << sysVar.get_usic()
                << "\n";
            }
            /** Every 10 waitTime's, show how many files found **/
            if ((display%10)==0) {
                cout << "(" << counter << "/" << n_total << ") files found;\n";
            } system_wait(waitsec);
        } ++display;
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
    
    ofstream progress(o.c_str());
    if (progress.is_open())
    {        
        progress
        << "Start       " << sysVar.get_startdatetime()<< "\n"
        << "Final_Write " << return_datetime()<< "\n\n";
        
        progress
        << "tau_alpha Cutoffs ";
        if (sysVar.get_systemUnit()=="real") {
            progress << "(timeUnit: fs)" << "\n";
        } else if (sysVar.get_systemUnit()=="metal") {
            progress << "(timeUnit: ps)" << "\n";
        } else if (sysVar.get_systemUnit()=="lj") {
            progress << "(timeUnit: tau)" << "\n";
        }
        progress
        << "------------------------------------------"
        << "\n";
        for (int i=0; i<sysVar.get_n_regime(); ++i) {
            double previous_teq,current_teq,teq;
            if (true) {
                teq = sysVar.get_equilibration_times()[i];
            } else {
                if (i>0) {
                    previous_teq = sysVar.get_equilibration_times()[i-1];
                    current_teq  = sysVar.get_equilibration_times()[i];
                    teq = max(previous_teq,current_teq);
                } else {
                    teq = sysVar.get_equilibration_times()[i];
                }
            }
            if (sysVar.get_n_equ_blocks_plus().size()>0) {
                progress << "Regime_"<<i<<"  "
                <<teq/sysVar.get_n_equ_blocks_plus().at(i)<< "\n";
            } else {
                progress << "Regime_"<<i<<"  "
                <<teq/sysVar.get_n_equ_blocks()<< "\n";
            }
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
            } progress << "\n";
        }
        if (sysVar.get_is_node1_info()) {
            progress
            << "Compute_1: ";
            for (size_t i=0; i<sysVar.get_node1_info().size(); ++i) {
                progress
                << sysVar.get_node1_info()[i] << " ";
            } progress << "\n";
        }
        if (sysVar.get_is_node2_info()) {
            progress
            << "Compute_2: ";
            for (size_t i=0; i<sysVar.get_node2_info().size(); ++i) {
                progress
                << sysVar.get_node2_info()[i] << " ";
            } progress << "\n";
        } progress << "\n";
        
        vector<double> teq=sysVar.get_equilibration_times();
        int current_regime=sysVar.get_current_regime();
        
        progress
        << "System Unit:    " << sysVar.get_systemUnit() << "\n"
        << "Quench Rate:    " << sysVar.get_quenchRate() << " units\n";
        
        double n_equ_blocks=sysVar.get_n_equ_blocks();
        double n_prd_blocks=sysVar.get_n_prd_blocks();
        if (sysVar.get_n_equ_blocks_plus().size()>0) {
            progress << "Equilibration:  ";
            for (int i=0;i<(int)sysVar.get_n_equ_blocks_plus().size();++i) {
                progress << sysVar.get_n_equ_blocks_plus().at(i) << " ";
            } progress << "tau_alpha(s)\n";
        } else {
            progress << "Equilibration:  "<<n_equ_blocks<<" tau_alpha(s)\n";
        }
        if (sysVar.get_n_prd_blocks_plus().size()>0) {
            progress << "Production:     ";
            for (int i=0;i<(int)sysVar.get_n_prd_blocks_plus().size();++i) {
                progress << sysVar.get_n_prd_blocks_plus().at(i) << " ";
            } progress << "block(s)\n\n";
        } else {
            progress << "Production:     "<<n_prd_blocks<<" block(s)\n\n";
        }
        
        progress
        << "Time Elapsed up to Current Regime         " << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << "  "
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
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_Tg_extrp()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_Tg_extrp()[i]);
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
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_m_extrp()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_m_extrp()[i]);
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
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_Tg_compu()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_Tg_compu()[i]);
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
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_m_compu()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_m_compu()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            progress
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        
        progress
        << "\n"
        << "Tg_VFT(trials) avg(Tg) stdev(Tg) stdev/avg   " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_Tg_VFT()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_Tg_VFT()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            progress
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        progress
        << "\n"
        << "m_VFT(trials) avg(m) stdev(m) stdev/avg      " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_m_VFT()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_m_VFT()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            progress
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
        
        progress
        << "\n"
        << "Tg_COOP(trials) avg(Tg) stdev(Tg) stdev/avg   " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_Tg_COOP()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_Tg_COOP()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            progress
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        }
        progress
        << "\n"
        << "m_COOP(trials) avg(m) stdev(m) stdev/avg      " << "\n"
        << "---------------------------------------------" << "\n";
        for (int i=0; i<=current_regime; ++i) {
            progress
            << "Regime_" << i << "  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                progress << sysVar.get_m_COOP()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_m_COOP()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            progress
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        }
    }
}





void autoWork::write_report(StructureClass& sysVar)
{
    const int n_regime  = sysVar.get_n_regime();
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar)+"/..");
    o.append("/report_");
    o.append(sysVar.get_usic());
    o.append(".txt");
    
    sysVar.set_sys_targets(o);
    
    cout
    << "Total Time Elapsed: "
    << sysVar.get_timediff_overall() << " seconds."
    << "\n\n";
    
    ofstream report(o.c_str());
    if (report.is_open())
    {
        report
        << "Start  " << sysVar.get_startdatetime()<< "\n"
        << "Finish " << return_datetime()<< "\n\n";
        
        report
        << "tau_alpha Cutoffs ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        } else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        } else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        }
        report << "------------------------------------------\n";
        for (int i=0; i<sysVar.get_n_regime(); ++i) {
            double previous_teq,current_teq,teq;
            if (true) {
                teq = sysVar.get_equilibration_times()[i];
            } else {
                if (i>0) {
                    previous_teq = sysVar.get_equilibration_times()[i-1];
                    current_teq  = sysVar.get_equilibration_times()[i];
                    teq = max(previous_teq,current_teq);
                } else {
                    teq = sysVar.get_equilibration_times()[i];
                }
            }
            if (sysVar.get_n_equ_blocks_plus().size()>0) {
                report << "Regime_"<<i<<"  "
                <<teq/sysVar.get_n_equ_blocks_plus().at(i)<< "\n";
            } else {
                report << "Regime_"<<i<<"  "
                <<teq/sysVar.get_n_equ_blocks()<< "\n";
            }
        }
        report << "------------------------------------------\n\n";
        
        /** report general information **/
        //----------------------------------------------------------------------
        report
        << "General Info" << "\n"
        << "------------------------------------------"
        << "\n"
        << "System_Unit               " << sysVar.get_systemUnit()
        << "\n\n"
        
        << "<Simulation>"
        << "\n"
        << "QuenchRate                " << sysVar.get_quenchRate()
        << "\n";
        
        double n_equ_blocks=sysVar.get_n_equ_blocks();
        double n_prd_blocks=sysVar.get_n_prd_blocks();
        if (sysVar.get_n_equ_blocks_plus().size()>0) {
            report << "Equilibration(tau_alpha)  ";
            for (int i=0;i<(int)sysVar.get_n_equ_blocks_plus().size();++i) {
                report << sysVar.get_n_equ_blocks_plus().at(i) << " ";
            } report << "\n";
        } else {
            report << "Equilibration(tau_alpha)  "<<n_equ_blocks<<"\n";
        }
        if (sysVar.get_n_prd_blocks_plus().size()>0) {
            report << "Production(block)         ";
            for (int i=0;i<(int)sysVar.get_n_prd_blocks_plus().size();++i) {
                report << sysVar.get_n_prd_blocks_plus().at(i) << " ";
            } report << "\n\n";
        } else {
            report << "Production(block)         "<<n_prd_blocks<<"\n\n";
        }
        report
        << "<Analysis>"
        << "\n"
        << "Relaxation_Model          ";
        //for (int i=0; i<sysVar.get_extrp_model().size();++i) {
        //report << sysVar.get_extrp_model().at(i)<<" ";
        report << sysVar.get_extrp_model().at(0);
        //}
        report
        << "\n"
        << "Extrapolation_Model       ";
        for (int i=0; i<sysVar.get_presq_model().size()-1;++i) {
            report << sysVar.get_presq_model().at(i) <<" ";
        }
        report
        << "\n"
        << "Wave_Number(q)            " << sysVar.get_wavenumber()
        << "\n"
        << "------------------------------------------"
        << "\n\n\n";
        
        vector<double> teq=sysVar.get_equilibration_times();
        
        /** report Extrapolated Tg **/
        //----------------------------------------------------------------------
        /** get_Tg_extrp() is a 2D vector;
         ** rows store per regime data, columns store per trial data **/
        report
        << "Extrapolated Tg [@("<<sysVar.get_extrp_time()<<")timeUnit]" << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_Tg_extrp()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_extrp()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_Tg_extrp()[n_regime-1]);
        double Tg_mean=sysVar.calc_mean();
        double Tg_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(Tg)          " << Tg_mean              << "\n"
        << "stdev(Tg)         " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------\n\n";
        
        
        report
        << "Tg(trials) avg(Tg) stdev(Tg) stdev/avg"     << "\n"
        << "------------------------------------------\n";
        for (int i=0; i<n_regime; ++i)
        {
            report << "Regime_"<<i<<"  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_Tg_extrp()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_Tg_extrp()[i]);
            double Tg_mean=sysVar.calc_mean();
            double Tg_stdev=sysVar.get_sample_stdev();
            report
            << Tg_mean  << " "
            << Tg_stdev << " "
            << Tg_stdev/Tg_mean << "\n";
        } report << "------------------------------------------\n\n";
        
        
        report
        << "Tg(VFT) [@("<<sysVar.get_extrp_time()<<")timeUnit]" << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_Tg_VFT()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_VFT()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_Tg_VFT()[n_regime-1]);
        Tg_mean=sysVar.calc_mean();
        Tg_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(Tg_VFT)      " << Tg_mean              << "\n"
        << "stdev(Tg_VFT)     " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------\n\n";
        
        
        report
        << "Tg(COOP) [@("<<sysVar.get_extrp_time()<<")timeUnit]" << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_Tg_COOP()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_COOP()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_Tg_COOP()[n_regime-1]);
        Tg_mean=sysVar.calc_mean();
        Tg_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(Tg_COOP)     " << Tg_mean              << "\n"
        << "stdev(Tg_COOP)    " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------\n\n\n\n";
        
        
        
        
        /** report fragility @ Tg_extrp **/
        //----------------------------------------------------------------------
        /** get_m_extrp() is a 2D vector;
         ** rows store per regime data, columns store per trial data **/
        report
        << "Fragility @ Tg_extrp                      " << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_m_extrp()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_extrp()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_m_extrp()[n_regime-1]);
        double m_mean=sysVar.calc_mean();
        double m_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(m_extrp)   " << m_mean                 << "\n"
        << "stdev(m_extrp)  " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------\n\n";
        
        
        report
        << "m(trials) avg(m) stdev(m) stdev/avg"        << "\n"
        << "------------------------------------------\n";
        for (int i=0; i<n_regime; ++i)
        {
            report << "Regime_"<<i<<"  ";
            for (int ii=0; ii<sysVar.get_n_trial(); ++ii) {
                report << sysVar.get_m_extrp()[i][ii] << " ";
            } sysVar.set_calcvector(sysVar.get_m_extrp()[i]);
            double m_mean=sysVar.calc_mean();
            double m_stdev=sysVar.get_sample_stdev();
            report
            << m_mean  << " "
            << m_stdev << " "
            << m_stdev/m_mean << "\n";
        } report << "------------------------------------------\n\n";
        
        
        report
        << "Fragility(VFT) @ Tg_VFT                   " << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_m_VFT()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_VFT()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_m_VFT()[n_regime-1]);
        m_mean=sysVar.calc_mean();
        m_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(m_VFT)     " << m_mean                 << "\n"
        << "stdev(m_VFT)    " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------\n\n";
        
        
        report
        << "Fragility(COOP) @ Tg_COOP                 " << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_m_COOP()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_COOP()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_m_COOP()[n_regime-1]);
        m_mean=sysVar.calc_mean();
        m_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(m_COOP)    " << m_mean                 << "\n"
        << "stdev(m_COOP)   " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------\n\n\n\n";
        
        
        
        
        /** report Computational Tg **/
        //----------------------------------------------------------------------
        /** get_Tg_compu() is a 2D vector;
         ** rows store per regime data, columns store per trial data **/
        report
        << "Computational Tg [@("<<sysVar.get_compu_time()<<")timeUnit]" << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_Tg_compu()[n_regime-1].size(); ++i) {
            report << sysVar.get_Tg_compu()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_Tg_compu()[n_regime-1]);
        Tg_mean=sysVar.calc_mean();
        Tg_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(Tg_compu)    " << Tg_mean              << "\n"
        << "stdev(Tg_compu)   " << Tg_stdev             << "\n"
        << "stdev/avg         " << Tg_stdev/Tg_mean     << "\n"
        << "------------------------------------------\n\n";
        
        
        /** report fragility @ Tg_compu **/
        //----------------------------------------------------------------------
        /** get_m_compu() is a 2D vector;
         ** rows store per regime data, columns store per trial data **/
        report
        << "Fragility @ Tg_compu                      " << "\n"
        << "------------------------------------------\n";
        for (size_t i=0; i<sysVar.get_m_compu()[n_regime-1].size(); ++i) {
            report << sysVar.get_m_compu()[n_regime-1][i] << "\n";
        } sysVar.set_calcvector(sysVar.get_m_compu()[n_regime-1]);
        m_mean=sysVar.calc_mean();
        m_stdev=sysVar.get_sample_stdev();
        report
        << "\n"
        << "mean(m_compu)   " << m_mean                 << "\n"
        << "stdev(m_compu)  " << m_stdev                << "\n"
        << "stdev/avg       " << m_stdev/m_mean         << "\n"
        << "------------------------------------------\n\n";
        
        
        /** report number of T's simulated in each regime **/
        //----------------------------------------------------------------------
        report
        << "Number of T's Simulated in Each Regime    " << "\n"
        << "------------------------------------------" << "\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "Regime_"<<i<<"  "
            <<sysVar.get_n_regime_temps()[i]<<"\n";
        } report << "------------------------------------------\n\n";
        
        
        /** report timestep size used in each regime **/
        //----------------------------------------------------------------------
        report
        << "Timestep size in each regime ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        } else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        } else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        } report << "------------------------------------------\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "ts[" << i << "]=  "
            << sysVar.get_ts_regime()[i]
            << "\n";
        } report << "------------------------------------------\n\n";
        
        
        /** report equilibration times **/
        //----------------------------------------------------------------------
        report
        << "Equilibration Times ";
        if (sysVar.get_systemUnit()=="real") {
            report << "(timeUnit: fs)" << "\n";
        } else if (sysVar.get_systemUnit()=="metal") {
            report << "(timeUnit: ps)" << "\n";
        } else if (sysVar.get_systemUnit()=="lj") {
            report << "(timeUnit: tau)" << "\n";
        } report << "------------------------------------------\n";
        for (int i=0; i<n_regime; ++i) {
            double teq = sysVar.get_equilibration_times()[i];
            report
            << "Regime_"<<i<<"  "<< teq << "\n";
        } report << "------------------------------------------\n\n";
        
        
        /** report if simulations were retried **/
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
            << "------------------------------------------\n";
            for (int i=0; i<n_regime; ++i) {
                report
                << "Regime_"<< sysVar.get_simulation_retry()[i][0]<<"  "
                << sysVar.get_simulation_retry()[i][1] << "\n";
            } report << "------------------------------------------\n\n";
        }
        
        
        /** report elapased time **/
        //----------------------------------------------------------------------
        report
        << "Time Elapsed for Each Temperature Regime  " << "\n"
        << "------------------------------------------\n";
        for (int i=0; i<n_regime; ++i) {
            report
            << "Regime_"<<i<<"  "
            << sysVar.get_elapsedtime()[i]
            << " seconds"
            << "\n";
        } report << "------------------------------------------\n\n";
        
        report
        << "Total Elapsed Time: "
        << sysVar.get_elapsedtime()[n_regime]
        << " seconds"
        << "\n";
        
        report.close();
    }
}





void autoWork::alglibfit_routine(const StructureClass& sysVar)
{
    FitData fd;
    fd.alglibfit_routine(sysVar);
}





void autoWork::call_system_bash(const string& strVar)
{
    string str0,str1;
    str0="bash ./";
    str1=str0+strVar;
    int ret=system(str1.c_str());
    vector<string> strVec;
    /** This is the list of files that should not be deleted **/
    //strVec.push_back("pkm.inp");
    int count=0;
    for (size_t i=0; i<strVec.size(); ++i) if(strVar==strVec[i]) ++count;
    if (count==0) {
        /** If not in the above list,
         ** the temporary files are deleted after use **/
        str0="rm ./";
        str1=str0+strVar;
        if (ret!=0) exit(EXIT_FAILURE);
        else system(str1.c_str());
    }
}





vector<vector<double>> autoWork::linearTemperatures(const SysInfo& sysVar)
{
    vector<double> hiT,loT,tinfo;
    vector<vector<double>> sortedTemps;
    
    if (sysVar.get_n_Temp()==0) { // use the upper and lower temperature bounds
        /* manually set above crossTemp */
        int n_hiT=
        (sysVar.get_startTemp()-sysVar.get_crossTemp())/sysVar.get_hTres();
        for (int i=0; i<n_hiT; ++i)	{
            tinfo.push_back(sysVar.get_startTemp()-sysVar.get_hTres()*i);
            hiT.push_back(sysVar.get_startTemp()-sysVar.get_hTres()*i);
        }
        /* manually set below crossTemp */
        int n_loT=
        ((sysVar.get_crossTemp()-sysVar.get_finalTemp())/sysVar.get_lTres())+1;
        for (int i=0; i<n_loT; ++i){
            tinfo.push_back(sysVar.get_crossTemp()-sysVar.get_lTres()*i);
            loT.push_back(sysVar.get_crossTemp()-sysVar.get_lTres()*i);
        }
        sortedTemps.push_back(tinfo); // index 0: tinfo
        sortedTemps.push_back(hiT);   // index 1: hiT
        sortedTemps.push_back(loT);   // index 2: loT
    } else { // use N starting temperatures
        double T=sysVar.get_startTemp();
        tinfo.push_back(T);
        if (T>sysVar.get_crossTemp()) {
            hiT.push_back(T);
            T -= sysVar.get_hTres();
        } else {
            loT.push_back(T);
            T -= sysVar.get_lTres();
        }
        for (int i=1; i<sysVar.get_n_Temp(); ++i) {
            if (T>sysVar.get_crossTemp()) {
                /* manually set above crossTemp */
                tinfo.push_back(T);
                hiT.push_back(T);
                T -= sysVar.get_hTres();
            } else if (T<=sysVar.get_crossTemp()) {
                /* manually set below crossTemp */
                tinfo.push_back(T);
                loT.push_back(T);
                T -= sysVar.get_lTres();
            }
        }
        sortedTemps.push_back(tinfo); // index 0: tinfo
        sortedTemps.push_back(hiT);   // index 1: hiT
        sortedTemps.push_back(loT);   // index 2: loT
    } return sortedTemps;
}





const string autoWork::return_year()
{
    time_t t=time(NULL);tm* timePtr=localtime(&t);
    stringstream ss; ss << timePtr->tm_year+1900;
    string year = ss.str();
    return year;
}





const string autoWork::return_date()
{
    /** format: MM/DD/YY **/
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
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}





void autoWork::make_SimulationsFolders(const StructureClass& sysVar)
{
    const string path_d   =sysVar.get_Path();
    const string simType_d=sysVar.get_simType();
    const string year_d   =sysVar.get_year();
    const string usicID_d =sysVar.get_usicID();
    const string usic_d   =sysVar.get_usic();
    
    string makeSimFolder=
    "Master="+path_d+";";
    /* make usic directory */
    makeSimFolder+=
    "usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";";
    /* check if 'simulations' folder exists */
    makeSimFolder+=
    "test -e ${usic_dir};"
    "if [ $? -ne 0 ]; then "
    "   mkdir -p ${usic_dir};"
    "else "
    "   read -t10 -p \""+usic_d+" folder exists, delete and make new? (y/n)\" VAR;"
    "	if [ ${VAR} == y ]; then "
    "		rm -r ${usic_dir};"
    "       echo \"removing $(basename \"${usic_dir}\")\"...;sleep 1;"
    "		wait;"
    "		mkdir -p ${usic_dir};"
    "       echo \"making   $(basename \"${usic_dir}\")\"...;sleep 1;"
    "	else"
    //"       read -t10 -p \"continue? (y/n)\" VAR;"
    //"       if [ ${VAR} == n ]; then "
    "           echo EXIT!;"
    "           exit 3;"
    //"       fi;"
    "	fi;"
    "fi;";
    /* make 'simulations' folder */
    makeSimFolder+=
    "simulations_dir=${usic_dir}/simulations;"
    "mkdir -p ${simulations_dir};";
    /* make 'packmol' folder */
    makeSimFolder+=
    "packmol=${simulations_dir}/packmol;"
    "mkdir -p ${packmol};";
    /* make 'lammps_inputs' folders */
    makeSimFolder+=
    "lammps_inputs=${simulations_dir}/lammps_inputs;"
    "mkdir -p ${lammps_inputs}/start_data;"
    "mkdir -p ${lammps_inputs}/generation;"
    "mkdir -p ${lammps_inputs}/quench;"
    "mkdir -p ${lammps_inputs}/equilibration;"
    "mkdir -p ${lammps_inputs}/tampering;"
    "mkdir -p ${lammps_inputs}/resize;"
    "mkdir -p ${lammps_inputs}/production;";
    /* packet folders to be put in each simulation phase */
    makeSimFolder+=
    "mkdir -p ${usic_dir}/packet/log;"
    "mkdir -p ${usic_dir}/packet/screen;"
    "mkdir -p ${usic_dir}/packet/trajectory;"
    "mkdir -p ${usic_dir}/packet/restart;"
    "mkdir -p ${usic_dir}/packet/submission_files;"
    "mkdir -p ${usic_dir}/packet/submission_files/cluster_out;";
    /* make simualtion (phase) folders */
    makeSimFolder+=
    "mkdir -p ${simulations_dir}/generation;"
    "mkdir -p ${simulations_dir}/quench;"
    "mkdir -p ${simulations_dir}/equilibration;"
    "mkdir -p ${simulations_dir}/tampering;"
    "mkdir -p ${simulations_dir}/resize;"
    "mkdir -p ${simulations_dir}/production;"
    "mkdir -p ${simulations_dir}/submission_scripts;";
    /* copy content of packet folder to each phase */
    makeSimFolder+=
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/generation/;"
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/quench/;"
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/equilibration/;"
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/tampering/;"
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/resize/;"
    "cp -r ${usic_dir}/packet/* ${simulations_dir}/production/;";
    /* make 'analysis' folder */
    makeSimFolder+=
    "analysis_dir=${usic_dir}/analysis;"
    "mkdir -p ${analysis_dir};";
    /* make corresponding directories in the Analysis folder */
    makeSimFolder+=
    "mkdir -p ${analysis_dir}/AMDAT_inputs;"
    "mkdir -p ${analysis_dir}/AMDAT_submission_files;"
    "mkdir -p ${analysis_dir}/AMDAT_submission_files/cluster_out;"
    "mkdir -p ${analysis_dir}/submission_scripts;"
    "mkdir -p ${analysis_dir}/screen;"
    "mkdir -p ${analysis_dir}/statistics;"
    "mkdir -p ${analysis_dir}/fit_data;";
    /* source files backup */
    makeSimFolder+=
    "src_backup=${usic_dir}/src_backup; mkdir -p ${src_backup};"
    "cp *.cpp *.h *.inp ${src_backup}/;"
    "cp ./packmol ${src_backup}/;";
    /* remove packet folder */
    makeSimFolder+="rm -rf ${usic_dir}/packet;";
    /* remove archive folder */
    makeSimFolder+="rm -rf ${Master}/archive;";
    
    int ret=system(makeSimFolder.c_str());
    if (ret!=0) {
        //cout
        //<< "in autoWork::make_SimulationsFolders():\n"
        //<< "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    } else {
        cout << "\nfileSystem made!\n";system_wait(1);
    }
}





void autoWork::makeNewAnanlysisFolder(const StructureClass& sysVar)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    const string mkdir=
    "read -t10 -p \"make new analysis folder? (y/n)\" VAR;"
    "if [ ${VAR} == y ]; then "
    "   Master="+path_d+";"
    "   usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "   cd ${usic_dir};"
    "   mv ./analysis ./analysis_"+return_datetime()+";"
    "   analysis_dir=${usic_dir}/analysis;"
    "   mkdir -p ${analysis_dir};"
    "   mkdir -p ${analysis_dir}/AMDAT_inputs;"
    "   mkdir -p ${analysis_dir}/AMDAT_submission_files;"
    "   mkdir -p ${analysis_dir}/AMDAT_submission_files/cluster_out;"
    "   mkdir -p ${analysis_dir}/submission_scripts;"
    "   mkdir -p ${analysis_dir}/screen;"
    "   mkdir -p ${analysis_dir}/statistics;"
    "   mkdir -p ${analysis_dir}/fit_data;"
    "   cd ${Master};"
    "   cp functions.cpp amdatanalysis.cpp fitdata.cpp theorytest.cpp "
    "   ${analysis_dir};"
    "fi;";
    int ret=system(mkdir.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::makeNewAnanlysisFolder():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::backupfitdata(const StructureClass& sysVar)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    const string mkdir=
    //"read -t10 -p \"backup fit_data folder? (y/n)\" VAR;"
    //"if [ ${VAR} == y ]; then "
    "   Master="+path_d+";"
    "   usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "   cd ${usic_dir}/analysis;"
    "   mv fit_data fit_data_"+return_datetime()+";"
    "   mkdir fit_data;";
    //"fi;";
    int ret=system(mkdir.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::backupfitdata():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::copy_from_to(const StructureClass& sysVar)
{
    const string cp=
    "cp "+sysVar.get_copy_source()+" "+sysVar.get_copy_destin()+";";
    int ret=system(cp.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::copy_from_to():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::test_AnalysisFolder(const StructureClass& sysVar,
                                   const AmdatAnalysis& aa)
{
    const string path_d=return_AnalysisFolderPath(sysVar);
    
    string testAnaFolder="test -e "+path_d;
    int test=system(testAnaFolder.c_str());
    if (test!=0) makeNewAnanlysisFolder(sysVar);
    
    string af;
    /* test "analysis" folder */
    //"test -e "+path_d+";"
    //"if [ $? -ne 0 ]; then "
    //"   mkdir "+path_d+";"
    //"  echo -e "+path_d+" does NOT exist!;"
    //"  exit 3;"
    //"fi;";
    /* create "screen" folder */
    //"test -e "+path_d+"/screen;"
    //"if [ $? -ne 0 ]; then "
    //"  mkdir "+path_d+"/screen;"
    //"fi;";
    /* create "submission_script" folder */
    //af+=
    //"test -e "+path_d+"/submission_scripts;"
    //"if [ $? -ne 0 ]; then "
    //"  mkdir "+path_d+"/submission_scripts;"
    //"fi;";
    /* create "fit_data" folder */
    //af+=
    //"test -e "+path_d+"/fit_data;"
    //"if [ $? -ne 0 ]; then "
    //"  mkdir "+path_d+"/fit_data;"
    //"fi;";
    /* create individual analysis folders */
    if (aa.get_is_strFac())
    {
        af+=
        "strFac_all="+path_d+"/statistics/strFac;"
        "test -e ${strFac_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${strFac_all};"
        "fi;";
    }
    if (aa.get_is_rdf())
    {
        af+=
        "rdf_all="+path_d+"/statistics/rdf;"
        "test -e ${rdf_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${rdf_all};"
        "fi;";
    }
    if (aa.get_is_msd())
    {
        af+=
        "msd_all="+path_d+"/statistics/msd;"
        "test -e ${msd_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${msd_all};"
        "fi;";
    }
    if (aa.get_is_ngp())
    {
        af+=
        "ngp_all="+path_d+"/statistics/ngp;"
        "test -e ${ngp_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${ngp_all};"
        "fi;";
    }
    if (aa.get_is_isfs())
    {
        af+=
        "isfs_all="+path_d+"/statistics/isfs;"
        "test -e ${isfs_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${isfs_all};"
        "fi;";
    }
    if (aa.get_is_baf())
    {
        af+=
        "baf_all="+path_d+"/statistics/baf;"
        "test -e ${baf_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${baf_all};"
        "fi;";
    }
    if (aa.get_is_composition())
    {
        af+=
        "comp_all="+path_d+"/statistics/composition;"
        "test -e ${comp_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${comp_all};"
        "fi;";
    }
    if (aa.get_is_u2dist())
    {
        af+=
        "u2dist_all="+path_d+"/statistics/u2dist;"
        "test -e ${u2dist_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${u2dist_all};"
        "fi;";
    }
    if (aa.get_is_stiffness_dist())
    {
        af+=
        "stiffness_dist_all="+path_d+"/statistics/stiffness_dist;"
        "test -e ${stiffness_dist_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${stiffness_dist_all};"
        "fi;";
    }
    if (aa.get_is_isf())
    {
        af+=
        "isf_all="+path_d+"/statistics/isf;"
        "test -e ${isf_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${isf_all};"
        "fi;";
    }
    if (aa.get_is_strings())
    {
        af+=
        "strings_all="+path_d+"/statistics/strings;"
        "test -e ${strings_all};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${strings_all};"
        "fi;";
    }
    if (aa.get_is_mbodies())
    {
        af+=
        "mbTraj="+path_d+"/statistics/mbTraj;"
        "test -e ${mbTraj};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${mbTraj};"
        "fi;";
    }
    if (aa.get_is_use_voroNeighbors())
    {
        af+=
        "voro="+path_d+"/statistics/voro;"
        "test -e ${voro};"
        "if [ $? -ne 0 ]; then "
        "  mkdir -p ${voro};"
        "fi;";
    }
    int ret=system(af.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::test_AnalysisFolder():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::move_InnerLoopFiles(const StructureClass& sysVar,
                                   const int n_sys,
                                   const int n_trl)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    
    /* folders */
    string mv=
    "Master="+path_d+";cd ${Master};"
    "simulations_dir="
    "${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+"/simulations;"
    "packmolDir=${simulations_dir}/packmol;"
    "lammps_inputs=${simulations_dir}/lammps_inputs;";
    /* move files to corresponding folders */
    mv+=
    "mv pkm_output.xyz "
    "pkm_00"+to_string((long long int)n_trl)+"_"+sysVar.get_nameString(n_sys)+".xyz;"
    "mv start.data input_00"+to_string((long long int)n_trl)+".data;"
    "mv pkm.inp "
    "pkmInput_00"+to_string((long long int)n_trl)+"_"+sysVar.get_nameString(n_sys)+".inp;"
    "mv pkm_* pkmInput_* ${packmolDir}/;"
    "mv input_* ${lammps_inputs}/start_data/;";
    int ret=system(mv.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::move_InnerLoopFiles():\n"
        << "system() returns "<<ret<<"\n\n";
        //exit(EXIT_FAILURE);
    }
}





void autoWork::move_OuterLoopFiles(const StructureClass& sysVar,
                                   const int n_sys)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    
    /* make folders */
    string mkdir=
    "Master="+path_d+";cd ${Master};"
    "usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "singleChain=${usic_dir}/singleChain; mkdir -p ${singleChain};";
    int ret=system(mkdir.c_str());
    if (ret!=0) exit(EXIT_FAILURE);
    /* move files to corresponding folders */
    string mv=
    "Master="+path_d+";cd ${Master};"
    "usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "singleChain=${usic_dir}/singleChain;"
    "mv singleChain.xyz single_"+sysVar.get_nameString(n_sys)+".xyz;"
    "mv Bonds.txt Bonds_"+sysVar.get_nameString(n_sys)+".txt;"
    "mv Angles.txt Angles_"+sysVar.get_nameString(n_sys)+".txt;"
    "mv sys_Params.txt sys_Params_"+sysVar.get_nameString(n_sys)+".txt;"
    "mv single_* Bonds_* Angles_* sys_Params_* ${singleChain}/;";
    ret=system(mv.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::move_OuterLoopFiles():\n"
        << "system() returns "<<ret<<"\n\n";
        //exit(EXIT_FAILURE);
    }
}





void autoWork::backup_simulation_traj(const StructureClass& sysVar,
                                      const string& phase)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    const string trajfolder=
    "Master="+path_d+";"
    "dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+"/simulations/;"
    "trajfolder=${dir}/"+phase+"/trajectory;"
    "cd ${trajfolder};"
    "mkdir Regime_"+to_string(sysVar.get_current_regime())+";"
    "mv trajectory_* Regime_"+to_string(sysVar.get_current_regime())+";";
    int ret=system(trajfolder.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::backup_simulation_traj():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::backup_analysis_stats(const StructureClass& sysVar)
{
    const string path_d    = sysVar.get_Path();
    const string simType_d = sysVar.get_simType();
    const string year_d    = sysVar.get_year();
    const string usicID_d  = sysVar.get_usicID();
    const string usic_d    = sysVar.get_usic();
    const string anafolder =
    "Master="+path_d+";"
    "usic_dir=${Master}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "cd ${usic_dir}/analysis;"
    "mv statistics statistics_Regime_"+to_string(sysVar.get_current_regime())+";";
    int ret=system(anafolder.c_str());
    if (ret!=0) {
        cout
        << "in autoWork::backup_analysis_stats():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





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





const string autoWork::extract_amdat_svn(const std::string& lineContent)
{
    string str,amdat_svn;
    istringstream iss(lineContent);
    while (iss>>str) {
        /** ex. amdat_svn = v.r117 **/
        if (str.at(0)=='v'&&str.at(1)=='.'&&str.at(2)=='r') {
            for (size_t i=3;i<str.size();++i) {
                amdat_svn+=str.at(i);
            }
        }
    } return amdat_svn;
}





void autoWork::system_initialization(StructureClass& sysVar,
                                     const UserInterface * const inpVar)
{
    sysVar.set_startdatetime(return_datetime());
    
    if (inpVar!=nullptr)
    {
        sysVar.set_startdatetime(inpVar->get_sysVar().get_startdatetime());
        sysVar.set_userName(inpVar->get_sysVar().get_userName());
        sysVar.set_masterFolder(inpVar->get_sysVar().get_masterFolder());
        if (inpVar->get_sysVar().get_computeCluster()=="Mycroft") {
            sysVar.set_Path("/home/"+sysVar.get_userName()+"/"+sysVar.get_masterFolder());
        } else if (inpVar->get_sysVar().get_computeCluster()=="myMac") {
            sysVar.set_Path("/Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork");
        } else if (inpVar->get_sysVar().get_computeCluster()=="archive") {
            sysVar.set_Path("/home/jh148/archive/Sean/AUTOWORK");
        } else if (inpVar->get_sysVar().get_computeCluster()=="storage2") {
            sysVar.set_Path("/home/jh148/storage2/Sean/AUTOWORK");
        } else if (inpVar->get_sysVar().get_computeCluster()=="") {
            cout
            << "Please specify the computeCluster. ex. Mycroft or "
            << "an absolute path to your master folder.\n";
            exit(EXIT_FAILURE);
        } else { /* set by the specified path */
            sysVar.set_Path(inpVar->get_sysVar().get_computeCluster());
        }
        sysVar.set_simType(inpVar->get_sysVar().get_simType());
        sysVar.set_year(inpVar->get_sysVar().get_year());
        sysVar.set_usicID(inpVar->get_sysVar().get_usicID());
        if (inpVar->get_sysVar().get_is_use_named_nameString())
        {
            sysVar.set_nameString(inpVar->get_sysVar().get_nameString());
        }
        if (inpVar->get_sysVar().get_is_specify_GPU())
        {
            sysVar.set_which_GPU(inpVar->get_sysVar().get_which_GPU());
        }
        sysVar.set_n_trial(inpVar->get_sysVar().get_n_trial());
        sysVar.set_systemUnit(inpVar->get_sysVar().get_systemUnit());
        sysVar.set_startTemp(inpVar->get_sysVar().get_startTemp());
        sysVar.set_crossTemp(inpVar->get_sysVar().get_crossTemp());
        sysVar.set_hTres(inpVar->get_sysVar().get_hTres());
        sysVar.set_lTres(inpVar->get_sysVar().get_lTres());
        //***
        sysVar.set_n_equ_blocks(inpVar->get_sysVar().get_n_equ_blocks());
        sysVar.set_n_prd_blocks(inpVar->get_sysVar().get_n_prd_blocks());
        sysVar.set_n_relaxation(inpVar->get_sysVar().get_n_relaxation());
        sysVar.set_equilibration_times(inpVar->get_sysVar().get_equilibration_times());
        //***
        sysVar.set_n_regime_temps(inpVar->get_sysVar().get_n_regime_temps());
        sysVar.set_ts_regime(inpVar->get_sysVar().get_ts_regime());
        sysVar.set_regime_beg(inpVar->get_sysVar().get_regime_beg());
        sysVar.set_regime_end(inpVar->get_sysVar().get_regime_end());
        sysVar.set_is_makeFileSystem(inpVar->get_sysVar().get_is_makeFileSystem());
        sysVar.set_is_defaultLmpData(inpVar->get_sysVar().get_is_defaultLmpData());
        sysVar.set_is_moltempLmpData(inpVar->get_sysVar().get_is_moltempLmpData());
        sysVar.set_is_LammpsPhase(inpVar->get_sysVar().get_is_LammpsPhase());
        sysVar.set_is_use_prepScripts(inpVar->get_sysVar().get_is_use_prepScripts());
        sysVar.set_is_Simulations(inpVar->get_sysVar().get_is_Simulations());
        sysVar.set_is_aging(inpVar->get_sysVar().get_is_aging());
        sysVar.set_is_AMDAT(inpVar->get_sysVar().get_is_AMDAT());
        sysVar.set_is_fitData(inpVar->get_sysVar().get_is_fitData());
        sysVar.set_is_fit_dielectric(inpVar->get_sysVar().get_is_fit_dielectric());
        sysVar.set_is_theorytest(inpVar->get_sysVar().get_is_theorytest());
        sysVar.set_is_cutoff_traject(inpVar->get_sysVar().get_is_cutoff_traject());
        sysVar.set_is_alglibfit(inpVar->get_sysVar().get_is_alglibfit());
        sysVar.set_is_directSub(inpVar->get_sysVar().get_is_directSub());
        sysVar.set_is_watch_hold(inpVar->get_sysVar().get_is_watch_hold());
        sysVar.set_is_GPU(inpVar->get_sysVar().get_is_GPU());
        sysVar.set_is_makeNewAnanlysisFolder(inpVar->get_sysVar().get_is_makeNewAnanlysisFolder());
        sysVar.set_is_backupfitdata(inpVar->get_sysVar().get_is_backupfitdata());
        sysVar.set_is_fullquench(inpVar->get_sysVar().get_is_fullquench());
        sysVar.set_is_singleTempTest(inpVar->get_sysVar().get_is_singleTempTest());
        //sysVar.set_is_makeNewAnanlysisFolder(inpVar->get_sysVar().get_is_makeNewAnanlysisFolder());
        //sysVar.set_is_backupfitdata(inpVar->get_sysVar().get_is_backupfitdata());
        sysVar.set_n_digits(inpVar->get_sysVar().get_n_digits());
        sysVar.set_precision(pow(10,sysVar.get_n_digits()));
        sysVar.set_is_amdatinp(inpVar->get_sysVar().get_is_amdatinp());
        sysVar.set_is_copy_from_to(inpVar->get_sysVar().get_is_copy_from_to());
        sysVar.set_copy_source(inpVar->get_sysVar().get_copy_source());
        sysVar.set_copy_destin(inpVar->get_sysVar().get_copy_destin());
    }
    
    //--------------------------------------------------------------------------
    /** define system unit related proprties **/
    /** NOTE: wrote on 20160609
     ** to use viscosity data, now convert it to relaxation time unit first (fs),
     ** fit it and then convert it back to viscosity unit (PaS)
     ** glass transition point may be 10^12 or 10^13 PaS for viscosity,
     ** and 10^17 or 10^18 fs for relaxation time **/
    //--------------------------------------------------------------------------
    if (sysVar.get_systemUnit()=="real") {
        /* time unit: fs */
        sysVar.set_corrFacforT(sysVar.get_precision());
        sysVar.set_extrp_time(1e+17);       // 10^17 fs = 100 seconds
        sysVar.set_is_use_viscosity(false); // true: using viscosity data
    } else if (sysVar.get_systemUnit()=="metal") {
        /* time unit: ps */
        sysVar.set_corrFacforT(sysVar.get_precision());
        sysVar.set_extrp_time(1e+14);       // 10^14 ps = 100 seconds
        sysVar.set_is_use_viscosity(false); // true: using viscosity data
    } else if (sysVar.get_systemUnit()=="lj") {
        /* time unit: tau ~ ps */
        sysVar.set_corrFacforT(sysVar.get_precision()*1000.0);
        sysVar.set_extrp_time(1e+14);       // 10^14 tau ~ 100 seconds
        sysVar.set_is_use_viscosity(false);
    }
    
    /** error checking **/
    error_checking(sysVar);
    
    /** store a cpoy of original equilibration time **/
    sysVar.set_equilibration_times_org(sysVar.get_equilibration_times());
    
    /** set number of executed regimes **/
    sysVar.set_n_regime(abs(sysVar.get_regime_end()-sysVar.get_regime_beg())+1);
    
    /** setup force fields **/
    sysVar.set_structureParam();
    
    if (sysVar.get_is_use_prepScripts()||sysVar.get_is_moltempLmpData())
    {
        sysVar.set_is_amdatinp(false);
    }
    
    /** Temperature containers **/
    //--------------------------------------------------------------------------
    // linearTemperatures() returns:
    // [0] -- overall Ts
    // [1] -- maually set hiTs
    // [2] -- maually set loTs
    //--------------------------------------------------------------------------
    sysVar.set_n_Temp(sysVar.get_n_regime_temps().at(0));
    vector<double> tinfo=linearTemperatures(sysVar).at(0);
    for (int i=0; i<tinfo.size(); ++i) {tinfo.at(i)*=sysVar.get_precision();}
    
    /** set the cutoff T for Arrhenius fit **/
    //**************************************************************************
    double cutT=0;
    if (true) {
        /** lowest T in first regime **/
        size_t s=tinfo.size();
        cutT=floor(tinfo[s-1]); // already expanded
    } else {
        /** arbitrary T based on experience **/
        cutT=1000.0;
        cutT=floor(cutT*sysVar.get_precision()); // use expanded form
    } sysVar.set_cutTforArrhenius(cutT/sysVar.get_corrFacforT());
    //**************************************************************************
    
    vector<double> tmpvd;
    vector<vector<double>> initialTs;
    vector<vector<double>> tinfo2D={tinfo,tinfo};
    //--------------------------------------------------------------------------
    // row:    trials and systems
    // column: elemnts per trial ans system
    //--------------------------------------------------------------------------
    int n_trl=0,n_sys=0;
    for (n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (n_sys=sysVar.get_n_sys_beg(); n_sys<=sysVar.get_n_sys_end(); ++n_sys) {
            initialTs.push_back(tinfo);
            sysVar.get_temperaturesInfo().push_back(tinfo); // 2D double
            sysVar.get_equilibratedTs().push_back(tinfo);   // 2D double
            sysVar.get_quenchingTs().push_back(tinfo);      // 2D double
            sysVar.get_extrpTmax().push_back(tinfo);        // 2D double
            sysVar.get_qchTdiff().push_back(0);             // 1D double
            sysVar.get_tauFit().push_back(tinfo2D);         // 3D double
            sysVar.get_tauEqu().push_back(tinfo2D);         // 3D double
        }
    } sysVar.set_initialtemps(initialTs);
    //--------------------------------------------------------------------------
    // row:    regime
    // column: trials and systems
    //--------------------------------------------------------------------------
    for (int i=0; i<sysVar.get_n_trial(); ++i) {
        for (int ii=sysVar.get_n_sys_beg(); ii<=sysVar.get_n_sys_end(); ++ii) {
            tmpvd.push_back(0.0);
        }
    }
    /** Initialize containers with dummy values to avoid segmentation fault **/
    //--------------------------------------------------------------------------
    for (int i=0; i<sysVar.get_equilibration_times().size(); ++i) {
        /** container has current regime data of all trials of systems **/
        sysVar.get_Tg_extrp().push_back(tmpvd);
        sysVar.get_Tg_compu().push_back(tmpvd);
        sysVar.get_m_extrp().push_back(tmpvd);
        sysVar.get_m_compu().push_back(tmpvd);
        sysVar.get_Tg_VFT().push_back(tmpvd);
        sysVar.get_Tg_COOP().push_back(tmpvd);
        sysVar.get_m_VFT().push_back(tmpvd);
        sysVar.get_m_COOP().push_back(tmpvd);
        sysVar.get_extrp_model().push_back("");
        sysVar.get_presq_model().push_back("");
        /** container has elapsed time of current regime **/
        sysVar.get_elapsedtime().push_back(0);
    }
}





void autoWork::moltemplate_initialization(StructureClass& sysVar,
                                          MoltemplateLmpData& moldata,
                                          const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        moldata.set_merSet(inpVar->get_moldata().get_merSet());
        moldata.set_merComp(inpVar->get_moldata().get_merComp());
        moldata.set_copolymerType(inpVar->get_moldata().get_copolymerType());
        moldata.set_tacticity(inpVar->get_moldata().get_tacticity());
        moldata.set_dop(inpVar->get_moldata().get_dop());
        moldata.set_is_restrictbynAtoms(inpVar->get_moldata().get_is_restrictbynAtoms());
        moldata.set_n_total_atoms(inpVar->get_moldata().get_n_total_atoms());
        moldata.set_sequenceNum(inpVar->get_moldata().get_sequenceNum());
        
        sysVar.set_chainLen(inpVar->get_moldata().get_dop());
    }
}





void autoWork::fitDielectric_initialization(StructureClass& sysVar,
                                            FitDielectric& ds,
                                            const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        if (sysVar.get_is_fit_dielectric()) {
            if (inpVar->get_ds().get_form()=="") {
                cout
                << "To run fitDielectric, form needs to be determined. "
                << "Available options are NPIC, NIST, or manual.     \n";
                exit(EXIT_FAILURE);
            }
        }
        ds.set_form(inpVar->get_ds().get_form());
        ds.set_func(inpVar->get_ds().get_func());
        ds.set_model(inpVar->get_ds().get_model());
        ds.set_is_use_dc_free(inpVar->get_ds().get_is_use_dc_free());
        ds.set_is_manualFitSingleDS(inpVar->get_ds().get_is_manualFitSingleDS());
        ds.set_sfreq_range(inpVar->get_ds().get_sfreq_range());        
    }
}





void autoWork::simulation_initialization(StructureClass& sysVar,
                                         WorkScripts& ws,
                                         const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        /* bool vars */
        ws.set_is_node0(inpVar->get_ws().get_is_node0());
        ws.set_is_node1(inpVar->get_ws().get_is_node1());
        ws.set_is_node2(inpVar->get_ws().get_is_node2());
        
        /* vars already assigned value in UserInterface constructor */
        ws.set_numCores(inpVar->get_ws().get_numCores());
        ws.set_run_cores(inpVar->get_ws().get_run_cores());
        
        //NOTE:
        //ws.cores (the c variable) does not have to reset the value here;
        //it is determined based on the number of configured cores (6 or 8) and
        //the number of temperatures simulated at the particular regime
        
        /* vars specifically specified in input script */
        if (inpVar->get_ws().get_is_node0()) {
            ws.set_node0(inpVar->get_ws().get_node0());
        }
        if (inpVar->get_ws().get_is_node1()) {
            ws.set_node1(inpVar->get_ws().get_node1());
        }
        if (inpVar->get_ws().get_is_node2()) {
            ws.set_node2(inpVar->get_ws().get_node2());
        }
        if (inpVar->get_ws().get_steps_gen()!=0) {
            ws.set_steps_gen(inpVar->get_ws().get_steps_gen());
        }
        if (inpVar->get_ws().get_quenchRate()!=0) {
            ws.set_quenchRate(inpVar->get_ws().get_quenchRate());
        }
    }
    sysVar.set_quenchRate(ws.get_quenchRate());
    sysVar.set_is_node0_info(ws.get_is_node0());
    sysVar.set_is_node1_info(ws.get_is_node1());
    sysVar.set_is_node2_info(ws.get_is_node2());
    sysVar.set_node0_info(ws.get_node0());
    sysVar.set_node1_info(ws.get_node1());
    sysVar.set_node2_info(ws.get_node2());
    
    if ((ws.get_is_node0()||ws.get_is_node1())&&(ws.get_is_node2()))
    {
        cout
        << "\n"
        << "Node-2 can NOT be used together with Node-0 or Node-1.\n"
        << "Please re-specify the compute nodes.\n\n";
        exit(EXIT_FAILURE);
    }
    
    /** Enfore configuration number for # of Ts for simulation when using GPU **/
    if (sysVar.get_is_GPU())
    {
        if (ws.get_is_node0()||ws.get_is_node1())
        {
            for (int i=0;i<sysVar.get_n_regime_temps().size();++i) {
                int n_temps=sysVar.get_n_regime_temps().at(i);
                if (n_temps>8) {
                    cout
                    << "\n"
                    << "At regime_"<<sysVar.get_current_regime()<<", "
                    << "number of Ts > 8, which exceeds the configuration number.\n"
                    << "For using GPU, Please choose from the suggested set {2,4,8}.\n"
                    << "Or, turn off GPU by setting is_GPU to no.\n";
                    exit(EXIT_FAILURE);
                }
                if (n_temps%2!=0&&n_temps<8) {
                    cout
                    << "\n"
                    << "At regime_"<<sysVar.get_current_regime()<<", "
                    << "number of Ts is "<<n_temps<<", which does not conform to the "
                    << "optimum configuration number.\n"
                    << "For using GPU, Please choose from the suggested set {2,4,8}.\n"
                    << "Or, turn off GPU by setting is_GPU to no.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
        if (ws.get_is_node2())
        {
            for (int i=0;i<sysVar.get_n_regime_temps().size();++i) {
                int n_temps=sysVar.get_n_regime_temps().at(i);
                if (n_temps>6) {
                    cout
                    << "\n"
                    << "At regime_"<<i<<", number of Ts > 6, "
                    << "which exceeds the configuration number for compute node-2.\n"
                    << "For using GPU, Please choose from the suggested set {2,3,6}.\n"
                    << "Or, turn off GPU by setting is_GPU to no.\n";
                    exit(EXIT_FAILURE);
                }
                if (n_temps%3!=0&&n_temps<6) {
                    cout
                    << "\n"
                    << "At regime_"<<sysVar.get_current_regime()<<", "
                    << "number of Ts is "<<n_temps<<", which does not conform to the "
                    << "optimum configuration number.\n"
                    << "For using GPU, Please choose from the suggested set {2,3,6}.\n"
                    << "Or, turn off GPU by setting is_GPU to no.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}





void autoWork::analysis_initialization(StructureClass& sysVar,
                                       AmdatAnalysis& aa,
                                       const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        /* bool vars */
        aa.set_is_mbodies(inpVar->get_aa().get_is_mbodies());
        aa.set_is_monoStruct(inpVar->get_aa().get_is_monoStruct());
        aa.set_is_strings(inpVar->get_aa().get_is_strings());
        aa.set_is_changeAtomType(inpVar->get_aa().get_is_changeAtomType());
        aa.set_is_use_voroNeighbors(inpVar->get_aa().get_is_use_voroNeighbors());
        
        /* vars specifically specified in input script */
        if (inpVar->get_aa().get_analysispart()!="") {
            aa.set_analysispart(inpVar->get_aa().get_analysispart());
        }
        if (inpVar->get_aa().get_relaxation_target()!="") {
            aa.set_relaxation_target(inpVar->get_aa().get_relaxation_target());
        }
        if (inpVar->get_aa().get_amdat_numCores()!=0) {
            aa.set_amdat_numCores(inpVar->get_aa().get_amdat_numCores());
        }
        if (inpVar->get_aa().get_logtaubeta()!=0) {
            aa.set_logtaubeta(inpVar->get_aa().get_logtaubeta());
        }
    }
}





void autoWork::fd_initialization(StructureClass& sysVar,
                                 FitData& fd,
                                 const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        /* bool vars */
        fd.set_is_use_FG(inpVar->get_fd().get_is_use_FG());
        fd.set_is_fit_Fs_by_spline(inpVar->get_fd().get_is_fit_Fs_by_spline());
        fd.set_is_fit_full_alpha(inpVar->get_fd().get_is_fit_full_alpha());
        fd.set_is_use_gammafunc(inpVar->get_fd().get_is_use_gammafunc());
        fd.set_is_fit_sExp(inpVar->get_fd().get_is_fit_sExp());
        fd.set_is_fit_lnFs(inpVar->get_fd().get_is_fit_lnFs());
        fd.set_is_1stmoment_tau(inpVar->get_fd().get_is_1stmoment_tau());
        fd.set_is_fit_Arrhenius(inpVar->get_fd().get_is_fit_Arrhenius());
        fd.set_is_fit_tauFit(inpVar->get_fd().get_is_fit_tauFit());
        fd.set_is_applytauFitcut(inpVar->get_fd().get_is_applytauFitcut());
        fd.set_is_use_KWWassist(inpVar->get_fd().get_is_use_KWWassist());
        fd.set_is_find_DWF(inpVar->get_fd().get_is_find_DWF());
        fd.set_is_find_NGP(inpVar->get_fd().get_is_find_NGP());
        fd.set_is_calc_thermoData(inpVar->get_fd().get_is_calc_thermoData());
        fd.set_is_normalModes(inpVar->get_fd().get_is_normalModes());
        fd.set_is_qspectrum(inpVar->get_fd().get_is_qspectrum());
        
        /* vars specifically specified in input script */
        if (inpVar->get_fd().get_definitionofTA()!="") {
            fd.set_definitionofTA(inpVar->get_fd().get_definitionofTA());
        }
        if (inpVar->get_fd().get_fcorr_model()!="") {
            fd.set_fcorr_model(inpVar->get_fd().get_fcorr_model());
        }
        if (inpVar->get_fd().get_extrp_model()!="") {
            fd.set_extrp_model(inpVar->get_fd().get_extrp_model());
        }
        if (inpVar->get_fd().get_presq_model()!="") {
            fd.set_presq_model(inpVar->get_fd().get_presq_model());
        }
        if (inpVar->get_fd().get_xtauA()!=0) {
            fd.set_xtauA(inpVar->get_fd().get_xtauA());
        }
        if (inpVar->get_fd().get_xTA()!=0) {
            fd.set_xTA(inpVar->get_fd().get_xTA());
        }
        if (inpVar->get_fd().get_waveindices().size()!=0) {
            fd.set_waveindices(inpVar->get_fd().get_waveindices());
        }
        if (inpVar->get_fd().get_logtselect().size()!=0) {
            fd.set_logtselect(inpVar->get_fd().get_logtselect());
        }
    }
    sysVar.get_extrp_model().at(sysVar.get_current_regime())=fd.get_extrp_model();
    sysVar.get_presq_model().at(sysVar.get_current_regime())=fd.get_presq_model();
}





void autoWork::tt_initialization(StructureClass& sysVar,
                                 TheoryTest& tt,
                                 const UserInterface * const inpVar)
{
    if (inpVar!=nullptr)
    {
        /* bool vars */
        tt.set_is_AGtest(inpVar->get_tt().get_is_AGtest());
        tt.set_is_GLMtest(inpVar->get_tt().get_is_GLMtest());
        tt.set_is_Leporinitest(inpVar->get_tt().get_is_Leporinitest());
        tt.set_is_HWtest(inpVar->get_tt().get_is_HWtest());
        tt.set_is_reltaur(inpVar->get_tt().get_is_reltaur());
    }
}





void autoWork::use_prepScripts(const StructureClass& sysVar)
{
    string path_scripts = sysVar.get_Path()+"/scripts/";
    string simfolder    = return_SimulationFolderPath(sysVar);
    string anafolder    = return_AnalysisFolderPath(sysVar);
    
    /** cp input_*.data **/
    string cp=
    "n_data=0;"
    "for i in "+simfolder+"/lammps_inputs/start_data/input_*;"
    "do test -e $i; "
    "   if [ $? -eq 0 ]; then "
    "       n_data=$[${n_data}+1];"
    "   fi; "
    "done;"
    "n=0;"
    "if [ ${n_data} -eq 0 ]; then "
    "   for i in "+path_scripts+"input_*;"
    "   do test -e $i; "
    "       if [ $? -eq 0 ]; then "
    "           n=$[$n+1];"
    "           cp $i "+simfolder+"/lammps_inputs/start_data/;"
    "           echo \"$(basename \"$i\") copied.\";"//" sleep 1;"
    "       fi; "
    "   done;"
    "   if [ $n -eq 0 ]; then "
    "       echo \"No data file found in 'scripts' folder\";"
    "       sleep 2;"
    "       exit 0;"
    "   fi;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    /** cp LAMMPS scripts **/
    cp=
    "n=0;"
    "geninp="+path_scripts+"generation.inp; test -e ${geninp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"gen* "+simfolder+"/lammps_inputs/generation/;"
    "   echo \"generation.inp copied.\";"//" sleep 1;"
    "fi;"
    "if [ $n -eq 0 ]; then "
    "   echo \"No generation.inp in the scripts folder!\";"
    "   exit 2;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    cp=
    "n=0;"
    "qchinp="+path_scripts+"quench.inp; test -e ${qchinp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"quench.inp "+simfolder+"/lammps_inputs/quench/;"
    "   echo \"quench.inp copied.\";"//" sleep 1;"
    "fi;"
    "if [ $n -eq 0 ]; then "
    "   echo \"No quench.inp in the scripts folder!\";"
    "   exit 3;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    cp=
    "n=0;"
    "equinp="+path_scripts+"equilibration.inp; test -e ${equinp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"equilibration.inp "+simfolder+"/lammps_inputs/equilibration/;"
    "   echo \"equilibration.inp copied.\";"//" sleep 1;"
    "fi;"
    "if [ $n -eq 0 ]; then "
    "   echo \"No equilibration.inp in the scripts folder!\";"
    "   exit 4;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    cp=
    "n=0;"
    "resinp="+path_scripts+"resize.inp; test -e ${resinp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"resize.inp "+simfolder+"/lammps_inputs/resize/;"
    "   echo \"resize.inp copied.\";"//" sleep 1;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    cp=
    "n=0;"
    "prdinp="+path_scripts+"production.inp; test -e ${prdinp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"production.inp "+simfolder+"/lammps_inputs/production/;"
    "   echo \"production.inp copied.\";"//" sleep 1;"
    "fi;"
    "if [ $n -eq 0 ]; then "
    "   echo \"No production.inp in the scripts folder!\";"
    "   exit 5;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
    
    /** copy amdat.inp scripts **/
    cp=
    "n=0;"
    "amdinp="+path_scripts+"amdat.inp; test -e ${amdinp}; "
    "if [ $? -eq 0 ]; then "
    "   n=$[$n+1];"
    "   cp "+path_scripts+"amdat.inp "+anafolder+"/AMDAT_inputs/;"
    "   echo \"amdat.inp copied.\";"//" sleep 1;"
    "fi;"
    "if [ $n -eq 0 ]; then "
    "   echo \"No amdat.inp in the scripts folder!\";"
    "   exit 6;"
    "fi;"; if (system(cp.c_str())!=0) exit(EXIT_FAILURE);
}





void autoWork::error_checking(const StructureClass& sysVar)
{
    if (sysVar.get_is_moltempLmpData()==true &&
        sysVar.get_is_defaultLmpData()==true) {
        cout << "\n"
        << "'Moltemplate' and 'Default' datafile gen switches are both turned on, "
        << "please just turn on one of the two.\n\n";
        exit(EXIT_FAILURE);
    }
    if (sysVar.get_is_moltempLmpData()==true &&
        sysVar.get_is_use_prepScripts()==true) {
        cout << "\n"
        << "'Moltemplate' and 'use_prepScripts' switches are both turned on, "
        << "please just turn on one of the two.\n\n";
        exit(EXIT_FAILURE);
    }
    if (sysVar.get_is_LammpsPhase()==true &&
        sysVar.get_is_use_prepScripts()==true) {
        cout << "\n"
        << "'LammpsPhase' and 'use_prepScripts' switches are both turned on, "
        << "please just turn on one of the two.\n\n";
        exit(EXIT_FAILURE);
    }
    if (sysVar.get_typeB()==0 && sysVar.get_typeS()==0) {
        if (sysVar.get_systemUnit()!="real" && sysVar.get_systemUnit()!="metal") {
            cout << "\n"
            << "species(0,0): system unit should be real or metal\n\n";
            exit(EXIT_FAILURE);
        }
    }
    if (sysVar.get_typeB()==0 && sysVar.get_typeS()==5) {
        if (sysVar.get_systemUnit()!="lj") {
            cout << "\n"
            << "species(0,5): system unit should be lj\n\n";
            exit(EXIT_FAILURE);
        }
    }
    if (sysVar.get_typeB()==5 && sysVar.get_typeS()==6) {
        if (sysVar.get_systemUnit()!="lj") {
            cout << "\n"
            << "species(5,6): system unit should be lj\n\n";
            exit(EXIT_FAILURE);
        }
    }
    if (sysVar.get_is_specify_GPU()) {
        if (sysVar.get_n_trial()!=1) {
            cout << "\n"
            << "To use specific GPU, n_trial should be 1. Please check.\n\n";
            exit(EXIT_FAILURE);
        }
    }
    if (sysVar.get_regime_beg()>sysVar.get_regime_end()) {
        cout
        <<"\nWarning: "
        <<"regime_beg("<<sysVar.get_regime_beg()<<") > "
        <<"regime_end("<<sysVar.get_regime_end()<<")\n\n";
        exit(EXIT_FAILURE);
    }
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
    for (size_t i=0; i<filesToClean.size(); ++i)
    {
        string tmpStr="test -e "+filesToClean.at(i);
        if (system(tmpStr.c_str())==0) {
            tmpStr=rm+filesToClean.at(i);
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





void autoWork::system_wait(const int waitsec)
{
    string wait="sleep "+to_string((long long int)waitsec);
    system(wait.c_str());
}





void autoWork::system_countdown(const int waitsec)
{
    string str=
    "time="+to_string((long long int)waitsec)+";"
    //"printf \"start in \";"
    "for ((i=0;i<time;++i)); do"
    "   ret=$[${time}-$i];"
    "   printf \"${ret} \"; sleep 1;"
    "done; echo;"; system(str.c_str());
}





bool autoWork::is_path_exist(const std::string& path)
{
    bool is_exist=false;
    string str="test -e "+path;
    int ret=system(str.c_str());
    if (ret==0) is_exist=true;
    else is_exist=false;
    return is_exist;
}





double autoWork::error_propagation(const std::vector<double>& dFi,
                                   const std::vector<double>& errori)
{
    if (dFi.size()==errori.size()) {
        double errorF=0;
        for (size_t i=0; i<dFi.size(); ++i) {
            errorF += pow(dFi.at(i),2)*pow(errori.at(i),2);
        } errorF=sqrt(errorF);
        return errorF;
    } else {
        // number of elements doesn't match
        return (-1.0);
    }
}





void autoWork::coreAllocation(const StructureClass& sysVar,
                              const int n,
                              int& c)
{
    if (sysVar.get_n_regime_temps().at(sysVar.get_current_regime())>n) {
        c=1;
    } else {
        c=(ceil)(n/sysVar.get_n_regime_temps().at(sysVar.get_current_regime()));
    }
    //cout << "n = " << n << "\n";
    //cout << "c = " << c << "\n";
}





void autoWork::makeltfilefrompdb()
{
    /** NOTE: RULEs for Building Molecules Using AVOGADRO
     ---------------------------------------------------------------------------
     ** Polymerizing Atoms:
     -- left-end  atom index = 1
     -- right-end atom index = 2
     
     ** Backbone Atoms Other Than the Polymerizing Atoms:
     -- From left to right 3,4,5,...
     -- ex. -(1-3-4-5-...-2)n-
     ------------------------------------------------------------------------**/
    
    string pdbfolderpath="/Users/SJH/Dropbox/UA_Research/Codes/cppWork/forcefieldgen/forcefieldgen/forcefieldgen/AA_Polybutadiene/";
    string monoBankpath="/Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork/Monomer_bank/";
    string newSpeciesName="AAPB_i";
    string newltfileID="S011";
    
    string newltfile=monoBankpath+newltfileID+".lt";
    ofstream writeFile(newltfile.c_str());
    
    string filepath=pdbfolderpath+newSpeciesName+".pdb";
    ifstream readFile(filepath.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        getline(readFile,lineContent);
        getline(readFile,lineContent);
        vector<double> xyz_one;
        vector<string> atomtypes;
        vector<vector<double>> xyz_all;
        vector<vector<int>> bonds_all;
        int int1=0,int2=0;
        double dubVar=0;
        string strVar;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strVar;
            if (strVar=="HETATM"||strVar=="ATOM") { //atom coordinates
                iss.clear();iss.str(lineContent);
                for (int i=0;i<5;++i) {
                    iss >> strVar;
                } xyz_one.clear();
                for (int i=0;i<3;++i) {
                    iss >> dubVar;
                    xyz_one.push_back(dubVar);
                } xyz_all.push_back(xyz_one);
                /** atom type (manually added) **/
                for (int i=0;i<3;++i) {
                    iss >> strVar;
                } iss >> strVar;
                atomtypes.push_back(strVar);
            }
            if (strVar=="CONECT") { //bonds
                iss.clear();iss.str(lineContent);
                iss >> strVar;
                iss >> int1;
                while (iss>>int2) {
                    if (int2>int1) bonds_all.push_back({int1,int2});
                }
            }
        } readFile.close();
        /* sanity check */
        if (xyz_all.size()!=atomtypes.size()) {
            cout
            << "in makeltfilefrompdb():\n"
            << "number of atoms and types don't match. Please check.\n\n";
            exit(EXIT_FAILURE);
        }
        /** convert pdb info into lt file **/
        /** atom coordinates **/
        writeFile
        << "import \"oplsaa.lt\"\n\n"
        << "## "<< newSpeciesName << "\n\n"
        << newltfileID<<" "<<"inherits OPLSAA {\n\n"
        << "    "<<"## atom-id  mol-id  atom-type  charge  x  y  z\n\n"
        << "    "<<"write(\"Data Atoms\") {\n";
        for (size_t i=0;i<xyz_all.size();++i) {
            writeFile
            << "        "
            << "$atom:C"<<i+1<<" "<<"$mol:..."<<" "
            <<"@atom:"<<atomtypes.at(i)<<" "<<"0"<<" "
            << xyz_all.at(i).at(0)<<" "
            << xyz_all.at(i).at(1)<<" "
            << xyz_all.at(i).at(2)<<"\n";
        } writeFile<< "    "<<"}\n";
        /** bonds **/
        writeFile
        << "    "<<"write('Data Bond List') {\n";
        for (size_t i=0;i<bonds_all.size();++i) {
            writeFile
            << "        "
            << "$bond:Bond"<<i+1<<" "
            << "$atom:C"<<bonds_all.at(i).at(0)<<" "
            << "$atom:C"<<bonds_all.at(i).at(1)<<"\n";
        } writeFile<< "    "<<"}\n";
        writeFile << "} ## "<<newSpeciesName<<"\n\n";
        writeFile.close();
    } else {
        cout
        << "in autoWork::makeltfilefrompdb():\n"
        << filepath+" cannot open. Please check.\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::makeTacltfile()
{
    /** number of chiral centers per mer **/
    int n_backboneAtoms=2;
    int n_chiralCentersPerMer=1;
    
    /** position indices of the chiral centers in a mer, from small to large **/
    vector<int> pos_chiralCenters={2};
    
    /** file path **/
    string monoBankpath="/Users/SJH/Dropbox/UA_Research/Codes/cppWork/cppWork/Monomer_bank/";
    string SpeciesName="AAPVC_i";
    string ltfileID="S009i";
    
    /** sanity check **/
    if (n_chiralCentersPerMer>n_backboneAtoms) {
        cout
        << "in autoWork::makeTacltfile():\n"
        << "n_chiralCentersPerMer("<<n_chiralCentersPerMer<<") should <= "
        << "n_backboneAtoms("<<n_backboneAtoms<<") \n\n";
        exit(EXIT_FAILURE);
    }
    string filepath=monoBankpath+ltfileID+".lt";
    ifstream readFile(filepath.c_str());
    if (readFile.is_open())
    {
        string lineContent,strVar;
        vector<vector<int>> linkAtomsSet;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strVar;
            /** navigate to chiral centers and reconize the atoms linked to the
             ** chiral atoms **/
            if (strVar=="write('Data") {
                iss >> strVar;
                if (strVar=="Bond") {
                    while (getline(readFile,lineContent)) {
                        iss.clear();iss.str(lineContent);
                        iss >> strVar;
                        if (strVar=="}") break;//break reading Bond List
                        iss >> strVar;//ex.$atom:C1
                        
                        iss >> strVar;//ex.$atom:C2
                        
                    }
                }
            }
        }
        /** for the linked atoms of a particular chiral atom, it's excluding
         ** other chiral atoms if there's any around the current one **/
        
        
        /** swap coordinates of two of the linked atoms **/
        
        
        /** store each swapped configuration to a new Tx lt file **/
        
        
        
    } else {
        cout
        << "in autoWork::makeTacltfile():\n"
        << filepath+" cannot open. Please check.\n\n";
        exit(EXIT_FAILURE);
    }
}





void autoWork::check_finished_blocks(const StructureClass& sysVar,
                                     const WorkScripts& ws,
                                     const vector<vector<double>>& tinfo)
{
    /** the block size of current regime **/
    int prd_blocksize=ws.get_prd_blocksize();
    
    /** if there are excess 0th frames **/
    //bool is_zerothframes=false;
    
    /** if total step number > 2e+9 steps **/
    bool is_divide=false;
    if (ws.get_divide_prd()>1) is_divide=true;
    
    for (int n_trl=0; n_trl<(int)tinfo.size(); ++n_trl) {
        for (int n_sys=sysVar.get_n_sys_beg(); n_sys<=sysVar.get_n_sys_end(); ++n_sys) {
            for (int T=0; T<(int)tinfo.at(n_trl).size(); ++T) {
                
                int total_frames=0;
                int finished_blocks=0;
                
                string file;
                file.append(return_SimulationFolderPath(sysVar));
                file.append("/production/trajectory/");
                file.append("trajectory_");
                file.append(sysVar.get_usic());
                file.append("_00"+to_string((long long int)n_trl));
                file.append("_"+sysVar.get_nameString(n_sys));
                file.append("_T"+to_string((long long int)tinfo.at(n_trl).at(T)));
                file.append(".prd.custom");
                ifstream readFile(file.c_str());
                
                if (readFile.is_open())
                {
                    bool is_readTimestep=false;
                    
                    int countStepkeyword=0;
                    int line_index=0;
                    int block_index=1;//1-based
                    int frame_index=0;
                    int countaccess=0;
                    
                    string lineContent;
                    string strData0,strData1;
                    
                    while (getline(readFile,lineContent))
                    {
                        ++line_index;
                        ++countStepkeyword;
                        istringstream iss(lineContent);
                        iss >> strData0;
                        
                        if (strData0=="ITEM:") {
                            iss >> strData1;
                            if (strData1=="TIMESTEP") {
                                ++countaccess;
                                countStepkeyword=0;
                            }
                        } is_readTimestep=(countStepkeyword==1)&&(countaccess>0);
                        
                        /** timestep number **/
                        //------------------------------------------------------
                        if (is_readTimestep) {
                            
                            ++frame_index;//1-based
                            ++total_frames;
                            
                            if (frame_index==prd_blocksize) {
                                ++block_index;//1-based
                                ++finished_blocks;
                                frame_index=0;
                            }
                        }
                    } readFile.close();
                    
                } else {
                    cout
                    << "in autoWork::check_finished_blocks():\n"
                    << file << " cannot open.\n";exit(EXIT_FAILURE);
                }
                
                string species;
                species.append(sysVar.get_usic());
                species.append("_00"+to_string((long long int)n_trl));
                species.append("_"+sysVar.get_nameString(n_sys));
                species.append("_T"+to_string((long long int)tinfo.at(n_trl).at(T)));
                cout << species << "\n"
                << "is_divide       " << is_divide       << "\n"
                << "prd_blocksize   " << prd_blocksize   << "\n"
                << "finished_blocks " << finished_blocks << "\n"
                << "total_frames    " << total_frames    << "\n\n";
            }
        }
    }
}





void autoWork::cutoff_trajectory(const StructureClass& sysVar,
                                 const WorkScripts& ws,
                                 const vector<vector<double>>& tinfo,
                                 const int cutoffBlock)
{
    /** the block size of current regime **/
    int prd_blocksize=ws.get_prd_blocksize();
    
    /** if there are excess 0th frames **/
    bool is_zerothframes=false;
    if (is_zerothframes) prd_blocksize+=1;
    
    /** if total step number > 2e+9 steps **/
    bool is_divide=false;
    if (ws.get_divide_prd()>1) is_divide=true;
    
    for (int n_trl=0; n_trl<(int)tinfo.size(); ++n_trl) {
        for (int n_sys=sysVar.get_n_sys_beg(); n_sys<=sysVar.get_n_sys_end(); ++n_sys) {
            for (int T=0; T<(int)tinfo.at(n_trl).size(); ++T) {
                
                string dir,file;
                dir.append(return_SimulationFolderPath(sysVar));
                dir.append("/production/trajectory/");
                file.append("trajectory_");
                file.append(sysVar.get_usic());
                file.append("_00"+to_string((long long int)n_trl));
                file.append("_"+sysVar.get_nameString(n_sys));
                file.append("_T"+to_string((long long int)tinfo.at(n_trl).at(T)));
                file.append(".prd.custom");
                
                /** backup old file **/
                string org=dir+file;
                string old=dir+"old_"+file;
                string mvname="mv "+org+" "+old;
                system(mvname.c_str());
                
                ofstream writeFile(org.c_str());
                ifstream readFile(old.c_str());
                
                if (readFile.is_open())
                {
                    bool is_final_frame=false;
                    bool is_readTimestep=false;
                    
                    int steps_zero=0;
                    int steps_corr=0;
                    int countStepkeyword=0;
                    int line_index=0;
                    int block_index=1;//1-based
                    int frame_index=0;
                    int countaccess=0;
                    
                    string lineContent;
                    string strData0,strData1;
                    
                    while (getline(readFile,lineContent))
                    {
                        ++line_index;
                        ++countStepkeyword;
                        istringstream iss(lineContent);
                        iss >> strData0;
                        
                        if (strData0=="ITEM:") {
                            iss >> strData1;
                            if (strData1=="TIMESTEP") {
                                ++countaccess;
                                countStepkeyword=0;
                            }
                        }
                        if (is_zerothframes) {
                            /** all blocks have equal number of frames **/
                            is_readTimestep=(countStepkeyword==1)&&(countaccess>0);
                        } else {
                            /** first block has one more frame than others **/
                            is_readTimestep=(countStepkeyword==1)&&(countaccess>1);
                        }
                        
                        /** timestep number **/
                        //------------------------------------------------------
                        if (is_readTimestep) {
                            
                            ++frame_index;//1-based
                            
                            if (frame_index==prd_blocksize) {
                                ++block_index;//1-based
                                frame_index=0;
                                is_final_frame=true;
                            } else {
                                is_final_frame=false;
                            }
                            
                            if (block_index>1) {
                                if (is_divide) {
                                    /** if is_divide is true, step is zeroed at new blocks **/
                                    writeFile << lineContent << "\n";
                                } else {
                                    writeFile << lineContent << "\n";
                                    /** reduce step number to per block basis **/
                                    if (false) {
                                        if (is_final_frame) {
                                            steps_zero=stoi(strData0);
                                            /* NOTE: block_index has increased */
                                            steps_corr=steps_zero-
                                            (block_index-2)*(int)pow(1.2,prd_blocksize-1);
                                            writeFile << steps_corr << "\n";
                                        } else {
                                            steps_zero=stoi(strData0);
                                            steps_corr=steps_zero-
                                            (block_index-1)*(int)pow(1.2,prd_blocksize-1);
                                            writeFile << steps_corr << "\n";
                                        }
                                    }
                                }
                            } else {
                                writeFile << lineContent << "\n";
                            }
                        }
                        /** atom coordinates & other keyWords **/
                        //------------------------------------------------------
                        else {
                            /** break while loop at the cutoff block **/
                            if ((block_index==cutoffBlock)&&(countStepkeyword==0)) break;
                            writeFile << lineContent << "\n";
                        }
                    } readFile.close(); writeFile.close();
                    
                } else {
                    cout
                    << "in autoWork::cutoff_trajectory():\n"
                    << file << " cannot open.\n";exit(EXIT_FAILURE);
                }
                
                string species;
                species.append(sysVar.get_usic());
                species.append("_00"+to_string((long long int)n_trl));
                species.append("_"+sysVar.get_nameString(n_sys));
                species.append("_T"+to_string((long long int)tinfo.at(n_trl).at(T)));
                cout << species << " finished!" << "\n";
            }
        }
    }
}





bool autoWork::randbool(double p)
{
    /** modified from randint;
     ** p has to be from 0 to 1,
     ** with 0 representing always returning false, 1 always true **/
    if (p<0) p=0;
    else if (p>1) p=1;
    int min=(int)(-2*1e7+p*2*1e7);
    int max=(int)(min+2*1e7);
    static time_t seed=time(NULL);
    srand((unsigned int)seed);
    seed+=3333;
    int dis=abs(max-min)+1;
    int res=abs((333333*rand())%dis);
    res += std::min(min,max);
    
    if (res<0) return false;
    else return true;
}





int autoWork::randint(const int min, const int max)
{
    /** returned results lie in the inclusive range of the two input integers,
     ** including crossover to negative integers **/
    static time_t seed=time(NULL);
    srand((unsigned int)seed);
    seed+=3333;
    int dis=abs(max-min)+1;
    int res=abs((333333*rand())%dis);
    res += std::min(min,max);
    return res;
}





#ifdef NEVER
const string autoWork::echo_NodesUsage(const vector<int> &node0_d,
                                       const vector<int> &node1_d)
{
    ofstream sh("num.sh");
    sh
    << "#!/bin/bash \n"
    << "node0=(";
    for (size_t i=0; i<node0_d.size(); ++i) {
        sh << node0_d[i];
        if(i!=(node0_d.size()-1)) sh << " ";
    } sh << ")";
    sh << "node1=(";
    for (size_t i=0; i<node1_d.size(); ++i) {
        sh << node1_d[i];
        if(i!=(node1_d.size()-1)) sh << " ";
    } sh << ")";
    sh
    << "numNode0=${#node0[@]} \n"
    << "numNode1=${#node1[@]} \n"
    << "sumNode=$[${numNode0}+${numNode1}] \n"
    << "hlfNode=$[${sumNode}/2] \n"
    << "numGPU=$[${sumNode}*2] \n"
    << "numCPU=$[${sumNode}*16] \n"
    << "echo 'numNode0 = '\"${#node0[@]}\" \n"
    << "echo 'node0    = ('\"${node0[@]}\"')' \n"
    << "echo 'numNode1 = '\"${#node1[@]}\" \n"
    << "echo 'node1    = ('\"${node1[@]}\"')' \n"
    << "echo 'sumNode  = '\"${sumNode}\" \n"
    << "echo 'hlfNode  = '\"${hlfNode}\" \n"
    << "echo 'numGPU   = '\"${numGPU}\" \n"
    << "echo 'numCPU   = '\"${numCPU}\" \n\n";
    sh.close(); return ("num.sh");
}
#endif




