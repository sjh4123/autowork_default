//
//  Created by Sean Jh H. on 7/4/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "workscripts.h"
#include "functions.h"
#include "structureclass.h"
#include "lmpscripts.h"

using namespace std;
using namespace autoWork;

WorkScripts::WorkScripts(const StructureClass& sysVar):

/* string */
lmp_exe(),

/* bool */
is_node0(),
is_node1(),

/* int */
typeB(sysVar.get_typeB()),
typeS(sysVar.get_typeS()),
CoresPerNode(16),
prd_blocksize_hiT(),
prd_blocksize_loT(),

/* double */
n_equ_blocks(),
n_prd_blocks(),
prd_exp_base(1.2),
timestep_size(),
time_equ(),
quenchRate(),
steps_gen(),
steps_qch(),
steps_equ(),
steps_res(5e+4),

/* vector<int> */
node0(),
node1(),
numCores(),
cores(),
run_cores(),
priority({0,0,0,0})

{
    lmp_exe="/home/jh148/lammps-26Jan2014-GPU2/bin/lmp_openmpi";
    
    const int current_regime=sysVar.get_current_regime();
    
    /* Set timestep size of current regime */
    //--------------------------------------------------------------------------
    set_timestep_size(sysVar.get_ts_regime()[current_regime]);
    //--------------------------------------------------------------------------
    
    
    /* Generation (@ highest specified T) */
    //--------------------------------------------------------------------------
    set_steps_gen(1e+6);
    //--------------------------------------------------------------------------
    
    
    /* Quenching */
    //--------------------------------------------------------------------------
    // NOTE:
    // unit of quench rate:
    // for real unit: 'kelvins per nano second'
    // for lj unit:   'T per tau'
    //--------------------------------------------------------------------------
    if (sysVar.get_systemUnit()=="real")
    {
        set_quenchRate(10.0);
    }
    else if (sysVar.get_systemUnit()=="lj")
    {
        set_quenchRate(1e-4);
    }
    else if (sysVar.get_systemUnit()=="metal")
    {
        set_quenchRate(10.0);
    }
    //--------------------------------------------------------------------------
    
    /* Equilibration */
    //--------------------------------------------------------------------------
    double previous_teq,current_teq,teq;
    if (current_regime>0) {
        previous_teq = sysVar.get_equilibration_times()[current_regime-1];
        current_teq  = sysVar.get_equilibration_times()[current_regime];
        teq=max(previous_teq,current_teq);
    }
    else {
        teq=sysVar.get_equilibration_times()[current_regime];
    }
    set_time_equ(teq);
    set_n_equ_blocks(sysVar.get_n_equ_blocks());
    //--------------------------------------------------------------------------
    
    
    /* Resize */
    //--------------------------------------------------------------------------
    set_steps_res(5e+4);
    //--------------------------------------------------------------------------
    
    
    /* Production */
    // NOTE:
    // # of production blocks should not be less than 1
    //--------------------------------------------------------------------------
    double n_prd_block=sysVar.get_n_prd_blocks();
    if (sysVar.get_n_prd_blocks()<1) n_prd_block=1;
    set_n_prd_blocks(n_prd_block);
    //--------------------------------------------------------------------------
}





void WorkScripts::make_GenerationSubFiles(const StructureClass& sysVar,
                                          const LmpScripts& lmp,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d)
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/generation/submission_files");
    o.append("/gen_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".qsub");
    //cout << o << "\n";
    
    ofstream gen0(o.c_str());
    if (gen0.is_open())	{
        
        gen0
        << "#!/bin/bash"                                                << "\n"
        << "#$ -V"                                                      << "\n"
        << "#$ -cwd"                                                    << "\n"
        << "#$ -j y"                                                    << "\n"
        << "#$ -pe orte " << numCores[0]                                << "\n"
        << "#$ -p " << priority[0]                                      << "\n";
        if (sysVar.get_is_GPU()) gen0 << "#$ -R y"                      << "\n";
        
        streamsize ss=gen0.precision();
        gen0 << fixed << setprecision(0);
        
        gen0
        << "#$ -N gen_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        gen0
        << "#$ -o ./generation/submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".gen.o"
        << "\n";
        
        gen0.precision(ss);
        gen0 << resetiosflags(ios::fixed|ios::showpoint);
        
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
        if (get_is_node0() && get_is_node1()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = (int)(get_node1().size());
            
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = 0;
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            numNode0 = 0;
            numNode1 = (int)(get_node1().size());
        }
        
        int sumNode = numNode0 + numNode1;
        
        static int counter=0;
        
        /*
         nonimal number of job allocation on a node; when this nominal number
         increases, it means change to the next node because the current node
         is filled up completely */
        int noml=counter/(get_CoresPerNode()/get_cores()[0]);
        
        /*
         index of the node */
        int tmpNode=noml%sumNode;
        int tmpNode1=0;
        
        /*
         using both node0 and node1 */
        if (get_is_node0() && get_is_node1()) {
            if (tmpNode<numNode0) {
                gen0
                << "#$ -q all.q@compute-0-"
                << get_node0()[tmpNode] << ".local"
                << "\n";
            }
            else {
                tmpNode1=(tmpNode-numNode0);
                gen0
                << "#$ -q all.q@compute-1-"
                << get_node1()[tmpNode1] << ".local"
                << "\n";
            }
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            gen0
            << "#$ -q all.q@compute-0-"
            << get_node0()[tmpNode] << ".local"
            << "\n";
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            gen0
            << "#$ -q all.q@compute-1-"
            << get_node1()[tmpNode] << ".local"
            << "\n";
        }
        //======================================================================
        // GPU Allocation
        //======================================================================
        /*
         gpuNum is switched between 0 and 1 */
        int gpuNum=0;
        /*
         number of half node cores */
        int hNode=get_CoresPerNode()/2;
        /*
         individual jobs on a node */
        int indv=counter%(get_CoresPerNode()/get_cores()[0]);
        /*
         assign GPU to individual jobs */
        if(sysVar.get_is_GPU()) {
            if (indv<(hNode/get_cores()[0])) gpuNum=0; else gpuNum=1;
        }
        /*
         Mycroft compute-0-16's GPU1 is broken */
        if((get_is_node0())||
           (get_is_node0() && get_is_node1()))
            if(get_node0()[tmpNode]==16) gpuNum=0;
        
        ++counter;
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        gen0 << "\n";
        
        
        /* Simulator Executable */
        //----------------------------------------------------------------------
        gen0
        << "mpirun -np " << get_run_cores()[0] << " "
        << get_lmp_exe();
        //----------------------------------------------------------------------
        
        
        /* Generation Inputfile */
        //----------------------------------------------------------------------
        gen0
        << " -in ./lammps_inputs/generation/generation.inp";
        //----------------------------------------------------------------------
        
        
        if (sysVar.get_is_GPU()) gen0 << " -sf gpu";
        
        
        /* random number generator for vseed */
        ////////////////////////////////////////////////////////////////////////
        //======================================================================
        static time_t seed=time(NULL);
        srand((unsigned int)seed);
        seed+=3333;
        
        int vseed=0; // rand number in (1e5,1e6)
        do {
            vseed=abs((333333*rand())%1000000);
        } while (vseed<1e+5);
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        
        
        /* input varaibles */
        //======================================================================
        gen0
        << " -var GPU "        << gpuNum
        << " -var usic "       << sysVar.get_usic()
        << " -var trial 00"    << n_trl
        << " -var namestring " << sysVar.get_nameString(n_sys)
        << " -var vseed "      << vseed
        << " -var ts "         << get_timestep_size()
        << " -var steps_gen "  << (int)ceil(get_steps_gen())
        << fixed << setprecision(0)
        << " -var Temp "       << Temp_d;
        //======================================================================
        gen0.precision(ss);
        gen0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /* log file */
        gen0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        gen0
        << " -log ./generation/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".gen.log";
        //----------------------------------------------------------------------
        
        /* screen file */
        //----------------------------------------------------------------------
        gen0
        << " > ./generation/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".gen.screen"
        << "\n";
        //----------------------------------------------------------------------
        gen0.precision(ss);
        gen0 << resetiosflags(ios::fixed|ios::showpoint);
        
        gen0.close();
    }
    else cout << "GenFile: 'gen.qsub' cannot open." << "\n";
    
}





void WorkScripts::make_QuenchSubFiles(StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys)
{
    const int n_system  = sysVar.get_n_system();
    const int n_sys_beg = sysVar.get_n_sys_beg();
    
    int index=n_trl*n_system+(n_sys-n_sys_beg);
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/quench/submission_files");
    o.append("/qch_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    // NOTE:
    o.append("_Regime"+to_string((long long int)sysVar.get_current_regime()));
    o.append(".qsub");
    //cout << o << "\n";
    
    ofstream qch0(o.c_str());
    if (qch0.is_open())	{
        
        qch0
        << "#!/bin/bash"                                                << "\n"
        << "#$ -V"                                                      << "\n"
        << "#$ -cwd"                                                    << "\n"
        << "#$ -j y"                                                    << "\n"
        << "#$ -pe orte " << numCores[1]                                << "\n"
        << "#$ -p " << priority[1]                                      << "\n";
        if (sysVar.get_is_GPU()) qch0 << "#$ -R y"                      << "\n";
        
        streamsize ss=qch0.precision();
        qch0 << fixed << setprecision(0);
        
        if (sysVar.get_current_regime()==0) {
            qch0
            << "#$ -hold_jid gen_"
            << sysVar.get_usic()
            << "_00" << n_trl << "_"
            << sysVar.get_nameString(n_sys)
            << "_T"	<< sysVar.get_initialtemps()[index][0]
            << "\n";
        }
        
        qch0
        << "#$ -N qch_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "\n";
        
        qch0
        << "#$ -o ./quench/submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.o"
        << "\n";
        
        qch0.precision(ss);
        qch0 << resetiosflags(ios::fixed|ios::showpoint);
        
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
        if (get_is_node0() && get_is_node1()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = (int)(get_node1().size());
            
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = 0;
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            numNode0 = 0;
            numNode1 = (int)(get_node1().size());
        }
        
        int sumNode = numNode0 + numNode1;
        
        static int counter=0;
        
        /*
         nonimal number of job allocation on a node; when this nominal number
         increases, it means change to the next node because the current node
         is filled up completely */
        int noml=counter/(get_CoresPerNode()/get_cores()[1]);
        
        /*
         index of the node */
        int tmpNode=noml%sumNode;
        int tmpNode1=0;
        
        /*
         using both node0 and node1 */
        if (get_is_node0() && get_is_node1()) {
            if (tmpNode<numNode0) {
                qch0
                << "#$ -q all.q@compute-0-"
                << get_node0()[tmpNode] << ".local"
                << "\n";
            }
            else {
                tmpNode1=(tmpNode-numNode0);
                qch0
                << "#$ -q all.q@compute-1-"
                << get_node1()[tmpNode1] << ".local"
                << "\n";
            }
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            qch0
            << "#$ -q all.q@compute-0-"
            << get_node0()[tmpNode] << ".local"
            << "\n";
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            qch0
            << "#$ -q all.q@compute-1-"
            << get_node1()[tmpNode] << ".local"
            << "\n";
        }
        //======================================================================
        // GPU Allocation
        //======================================================================
        /*
         gpuNum is switched between 0 and 1 */
        int gpuNum=0;
        /*
         number of half node cores */
        int hNode=get_CoresPerNode()/2;
        /*
         individual jobs on a node */
        int indv=counter%(get_CoresPerNode()/get_cores()[1]);
        /*
         assign GPU to individual jobs */
        if(sysVar.get_is_GPU()) {
            if (indv<(hNode/get_cores()[1])) gpuNum=0; else gpuNum=1;
        }
        /*
         Mycroft compute-0-16's GPU1 is broken */
        if((get_is_node0())||
           (get_is_node0() && get_is_node1()))
            if(get_node0()[tmpNode]==16) gpuNum=0;
        
        ++counter;
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        qch0 << "\n";
        
        
        /* Simulator Executable */
        //----------------------------------------------------------------------
        qch0
        << "mpirun -np " << get_run_cores()[1] << " "
        << get_lmp_exe();
        //----------------------------------------------------------------------
        
        
        /* Quench Inputfile */
        //----------------------------------------------------------------------
        qch0
        << " -in ./lammps_inputs/quench/quench.inp";
        //----------------------------------------------------------------------
        
        
        if (sysVar.get_is_GPU()) qch0 << " -sf gpu";
        
        
        ////////////////////////////////////////////////////////////////////////
        //======================================================================
        // Deciding Quenching Temperatures
        //======================================================================
        int n_qchTemps=0;       // # of quenching T's
        
        vector<double> vecA;    // stores higher T's for quenching
        vector<double> vecB;    // stores lower T's for quenching
        
        double steps_qch  = 0;     // # of quench steps
        
        double qchTdiff   = 0;  // T difference of highest/lowest qch T's
        double quenchRate = get_quenchRate(); // qch rate: kelvins per ns
        double ts         = get_timestep_size();
        
        
        //----------------------------------------------------------------------
        // Full Quench Mode:
        // Quench straight from highest T to lowest T
        //----------------------------------------------------------------------
        if (sysVar.get_is_fullquench()) {
            
            double hiT=sysVar.get_initialtemps()[index][0];
            double loT=sysVar.get_finalTemp();
            
            qchTdiff=fabs(hiT-loT);
            
            double currentT=hiT;
            while (currentT>loT) {
                vecA.push_back(currentT);
                ++n_qchTemps;
                currentT -= sysVar.get_precision();
            }
            ++n_qchTemps; // total number of quenching temperatures
            
            currentT=hiT-sysVar.get_precision();
            while (currentT>=loT) {
                vecB.push_back(currentT);
                currentT -= sysVar.get_precision();
            }
        }
        //----------------------------------------------------------------------
        
        
        //----------------------------------------------------------------------
        // Stepwise Quench Mode:
        // Quench from highest T to lowest T of "current regime" only
        // and repeat the same quenching pattern until all regimes finish
        //----------------------------------------------------------------------
        else {
            
            vector<double> quenchingTs;
            
            if (sysVar.get_current_regime()==0) {
                
                quenchingTs=sysVar.get_temperaturesInfo()[index];
                
                /* set T difference in quenching */
                size_t qchTsize=quenchingTs.size();
                qchTdiff=fabs(quenchingTs[0]-quenchingTs[qchTsize-1]);
                sysVar.get_qchTdiff()[index]=qchTdiff;
                
            }
            else {
                
                /* equlibrated temperatures of previous regime */
                size_t equTSize = sysVar.get_equilibratedTs()[index].size();
                double equT0 = sysVar.get_equilibratedTs()[index][equTSize-1];
                
                /* lowest equilibrated T of previous regime */
                quenchingTs.push_back(equT0);
                
                /* put extrapolated T's into the temperature container */
                for (size_t i=0; i<sysVar.get_quenchingTs()[index].size(); ++i) {
                    quenchingTs.push_back(sysVar.get_quenchingTs()[index][i]);
                }
                
                /* set T difference for quenching */
                size_t qchTsize=quenchingTs.size();
                qchTdiff=fabs(quenchingTs[0]-quenchingTs[qchTsize-1]);
                sysVar.get_qchTdiff()[index]=qchTdiff;
                
            }
            
            // # of quenching temperatures
            n_qchTemps=int(quenchingTs.size());
            
            for (int i=0; i<(n_qchTemps-1); ++i) {
                vecA.push_back(quenchingTs[i]);
            }
            for (int i=1; i<n_qchTemps; ++i) {
                vecB.push_back(quenchingTs[i]);
            }
        }
        //----------------------------------------------------------------------
        
        
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        
        // qchTdiff is determined in the above block by either of quench mode
        qchTdiff *= pow(sysVar.get_corrFacforT(),-1);
        
        if (sysVar.get_systemUnit()=="real") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*ts*1e-6));
        }
        else if (sysVar.get_systemUnit()=="metal") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*ts*1e-3));
        }
        else if (sysVar.get_systemUnit()=="lj") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*ts));
        }
        
        if (steps_qch<1e+3) steps_qch=1e+3;
        
        set_steps_qch(steps_qch);
        
        /* input varaibles */
        //======================================================================
        qch0
        << " -var GPU "         << gpuNum
        << " -var usic "        << sysVar.get_usic()
        << " -var trial 00"     << n_trl
        << " -var namestring "  << sysVar.get_nameString(n_sys)
        << " -var Regime "      << sysVar.get_current_regime()
        << " -var n_intervals " << (n_qchTemps-1)
        << " -var ts "          << get_timestep_size()
        << " -var steps_qch "   << (int)(ceil)(steps_qch);
        
        qch0 << fixed << setprecision(0);
        
        qch0
        << " -var TempA ";
        for (size_t i=0; i<vecA.size(); ++i) {
            if (i!=(vecA.size()-1)) {
                qch0 << vecA[i] << " ";}
            else {
                qch0 << vecA[i];}
        }
        qch0
        << " -var TempB ";
        for (size_t i=0; i<vecB.size(); ++i) {
            if (i!=(vecB.size()-1)) {
                qch0 << vecB[i] << " ";}
            else {
                qch0 << vecB[i];}
        }
        
        if (sysVar.get_current_regime()==0) {
            qch0
            << " -var readRestartPath "
            << "./generation/restart/restart_"
            << sysVar.get_usic()
            << "_00" << n_trl
            << "_"   << sysVar.get_nameString(n_sys)
            << "_T"  << sysVar.get_initialtemps()[index][0]
            << ".gen.restart";
        }
        else {
            
            size_t vecSize = sysVar.get_equilibratedTs()[index].size();
            double restartT= sysVar.get_equilibratedTs()[index][vecSize-1];
            qch0
            << " -var readRestartPath "
            << "./production/restart/restart_"
            << sysVar.get_usic()
            << "_00" << n_trl
            << "_"   << sysVar.get_nameString(n_sys)
            << "_T"  << restartT
            << ".prd.restart";
        }
        //======================================================================
        qch0.precision(ss);
        qch0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /* log file */
        qch0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        qch0
        << " -log ./quench/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.log";
        //----------------------------------------------------------------------
        
        /* screen file */
        //----------------------------------------------------------------------
        qch0
        << " > ./quench/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.screen"
        << "\n";
        //----------------------------------------------------------------------
        qch0.precision(ss);
        qch0 << resetiosflags(ios::fixed|ios::showpoint);
        
        qch0.close();
    }
    else cout << "QchFile: 'qch.qsub' cannot open." << "\n";
    
}





void WorkScripts::make_EquilibrationSubfiles(const StructureClass& sysVar,
                                             const LmpScripts& lmp,
                                             const int n_trl,
                                             const int n_sys,
                                             const double Temp_d)
{
    // NOTE:
    // the upper limit for a signed 32-bit interger to store is 2^31
    // No larger than this number should be used in a single run
    double run_steps=get_steps_equ();
    if (sysVar.get_systemUnit()=="real")  if (run_steps>1e+9) run_steps=1e+9;
    if (sysVar.get_systemUnit()=="metal") if (run_steps>1e+9) run_steps=1e+9;
    if (sysVar.get_systemUnit()=="lj")    if (run_steps>2e+9) run_steps=2e+9;
    
    divide_equ=1;
    if (run_steps>2e+9) {
        do {
            ++divide_equ;
        }
        while ((run_steps/(double)divide_equ)>2e+9);
        run_steps /= (double)divide_equ;
        run_steps = (ceil)(run_steps);
    }
    
    bool   is_set_CheckPoint=false;
    double restartf=0;
    if (sysVar.get_systemUnit()=="real")  restartf=1e+7;
    if (sysVar.get_systemUnit()=="metal") restartf=1e+7;
    if (sysVar.get_systemUnit()=="lj")    restartf=1e+8;
    if (run_steps>restartf) is_set_CheckPoint=true;
    
    string o;
    
    for (int run_phase=1; run_phase<=divide_equ; ++run_phase) {
        
        if (divide_equ>1) {
            
            o.clear();
            o.append(return_SimulationFolderPath(sysVar));
            o.append("/equilibration/submission_files");
            o.append("/equ_");
            o.append(to_string((long long int)(run_phase))+"_");
            o.append(to_string((long long int)divide_equ)+"_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append("_T"+to_string((long long int)Temp_d));
            o.append(".qsub");
            //cout << o << "\n";
        }
        else {
            
            o.clear();
            o.append(return_SimulationFolderPath(sysVar));
            o.append("/equilibration/submission_files");
            o.append("/equ_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append("_T"+to_string((long long int)Temp_d));
            o.append(".qsub");
            //cout << o << "\n";
        }
        
        ofstream equ0(o.c_str());
        if (equ0.is_open())	{
            
            equ0
            << "#!/bin/bash"                                            << "\n"
            << "#$ -V"                                                  << "\n"
            << "#$ -cwd"                                                << "\n"
            << "#$ -j y"                                                << "\n"
            << "#$ -pe orte " << numCores[2]                            << "\n"
            << "#$ -p " << priority[2]                                  << "\n";
            if (sysVar.get_is_GPU()) equ0 << "#$ -R y"                  << "\n";
            
            streamsize ss=equ0.precision();
            equ0 << fixed << setprecision(0);
            
            if (divide_equ>1) {
                
                if (run_phase==1) {
                    
                    equ0
                    << "#$ -hold_jid qch_"
                    << sysVar.get_usic()
                    << "_00" << n_trl << "_"
                    << sysVar.get_nameString(n_sys)
                    << "\n";
                }
                else {
                    
                    equ0
                    << "#$ -hold_jid equ_"
                    << (run_phase-1) << "_"
                    << divide_equ << "_"
                    << sysVar.get_usic()
                    << "_00" << n_trl << "_"
                    << sysVar.get_nameString(n_sys)
                    << "_T"	<< Temp_d
                    << "\n";
                }
                
                equ0
                << "#$ -N equ_"
                << run_phase << "_"
                << divide_equ << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
                
                equ0
                << "#$ -o ./equilibration/submission_files/cluster_out/out_"
                << run_phase << "_"
                << divide_equ << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".equ.o"
                << "\n";
            }
            
            else {
                
                equ0
                << "#$ -hold_jid qch_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "\n";
                
                equ0
                << "#$ -N equ_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
                
                equ0
                << "#$ -o ./equilibration/submission_files/cluster_out/out_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".equ.o"
                << "\n";
            }
            
            equ0.precision(ss);
            equ0 << resetiosflags(ios::fixed|ios::showpoint);
            
            ////////////////////////////////////////////////////////////////////
            //==================================================================
            // Nodes Allocation
            //==================================================================
            /*
             For compute nodes alloccation: within a node, jobs are filled on a half
             node with the same GPU, and then to the next half with the other GPU;
             the same procedure continues in the next requested node. */
            
            int numNode0 = 0; // number of nodes on node0
            int numNode1 = 0; // number of nodes on node1
            
            /*
             using both node0 and node1 */
            if (get_is_node0() && get_is_node1()) {
                numNode0 = (int)(get_node0().size());
                numNode1 = (int)(get_node1().size());
                
            }
            /*
             only using node0 */
            else if (get_is_node0()) {
                numNode0 = (int)(get_node0().size());
                numNode1 = 0;
            }
            /*
             only using node1 */
            else if (get_is_node1()) {
                numNode0 = 0;
                numNode1 = (int)(get_node1().size());
            }
            
            int sumNode = numNode0 + numNode1;
            
            static int counter=0;
            
            /*
             nonimal number of job allocation on a node; when this nominal number
             increases, it means change to the next node because the current node
             is filled up completely */
            int noml=counter/(get_CoresPerNode()/get_cores()[2]);
            
            /*
             index of the node */
            int tmpNode=noml%sumNode;
            int tmpNode1=0;
            
            /*
             using both node0 and node1 */
            if (get_is_node0() && get_is_node1()) {
                if (tmpNode<numNode0) {
                    equ0
                    << "#$ -q all.q@compute-0-"
                    << get_node0()[tmpNode] << ".local"
                    << "\n";
                }
                else {
                    tmpNode1=(tmpNode-numNode0);
                    equ0
                    << "#$ -q all.q@compute-1-"
                    << get_node1()[tmpNode1] << ".local"
                    << "\n";
                }
            }
            /*
             only using node0 */
            else if (get_is_node0()) {
                equ0
                << "#$ -q all.q@compute-0-"
                << get_node0()[tmpNode] << ".local"
                << "\n";
            }
            /*
             only using node1 */
            else if (get_is_node1()) {
                equ0
                << "#$ -q all.q@compute-1-"
                << get_node1()[tmpNode] << ".local"
                << "\n";
            }
            //==================================================================
            // GPU Allocation
            //==================================================================
            /*
             gpuNum is switched between 0 and 1 */
            int gpuNum=0;
            /*
             number of half node cores */
            int hNode=get_CoresPerNode()/2;
            /*
             individual jobs on a node */
            int indv=counter%(get_CoresPerNode()/get_cores()[2]);
            /*
             assign GPU to individual jobs */
            if(sysVar.get_is_GPU()) {
                if (indv<(hNode/get_cores()[2])) gpuNum=0; else gpuNum=1;
            }
            /*
             Mycroft compute-0-16's GPU1 is broken */
            if((get_is_node0())||
               (get_is_node0() && get_is_node1()))
                if(get_node0()[tmpNode]==16) gpuNum=0;
            
            ++counter;
            //==================================================================
            ////////////////////////////////////////////////////////////////////
            equ0 << "\n";
            
            
            /* Simulator Executable */
            //------------------------------------------------------------------
            equ0
            << "mpirun -np "
            << get_run_cores()[2] << " "
            << get_lmp_exe();
            //------------------------------------------------------------------
            
            
            /* Equilibration Inputfile */
            //------------------------------------------------------------------
            equ0
            << " -in ./lammps_inputs/equilibration/equilibration.inp";
            //------------------------------------------------------------------
            
            
            if (sysVar.get_is_GPU()) equ0 << " -sf gpu";
            
            
            /* input varaibles */
            //==================================================================
            equ0
            << " -var GPU "            << gpuNum
            << " -var usic "           << sysVar.get_usic()
            << " -var trial 00"        << n_trl
            << " -var namestring "     << sysVar.get_nameString(n_sys)
            << " -var Regime "         << sysVar.get_current_regime()
            << " -var ts "             << get_timestep_size()
            << " -var steps_equ "      << (int)run_steps
            << " -var run_phase "      << run_phase
            << " -var divide_equ "     << divide_equ
            << " -var set_CheckPoint " << is_set_CheckPoint
            << " -var restartf "       << (int)restartf
            << fixed << setprecision(0)
            << " -var Temp "           << Temp_d;
            
            equ0 << " -var read_res ";
            string target;
            if (run_phase==1) {
                
                target.append(return_SimulationFolderPath(sysVar));
                target.append("/quench/restart");
                target.append("/restart_");
                target.append(sysVar.get_usic());
                target.append("_00"+to_string((long long int)n_trl));
                target.append("_"+sysVar.get_nameString(n_sys));
                target.append("_T"+to_string((long long int)Temp_d));
                target.append(".qch.restart");
            }
            else {
                
                target.append(return_SimulationFolderPath(sysVar));
                target.append("/equilibration/restart");
                target.append("/restart_");
                target.append(sysVar.get_usic());
                target.append("_00"+to_string((long long int)n_trl));
                target.append("_"+sysVar.get_nameString(n_sys));
                target.append("_T"+to_string((long long int)Temp_d));
                target.append(".equ");
                target.append("."+to_string((long long int)(run_phase-1)));
                target.append("."+to_string((long long int)divide_equ));
                target.append(".restart");
            }
            equ0 << target;
            
            
            equ0 << " -var write_res ";
            target.clear();
            if (run_phase<divide_equ) {
                
                target.append(return_SimulationFolderPath(sysVar));
                target.append("/equilibration/restart");
                target.append("/restart_");
                target.append(sysVar.get_usic());
                target.append("_00"+to_string((long long int)n_trl));
                target.append("_"+sysVar.get_nameString(n_sys));
                target.append("_T"+to_string((long long int)Temp_d));
                target.append(".equ");
                target.append("."+to_string((long long int)run_phase));
                target.append("."+to_string((long long int)divide_equ));
                target.append(".restart");
            }
            else {
                
                target.append(return_SimulationFolderPath(sysVar));
                target.append("/equilibration/restart");
                target.append("/restart_");
                target.append(sysVar.get_usic());
                target.append("_00"+to_string((long long int)n_trl));
                target.append("_"+sysVar.get_nameString(n_sys));
                target.append("_T"+to_string((long long int)Temp_d));
                target.append(".equ.restart");
            }
            equ0 << target;
            
            // NOTE:
            // In LAMMPS, the timesteps N must fit in a signed 32-bit integer
            // so limited to ~ 2 billion steps (2^31) in a "single run"
            //==================================================================
            equ0.precision(ss);
            equ0 << resetiosflags(ios::fixed|ios::showpoint);
            
            /* log file */
            equ0 << fixed << setprecision(0);
            //------------------------------------------------------------------
            equ0
            << " -log ./equilibration/log/"
            << sysVar.get_usic()
            << "_00" << n_trl << "_"
            << sysVar.get_nameString(n_sys)
            << "_T" << Temp_d
            << ".equ.log";
            //------------------------------------------------------------------
            
            /* screen file */
            //------------------------------------------------------------------
            if (divide_equ>1) {
                equ0
                << " > ./equilibration/screen/"
                << "equ_"
                << run_phase << "_"
                << divide_equ << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".screen" << "\n";
            }
            else {
                equ0
                << " > ./equilibration/screen/"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".equ.screen" << "\n";
            }
            //------------------------------------------------------------------
            equ0.precision(ss);
            equ0 << resetiosflags(ios::fixed|ios::showpoint);
            
            equ0.close();
        }
        else cout << "EquFile: 'equ.qsub' cannot open." << "\n";
        
    }
}





void WorkScripts::make_ResizeSubfiles(const StructureClass& sysVar,
                                      const LmpScripts& lmp,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/resize/submission_files");
    o.append("/res_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".qsub");
    //cout << o << "\n";
    
    ofstream res0(o.c_str());
    if (res0.is_open())	{
        
        res0
        << "#!/bin/bash"                                                << "\n"
        << "#$ -V"                                                      << "\n"
        << "#$ -cwd"                                                    << "\n"
        << "#$ -j y"                                                    << "\n"
        << "#$ -pe orte " << numCores[3]                                << "\n"
        << "#$ -p " << priority[3]                                      << "\n";
        if (sysVar.get_is_GPU()) res0 << "#$ -R y"                      << "\n";
        
        streamsize ss=res0.precision();
        res0 << fixed << setprecision(0);
        
        res0
        << "#$ -hold_jid equ_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        res0
        << "#$ -N res_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        res0
        << "#$ -o ./resize/submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".res.o"
        << "\n";
        
        res0.precision(ss);
        res0 << resetiosflags(ios::fixed|ios::showpoint);
        
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
        if (get_is_node0() && get_is_node1()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = (int)(get_node1().size());
            
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = 0;
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            numNode0 = 0;
            numNode1 = (int)(get_node1().size());
        }
        
        int sumNode = numNode0 + numNode1;
        
        static int counter=0;
        
        /*
         nonimal number of job allocation on a node; when this nominal number
         increases, it means change to the next node because the current node
         is filled up completely */
        int noml=counter/(get_CoresPerNode()/get_cores()[3]);
        
        /*
         index of the node */
        int tmpNode=noml%sumNode;
        int tmpNode1=0;
        
        /*
         using both node0 and node1 */
        if (get_is_node0() && get_is_node1()) {
            if (tmpNode<numNode0) {
                res0
                << "#$ -q all.q@compute-0-"
                << get_node0()[tmpNode] << ".local"
                << "\n";
            }
            else {
                tmpNode1=(tmpNode-numNode0);
                res0
                << "#$ -q all.q@compute-1-"
                << get_node1()[tmpNode1] << ".local"
                << "\n";
            }
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            res0
            << "#$ -q all.q@compute-0-"
            << get_node0()[tmpNode] << ".local"
            << "\n";
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            res0
            << "#$ -q all.q@compute-1-"
            << get_node1()[tmpNode] << ".local"
            << "\n";
        }
        //======================================================================
        // GPU Allocation
        //======================================================================
        /*
         gpuNum is switched between 0 and 1 */
        int gpuNum=0;
        /*
         number of half node cores */
        int hNode=get_CoresPerNode()/2;
        /*
         individual jobs on a node */
        int indv=counter%(get_CoresPerNode()/get_cores()[3]);
        /*
         assign GPU to individual jobs */
        if(sysVar.get_is_GPU()) {
            if (indv<(hNode/get_cores()[3])) gpuNum=0; else gpuNum=1;
        }
        /*
         Mycroft compute-0-16's GPU1 is broken */
        if((get_is_node0())||
           (get_is_node0() && get_is_node1()))
            if(get_node0()[tmpNode]==16) gpuNum=0;
        
        ++counter;
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        res0 << "\n";
        
        
        /* Simulator Executable */
        //----------------------------------------------------------------------
        res0
        << "mpirun -np "
        << get_run_cores()[3] << " "
        << get_lmp_exe();
        //----------------------------------------------------------------------
        
        
        /* Resize Inputfile */
        //----------------------------------------------------------------------
        res0
        << " -in ./lammps_inputs/resize/resize.inp";
        //----------------------------------------------------------------------
        
        
        if (sysVar.get_is_GPU()) res0 << " -sf gpu";
        
        
        /* input varaibles */
        //======================================================================
        res0
        << " -var GPU "         << gpuNum
        << " -var usic "        << sysVar.get_usic()
        << " -var trial 00"     << n_trl
        << " -var namestring "  << sysVar.get_nameString(n_sys)
        << " -var Regime "      << sysVar.get_current_regime();
        
        vector<double> avgboxL=
        autoWork::thermocalc_equ_avgL(sysVar,n_sys,Temp_d);
        
        res0
        << " -var halfX "     << avgboxL[0]/2.0
        << " -var halfY "     << avgboxL[1]/2.0
        << " -var halfZ "     << avgboxL[2]/2.0
        << " -var ts "        << get_timestep_size()
        << " -var steps_res " << (int)(ceil)(get_steps_res())
        << fixed << setprecision(0)
        << " -var Temp "      << Temp_d;
        //======================================================================
        res0.precision(ss);
        res0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /* log file */
        res0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        res0
        << " -log ./resize/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".res.log";
        //----------------------------------------------------------------------
        
        /* screen file */
        //----------------------------------------------------------------------
        res0
        << " > ./resize/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".res.screen" << "\n";
        //----------------------------------------------------------------------
        res0.precision(ss);
        res0 << resetiosflags(ios::fixed|ios::showpoint);
        
        res0.close();
    }
    else cout << "ResFile: 'res.qsub' cannot open." << "\n";
}





void WorkScripts::make_ProductionSubFiles(StructureClass& sysVar,
                                          const LmpScripts& lmp,
                                          const int n_trl,
                                          const int n_sys,
                                          const double Temp_d)
{
    // NOTE:
    // the upper limit for a signed 32-bit interger to store is 2^31
    // No larger than this number should be used in a single run
    double timePerBlock=get_time_equ()/get_n_equ_blocks(); // one tau_alpha
    double run_steps_block=(ceil)(timePerBlock/get_timestep_size());
    double run_steps=get_n_prd_blocks()*run_steps_block;
    
    divide_prd=1;
    if (run_steps>2e+9) {
        
        do {
            ++divide_prd;
        }
        while ((run_steps/(double)divide_prd)>2e+9);
        run_steps /= (double)divide_prd;
        run_steps = (ceil)(run_steps);
    }
    
    bool   is_set_CheckPoint=false;
    double restartf=0;
    if (sysVar.get_systemUnit()=="real")  restartf=1e+6;
    if (sysVar.get_systemUnit()=="metal") restartf=1e+6;
    if (sysVar.get_systemUnit()=="lj")    restartf=1e+7;
    if (run_steps_block>restartf) is_set_CheckPoint=true;
    sysVar.set_is_set_CheckPoint(is_set_CheckPoint);
    
    string o;
    
    //==========================================================================
    /* Exponential Timestepping */
    //==========================================================================
    
    o.clear();
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/production/submission_files");
    o.append("/prd_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".qsub");
    //cout << o << "\n";
    
    ofstream prd0(o.c_str());
    if (prd0.is_open()){
        
        prd0
        << "#!/bin/bash"                                            << "\n"
        << "#$ -V"                                                  << "\n"
        << "#$ -cwd"                                                << "\n"
        << "#$ -j y"                                                << "\n"
        << "#$ -pe orte " << numCores[4]                            << "\n"
        << "#$ -p " << priority[4]                                  << "\n";
        if (sysVar.get_is_GPU()) prd0 << "#$ -R y"                  << "\n";
        
        streamsize ss=prd0.precision();
        prd0 << fixed << setprecision(0);
        
        if (sysVar.get_is_singleTempTest()) {
            
            prd0
            << "#$ -hold_jid gen_"
            << sysVar.get_usic()
            << "_00" << n_trl << "_"
            << sysVar.get_nameString(n_sys)
            << "_T"	<< Temp_d
            << "\n";
        }
        else {
            
            if (lmp.get_is_resize()) {
                
                prd0
                << "#$ -hold_jid res_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
            }
            else {
                
                prd0
                << "#$ -hold_jid equ_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
            }
        }
        
        prd0
        << "#$ -N prd_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        prd0
        << "#$ -o ./production/submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".prd.o"
        << "\n";
        
        prd0.precision(ss);
        prd0 << resetiosflags(ios::fixed|ios::showpoint);
        
        ////////////////////////////////////////////////////////////////////
        //==================================================================
        // Nodes Allocation
        //==================================================================
        /*
         For compute nodes alloccation: within a node, jobs are filled on a half
         node with the same GPU, and then to the next half with the other GPU;
         the same procedure continues in the next requested node. */
        
        int numNode0 = 0; // number of nodes on node0
        int numNode1 = 0; // number of nodes on node1
        
        /*
         using both node0 and node1 */
        if (get_is_node0() && get_is_node1()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = (int)(get_node1().size());
            
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            numNode0 = (int)(get_node0().size());
            numNode1 = 0;
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            numNode0 = 0;
            numNode1 = (int)(get_node1().size());
        }
        
        int sumNode = numNode0 + numNode1;
        
        static int counter=0;
        
        /*
         nonimal number of job allocation on a node; when this nominal number
         increases, it means change to the next node because the current node
         is filled up completely */
        int noml=counter/(get_CoresPerNode()/get_cores()[4]);
        
        /*
         index of the node */
        int tmpNode=noml%sumNode;
        int tmpNode1=0;
        
        /*
         using both node0 and node1 */
        if (get_is_node0() && get_is_node1()) {
            if (tmpNode<numNode0) {
                prd0
                << "#$ -q all.q@compute-0-"
                << get_node0()[tmpNode] << ".local"
                << "\n";
            }
            else {
                tmpNode1=(tmpNode-numNode0);
                prd0
                << "#$ -q all.q@compute-1-"
                << get_node1()[tmpNode1] << ".local"
                << "\n";
            }
        }
        /*
         only using node0 */
        else if (get_is_node0()) {
            prd0
            << "#$ -q all.q@compute-0-"
            << get_node0()[tmpNode] << ".local"
            << "\n";
        }
        /*
         only using node1 */
        else if (get_is_node1()) {
            prd0
            << "#$ -q all.q@compute-1-"
            << get_node1()[tmpNode] << ".local"
            << "\n";
        }
        //==================================================================
        // GPU Allocation
        //==================================================================
        /*
         gpuNum is switched between 0 and 1 */
        int gpuNum=0;
        /*
         number of half node cores */
        int hNode=get_CoresPerNode()/2;
        /*
         individual jobs on a node */
        int indv=counter%(get_CoresPerNode()/get_cores()[4]);
        /*
         assign GPU to individual jobs */
        if(sysVar.get_is_GPU()) {
            if (indv<(hNode/get_cores()[4])) gpuNum=0; else gpuNum=1;
        }
        /*
         Mycroft compute-0-16's GPU1 is broken */
        if((get_is_node0())||
           (get_is_node0() && get_is_node1()))
            if(get_node0()[tmpNode]==16) gpuNum=0;
        
        ++counter;
        //==================================================================
        ////////////////////////////////////////////////////////////////////
        prd0 << "\n";
        
        
        /* Simulator Executable */
        //------------------------------------------------------------------
        prd0
        << "mpirun -np "
        << get_run_cores()[4] << " "
        << get_lmp_exe();
        //------------------------------------------------------------------
        
        
        /* Production Inputfile */
        //------------------------------------------------------------------
        prd0
        << " -in ./lammps_inputs/production/production.inp";
        //------------------------------------------------------------------
        
        
        if (sysVar.get_is_GPU()) prd0 << " -sf gpu";
        
        
        /* input varaibles */
        //==================================================================
        prd0
        << " -var GPU "            << gpuNum
        << " -var usic "           << sysVar.get_usic()
        << " -var trial 00"        << n_trl
        << " -var namestring "     << sysVar.get_nameString(n_sys)
        << " -var Regime "         << sysVar.get_current_regime()
        << " -var ts "             << get_timestep_size()
        << " -var exp_base "       << get_prd_exp_base()
        << " -var n_prd_blocks "   << (int)get_n_prd_blocks()
        << " -var blocksize "      << get_prd_blocksize()
        << " -var divide_prd "     << divide_prd
        << " -var set_CheckPoint " << is_set_CheckPoint
        << fixed << setprecision(0)
        << " -var Temp "           << Temp_d;
        //==================================================================
        prd0.precision(ss);
        prd0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /* log file */
        prd0 << fixed << setprecision(0);
        //------------------------------------------------------------------
        prd0
        << " -log ./production/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".prd.log";
        //------------------------------------------------------------------
        
        /* screen file */
        //------------------------------------------------------------------
        prd0
        << " > ./production/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".prd.screen" << "\n";
        //------------------------------------------------------------------
        prd0.precision(ss);
        prd0 << resetiosflags(ios::fixed|ios::showpoint);
        
        prd0.close();
    }
    else cout << "PrdFile_Exponential: 'prd.qsub' cannot open." << "\n";
    
    
    if (0) {
        
        //======================================================================
        /* Linear Timestepping */
        //======================================================================
        
        for (int run_phase=1; run_phase<=divide_prd; ++run_phase) {
            
            o.clear();
            o.append(return_SimulationFolderPath(sysVar));
            o.append("/production/submission_files");
            o.append("/prd_");
            o.append(to_string((long long int)(run_phase))+"_");
            o.append(to_string((long long int)divide_prd)+"_");
            o.append(sysVar.get_usic());
            o.append("_00"+to_string((long long int)n_trl));
            o.append("_"+sysVar.get_nameString(n_sys));
            o.append("_T"+to_string((long long int)Temp_d));
            o.append(".qsub");
            //cout << o << "\n";
            
            ofstream prd0(o.c_str());
            if (prd0.is_open())	{
                
                prd0
                << "#!/bin/bash"                                        << "\n"
                << "#$ -V"                                              << "\n"
                << "#$ -cwd"                                            << "\n"
                << "#$ -j y"                                            << "\n"
                << "#$ -pe orte " << numCores[4]                        << "\n"
                << "#$ -p " << priority[4]                              << "\n";
                if (sysVar.get_is_GPU()) prd0 << "#$ -R y"              << "\n";
                
                prd0 << fixed << setprecision(0);
                
                if (run_phase==1) {
                    
                    prd0
                    << "#$ -hold_jid equ_"
                    << divide_equ << "_"
                    << divide_equ << "_"
                    << sysVar.get_usic()
                    << "_00" << n_trl << "_"
                    << sysVar.get_nameString(n_sys)
                    << "\n";
                }
                else {
                    
                    prd0
                    << "#$ -hold_jid prd_"
                    << (run_phase-1) << "_"
                    << divide_prd << "_"
                    << sysVar.get_usic()
                    << "_00" << n_trl << "_"
                    << sysVar.get_nameString(n_sys)
                    << "_T"	<< Temp_d
                    << "\n";
                }
                
                prd0
                << "#$ -N prd_"
                << run_phase << "_"
                << divide_prd << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
                
                prd0
                << "#$ -o ./production/submission_files/cluster_out/out_"
                << run_phase << "_"
                << divide_prd << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".prd.o"
                << "\n";
                
                prd0.unsetf(ios_base::floatfield);
                
                ////////////////////////////////////////////////////////////////
                //==============================================================
                // Nodes Allocation
                //==============================================================
                /*
                 For compute nodes alloccation: within a node, jobs are filled on a half
                 node with the same GPU, and then to the next half with the other GPU;
                 the same procedure continues in the next requested node. */
                
                int numNode0 = 0; // number of nodes on node0
                int numNode1 = 0; // number of nodes on node1
                
                /*
                 using both node0 and node1 */
                if (get_is_node0() && get_is_node1()) {
                    numNode0 = (int)(get_node0().size());
                    numNode1 = (int)(get_node1().size());
                    
                }
                /*
                 only using node0 */
                else if (get_is_node0()) {
                    numNode0 = (int)(get_node0().size());
                    numNode1 = 0;
                }
                /*
                 only using node1 */
                else if (get_is_node1()) {
                    numNode0 = 0;
                    numNode1 = (int)(get_node1().size());
                }
                
                int sumNode = numNode0 + numNode1;
                
                static int counter=0;
                
                /*
                 nonimal number of job allocation on a node; when this nominal number
                 increases, it means change to the next node because the current node
                 is filled up completely */
                int noml=counter/(get_CoresPerNode()/get_cores()[4]);
                
                /*
                 index of the node */
                int tmpNode=noml%sumNode;
                int tmpNode1=0;
                
                /*
                 using both node0 and node1 */
                if (get_is_node0() && get_is_node1()) {
                    if (tmpNode<numNode0) {
                        prd0
                        << "#$ -q all.q@compute-0-"
                        << get_node0()[tmpNode] << ".local"
                        << "\n";
                    }
                    else {
                        tmpNode1=(tmpNode-numNode0);
                        prd0
                        << "#$ -q all.q@compute-1-"
                        << get_node1()[tmpNode1] << ".local"
                        << "\n";
                    }
                }
                /*
                 only using node0 */
                else if (get_is_node0()) {
                    prd0
                    << "#$ -q all.q@compute-0-"
                    << get_node0()[tmpNode] << ".local"
                    << "\n";
                }
                /*
                 only using node1 */
                else if (get_is_node1()) {
                    prd0
                    << "#$ -q all.q@compute-1-"
                    << get_node1()[tmpNode] << ".local"
                    << "\n";
                }
                //==============================================================
                // GPU Allocation
                //==============================================================
                /*
                 gpuNum is switched between 0 and 1 */
                int gpuNum=0;
                /*
                 number of half node cores */
                int hNode=get_CoresPerNode()/2;
                /*
                 individual jobs on a node */
                int indv=counter%(get_CoresPerNode()/get_cores()[4]);
                /*
                 assign GPU to individual jobs */
                if(sysVar.get_is_GPU()) {
                    if (indv<(hNode/get_cores()[4])) gpuNum=0; else gpuNum=1;
                }
                /*
                 Mycroft compute-0-16's GPU1 is broken */
                if((get_is_node0())||
                   (get_is_node0() && get_is_node1()))
                    if(get_node0()[tmpNode]==16) gpuNum=0;
                
                ++counter;
                //==============================================================
                ////////////////////////////////////////////////////////////////
                prd0 << "\n";
                
                
                /* Simulator Executable */
                //--------------------------------------------------------------
                prd0
                << "mpirun -np "
                << get_run_cores()[4] << " "
                << get_lmp_exe();
                //--------------------------------------------------------------
                
                
                /* production Inputfile */
                //--------------------------------------------------------------
                prd0
                << " -in ./lammps_inputs/production/production_linear.inp";
                //--------------------------------------------------------------
                
                
                if (sysVar.get_is_GPU()) prd0 << " -sf gpu";
                
                
                /* input varaibles */
                //==============================================================
                prd0
                << " -var GPU "          << gpuNum
                << " -var usic "         << sysVar.get_usic()
                << " -var trial 00"      << n_trl
                << " -var namestring "   << sysVar.get_nameString(n_sys)
                << " -var Regime "       << sysVar.get_current_regime()
                << " -var ts "           << get_timestep_size()
                << " -var steps_prd "    << (int)run_steps
                << " -var run_phase "    << run_phase
                << " -var divide_prd "   << divide_prd
                << fixed << setprecision(0)
                << " -var Temp "         << Temp_d;
                
                prd0 << " -var read_res ";
                string target;
                if (run_phase==1) {
                    
                    target.append(return_SimulationFolderPath(sysVar));
                    target.append("/equilibration/restart");
                    target.append("/restart_");
                    target.append(sysVar.get_usic());
                    target.append("_00"+to_string((long long int)n_trl));
                    target.append("_"+sysVar.get_nameString(n_sys));
                    target.append("_T"+to_string((long long int)Temp_d));
                    target.append(".equ.restart");
                }
                else {
                    
                    target.append(return_SimulationFolderPath(sysVar));
                    target.append("/production/restart");
                    target.append("/restart_");
                    target.append(sysVar.get_usic());
                    target.append("_00"+to_string((long long int)n_trl));
                    target.append("_"+sysVar.get_nameString(n_sys));
                    target.append("_T"+to_string((long long int)Temp_d));
                    target.append(".prd");
                    target.append("."+to_string((long long int)(run_phase-1)));
                    target.append("."+to_string((long long int)divide_prd));
                    target.append(".restart");
                }
                prd0 << target;
                
                
                prd0 << " -var write_res ";
                target.clear();
                if (run_phase<divide_prd) {
                    
                    target.append(return_SimulationFolderPath(sysVar));
                    target.append("/production/restart");
                    target.append("/restart_");
                    target.append(sysVar.get_usic());
                    target.append("_00"+to_string((long long int)n_trl));
                    target.append("_"+sysVar.get_nameString(n_sys));
                    target.append("_T"+to_string((long long int)Temp_d));
                    target.append(".prd");
                    target.append("."+to_string((long long int)run_phase));
                    target.append("."+to_string((long long int)divide_prd));
                    target.append(".restart");
                }
                else {
                    
                    target.append(return_SimulationFolderPath(sysVar));
                    target.append("/production/restart");
                    target.append("/restart_");
                    target.append(sysVar.get_usic());
                    target.append("_00"+to_string((long long int)n_trl));
                    target.append("_"+sysVar.get_nameString(n_sys));
                    target.append("_T"+to_string((long long int)Temp_d));
                    target.append(".prd.restart");
                }
                prd0 << target;
                //==============================================================
                prd0.unsetf(ios_base::floatfield);
                
                
                /* log file */
                prd0 << fixed << setprecision(0);
                //--------------------------------------------------------------
                prd0
                // reset precision to zero for file naming
                << " -log ./production/log/"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".prd.log";
                //--------------------------------------------------------------
                
                /* screen file */
                //--------------------------------------------------------------
                prd0
                << " > ./production/screen/"
                << "prd_"
                << run_phase << "_"
                << divide_prd << "_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".screen" << "\n";
                //--------------------------------------------------------------
                prd0.unsetf(ios_base::floatfield);
                
                prd0.close();
            }
            else cout << "PrdFile_Linear: 'prd.qsub' cannot open." << "\n";
        }
    }
    
}





void WorkScripts::make_GenerationSubScripts(const StructureClass& sysVar) const
{
    const int n_trial   = sysVar.get_n_trial();
    const int n_system  = sysVar.get_n_system();
    const int n_sys_beg = sysVar.get_n_sys_beg();
    const int n_sys_end = sysVar.get_n_sys_end();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/genSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    ofstream gen1(o.c_str());
    if (gen1.is_open())	{
        
        gen1 << "#!/bin/bash" << "\n";
        gen1 << "\n";
        
        /* Trial */
        //----------------------------------------------------------------------
        gen1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            gen1 << i;
            if (i!=(n_trial-1)) gen1 << " ";
        }
        gen1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Namestring */
        //----------------------------------------------------------------------
        gen1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            gen1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) gen1 << " ";
        }
        gen1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Temperature */
        streamsize ss=gen1.precision();
        gen1 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        gen1 << "Temp=(";
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                gen1
                << "\""
                << sysVar.get_initialtemps()[index][0]
                << "\"";
                if (n_trl!=(n_trial-1)) gen1 << " ";
            }
        }
        gen1 << ")" << "\n";
        //----------------------------------------------------------------------
        gen1.precision(ss);
        gen1 << resetiosflags(ios::fixed|ios::showpoint);
        
        gen1 << "\n";
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        gen1
        << "cd "
        << return_SimulationFolderPath(sysVar) << "\n";
        
        gen1 << "\n";
        
        gen1
        << "for i in ${trial[@]}; do"                               << "\n"
        << "  for ii in ${namestring[@]}; do"                       << "\n"
        << "    for iii in ${Temp[$i]}; do"                         << "\n"
        << "      qsub_File=./generation/submission_files/"
        
        << "gen_"
        << sysVar.get_usic()
        << "_00${i}"
        << "_${ii}"
        << "_T${iii}"
        << ".qsub"                                                  << "\n"
        
        << "      test -e ${qsub_File}"                             << "\n"
        << "      if [ \"$?\" -eq \"0\" ]; then"                    << "\n"
        << "        qsub ${qsub_File}"                              << "\n"
        << "      else"                                             << "\n"
        << "        continue"                                       << "\n"
        << "      fi"                                               << "\n"
        << "    done"                                               << "\n"
        << "  done"                                                 << "\n"
        << "done"                                                   << "\n";
        //----------------------------------------------------------------------
        
        gen1.close();
    }
    else cout << "GenScript: 'gen.sh' cannot open." << "\n";
    
    
    //======================================================
    // Direct submission of all generation jobs
    //======================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/genSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void WorkScripts::make_QuenchSubScripts(const StructureClass& sysVar) const
{
    const int n_trial  = sysVar.get_n_trial();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/qchSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    ofstream qch1(o.c_str());
    if (qch1.is_open())	{
        
        qch1 << "#!/bin/bash" << "\n";
        qch1 << "\n";
        
        /* Trial */
        //----------------------------------------------------------------------
        qch1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            qch1 << i;
            if (i!=(n_trial-1)) qch1 << " ";
        }
        qch1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Namestring */
        //----------------------------------------------------------------------
        qch1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            qch1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) qch1 << " ";
        }
        qch1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        qch1 << "\n";
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        qch1
        << "cd "
        << return_SimulationFolderPath(sysVar) << "\n";
        
        qch1 << "\n";
        
        qch1
        << "for i in ${trial[@]}; do"                               << "\n"
        << "  for ii in ${namestring[@]}; do"                       << "\n"
        << "      qsub_File=./quench/submission_files/"
        
        << "qch_"
        << sysVar.get_usic()
        << "_00${i}"
        << "_${ii}"
        << "_Regime" << sysVar.get_current_regime()
        << ".qsub"                                                  << "\n"
        
        << "      test -e ${qsub_File}"                             << "\n"
        << "      if [ \"$?\" -eq \"0\" ]; then"                    << "\n"
        << "        qsub ${qsub_File}"                              << "\n"
        << "      else"                                             << "\n"
        << "        continue"                                       << "\n"
        << "      fi"                                               << "\n"
        << "  done"                                                 << "\n"
        << "done"                                                   << "\n";
        //----------------------------------------------------------------------
        
        
        qch1.close();
    }
    else cout << "QchScript: 'qch.sh' cannot open." << "\n";
    
    
    //======================================================
    // Direct submission of all quench jobs
    //======================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/qchSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void WorkScripts::make_EquilibraitionSubScripts(const StructureClass& sysVar,
                                                const std::vector<std::vector<double>> tinfo) const
{
    const int n_trial   = sysVar.get_n_trial();
    const int n_system  = sysVar.get_n_system();
    const int n_sys_beg = sysVar.get_n_sys_beg();
    const int n_sys_end = sysVar.get_n_sys_end();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/equSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    ofstream equ1(o.c_str());
    if (equ1.is_open())	{
        
        equ1 << "#!/bin/bash" << "\n";
        equ1 << "\n";
        
        /* Run Phase (if single run > 2e+9 steps) */
        //----------------------------------------------------------------------
        if (divide_equ>1) {
            equ1 << "runPhase=(";
            for (int i=1; i<=divide_equ; ++i) {
                equ1 << i;
                if (i!=divide_equ) equ1 << " ";
            }
            equ1 << ")" << "\n";
        }
        //----------------------------------------------------------------------
        
        /* Trial */
        //----------------------------------------------------------------------
        equ1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            equ1 << i;
            if (i!=(n_trial-1)) equ1 << " ";
        }
        equ1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Namestring */
        //----------------------------------------------------------------------
        equ1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            equ1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) equ1 << " ";
        }
        equ1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Temperature */
        streamsize ss=equ1.precision();
        equ1 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        equ1 << "Temp=(";
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                int n_Temp=(int)tinfo[index].size();
                
                equ1 << "\"";
                for (int T=0; T<n_Temp; ++T) {
                    equ1
                    << tinfo[index][T];
                    if (T!=(n_Temp-1)) equ1 << " ";
                }
                if (n_trl!=(n_trial-1)) equ1 << "\"" << " ";
            }
        }
        equ1
        << "\""
        << ")"
        << "\n";
        //----------------------------------------------------------------------
        equ1.precision(ss);
        equ1 << resetiosflags(ios::fixed|ios::showpoint);
        
        equ1 << "\n";
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        equ1
        << "cd "
        << return_SimulationFolderPath(sysVar) << "\n";
        
        equ1 << "\n";
        
        if (divide_equ>1) {
            equ1
            << "for k in ${runPhase[@]}; do"                            << "\n"
            << "for i in ${trial[@]}; do"                               << "\n"
            << "  for ii in ${namestring[@]}; do"                       << "\n"
            << "    for iii in ${Temp[$i]}; do"                         << "\n"
            << "      qsub_File=./equilibration/submission_files/"
            
            << "equ_"
            << "${k}_"
            << divide_equ << "_"
            << sysVar.get_usic()
            << "_00${i}"
            << "_${ii}"
            << "_T${iii}"
            << ".qsub"
            << "\n"
            
            << "      test -e ${qsub_File}"                             << "\n"
            << "      if [ \"$?\" -eq \"0\" ]; then"                    << "\n"
            << "        qsub ${qsub_File}"                              << "\n"
            << "      else"                                             << "\n"
            << "        continue"                                       << "\n"
            << "      fi"                                               << "\n"
            << "    done"                                               << "\n"
            << "  done"                                                 << "\n"
            << "done"                                                   << "\n"
            << "done"                                                   << "\n";
        }
        else {
            equ1
            << "for i in ${trial[@]}; do"                               << "\n"
            << "  for ii in ${namestring[@]}; do"                       << "\n"
            << "    for iii in ${Temp[$i]}; do"                         << "\n"
            << "      qsub_File=./equilibration/submission_files/"
            
            << "equ_"
            << sysVar.get_usic()
            << "_00${i}"
            << "_${ii}"
            << "_T${iii}"
            << ".qsub"
            << "\n"
            
            << "      test -e ${qsub_File}"                             << "\n"
            << "      if [ \"$?\" -eq \"0\" ]; then"                    << "\n"
            << "        qsub ${qsub_File}"                              << "\n"
            << "      else"                                             << "\n"
            << "        continue"                                       << "\n"
            << "      fi"                                               << "\n"
            << "    done"                                               << "\n"
            << "  done"                                                 << "\n"
            << "done"                                                   << "\n";
        }
        //----------------------------------------------------------------------
        
        equ1.close();
    }
    else cout << "EquScript: 'equ.sh' cannot open." << "\n";
    
    
    //======================================================
    // Direct submission of all equilibration jobs
    //======================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/equSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void WorkScripts::make_ResizeSubScripts(const StructureClass& sysVar,
                                        const std::vector<std::vector<double>> tinfo) const
{
    const int n_trial   = sysVar.get_n_trial();
    const int n_system  = sysVar.get_n_system();
    const int n_sys_beg = sysVar.get_n_sys_beg();
    const int n_sys_end = sysVar.get_n_sys_end();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/resSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    ofstream res1(o.c_str());
    if (res1.is_open())	{
        
        res1 << "#!/bin/bash" << "\n";
        res1 << "\n";
        
        /* Trial */
        //----------------------------------------------------------------------
        res1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            res1 << i;
            if (i!=(n_trial-1)) res1 << " ";
        }
        res1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Namestring */
        //----------------------------------------------------------------------
        res1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            res1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) res1 << " ";
        }
        res1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Temperature */
        streamsize ss=res1.precision();
        res1 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        res1 << "Temp=(";
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                int n_Temp=(int)tinfo[index].size();
                
                res1 << "\"";
                for (int T=0; T<n_Temp; ++T) {
                    res1
                    << tinfo[index][T];
                    if (T!=(n_Temp-1)) res1 << " ";
                }
                if (n_trl!=(n_trial-1)) res1 << "\"" << " ";
            }
        }
        res1
        << "\""
        << ")"
        << "\n";
        //----------------------------------------------------------------------
        res1.precision(ss);
        res1 << resetiosflags(ios::fixed|ios::showpoint);
        
        res1 << "\n";
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        res1
        << "cd "
        << return_SimulationFolderPath(sysVar) << "\n";
        
        res1 << "\n";
        
        res1
        << "for i in ${trial[@]}; do"                                   << "\n"
        << "  for ii in ${namestring[@]}; do"                           << "\n"
        << "    for iii in ${Temp[$i]}; do"                             << "\n"
        << "      qsub_File=./resize/submission_files/"
        
        << "res_"
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
        //----------------------------------------------------------------------
        
        res1.close();
    }
    else cout << "ResScript: 'res.sh' cannot open." << "\n";
    
    
    //======================================================
    // Direct submission of all resize jobs
    //======================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/resSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void WorkScripts::make_ProductionSubScripts(const StructureClass& sysVar,
                                            const std::vector<std::vector<double>> tinfo) const
{
    const int n_trial   = sysVar.get_n_trial();
    const int n_system  = sysVar.get_n_system();
    const int n_sys_beg = sysVar.get_n_sys_beg();
    const int n_sys_end = sysVar.get_n_sys_end();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/prdSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    //cout << o << "\n";
    
    ofstream prd1(o.c_str());
    if (prd1.is_open())	{
        
        prd1 << "#!/bin/bash" << "\n";
        prd1 << "" << "\n";
        
        /* Trial */
        //----------------------------------------------------------------------
        prd1 << "trial=(";
        for (int i=0; i<n_trial; ++i) {
            prd1 << i;
            if (i!=(n_trial-1)) prd1 << " ";
        }
        prd1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Namestring */
        //----------------------------------------------------------------------
        prd1 << "namestring=(";
        for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
            prd1 << sysVar.get_nameString(i);
            if (i!=sysVar.get_n_sys_end()) prd1 << " ";
        }
        prd1 << ")" << "\n";
        //----------------------------------------------------------------------
        
        /* Temperature */
        streamsize ss=prd1.precision();
        prd1 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        prd1 << "Temp=(";
        for (int n_trl=0; n_trl<n_trial; ++n_trl) {
            for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
                
                int index=n_trl*n_system+(n_sys-n_sys_beg);
                
                int n_Temp=(int)tinfo[index].size();
                
                prd1 << "\"";
                for (int T=0; T<n_Temp; ++T) {
                    prd1
                    << tinfo[index][T];
                    if (T!=(n_Temp-1)) prd1 << " ";
                }
                if (n_trl!=(n_trial-1)) prd1 << "\"" << " ";
            }
        }
        prd1
        << "\""
        << ")"
        << "\n";
        //----------------------------------------------------------------------
        prd1.precision(ss);
        prd1 << resetiosflags(ios::fixed|ios::showpoint);
        
        prd1 << "\n";
        
        /* make 'simulations' folder the current working directory */
        //----------------------------------------------------------------------
        prd1
        << "cd "
        << return_SimulationFolderPath(sysVar) << "\n";
        
        prd1 << "\n";
        
        prd1
        << "for i in ${trial[@]}; do"                               << "\n"
        << "  for ii in ${namestring[@]}; do"                       << "\n"
        << "    for iii in ${Temp[$i]}; do"                         << "\n"
        << "      qsub_File=./production/submission_files/"
        
        << "prd_"
        << sysVar.get_usic()
        << "_00${i}"
        << "_${ii}"
        << "_T${iii}"
        << ".qsub"
        << "\n"
        
        << "      test -e ${qsub_File}"                             << "\n"
        << "      if [ \"$?\" -eq \"0\" ]; then"                    << "\n"
        << "        qsub ${qsub_File}"                              << "\n"
        << "      else"                                             << "\n"
        << "        continue"                                       << "\n"
        << "      fi"                                               << "\n"
        << "    done"                                               << "\n"
        << "  done"                                                 << "\n"
        << "done"                                                   << "\n";
        //----------------------------------------------------------------------
        
        prd1.close();
    }
    else cout << "PrdScript: 'prd.sh' cannot open." << "\n";
    
    
    //======================================================
    // Direct submission of all production jobs
    //======================================================
    if (sysVar.get_is_directSub()==true) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/prdSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) exit(EXIT_FAILURE);
    }
}





void WorkScripts::echoNodeUsage()
{
    cout << "==============================\n";
    string str=echo_NodesUsage(get_node0(),get_node1());
    call_system_bash(str);
    cout << "==============================\n";
}





/*==( public setters )==*/

/* string */
void WorkScripts::set_lmp_exe(const string& s){lmp_exe=s;}

/* bool */
void WorkScripts::set_is_node0(const bool b){is_node0=b;}
void WorkScripts::set_is_node1(const bool b){is_node1=b;}
/* double */
void WorkScripts::set_n_equ_blocks(const double d){n_equ_blocks=d;}
void WorkScripts::set_n_prd_blocks(const double d){n_prd_blocks=d;}
void WorkScripts::set_prd_exp_base(const double d){prd_exp_base=d;}
void WorkScripts::set_timestep_size(const double d){timestep_size=d;}
void WorkScripts::set_time_equ(const double d){time_equ=d;}
void WorkScripts::set_quenchRate(const double d){quenchRate=d;}
void WorkScripts::set_steps_gen(const double d){steps_gen=d;}
void WorkScripts::set_steps_qch(const double d){steps_qch=d;}
void WorkScripts::set_steps_res(const double d){steps_res=d;}
/* vector<int> */
void WorkScripts::set_node0(const vector<int>& vi){node0=vi;}
void WorkScripts::set_node1(const vector<int>& vi){node1=vi;}
void WorkScripts::set_numCores(const vector<int>& vi){numCores=vi;}
void WorkScripts::set_cores(const vector<int>& vi){cores=vi;}
void WorkScripts::set_run_cores(const std::vector<int>& vi){run_cores=vi;}
void WorkScripts::set_priority(const vector<int>& vi){priority=vi;}

/*==( public getters )==*/

double WorkScripts::get_steps_equ() const
{
    double steps_equ=(ceil)(get_time_equ()/get_timestep_size());
    
    if (steps_equ<1e+2) steps_equ=1e+2; // limit for extremely short equil.
    
    return steps_equ;
}

int WorkScripts::get_prd_blocksize() const
{
    double timePerBlock;
    int n_steps;
    int exponent;
    int blocksize;
    
    timePerBlock=get_time_equ()/get_n_equ_blocks(); // one tau_alpha length
    n_steps=(int)(ceil)(timePerBlock/get_timestep_size());
    
    if (n_steps>2e+9) n_steps=(int)2e+9;
    
    exponent=(int)(ceil)(log((double)n_steps)/log(get_prd_exp_base()));
    
    if (n_steps<1e+8) blocksize=exponent+1;
    else blocksize=exponent;
    
    
    
    static int count_access=0;
    
    if (blocksize==102) { // fix for Martini PS30
        ++count_access;
        if ((typeB==0)&&(typeS==0)) {
            if (count_access==1) {
                cout
                << "NOTE: blocksize==102" << "\n"
                << "fix for MartiniPS30, blocksize==103" << "\n";
            }
            blocksize=103;
        }
    }
    if (blocksize==105) { // fix for stdFENE KG
        ++count_access;
        if ((typeB==0)&&(typeS==5)) {
            if (count_access==1) {
                cout
                << "NOTE: blocksize==105" << "\n"
                << "fix for KGPolymer, blocksize==106" << "\n";
            }
            blocksize=106;
        }
    }
    
    return blocksize;
}
