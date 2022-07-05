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

#include "workscripts.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

WorkScripts::WorkScripts(const StructureClass& sysVar):
/** bool **/
is_node0(false),
is_node1(false),
is_node2(false),
/** int **/
indexi(0),
indexii(0),
typeB(sysVar.get_typeB()),
typeS(sysVar.get_typeS()),
current_regime(sysVar.get_current_regime()),
n_cores_node(16),
n_gpus_node(2),
prd_blocksize_hiT(1),
prd_blocksize_loT(1),
divide_equ(1),
divide_prd(1),
gpuNum(0),
/** double **/
n_equ_blocks(sysVar.get_n_equ_blocks()),
n_prd_blocks(sysVar.get_n_prd_blocks()),
n_relaxation(sysVar.get_n_relaxation()),
prd_exp_base(1.2),
timestep_size(sysVar.get_ts_regime().at(current_regime)),
time_equ(1.0),
quenchRate(1000.0),
steps_gen(1e+6),
steps_qch(0),
steps_equ(0),
steps_res(5e+4),
/** string **/
lmp_exe(""),
/** STL **/
node0(),
node1(),
node2(),
numCores(),
cores(),
run_cores(),
priority({0,0,0,0})
{
    if (sysVar.get_is_GPU()) {
        lmp_exe="/share/apps/lammps-26Jan2014-GPU2/bin/lmp_openmpi";//GPU LAMMPS
    } else {
        lmp_exe="/share/apps/lammps/bin/lmp_openmpi";//CPU LAMMPS
    }
    
    /** Quench **/
    //--------------------------------------------------------------------------
    if (sysVar.get_systemUnit()=="real") {
        quenchRate = 1000.0;// kelvins per nano second
    } else if (sysVar.get_systemUnit()=="lj") {
        quenchRate = 1e-4;  // T per tau
    } else if (sysVar.get_systemUnit()=="metal") {
        quenchRate = 1000.0;// kelvins per nano second
    }
    
    /** Equilibration **/
    //--------------------------------------------------------------------------
    double previous_teq,current_teq,teq;
    if (true) {
        teq=sysVar.get_equilibration_times().at(current_regime);
    } else {
        if (current_regime>0) {
            previous_teq = sysVar.get_equilibration_times().at(current_regime-1);
            current_teq  = sysVar.get_equilibration_times().at(current_regime);
            teq=max(previous_teq,current_teq);
        } else {
            teq=sysVar.get_equilibration_times().at(current_regime);
        }
    } time_equ=teq;
    
    if (n_equ_blocks<0) {
        cout
        << "WARNING: n_equ_blocks cannot be < 0. Please check.\n";
        exit(EXIT_FAILURE);
    }
    else if (n_equ_blocks<1) {
        if (!sysVar.get_is_aging()) {
            cout
            << "WARNING: n_equ_blocks < 1. Program resets n_equ_blocks to 1.\n";
            system_wait(5); n_equ_blocks=1;
        }
    }
    
    /** Production **/
    //--------------------------------------------------------------------------
    if (n_prd_blocks<0) {
        cout
        << "WARNING: n_prd_blocks cannot be < 0. Please check.\n";
        exit(EXIT_FAILURE);
    } else if (n_prd_blocks<1) {
        cout
        << "WARNING: n_prd_block < 1. Program resets n_prd_block to 1.\n";
        system_wait(5); n_prd_blocks=1;
    }
}





void WorkScripts::setup_sge_env(const StructureClass& sysVar,
                                const string& phase,
                                const int counter,
                                std::ofstream& writefile)
{
    string jobname;
    int cores_d=0;
    
    if      (phase=="generation")    cores_d=numCores.at(0);
    else if (phase=="quench")        cores_d=numCores.at(1);
    else if (phase=="equilibration") cores_d=numCores.at(2);
    else if (phase=="resize")        cores_d=numCores.at(3);
    else if (phase=="production")    cores_d=numCores.at(4);
    
    writefile
    << "#!/bin/bash"                                                << "\n"
    << "#$ -V"                                                      << "\n"
    << "#$ -cwd"                                                    << "\n"
    << "#$ -j y"                                                    << "\n"
    << "#$ -pe orte " << cores_d                                    << "\n"
    << "#$ -p " << priority[1]                                      << "\n";
    if (sysVar.get_is_GPU()) writefile << "#$ -R y"                 << "\n";
    
    writefile
    << "#$ -N "<< jobname << "\n";
    
    if (phase=="quench") {
        writefile
        << "#$ -o ./quench/submission_files/cluster_out/out_"
        << jobname
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.o"
        << "\n";
    }
    
}





void WorkScripts::node_allocation(const StructureClass& sysVar,
                                  const string& phase,
                                  const int counter,
                                  std::ofstream& writefile,
                                  const int n_trl)
{
    int cores_d=0;
    
    if      (phase=="generation")    cores_d=cores.at(0);
    else if (phase=="quench")        cores_d=cores.at(1);
    else if (phase=="equilibration") cores_d=cores.at(2);
    else if (phase=="resize")        cores_d=cores.at(3);
    else if (phase=="production")    cores_d=cores.at(4);
    
    /* at least one node needs to be used */
    if (!(is_node0||is_node1||is_node2)) {
        cout
        << "\n"
        << "All compute nodes are switched off! Please switch on at least one."
        << "\n\n";
        exit(EXIT_FAILURE);
    }
    /* node2 can NOT be mixed in use with node0 or node1 */
    if ((is_node0||is_node1)&&(is_node2)) {
        cout
        << "\n"
        << "Node-2 can NOT be used together with Node-0 or Node-1. "
        << "Please re-specify the compute nodes.\n\n";
        exit(EXIT_FAILURE);
    }
    
    /** Compute Nodes Allocation ("Global Node Level") **/
    //--------------------------------------------------------------------------
    /** The total number of available compute nodes is represented by 'sumNode';
     ** a nominal value (noml) is then used to represent the number of nodes
     ** that have been filled; say, if noml=2, it means two of the availabe
     ** nodes have been filled and the current job allocation is on the third
     ** available node (if there is one) **/
    
    int numNode0=0;//number of nodes used on node0
    int numNode1=0;//number of nodes used on node1
    int numNode2=0;//number of nodes used on node2
    
    if (is_node0) numNode0=(int)(node0.size());
    if (is_node1) numNode1=(int)(node1.size());
    if (is_node2) numNode2=(int)(node2.size());
    
    if (is_node2) {
        n_cores_node = 24;//24 CPU cores per node-2
        n_gpus_node  =  4;// 4 GPU cores per node-2
    } else {
        n_cores_node = 16;
        n_gpus_node  =  2;
    }
    
    // sumNodes is the total number of nodes available for job allocation
    int sumNode=numNode0+numNode1+numNode2;
    int noml=counter/(n_cores_node/cores_d);
    int tmpNode=noml%sumNode;
    int tmpNode1=0;
    
    /** Node-0 and Node-1 assignment **/
    if (is_node0&&is_node1) { /* using both node0 and node1 */
        if (tmpNode<numNode0) {
            writefile << "#$ -q all.q@compute-0-" << node0.at(tmpNode) << ".local\n";
        } else {
            tmpNode1=(tmpNode-numNode0);
            writefile << "#$ -q all.q@compute-1-" << node1.at(tmpNode1) << ".local\n";
        }
    } else if (is_node0) { /* only using node0 */
        writefile << "#$ -q all.q@compute-0-" << node0.at(tmpNode) << ".local\n";
    } else if (is_node1) { /* only using node1 */
        writefile << "#$ -q all.q@compute-1-" << node1.at(tmpNode) << ".local\n";
    }
    /** Node-2 assignment **/
    if (is_node2) { /* only using node2 */
        writefile << "#$ -q all.q@compute-2-" << node2.at(tmpNode) << ".local\n";
    } writefile << "\n";
    
    
    /** GPU# Allocation ("Local Node Level" or "Per Node Basis") **/
    //--------------------------------------------------------------------------
    /** GPU indices are assigned based on partial node;
     ** Original desgin (for Node-0 and Node-1) was based on 2 GPUs per node:
     ** jobs that fill the first half of a node are assigned to GPU0, the rest
     ** that fill the second half assigned to GPU1.
     ** For Node-2, the allocation rule appiles to a quarter of a node:
     ** first quarter jobs to GPU0, second quarter to GPU1, and so on **/
    
    /* number of partial node cores */
    int pNode=n_cores_node/n_gpus_node;//8 for Node-0 and Node-1, 6 for Node-2
    /* individual jobs on a node */
    int indv=counter%(n_cores_node/cores_d);
    /* assign GPU# on Node-0 and Node-1 */
    if (is_node0||is_node1) {
        if (indv<(pNode/cores_d)) {
            gpuNum=0;//GPU0
        } else {
            gpuNum=1;//GPU1
        }
    }
    /* assign GPU# on Node-2 */
    if (is_node2) {
        if (indv<(pNode/cores_d)) {
            gpuNum=0;//GPU0
        } else if (indv<2*(pNode/cores_d)) {
            gpuNum=1;//GPU1
        } else if (indv<3*(pNode/cores_d)) {
            gpuNum=2;//GPU2
        } else {
            gpuNum=3;//GPU3
        }
    }
    /* specify a particular GPU if it is requested (for the only trial) */
    if (sysVar.get_is_specify_GPU())
    {
        int gpu=sysVar.get_which_GPU();
        if (is_node0||is_node1)
        {
            if (gpu!=0&&gpu!=1) {
                cout << "GPU has to be 0 or 1.\n";
                exit(EXIT_FAILURE);
            }
        }
        if (is_node2)
        {
            if (gpu!=0&&gpu!=1&&gpu!=2&&gpu!=3) {
                cout << "GPU has to be 0, 1, 2, or 3.\n";
                exit(EXIT_FAILURE);
            }
        }
        gpuNum=sysVar.get_which_GPU();
    }
    /** Exception:
     ** Mycroft compute-0-16's GPU1 is broken (fixed) **/
    //if(!is_node1)if(is_node0)if(node0.at(tmpNode)==16) gpuNum=0;
}





void WorkScripts::append_adhoc_PATH(std::ofstream& subscript)
{
    subscript
    << "PATH_org=${PATH}\n"
    << "LD_LIBRARY_PATH_org=${LD_LIBRARY_PATH}\n"
    << "PATH=/usr/local/cuda/bin:${PATH}\n"
    << "LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:${LD_LIBRARY_PATH}\n"
    << "\n";
}
void WorkScripts::clean_adhoc_PATH(std::ofstream& subscript)
{
    subscript
    << "\n"
    << "PATH=${PATH_org}\n"
    << "LD_LIBRARY_PATH=${LD_LIBRARY_PATH_org}\n";
}




void WorkScripts::make_GenerationSubFiles(StructureClass& sysVar,
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
        
        /** Compute node and GPU allocation **/
        //----------------------------------------------------------------------
        static int counter=0;
        node_allocation(sysVar,"generation",counter,gen0,n_trl);
        ++counter;
        
        /** Simulation executable **/
        //----------------------------------------------------------------------
        gen0 << "mpirun -np " << run_cores.at(0) << " " << lmp_exe;
        
        /** Simulation inputfile **/
        //----------------------------------------------------------------------
        gen0 << " -in ./lammps_inputs/generation/generation.inp";
        
        /** LAMMPS switch for GPU package **/
        //----------------------------------------------------------------------
        if (sysVar.get_is_GPU()) {
            gen0 << " -sf gpu";
            gen0 << " -var is_GPU 1";
        } else {
            gen0 << " -var is_GPU 0";
        }
        
        /* Random number generator for vseed */
        ////////////////////////////////////////////////////////////////////////
        //======================================================================
        int vseed=randint(100000,1000000); // rand integer number in (1e5,1e6)
        //======================================================================
        ////////////////////////////////////////////////////////////////////////
        
        /** input varaibles **/
        //======================================================================
        gen0
        << " -var GPU "        << gpuNum
        << " -var usic "       << sysVar.get_usic()
        << " -var trial 00"    << n_trl
        << " -var namestring " << sysVar.get_nameString(n_sys)
        << " -var n_poly "     << sysVar.get_n_poly()
        << " -var vseed "      << vseed
        << " -var ts "         << timestep_size
        << " -var steps_gen "  << (int)(ceil)(steps_gen);
        
        if (sysVar.get_types_all().size()>0) {
            gen0 << " -var types_all";
            for (indexi=0; indexi<(int)sysVar.get_types_all().size(); ++indexi) {
                gen0 << " " << sysVar.get_types_all().at(indexi);
            }
        }
        if (sysVar.get_types_heavy().size()>0) {
            gen0 << " -var types_heavy";
            for (indexi=0; indexi<(int)sysVar.get_types_heavy().size(); ++indexi) {
                gen0 << " " << sysVar.get_types_heavy().at(indexi);
            }
        }
        if (sysVar.get_types_light().size()>0) {
            gen0 << " -var types_light";
            for (indexi=0; indexi<(int)sysVar.get_types_light().size(); ++indexi) {
                gen0 << " " << sysVar.get_types_light().at(indexi);
            }
        }
        
        gen0 << fixed << setprecision(0)
        << " -var Temp "       << Temp_d;
        //======================================================================
        gen0.precision(ss);
        gen0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /** log file path **/
        gen0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        gen0
        << " -log ./generation/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".gen.log";
        
        /** screen file path **/
        //----------------------------------------------------------------------
        gen0
        << " > ./generation/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".gen.screen"
        << "\n";
        
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
            << "_T"	<< sysVar.get_initialtemps().at(index)[0]
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
        
        /** Compute node and GPU allocation **/
        //----------------------------------------------------------------------
        static int counter=0;
        node_allocation(sysVar,"quench",counter,qch0,n_trl);
        ++counter;
        
        /** Simulation executable **/
        //----------------------------------------------------------------------
        qch0 << "mpirun -np " << run_cores.at(1) << " " << lmp_exe;
        
        /** Simulation inputfile **/
        //----------------------------------------------------------------------
        qch0 << " -in ./lammps_inputs/quench/quench.inp";
        
        /** LAMMPS switch for GPU package **/
        //----------------------------------------------------------------------
        if (sysVar.get_is_GPU()) {
            qch0 << " -sf gpu";
            qch0 << " -var is_GPU 1";
        } else {
            qch0 << " -var is_GPU 0";
        }
        
        //======================================================================
        /** Deciding Quenching Temperatures **/
        //======================================================================
        int n_qchTemps=0;       // # of quenching T's
        vector<double> vecA;    // stores higher T's for quenching
        vector<double> vecB;    // stores lower T's for quenching
        double steps_qch = 0;   // # of quench steps
        double qchTdiff  = 0;   // T difference of highest/lowest qch T's
        
        //----------------------------------------------------------------------
        // Full Quench
        // Direct quench from highest to lowest Ts (all previously set)
        //----------------------------------------------------------------------
        if (sysVar.get_is_fullquench()) {
            
            double hiT=sysVar.get_initialtemps().at(index)[0];
            double loT=sysVar.get_finalTemp();
            
            qchTdiff=fabs(hiT-loT);
            
            double currentT=hiT;
            while (currentT>loT) {
                vecA.push_back(currentT);
                ++n_qchTemps;
                currentT -= sysVar.get_precision();
            } ++n_qchTemps; // total number of quenching temperatures
            
            currentT=hiT-sysVar.get_precision();
            while (currentT>=loT) {
                vecB.push_back(currentT);
                currentT -= sysVar.get_precision();
            }
        }
        //----------------------------------------------------------------------
        // Stepwise Quench
        // Quench from highest to lowest Ts bound by the tau_alpha targets
        //----------------------------------------------------------------------
        else {
            
            vector<double> quenchingTs;
            if (sysVar.get_current_regime()==0) {
                
                quenchingTs=sysVar.get_temperaturesInfo().at(index);
                
                /** set T difference in quenching **/
                size_t qchTsize=quenchingTs.size();
                qchTdiff=fabs(quenchingTs[0]-quenchingTs.at(qchTsize-1));
                sysVar.get_qchTdiff().at(index)=qchTdiff;
                
            } else {
                
                /** in-equilibrium temperatures of previous regime **/
                size_t equSize = sysVar.get_equilibratedTs().at(index).size();
                double equTmin = sysVar.get_equilibratedTs().at(index).at(equSize-1);
                
                /** lowest in-equilibrium T of previous regime **/
                quenchingTs.push_back(equTmin);
                
                /** put quench T's into the temperature container **/
                for (size_t i=0; i<sysVar.get_quenchingTs().at(index).size(); ++i) {
                    quenchingTs.push_back(sysVar.get_quenchingTs().at(index).at(i));
                }
                
                /** difference of highest and lowest quench temperatures **/
                size_t qchTsize=quenchingTs.size();
                qchTdiff=fabs(quenchingTs[0]-quenchingTs.at(qchTsize-1));
                sysVar.get_qchTdiff().at(index)=qchTdiff;
            }
            
            /** setup quench temperature arrays **/
            n_qchTemps=int(quenchingTs.size());
            for (int i=0; i<(n_qchTemps-1); ++i) {
                vecA.push_back(quenchingTs.at(i));
            }
            for (int i=1; i<n_qchTemps; ++i) {
                vecB.push_back(quenchingTs.at(i));
            }
        }
        
        qchTdiff *= pow(sysVar.get_corrFacforT(),-1);
        
        if (sysVar.get_systemUnit()=="real") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*timestep_size*1e-6));
        } else if (sysVar.get_systemUnit()=="metal") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*timestep_size*1e-3));
        } else if (sysVar.get_systemUnit()=="lj") {
            steps_qch = (ceil)(qchTdiff/(quenchRate*timestep_size));
        }
        
        if (steps_qch<1e+3) steps_qch=1e+3;
        
        if (sysVar.get_is_aging()) {
            set_steps_qch(0);
        } else {
            set_steps_qch(steps_qch);
        }
        
        /** input varaibles **/
        //======================================================================
        qch0
        << " -var GPU "         << gpuNum
        << " -var usic "        << sysVar.get_usic()
        << " -var trial 00"     << n_trl
        << " -var namestring "  << sysVar.get_nameString(n_sys)
        << " -var Regime "      << sysVar.get_current_regime()
        << " -var n_intervals " << (n_qchTemps-1)
        << " -var ts "          << timestep_size
        << " -var steps_qch "   << (int)(ceil)(fabs(get_steps_qch()));
        
        qch0 << fixed << setprecision(0);
        
        qch0
        << " -var TempA ";
        for (size_t i=0; i<vecA.size(); ++i) {
            if (i!=(vecA.size()-1)) {
                qch0 << vecA.at(i) << " ";
            } else {
                qch0 << vecA.at(i);
            }
        }
        qch0
        << " -var TempB ";
        for (size_t i=0; i<vecB.size(); ++i) {
            if (i!=(vecB.size()-1)) {
                qch0 << vecB.at(i) << " ";
            } else {
                qch0 << vecB.at(i);
            }
        }
        
        if (sysVar.get_current_regime()==0) {
            qch0
            << " -var readRestartPath "
            << "./generation/restart/restart_"
            << sysVar.get_usic()
            << "_00" << n_trl
            << "_"   << sysVar.get_nameString(n_sys)
            << "_T"  << sysVar.get_initialtemps().at(index)[0]
            << ".gen.restart";
        } else {
            size_t vecSize = sysVar.get_equilibratedTs().at(index).size();
            double restartT= sysVar.get_equilibratedTs().at(index).at(vecSize-1);
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
        
        /** log file path **/
        qch0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        qch0
        << " -log ./quench/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.log";
        
        /** screen file **/
        //----------------------------------------------------------------------
        qch0
        << " > ./quench/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_Regime" << sysVar.get_current_regime()
        << ".qch.screen"
        << "\n";
        
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
    double run_steps=get_steps_equ(sysVar);
    if (run_steps>2e+9) run_steps=2e+9;
    //if (sysVar.get_systemUnit()=="real")  if (run_steps>1e+9) run_steps=1e+9;
    //if (sysVar.get_systemUnit()=="metal") if (run_steps>1e+9) run_steps=1e+9;
    //if (sysVar.get_systemUnit()=="lj")    if (run_steps>2e+9) run_steps=2e+9;
    
    //NOTE: below correction is not changing the run steps in equilibration
    divide_equ=1;
    if (run_steps>2e+9)
    {
        do {
            ++divide_equ;
        } while ((run_steps/(double)divide_equ)>2e+9);
        run_steps /= (double)divide_equ;
        run_steps  = (ceil)(run_steps);
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
        } else {
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
                } else {
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
                
            } else {
                
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
            
            /** Compute node and GPU allocation **/
            //------------------------------------------------------------------
            static int counter=0;
            node_allocation(sysVar,"equilibration",counter,equ0,n_trl);
            ++counter;
            
            /** Simulation executable **/
            //------------------------------------------------------------------
            equ0 << "mpirun -np " << run_cores.at(2) << " " << lmp_exe;
            
            /** Simulation inputfile **/
            //------------------------------------------------------------------
            equ0 << " -in ./lammps_inputs/equilibration/equilibration.inp";
            
            /** LAMMPS switch for GPU package **/
            //------------------------------------------------------------------
            if (sysVar.get_is_GPU()) {
                equ0 << " -sf gpu";
                equ0 << " -var is_GPU 1";
            } else {
                equ0 << " -var is_GPU 0";
            }
            
            /** input varaibles **/
            //==================================================================
            equ0
            << " -var GPU "            << gpuNum
            << " -var usic "           << sysVar.get_usic()
            << " -var trial 00"        << n_trl
            << " -var namestring "     << sysVar.get_nameString(n_sys)
            << " -var Regime "         << sysVar.get_current_regime()
            << " -var ts "             << timestep_size
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
            } else {
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
            } equ0 << target;
            
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
            } else {
                target.append(return_SimulationFolderPath(sysVar));
                target.append("/equilibration/restart");
                target.append("/restart_");
                target.append(sysVar.get_usic());
                target.append("_00"+to_string((long long int)n_trl));
                target.append("_"+sysVar.get_nameString(n_sys));
                target.append("_T"+to_string((long long int)Temp_d));
                target.append(".equ.restart");
            } equ0 << target;
            
            // NOTE:
            // In LAMMPS, the timesteps N must fit in a signed 32-bit integer
            // so limited to ~ 2 billion steps (2^31) in a "single run"
            //==================================================================
            equ0.precision(ss);
            equ0 << resetiosflags(ios::fixed|ios::showpoint);
            
            /** log file path **/
            equ0 << fixed << setprecision(0);
            //------------------------------------------------------------------
            equ0
            << " -log ./equilibration/log/"
            << sysVar.get_usic()
            << "_00" << n_trl << "_"
            << sysVar.get_nameString(n_sys)
            << "_T" << Temp_d
            << ".equ.log";
            
            /** screen file path **/
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
            } else {
                equ0
                << " > ./equilibration/screen/"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".equ.screen" << "\n";
            }
            
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
        
        /** Compute node and GPU allocation **/
        //----------------------------------------------------------------------
        static int counter=0;
        node_allocation(sysVar,"resize",counter,res0,n_trl);
        ++counter;
        
        /** Simulation executable **/
        //----------------------------------------------------------------------
        res0 << "mpirun -np " << run_cores.at(3) << " " << lmp_exe;
        
        /** Simulation inputfile **/
        //----------------------------------------------------------------------
        res0 << " -in ./lammps_inputs/resize/resize.inp";
        
        /** LAMMPS switch for GPU package **/
        //----------------------------------------------------------------------
        if (sysVar.get_is_GPU()) {
            res0 << " -sf gpu";
            res0 << " -var is_GPU 1";
        } else {
            res0 << " -var is_GPU 0";
        }
        
        /** input varaibles **/
        //======================================================================
        res0
        << " -var GPU "         << gpuNum
        << " -var usic "        << sysVar.get_usic()
        << " -var trial 00"     << n_trl
        << " -var namestring "  << sysVar.get_nameString(n_sys)
        << " -var Regime "      << sysVar.get_current_regime();
        
        double avgboxL=thermocalc_equ(sysVar,n_sys,Temp_d,"Lx");//assume cubic
        
        res0
        << " -var halfX "     << avgboxL/2.0
        << " -var halfY "     << avgboxL/2.0
        << " -var halfZ "     << avgboxL/2.0
        << " -var ts "        << timestep_size
        << " -var steps_res " << (int)(ceil)(steps_res)
        << fixed << setprecision(0)
        << " -var Temp "      << Temp_d;
        //======================================================================
        res0.precision(ss);
        res0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /** log file path **/
        res0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        res0
        << " -log ./resize/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".res.log";
        
        /** screen file path **/
        //----------------------------------------------------------------------
        res0
        << " > ./resize/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".res.screen" << "\n";
        
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
    double timePerBlock=time_equ/get_n_equ_blocks(); // one tau_alpha
    double run_steps_block=(ceil)(timePerBlock/timestep_size);
    double run_steps=get_n_prd_blocks()*run_steps_block;
    
    divide_prd=1;
    if (run_steps>2e+9)
    {
        do {
            ++divide_prd;
        } while ((run_steps/(double)divide_prd)>2e+9);
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
    
    /** Exponential Timestepping **/
    
    string o;
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
            
        } else {
            
            if (lmp.get_is_resize()) {
                prd0
                << "#$ -hold_jid res_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T"	<< Temp_d
                << "\n";
            } else {
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
        
        /** Compute node and GPU allocation **/
        //----------------------------------------------------------------------
        static int counter=0;
        node_allocation(sysVar,"production",counter,prd0,n_trl);
        ++counter;
        
        /** Simulation executable **/
        //----------------------------------------------------------------------
        prd0 << "mpirun -np " << run_cores.at(4) << " " << lmp_exe;
        
        /** Simulation inputfile **/
        //----------------------------------------------------------------------
        prd0 << " -in ./lammps_inputs/production/production.inp";
        
        /** LAMMPS switch for GPU package **/
        //----------------------------------------------------------------------
        if (sysVar.get_is_GPU()) {
            prd0 << " -sf gpu";
            prd0 << " -var is_GPU 1";
        } else {
            prd0 << " -var is_GPU 0";
        }
        
        /** input varaibles **/
        //======================================================================
        prd0
        << " -var GPU "            << gpuNum
        << " -var usic "           << sysVar.get_usic()
        << " -var trial 00"        << n_trl
        << " -var namestring "     << sysVar.get_nameString(n_sys)
        << " -var Regime "         << sysVar.get_current_regime()
        << " -var ts "             << timestep_size
        << " -var exp_base "       << prd_exp_base
        << " -var n_prd_blocks "   << (int)get_n_prd_blocks()
        << " -var blocksize "      << get_prd_blocksize()
        << " -var divide_prd "     << divide_prd
        << " -var set_CheckPoint " << is_set_CheckPoint
        << fixed << setprecision(0)
        << " -var Temp "           << Temp_d;
        //======================================================================
        prd0.precision(ss);
        prd0 << resetiosflags(ios::fixed|ios::showpoint);
        
        
        /** log file path **/
        prd0 << fixed << setprecision(0);
        //------------------------------------------------------------------
        prd0
        << " -log ./production/log/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".prd.log";
        
        /** screen file path **/
        //------------------------------------------------------------------
        prd0
        << " > ./production/screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".prd.screen" << "\n";
        
        prd0.precision(ss);
        prd0 << resetiosflags(ios::fixed|ios::showpoint);
        
        prd0.close();
    }
    else cout << "PrdFile_Exponential: 'prd.qsub' cannot open." << "\n";
    
    
    if (0) { /** Linear Timestepping **/
        
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
                } else {
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
                
                /** Compute node and GPU allocation **/
                //--------------------------------------------------------------
                static int counter=0;
                node_allocation(sysVar,"production",counter,prd0,n_trl);
                ++counter;
                
                /** Simulation executable **/
                //--------------------------------------------------------------
                prd0 << "mpirun -np " << run_cores.at(4) << " " << lmp_exe;
                
                /** Simulation inputfile **/
                //--------------------------------------------------------------
                prd0 << " -in ./lammps_inputs/production/production.inp";
                
                /** LAMMPS switch for GPU package **/
                //--------------------------------------------------------------
                if (sysVar.get_is_GPU()) {
                    prd0 << " -sf gpu";
                    prd0 << " -var is_GPU 1";
                } else {
                    prd0 << " -var is_GPU 0";
                }
                
                /** input varaibles **/
                //==============================================================
                prd0
                << " -var GPU "          << gpuNum
                << " -var usic "         << sysVar.get_usic()
                << " -var trial 00"      << n_trl
                << " -var namestring "   << sysVar.get_nameString(n_sys)
                << " -var Regime "       << sysVar.get_current_regime()
                << " -var ts "           << timestep_size
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
                } else {
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
                } prd0 << target;
                
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
                } else {
                    target.append(return_SimulationFolderPath(sysVar));
                    target.append("/production/restart");
                    target.append("/restart_");
                    target.append(sysVar.get_usic());
                    target.append("_00"+to_string((long long int)n_trl));
                    target.append("_"+sysVar.get_nameString(n_sys));
                    target.append("_T"+to_string((long long int)Temp_d));
                    target.append(".prd.restart");
                } prd0 << target;
                
                prd0.unsetf(ios_base::floatfield);
                
                /** log file path **/
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
                
                /** screen file path **/
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
                
                prd0.unsetf(ios_base::floatfield);
                
                prd0.close();
            }
            else cout << "PrdFile_Linear: 'prd.qsub' cannot open." << "\n";
        }
    }
    
}





void WorkScripts::make_GenerationSubScripts(const StructureClass& sysVar)
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
    ofstream gen1(o.c_str());
    
    gen1 << "#!/bin/bash" << "\n\n";
    
    /* Trial */
    //--------------------------------------------------------------------------
    gen1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        gen1 << i;
        if (i!=(n_trial-1)) gen1 << " ";
    } gen1 << ")\n";
    
    /* Namestring */
    //--------------------------------------------------------------------------
    gen1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        gen1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) gen1 << " ";
    } gen1 << ")\n";
    
    /* Temperature */
    //--------------------------------------------------------------------------
    streamsize ss=gen1.precision();
    gen1 << fixed << setprecision(0);
    gen1 << "Temp=(";
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            int index=n_trl*n_system+(n_sys-n_sys_beg);
            gen1
            << "\""
            << sysVar.get_initialtemps().at(index)[0]
            << "\"";
            if (n_trl!=(n_trial-1)) gen1 << " ";
        }
    } gen1 << ")\n\n";
    gen1.precision(ss);
    gen1 << resetiosflags(ios::fixed|ios::showpoint);
    
    //append_adhoc_PATH(gen1);
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //clean_adhoc_PATH(gen1);
    
    gen1.close();
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/genSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in WorkScripts::make_GenerationSubScripts():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void WorkScripts::make_QuenchSubScripts(const StructureClass& sysVar)
{
    const int n_trial  = sysVar.get_n_trial();
    
    string o;
    o.append(return_SimulationFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/qchSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    ofstream qch1(o.c_str());
    
    qch1 << "#!/bin/bash" << "\n\n";
    
    /* Trial */
    //--------------------------------------------------------------------------
    qch1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        qch1 << i;
        if (i!=(n_trial-1)) qch1 << " ";
    } qch1 << ")\n";
    
    /* Namestring */
    //--------------------------------------------------------------------------
    qch1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        qch1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) qch1 << " ";
    } qch1 << ")\n\n";
    
    //append_adhoc_PATH(qch1);
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //clean_adhoc_PATH(qch1);
    
    qch1.close();
    
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/qchSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        //cout << sub.c_str() << "\n";
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in WorkScripts::make_QuenchSubScripts():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void WorkScripts::make_EquilibraitionSubScripts(const StructureClass& sysVar,
                                                const std::vector<std::vector<double>> tinfo)
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
    ofstream equ1(o.c_str());
    
    equ1 << "#!/bin/bash" << "\n\n";
    
    /* Run Phase (if single run > 2e+9 steps) */
    //--------------------------------------------------------------------------
    if (divide_equ>1) {
        
        cout << "\n"
        << "in WorkScripts::make_EquilibraitionSubScripts():\n"
        << "divide_equ>1, # of run phase = " << divide_equ << "\n\n";
        
        equ1 << "runPhase=(";
        for (int i=1; i<=divide_equ; ++i) {
            equ1 << i;
            if (i!=divide_equ) equ1 << " ";
        } equ1 << ")\n";
    }
    
    /* Trial */
    //--------------------------------------------------------------------------
    equ1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        equ1 << i;
        if (i!=(n_trial-1)) equ1 << " ";
    } equ1 << ")\n";
    
    /* Namestring */
    //--------------------------------------------------------------------------
    equ1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        equ1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) equ1 << " ";
    } equ1 << ")\n";
    
    /* Temperature */
    //--------------------------------------------------------------------------
    streamsize ss=equ1.precision();
    equ1 << fixed << setprecision(0);
    equ1 << "Temp=(";
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            int index=n_trl*n_system+(n_sys-n_sys_beg);
            
            int n_Temp=(int)tinfo.at(index).size();
            
            equ1 << "\"";
            for (int T=0; T<n_Temp; ++T) {
                equ1
                << tinfo.at(index)[T];
                if (T!=(n_Temp-1)) equ1 << " ";
            }
            if (n_trl!=(n_trial-1)) equ1 << "\"" << " ";
        }
    } equ1 << "\")\n\n";
    equ1.precision(ss);
    equ1 << resetiosflags(ios::fixed|ios::showpoint);
    
    //append_adhoc_PATH(equ1);
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //clean_adhoc_PATH(equ1);
    
    equ1.close();
    
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/equSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in WorkScripts::make_EquilibraitionSubScripts():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void WorkScripts::make_ResizeSubScripts(const StructureClass& sysVar,
                                        const std::vector<std::vector<double>> tinfo)
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
    ofstream res1(o.c_str());
    
    res1 << "#!/bin/bash" << "\n\n";
    
    /* Trial */
    //--------------------------------------------------------------------------
    res1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        res1 << i;
        if (i!=(n_trial-1)) res1 << " ";
    } res1 << ")\n";
    
    /* Namestring */
    //--------------------------------------------------------------------------
    res1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        res1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) res1 << " ";
    } res1 << ")\n";
    
    /* Temperature */
    //--------------------------------------------------------------------------
    streamsize ss=res1.precision();
    res1 << fixed << setprecision(0);
    res1 << "Temp=(";
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            int index=n_trl*n_system+(n_sys-n_sys_beg);
            
            int n_Temp=(int)tinfo.at(index).size();
            
            res1 << "\"";
            for (int T=0; T<n_Temp; ++T) {
                res1
                << tinfo.at(index)[T];
                if (T!=(n_Temp-1)) res1 << " ";
            }
            if (n_trl!=(n_trial-1)) res1 << "\"" << " ";
        }
    } res1 << "\")\n\n";
    res1.precision(ss);
    res1 << resetiosflags(ios::fixed|ios::showpoint);
    
    //append_adhoc_PATH(res1);
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //clean_adhoc_PATH(res1);
    
    res1.close();
    
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/resSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in WorkScripts::make_ResizeSubScripts():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void WorkScripts::make_ProductionSubScripts(const StructureClass& sysVar,
                                            const std::vector<std::vector<double>> tinfo)
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
    ofstream prd1(o.c_str());
    
    prd1 << "#!/bin/bash" << "\n\n";
    
    /* Trial */
    //--------------------------------------------------------------------------
    prd1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        prd1 << i;
        if (i!=(n_trial-1)) prd1 << " ";
    } prd1 << ")\n";
    
    /* Namestring */
    //--------------------------------------------------------------------------
    prd1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        prd1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) prd1 << " ";
    } prd1 << ")\n";
    
    /* Temperature */
    //--------------------------------------------------------------------------
    streamsize ss=prd1.precision();
    prd1 << fixed << setprecision(0);
    prd1 << "Temp=(";
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=n_sys_beg; n_sys<=n_sys_end; ++n_sys) {
            
            int index=n_trl*n_system+(n_sys-n_sys_beg);
            
            int n_Temp=(int)tinfo.at(index).size();
            
            prd1 << "\"";
            for (int T=0; T<n_Temp; ++T) {
                prd1
                << tinfo.at(index)[T];
                if (T!=(n_Temp-1)) prd1 << " ";
            }
            if (n_trl!=(n_trial-1)) prd1 << "\"" << " ";
        }
    } prd1 << "\")\n\n";
    prd1.precision(ss);
    prd1 << resetiosflags(ios::fixed|ios::showpoint);
    
    //append_adhoc_PATH(prd1);
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //clean_adhoc_PATH(prd1);
    
    prd1.close();
    
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_SimulationFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/prdSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in WorkScripts::make_ProductionSubScripts():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





int WorkScripts::get_keyWordIndex(const std::string& phase,const std::string& keyWord)
{
    int index=0;
    if (phase=="equilibration_log")
    {
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
        if (keyWord=="Step")         index=0;
        else if (keyWord=="Temp")    index=1;
        else if (keyWord=="Press")   index=2;
        else if (keyWord=="PotEng")  index=3;
        else if (keyWord=="KinEng")  index=4;
        else if (keyWord=="TotEng")  index=5;
        else if (keyWord=="E_bond")  index=6;
        else if (keyWord=="E_angle") index=7;
        else if (keyWord=="E_pair")  index=8;
        else if (keyWord=="Lx")      index=9;
        else if (keyWord=="Ly")      index=10;
        else if (keyWord=="Lz")      index=11;
        else if (keyWord=="Volume")  index=12;
        else if (keyWord=="Density") index=13;
        else if (keyWord=="Dt")      index=14;
        else if (keyWord=="Time")    index=15;
        else if (keyWord=="CPU")     index=16;
        else if (keyWord=="T/CPU")   index=17;
        else if (keyWord=="S/CPU")   index=18;
        else if (keyWord=="CPULeft") index=19;
        else index=-1;
    }
    if (phase=="production_log")
    {
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
        if (keyWord=="Step")         index=0;
        else if (keyWord=="Temp")    index=1;
        else if (keyWord=="Press")   index=2;
        else if (keyWord=="PotEng")  index=3;
        else if (keyWord=="KinEng")  index=4;
        else if (keyWord=="TotEng")  index=5;
        else if (keyWord=="E_bond")  index=6;
        else if (keyWord=="E_angle") index=7;
        else if (keyWord=="E_pair")  index=8;
        else if (keyWord=="Lx")      index=9;
        else if (keyWord=="Ly")      index=10;
        else if (keyWord=="Lz")      index=11;
        else if (keyWord=="Volume")  index=12;
        else if (keyWord=="Density") index=13;
        else if (keyWord=="Dt")      index=14;
        else if (keyWord=="Time")    index=15;
        else if (keyWord=="CPU")     index=16;
        else if (keyWord=="T/CPU")   index=17;
        else if (keyWord=="S/CPU")   index=18;
        else index=-1;
    }
    return index;
}





double WorkScripts::thermocalc_equ(const StructureClass& sysVar,
                                   const int n_sys,
                                   const double Temp_d,
                                   const string& keyWord)
{
    int     countFiles=0; // count of files used for average
    double  value=0,avgValue=0;
    SysInfo calc_stats;
    string  in,out;
    vector<double> stats;
    vector<double> calc,calc_tmp;
    
    /* # of columns of data to be stored */
    int n_columns_of_data=14;
    int keyWordIndex=get_keyWordIndex("equilibration_log",keyWord);
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl)
    {
        /** equilibration log file **/
        in.clear();
        in.append(return_SimulationFolderPath(sysVar));
        in.append("/equilibration/log/");
        in.append(sysVar.get_usic());
        in.append("_00"+to_string((long long int)n_trl));
        in.append("_"+sysVar.get_nameString(n_sys));
        in.append("_T"+to_string((long long int)Temp_d));
        in.append(".equ.log");
        ifstream readFile(in.c_str());
        if (readFile.is_open())
        {
            ++countFiles; // if the file exists
            
            bool   is_inside_thermoSection=false;
            int    line_index=0;
            string lineContent;
            string strData0;
            
            vector<double> vecData; // used for storing line content
            while (getline(readFile,lineContent))
            {
                ++line_index;
                istringstream iss(lineContent);
                iss >> strData0; // the 1st deliminated string of line
                /** Loop is the keyWord right after the thermo section **/
                if (strData0=="Loop") is_inside_thermoSection=false;
                if (is_inside_thermoSection) {
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0); iss >> vecData.at(i);
                    }
                    /** offset by 1 because the first column is for strData0 **/
                    value=vecData.at(keyWordIndex-1);
                    calc_tmp.push_back(value);//store "all" values
                    vecData.clear();//clear container content for next use
                }
                /** Step is the keyWord right before the thermo section **/
                if (strData0=="Step") is_inside_thermoSection=true;
            } readFile.clear();
            /** get average from second half of equilibration **/
            for (indexi=(int)(calc_tmp.size()-1)/2;indexi<(int)calc_tmp.size();++indexi) {
                //for size=51, from index=25 to index=50 (n=26)
                calc.push_back(calc_tmp.at(indexi));
            } calc_stats.set_calcvector(calc);
            calc_tmp.clear(); calc.clear();
            double avg=calc_stats.calc_mean();
            avgValue+=avg;
        } readFile.close();
    }//looping over all trials
    
    if (countFiles==0) {
        return -1;
    } else {
        return avgValue/=(double)countFiles;
    }
}





void WorkScripts::thermodataprocess_prd(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys,
                                        const double Temp_d)
{
    /** Modified: SJH (20170429) **/
    
    string in,out;
    /** input **/
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/log/");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.log");
    ifstream readFile(in.c_str());
    /** output **/
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/log/");
    out.append("new_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.log");
    ofstream writeFile(out.c_str());
    if (readFile.is_open())
    {
        bool   is_thermo=false;
        int    now=0,pre=0;
        int    countaccess=0;
        string lineContent;
        string strData0;
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            
            /** write out thermo values **/
            if (is_thermo && strData0!="Loop") {
                now=stoi(strData0);//current step number
                if ((countaccess==1)||(now!=pre)) {
                    writeFile << lineContent << "\n";
                } pre=now;
            }
            
            /** turn on/ff thermo section switch **/
            if (strData0=="Step") {
                ++countaccess;
                // write out the thermo_style items only one time
                if (countaccess==1) {
                    writeFile << lineContent << "\n";
                } is_thermo=true;
            } else if (strData0=="Loop") {
                is_thermo=false;
            }
        }
        readFile.close();
        writeFile.close();
    } else {
        cout
        << "in WorkScripts::thermodataprocess_prd():\n"
        << in+" cannot open.\n";
        //exit(EXIT_FAILURE);
        return;
    }
}





void WorkScripts::thermocalc_prd(const StructureClass& sysVar,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d)
{
    double  value=0;
    string  keyWord;
    SysInfo calc_stats;
    vector<double> stats,calc;
    
    vector<string> allkeyWord;
    allkeyWord.push_back("Density");
    allkeyWord.push_back("Volume");
    allkeyWord.push_back("Volume2");
    allkeyWord.push_back("TotEng");
    allkeyWord.push_back("Lx");
    //allkeyWord.push_back("Ly");
    //allkeyWord.push_back("Lz");
    allkeyWord.push_back("Temp");
    
    /** # of columns of data to be stored **/
    int n_columns_of_data=14;
    int keyWordIndex=0;
    
    /** input **/
    string in,out;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/log/");
    in.append("new_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.log");
    ifstream readFile(in.c_str());
    auto start_position=readFile.tellg();//NOTE: mark start
    if (readFile.is_open())
    {
        for (indexi=0; indexi<allkeyWord.size(); ++indexi)
        {
            keyWord=allkeyWord[indexi];
            if (keyWord=="Volume2") {
                keyWordIndex=get_keyWordIndex("production_log","Volume");
            } else {
                keyWordIndex=get_keyWordIndex("production_log",keyWord);
            }
            /** output **/
            out.clear();
            out.append(return_AnalysisFolderPath(sysVar));
            out.append("/fit_data/");
            out.append("thermo_");
            out.append(keyWord+"_");
            out.append(sysVar.get_usic());
            out.append("_00"+to_string((long long int)n_trl));
            out.append("_"+sysVar.get_nameString(n_sys));
            out.append(".dat");
            ofstream writeFile(out.c_str(),ofstream::app); // 'append'
            
            bool   is_inside_thermoSection=false;
            int    line_index=0;
            string lineContent,strData0;
            
            readFile.clear();//NOTE: reset flag
            readFile.seekg(start_position);//NOTE: rewind
            vector<double> vecData;//used for storing line content
            while (getline(readFile,lineContent))
            {
                ++line_index;
                istringstream iss(lineContent);
                iss >> strData0; // 1st deliminated string of input
                if (is_inside_thermoSection) {
                    for (int i=0; i<n_columns_of_data; ++i) {
                        vecData.push_back(0); iss >> vecData.at(i);
                    }
                    // offset by 1 because the first column is for strData0
                    // so basically the first item is wasted, i.e. "Step"
                    value=vecData.at(keyWordIndex-1);
                    if (keyWord=="Volume2") {
                        calc.push_back(pow(value,2));
                    } else {
                        calc.push_back(value);//store per "keyWord" values
                    } vecData.clear();//clear container content for next use
                }
                /** Step is the keyWord right before the thermo section **/
                if (strData0=="Step") is_inside_thermoSection=true;
            }
            calc_stats.set_calcvector(calc);
            calc.clear();
            
            double avg=calc_stats.calc_mean();
            double stdev=calc_stats.get_sample_stdev();
            
            writeFile
            << Temp_d*pow(sysVar.get_corrFacforT(),-1) //actual
            << " "<<avg<<" "<<stdev<<"\n";
            writeFile.close();
            
        } readFile.close();
        
    } else {
        cout
        << "in WorkScripts::thermocalc_prd():\n"
        << in+" cannot open.\n";
        //exit(EXIT_FAILURE);
        return;
    }
    /** keep the original production log file **/
    bool is_keep_original_prd=true;
    if(!is_keep_original_prd){
        string rm="rm "+in; system(rm.c_str());
    }
}





/** public setters **/
//------------------------------------------------------------------------------
/* string */
void WorkScripts::set_lmp_exe(const string& s){lmp_exe=s;}

/* bool */
void WorkScripts::set_is_node0(const bool b){is_node0=b;}
void WorkScripts::set_is_node1(const bool b){is_node1=b;}
void WorkScripts::set_is_node2(const bool b){is_node2=b;}
/* int */
void WorkScripts::set_n_cores_node(const int i){n_cores_node=i;}
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
void WorkScripts::set_node2(const vector<int>& vi){node2=vi;}
void WorkScripts::set_numCores(const vector<int>& vi){numCores=vi;}
void WorkScripts::set_cores(const vector<int>& vi){cores=vi;}
void WorkScripts::set_run_cores(const std::vector<int>& vi){run_cores=vi;}
void WorkScripts::set_priority(const vector<int>& vi){priority=vi;}





/** public getters **/
//------------------------------------------------------------------------------
double WorkScripts::get_steps_equ(const SysInfo& sysVar) const
{
    double steps_equ=(ceil)(time_equ/timestep_size);
    if (sysVar.get_is_aging()) {
        steps_equ=1;//if is_aging, do instantaneous equilibration
    } else {
        if (steps_equ<100) steps_equ=100;       //in case of extremely short run
        else if (steps_equ>2e+9) steps_equ=2e+9;//in case of extremely long  run
    } return steps_equ;
}

int WorkScripts::get_prd_blocksize() const
{
    int    n_steps=0,exponent=0,blocksize=0;
    double timePerBlock=0;
    
    /* Each block is a relaxation unit long (tau_alpha or tau_rouse) */
    timePerBlock=time_equ/n_equ_blocks;
    
    n_steps=(int)(ceil)(timePerBlock/timestep_size);
    if (n_steps>2e+9) n_steps=(int)2e+9;
    exponent=(int)(ceil)(log((double)n_steps)/log(prd_exp_base));
    if (n_steps<1e+8) blocksize=exponent+1;
    else blocksize=exponent;
    
    if (false) //OLD style
    {
        static int count_access=0;
        if (blocksize==102) { // fix for Martini PS30
            ++count_access;
            if ((typeB==0)&&(typeS==0)) {
                if (count_access==1) {
                    cout
                    << "NOTE: blocksize==102" << "\n"
                    << "fix for MartiniPS30, blocksize==103" << "\n";
                } blocksize=103;
            }
        }
        if (blocksize==105) { // fix for stdFENE KG
            ++count_access;
            if ((typeB==0)&&(typeS==5)) {
                if (count_access==1) {
                    cout
                    << "NOTE: blocksize==105" << "\n"
                    << "fix for KGPolymer, blocksize==106" << "\n";
                } blocksize=106;
            }
        }
    }
    return blocksize;
}




