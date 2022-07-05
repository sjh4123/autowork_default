//
//  Created by Sean Jh H. on 6/21/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef WORKSCRIPTS_H
#define WORKSCRIPTS_H

#include "structureclass.h"
#include "lmpscripts.h"

namespace autoWork {
    
    class WorkScripts
    {
        std::string lmp_exe;
        
        bool is_node0;
        bool is_node1;
        
        int typeB;         // type of backbonedouble
        int typeS;         // type of side group
        
        int CoresPerNode; // default = 16
        
        int prd_blocksize_hiT;
        int prd_blocksize_loT;
        
        int divide_equ;
        int divide_prd;
        
        double n_equ_blocks;
        double n_prd_blocks;
        
        double prd_exp_base; // default = 1.2
        double timestep_size;
        double time_equ;
        double quenchRate;
        
        double steps_gen;
        double steps_qch;
        double steps_equ;
        double steps_res; // default = 5e+4
        
        std::vector<int> node0;
        std::vector<int> node1;
        std::vector<int> numCores;
        std::vector<int> cores;
        std::vector<int> run_cores;
        std::vector<int> priority;  // default = {0,0,0,0}
        
    public:
        
        /* Constructors */
        WorkScripts()=delete;
        WorkScripts(const StructureClass&);
        WorkScripts(const WorkScripts&)=default;
        
        /* Assignment */
        WorkScripts& operator= (const WorkScripts&)=default;
        
        /* Destructor */
        ~WorkScripts()=default;
        
        
        //======================================================
        // make submission files of each simulation phase
        //======================================================
        // submission file contains information of simulation
        
        void make_GenerationSubFiles(const StructureClass&,
                                     const LmpScripts&,
                                     const int n_trl,
                                     const int n_sys,
                                     const double T);
        
        void make_QuenchSubFiles(StructureClass&, // T changes w/ regime
                                 const int n_trl,
                                 const int n_sys);
        
        void make_EquilibrationSubfiles(const StructureClass&,
                                        const LmpScripts&,
                                        const int n_trl,
                                        const int n_sys,
                                        const double T);
        
        void make_ResizeSubfiles(const StructureClass&,
                                 const LmpScripts&,
                                 const int n_trl,
                                 const int n_sys,
                                 const double T);
        
        void make_ProductionSubFiles(StructureClass&, // is_set_CheckPoint
                                     const LmpScripts&,
                                     const int n_trl,
                                     const int n_sys,
                                     const double T);
        
        
        //======================================================
        // make job submisstion scripts of each simulation phase
        //======================================================
        // submission script is for submitting jobs
        
        void make_GenerationSubScripts(const StructureClass&) const;
        
        void make_QuenchSubScripts(const StructureClass&) const;
        
        void make_EquilibraitionSubScripts(const StructureClass&,
                                           const std::vector<std::vector<double>>) const;
        
        void make_ResizeSubScripts(const StructureClass&,
                                   const std::vector<std::vector<double>>) const;
        
        void make_ProductionSubScripts(const StructureClass&,
                                       const std::vector<std::vector<double>>) const;
        
        
        /* echo the simulation nodes information */
        void echoNodeUsage();
        
        
        
        /*==( public setters )==*/
        
        /* string */
        void set_lmp_exe(const std::string&);
        
        /* bool */
        void set_is_node0(const bool);
        void set_is_node1(const bool);
        
        /* double */
        void set_n_equ_blocks(const double);
        void set_n_prd_blocks(const double);
        void set_prd_exp_base(const double);
        void set_timestep_size(const double);
        void set_time_equ(const double);
        void set_quenchRate(const double);
        void set_steps_gen(const double);
        void set_steps_qch(const double);
        void set_steps_res(const double);
        
        /* vector<int> */
        void set_node0(const std::vector<int>&);
        void set_node1(const std::vector<int>&);
        void set_numCores(const std::vector<int>&);
        void set_cores(const std::vector<int>&);
        void set_run_cores(const std::vector<int>&);
        void set_priority(const std::vector<int>&);
        
        
        
        /*==( public getters )==*/
        
        /* string */
        const std::string get_lmp_exe() const {return lmp_exe;}
        
        /* bool */
        bool get_is_node0() const {return is_node0;}
        bool get_is_node1() const {return is_node1;}
        
        /* int */
        int get_CoresPerNode() const {return CoresPerNode;}
        int get_prd_blocksize() const;
        int get_divide_equ() const {return divide_equ;}
        int get_divide_prd() const {return divide_prd;}
        
        /* double */
        double get_n_equ_blocks() const {return n_equ_blocks;}
        double get_n_prd_blocks() const {return n_prd_blocks;}
        double get_prd_exp_base() const {return prd_exp_base;}
        double get_timestep_size() const {return timestep_size;}
        double get_time_equ() const {return time_equ;}
        double get_quenchRate() const {return quenchRate;}
        
        double get_steps_gen() const {return steps_gen;}
        double get_steps_qch() const {return steps_qch;}
        double get_steps_res() const {return steps_res;}
        double get_steps_equ() const;
        
        /* vector<int> */
        std::vector<int> get_node0() const {return node0;}
        std::vector<int> get_node1() const {return node1;}
        std::vector<int> get_numCores() const {return numCores;}
        std::vector<int> get_cores() const {return cores;}
        std::vector<int> get_run_cores() const {return run_cores;}
        std::vector<int> get_priority() const {return priority;}
    };
    
}
#endif