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

#ifndef WORKSCRIPTS_H
#define WORKSCRIPTS_H

#include "structureclass.h"
#include "lmpscripts.h"

namespace autoWork
{
    class WorkScripts
    {
        bool is_node0;
        bool is_node1;
        bool is_node2;
        bool is_fixResize;
        
        int indexi,indexii;
        int typeB;
        int typeS;
        int current_regime;
        int n_cores_node;
        int n_gpus_node;
        int prd_blocksize_hiT;
        int prd_blocksize_loT;
        int divide_equ;
        int divide_prd;
        int gpuNum;
        
        double n_equ_blocks;
        double n_prd_blocks;
        double n_relaxation;
        double prd_exp_base;
        double timestep_size;
        double time_equ;
        double quenchRate;
        double steps_gen;
        double steps_qch;
        double steps_equ;
        double steps_res;
        double aratio;
        
        std::string lmp_exe;
        std::string resizeShape;
        
        /** STL containers **/
        std::vector<int> node0;
        std::vector<int> node1;
        std::vector<int> node2;
        std::vector<int> numCores;
        std::vector<int> cores;
        std::vector<int> run_cores;
        std::vector<int> priority;
        
        /** private utility functions **/
        /* SGE environment settings */
        void setup_sge_env(const StructureClass&,
                           const std::string&,
                           const int counter,
                           std::ofstream&);
        
        /* Compute nodes and GPUs allocation */
        void node_allocation(const StructureClass&,
                             const std::string&,
                             const int counter,
                             std::ofstream&,
                             const int n_trl=0);
        
        /* make ad hoc modification to PATH and LD_LIBRARY_PATH */
        void append_adhoc_PATH(std::ofstream&);
        void clean_adhoc_PATH(std::ofstream&);
        
        
    public:
        
        WorkScripts()=default;
        WorkScripts(const StructureClass&);
        WorkScripts(const WorkScripts&)=default;
        WorkScripts& operator= (const WorkScripts&)=default;
        ~WorkScripts()=default;
        
        void make_GenerationSubFiles(StructureClass&,const LmpScripts&,const int,const int,const double);/** non-const: b/c nodeAlloc **/
        void make_QuenchSubFiles(StructureClass&,const int,const int);/** non-const: b/c T changes w/ regime **/
        void make_EquilibrationSubfiles(const StructureClass&,const LmpScripts&,const int,const int,const double);
        void make_ResizeSubfiles(const StructureClass&,const LmpScripts&,const int,const int,const double);
        void make_ProductionSubFiles(StructureClass&,const LmpScripts&,const int,const int,const double);/** non-const: b/c is_set_CheckPoint **/
        
        void make_GenerationSubScripts(const StructureClass&);
        void make_QuenchSubScripts(const StructureClass&);
        void make_EquilibraitionSubScripts(const StructureClass&,const std::vector<std::vector<double>>);
        void make_ResizeSubScripts(const StructureClass&,const std::vector<std::vector<double>>);
        void make_ProductionSubScripts(const StructureClass&,const std::vector<std::vector<double>>);
        
        void   thermodataprocess_prd(const StructureClass&,const int,const int,const double);
        void   thermocalc_prd(const StructureClass&,const int,const int,const double);
        int    get_keyWordIndex(const std::string&,const std::string&);
        std::vector<double> thermocalc_equ(const StructureClass&,const int,const double,const std::string&,const std::string& shape="cubic");
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_node0(const bool);
        void set_is_node1(const bool);
        void set_is_node2(const bool);
        void set_is_fixResize(const bool);
        /* int */
        void set_n_cores_node(const int);
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
        void set_aratio(const double);
        /* string */
        void set_lmp_exe(const std::string&);
        void set_resizeShape(const std::string&);
        /* vector<int> */
        void set_node0(const std::vector<int>&);
        void set_node1(const std::vector<int>&);
        void set_node2(const std::vector<int>&);
        void set_numCores(const std::vector<int>&);
        void set_cores(const std::vector<int>&);
        void set_run_cores(const std::vector<int>&);
        void set_priority(const std::vector<int>&);
        
        /** public getters **/
        //----------------------------------------------------------------------
        /* bool */
        bool get_is_node0() const {return is_node0;}
        bool get_is_node1() const {return is_node1;}
        bool get_is_node2() const {return is_node2;}
        bool get_is_fixResize() const {return is_fixResize;}
        /* int */
        int get_n_cores_node() const {return n_cores_node;}
        int get_prd_blocksize() const;
        int get_divide_equ() const {return divide_equ;}
        int get_divide_prd() const {return divide_prd;}
        /* double */
        double get_n_equ_blocks() const {return n_equ_blocks;}
        double get_n_prd_blocks() const {return n_prd_blocks;}
        double get_prd_exp_base() const {return prd_exp_base;}
        double get_steps_gen() const {return steps_gen;}
        double get_steps_qch() const {return steps_qch;}
        double get_steps_equ(const SysInfo&) const;
        double get_aratio() const {return aratio;}
        double get_quenchRate() const {return quenchRate;}
        double get_timestep_size() const {return timestep_size;}
        double get_time_equ() const {return time_equ;}
        /* string */
        const std::string get_lmp_exe() const {return lmp_exe;}
        const std::string get_resizeShape() const {return resizeShape;}
        /* vector<int> */
        const std::vector<int>& get_node0() const {return node0;}
        const std::vector<int>& get_node1() const {return node1;}
        const std::vector<int>& get_node2() const {return node2;}
        const std::vector<int>& get_numCores() const {return numCores;}
        const std::vector<int>& get_run_cores() const {return run_cores;}
        const std::vector<int>& get_priority() const {return priority;}
    };
}
#endif /* WORKSCRIPTS_H */




