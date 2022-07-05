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

#ifndef AMDATANALYSIS_H
#define AMDATANALYSIS_H

#include "structureclass.h"
#include "workscripts.h"

namespace autoWork
{
    class AmdatAnalysis
    {
        bool is_check_amdat_version;//AMDAT version consistency check
        bool is_changeAtomOrder;
        bool is_changeAtomType;
        bool is_recoverAtomType;
        bool is_changeFrameSteps;
        bool is_deleteZerothFrames;
        bool is_NPT;
        bool is_fixed_q;
        bool is_fixed_maxL;
        bool is_auto_qvalue;
        bool is_mbodies;
        bool is_keep_new_trajectory;
        bool is_use_displacement_list;
        bool is_use_find_fast;
        bool is_write_pdb_file;
        bool is_write_list_trajectory;
        bool is_use_strFac_frame;
        bool is_rouseModes;
        bool is_use_voroNeighbors;
        bool is_stringmb;
        bool is_streamlined_strings;
        bool is_monoStruct;
        bool is_mb_all_molecule;
        bool is_binning;
        bool is_cusstrfolder;
        
        /** typical AMDAT analyses **/
        bool is_strFac;
        bool is_rdf;
        bool is_msd;
        bool is_ngp;
        bool is_isfs;
        
        /** more AMDAT analyses **/
        bool is_baf;
        bool is_composition;
        bool is_u2dist;
        bool is_stiffness_dist;
        bool is_isf;
        bool is_strings;
        bool is_peak_frame;
        
        int indexi;
        int indexii;
        int amdat_numCores;
        int amdat_run_cores;
        int amdat_priority;
        int fullblock;
        int rdf_nbins;
        int vhs_nbins;
        int u2dist_nbins;
        int stiffness_nbins;
        int prd_blocksize;
        int waveindex;
        int n_poly;
        int chainLen;
        int n_moments;
        int strFac_frame;
        int n_species_types;
        int n_rouseModes;
        
        double prd_exp_base;
        double n_prd_blocks;
        double DWF_time;
        double maxLenScale;
        double vhs_rangebinned;
        double max_u2;
        double max_stiffness;
        double threshold_percentile_lo;
        double threshold_percentile_hi;
        double strings_threshold;
        double threshold_percentile;
        double logtaubeta;
        
        std::string amdat_svn;
        std::string amdat_exe;
        std::string threshold_keyword;
        std::string symmetry;
        std::string geometry;
        std::string relaxation_target;
        std::string centertype;
        std::string species;
        std::string segmode;
        std::string analysispart;
        std::string analysispart_bin;
        std::string cusstrfolder;
        
        std::vector<std::string> speciesName;
        std::vector<std::string> speciesType;
        std::vector<std::vector<int>> n_typeSet;
        std::vector<std::string> mbodyList;
        std::vector<std::string> mbodyType;
        std::vector<std::vector<double>> sigmaMatrix;
        
        /** allocate compute nodes to submitted jobs **/
        void node_allocation(const StructureClass&,const WorkScripts&,const int,std::ofstream&,const int n_trl=0);
        
        /* file strings */
        const std::string get_filenameString();
        
        /** find the peak frame in non-Gaussian parameter (in-block data) **/
        int find_peak_ngp_frame(const StructureClass&,const int,const int,const double);
        
        /** find waveindex of first peak from strFac **/
        int find_waveindex(const StructureClass&,const double);
        
        /** AMDAT string inputfile generation **/
        void systemBlock_definition(std::ofstream&);
        void systemBlock_multibody_allmolecule(const StructureClass&,std::ofstream&);
        void write_amdat_loopstring(std::ofstream&,const StructureClass&,const int rouse=0);//since amdat r128
        void write_amdat_typical(std::ofstream&,const StructureClass&,const int rouse=0);
        void write_amdat_stringInput(std::ofstream&,const StructureClass&,const int frame_index=0);
        
        /** strings **/
        void find_strings_syntax(std::ofstream&,const int);
        void stringmb_syntax(std::ofstream&,const int);
        
        void displacement_value_list(std::ofstream&,const int);
        
    public:
        
        AmdatAnalysis()=default;
        AmdatAnalysis(const StructureClass&,const WorkScripts&);
        AmdatAnalysis(const AmdatAnalysis&)=default;
        AmdatAnalysis& operator= (const AmdatAnalysis&)=default;
        ~AmdatAnalysis()=default;
        
        /** Check for folder existence **/
        void test_AnalysisFolder(const StructureClass&);
        
        /** AMDAT input scripts generation **/
        //----------------------------------------------------------------------
        void make_sigmamatrix(const StructureClass&);
        void make_amdatInputFile(const StructureClass&,const WorkScripts&);
        std::vector<int> make_amdatInputFile_strings_loop(const StructureClass&,const WorkScripts&,const int,const int,const double,const int,const int,const int beg=-1,const int end=-1);
        
        /** AMDAT submissoin generation **/
        //----------------------------------------------------------------------
        void make_amdatSubFile(const StructureClass&,const WorkScripts&,const int,const int,const double,const int frame=0);
        
        /** AMDAT submissoin bash scripts generation **/
        //----------------------------------------------------------------------
        void make_amdatSubScript(const StructureClass&,const std::vector<std::vector<double>>) const;
        
        /** AMDAT target paths **/
        const std::string amdat_target(const StructureClass&,const int,const int,const double,const std::string&) const;
        const std::string amdat_target(const StructureClass&,const int,const int,const double,const std::string&,const int) const;
        const std::string amdat_target(const StructureClass&,const int,const int,const double,const std::string&,const std::string&,const int) const;
        
        /** change the order of atom types in the custom trajectory file **/
        void change_atom_order(const StructureClass&,const int,const int,const double);
        
        /** change atom types in the custom trajectory file **/
        void change_productionAtomTypes(const StructureClass&,const int,const int,const double);
        void recover_productionAtomTypes(const StructureClass&,const int,const int,const double);
        
        /** change frame timesteps in the custom trajectory file **/
        void change_prdcustomfilesteps(const StructureClass&,const int,const int,const double);
        
        /** change frame timesteps in the custom trajectory file **/
        void delete_zerothFrames(const StructureClass&,const int,const int,const double);
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_changeAtomOrder(const bool);
        void set_is_changeAtomType(const bool);
        void set_is_recoverAtomType(const bool);
        void set_is_changeFrameSteps(const bool);
        void set_is_deleteZerothFrames(const bool);
        void set_is_NPT(const bool);
        void set_is_fixed_q(const bool);
        void set_is_fixed_maxL(const bool);
        void set_is_auto_qvalue(const bool);
        void set_is_strFac(const bool);
        void set_is_rdf(const bool);
        void set_is_msd(const bool);
        void set_is_ngp(const bool);
        void set_is_isfs(const bool);
        void set_is_baf(const bool);
        void set_is_composition(const bool);
        void set_is_u2dist(const bool);
        void set_is_stiffness_dist(const bool);
        void set_is_isf(const bool);
        void set_is_mbodies(const bool);
        void set_is_strings(const bool);
        void set_is_peak_frame(const bool);
        void set_is_keep_new_trajectory(const bool);
        void set_is_streamlined_strings(const bool);
        void set_is_monoStruct(const bool);
        void set_is_use_voroNeighbors(const bool);
        void set_is_binning(const bool);
        void set_is_cusstrfolder(const bool);
        /* int */
        void set_amdat_numCores(const int);
        void set_amdat_run_cores(const int);
        void set_amdat_priority(const int);
        void set_fullblock(const int);
        /* double */
        void set_logtaubeta(const double);
        void set_strings_threshold(const double);
        /* string */
        void set_amdat_exe(const std::string&);
        void set_symmetry(const std::string&);
        void set_geometry(const std::string&);
        void set_relaxation_target(const std::string&);
        void set_centertype(const std::string&);
        void set_analysispart(const std::string&);
        void set_analysispart_bin(const std::string&);
        void set_segmode(const std::string&);
        void set_cusstrfolder(const std::string&);
        void set_species(const std::string&);
        /* STL */
        void set_n_typeSet(const std::vector<std::vector<int>>&);
        void set_speciesName(const std::vector<std::string>&);
        void set_speciesType(const std::vector<std::string>&);
        void set_sigmaMatrix(const std::vector<std::vector<double>>&);
        
        /** public getters **/
        //----------------------------------------------------------------------
        /* bool */
        bool get_is_check_amdat_version() const {return is_check_amdat_version;}
        bool get_is_changeAtomOrder() const {return is_changeAtomOrder;}
        bool get_is_changeAtomType() const {return is_changeAtomType;}
        bool get_is_recoverAtomType() const {return is_recoverAtomType;}
        bool get_is_changeFrameSteps() const {return is_changeFrameSteps;}
        bool get_is_deleteZerothFrames() const {return is_deleteZerothFrames;}
        bool get_is_NPT() const {return is_NPT;}
        bool get_is_fixed_q() const {return is_fixed_q;}
        bool get_is_fixed_maxL() const {return is_fixed_maxL;}
        bool get_is_auto_qvalue() const {return is_auto_qvalue;}
        bool get_is_strFac() const {return is_strFac;}
        bool get_is_rdf() const {return is_rdf;}
        bool get_is_msd() const {return is_msd;}
        bool get_is_ngp() const {return is_ngp;}
        bool get_is_isfs() const {return is_isfs;}
        bool get_is_baf() const {return is_baf;}
        bool get_is_composition() const {return is_composition;}
        bool get_is_u2dist() const {return is_u2dist;}
        bool get_is_stiffness_dist() const {return is_stiffness_dist;}
        bool get_is_isf() const {return is_isf;}
        bool get_is_mbodies() const {return is_mbodies;}
        bool get_is_strings() const {return is_strings;}
        bool get_is_peak_frame() const {return is_peak_frame;}
        bool get_is_use_voroNeighbors() const {return is_use_voroNeighbors;}
        bool get_is_stringmb() const {return is_stringmb;}
        bool get_is_streamlined_strings() const {return is_streamlined_strings;}
        bool get_is_monoStruct() const {return is_monoStruct;}
        bool get_is_use_displacement_list() const {return is_use_displacement_list;}
        bool get_is_binning() const {return is_binning;}
        bool get_is_cusstrfolder() const {return is_cusstrfolder;}
        /* int */
        int get_amdat_numCores() const {return amdat_numCores;}
        int get_amdat_run_cores() const {return amdat_run_cores;}
        int get_amdat_priority() const {return amdat_priority;}
        int get_fullblock() const {return fullblock;}
        int get_vhs_nbins() const {return vhs_nbins;}
        int get_n_moments() const {return n_moments;}
        /* double */
        double get_vhs_rangebinned() const {return vhs_rangebinned;}
        double get_strings_threshold() const {return strings_threshold;}
        double get_DWF_time() const {return DWF_time;}
        double get_threshold_percentile() const {return threshold_percentile;}
        double get_logtaubeta() const {return logtaubeta;}
        /* string */
        const std::string get_amdat_svn() const {return amdat_svn;}
        const std::string get_amdat_exe() const {return amdat_exe;}
        const std::string get_symmetry() const {return symmetry;}
        const std::string get_geometry() const {return geometry;}
        const std::string get_segmode() const {return segmode;}
        const std::string get_cusstrfolder() const {return cusstrfolder;}
        const std::string get_analysispart() const {return analysispart;}
        const std::string get_analysispart_bin() const {return analysispart_bin;}
        const std::string get_relaxation_target() const {return relaxation_target;}
        const std::string get_centertype() const {return centertype;}
        const std::string get_species() const {return species;}
        /* STL */
        const std::vector<std::vector<double>>& get_sigmaMatrix() const {return sigmaMatrix;}
        const std::vector<std::string> get_mbodyList() const {return mbodyList;}
        const std::vector<std::string> get_speciesName() const {return speciesName;}
        const std::vector<std::string> get_speciesType() const {return speciesType;}
    };
}
#endif /* AMDATANALYSIS_H */




